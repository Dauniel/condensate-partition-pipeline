"""
pipeline.py — Spring 2026 Cellpose 3 condensate analysis pipeline.

Usage:
    python pipeline.py --cond <condensate_tif> --nuc <nuclei_tif> [options]

Options:
    --cond        Path to condensate channel Z-stack TIF (required)
    --nuc         Path to nuclei channel Z-stack TIF (required)
    --output      Output directory (default: spring_implementation/outputs)
    --voxel-xy    XY pixel size in µm, e.g. 0.065 (optional; volumes reported in voxels if omitted)
    --voxel-z     Z-slice spacing in µm, e.g. 0.3   (optional)
    --diameter      Cellpose condensate diameter in pixels, None = auto-detect (default: None)
    --nuc-diameter  Cellpose nuclei diameter in pixels, None = auto-detect (default: None)
    --nuc-cellprob  Nuclei cell probability threshold (default: -2, lower = more merging)
    --no-gpu        Disable GPU even if available

Outputs (all written to --output):
    cond_restored.tif               denoised condensate stack
    nuc_restored.tif                denoised nuclei stack
    condensate_masks.tif            3D condensate instance labels
    nuclei_masks.tif                3D nuclei instance labels
    condensate_measurements.csv     per-slice regionprops (label, area, centroid, mean_intensity, z)
    nuclei_measurements.csv         per-slice regionprops
    condensate_volumes.csv          per-object 3D volume (voxels + µm³ if voxel sizes given)
    nuclei_volumes.csv              per-object 3D volume
    summary.csv                     partition coefficient + metadata
    results.png                     summary figure
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile as tiff
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from skimage.measure import regionprops_table, label
from skimage.feature import blob_log
from scipy.ndimage import binary_fill_holes, binary_closing

import torch
from cellpose import models, core, denoise


# ── Per-construct calibration ────────────────────────────────────────────────
# Maps raw pipeline_pc -> Imaris-reference PC using a piecewise-linear
# (isotonic-style) curve per construct. Each entry has 'xs' and 'ys' learned
# from analyze_calibration.py running the trained model with cond_topx=100
# and cellprob=-2. If 'kind' == 'linear' we use slope/intercept instead.
# Empty / missing entries -> no calibration.
import json

_CALIB_PATH = Path(__file__).parent / "outputs" / "calibration_table.json"

def _load_calibration() -> dict:
    if _CALIB_PATH.exists():
        try:
            return json.loads(_CALIB_PATH.read_text())
        except Exception as e:
            print(f"  [warn] failed to load {_CALIB_PATH}: {e}")
    return {}

CALIBRATION = _load_calibration()


def apply_calibration(pipeline_pc: float, construct: str | None) -> tuple[float, bool]:
    """Return (calibrated_pc, was_calibrated)."""
    if construct is None or construct not in CALIBRATION:
        return pipeline_pc, False
    entry = CALIBRATION[construct]
    if entry.get("kind") == "linear":
        return entry["slope"] * pipeline_pc + entry["intercept"], True
    if entry.get("kind") == "isotonic":
        # Piecewise-linear interpolation, clamped at endpoints.
        xs, ys = np.asarray(entry["xs"]), np.asarray(entry["ys"])
        return float(np.interp(pipeline_pc, xs, ys)), True
    return pipeline_pc, False


# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description="Cellpose 3 condensate segmentation + partition coefficient pipeline."
    )
    p.add_argument("--roi",         default=None,   type=Path, help="Multi-channel ROI TIF (Z, 2, Y, X): ch0=nuclei, ch1=condensate")
    p.add_argument("--cond",       default=None,   type=Path, help="Condensate channel TIF (required if --roi not given)")
    p.add_argument("--nuc",        default=None,   type=Path, help="Nuclei channel TIF (required if --roi not given)")
    p.add_argument("--single",     default=None,   type=Path, help="Segment ONE single-channel TIF on its own (no PC). Use with --channel.")
    p.add_argument("--channel",    default="condensate", choices=["nuclei", "condensate"],
                                                              help="What the --single image is (default: condensate)")
    p.add_argument("--output",     default=None,   type=Path, help="Output directory")
    p.add_argument("--voxel-xy",   default=None,   type=float, help="XY pixel size in µm")
    p.add_argument("--voxel-z",    default=None,   type=float, help="Z-slice spacing in µm")
    p.add_argument("--diameter",      default=None,  type=float, help="Cellpose condensate diameter (px)")
    p.add_argument("--nuc-diameter",  default=None,  type=float, help="Nuclei diameter (px), None = auto-detect")
    p.add_argument("--nuc-cellprob",  default=-2.0,  type=float, help="Nuclei cellprob_threshold")
    p.add_argument("--nuc-close",      default=0,     type=int,   help="Per-slice morphological closing iterations on nuclei before hole-fill (0 = off). Seals open bays so condensate voids count as intra-nuclear.")
    p.add_argument("--cond-topx",     default=75.0,  type=float, help="Use mean of top-X%% brightest voxels for cond density (default 75)")
    p.add_argument("--cond-model",    default=None,  type=Path,  help="Fine-tuned Cellpose model path for condensates (default: cyto3)")
    p.add_argument("--cond-cellprob", default=0.0,   type=float, help="cellprob_threshold for condensate seg")
    p.add_argument("--construct",     default=None,  type=str,   help="Construct name for per-construct PC calibration")
    p.add_argument("--detector",      default="auto", choices=["auto", "cellpose", "blob_log"],
                                                                  help="Condensate detector. 'auto' picks blob_log for JABr, Cellpose otherwise (default: auto)")
    p.add_argument("--blob-threshold",default=0.03,  type=float, help="blob_log threshold (default 0.03)")
    p.add_argument("--min-sigma",     default=1.5,   type=float, help="blob_log min sigma")
    p.add_argument("--max-sigma",     default=6.0,   type=float, help="blob_log max sigma")
    p.add_argument("--num-sigma",     default=8,     type=int,   help="blob_log num sigma steps")
    p.add_argument("--no-gpu",        action="store_true",       help="Disable GPU")
    p.add_argument("--view",          action="store_true",       help="Also write an interactive zviewer.html (scroll Z-slices)")
    return p.parse_args()


# ── Step 1: Load ──────────────────────────────────────────────────────────────

def load_stacks(cond_path: Path | None, nuc_path: Path | None, roi_path: Path | None = None):
    """Load condensate and nuclei Z-stacks from TIF files.

    Accepts either:
      - A single multi-channel TIF via roi_path (OME or plain).
        Handles any axis order (CZYX, ZCYX, TCZYX, etc.) by reading OME
        metadata when available. Ch0 = nuclei, Ch1 = condensate.
      - Separate single-channel TIFs via cond_path + nuc_path.
    """
    if roi_path is not None:
        with tiff.TiffFile(roi_path) as tf:
            if tf.is_ome and tf.series:
                data = tf.series[0].asarray()
                axes = tf.series[0].axes.upper()
                print(f"  OME-TIFF detected  axes={axes}  shape={data.shape}")
            else:
                data = tf.asarray()
                axes = None

        # Normalise to (C, Z, Y, X)
        if axes is not None:
            if 'T' in axes:                        # drop time axis (first T)
                data = data.take(0, axis=axes.index('T'))
                axes = axes.replace('T', '')
            if 'C' not in axes:
                raise ValueError(f"No channel axis found in OME axes '{axes}'")
            if axes.index('C') != 0:               # bring C to front
                data = np.moveaxis(data, axes.index('C'), 0)
                axes = 'C' + axes.replace('C', '')
        else:
            # Plain TIF: infer channel axis from shape
            if data.ndim == 4:
                # (Z, C, Y, X) → move C to front
                if data.shape[1] == 2:
                    data = np.moveaxis(data, 1, 0)
                elif data.shape[0] == 2:
                    pass  # already (C, Z, Y, X)
                else:
                    raise ValueError(
                        f"Cannot find channel axis in shape {data.shape}. "
                        "Expected 2 channels."
                    )
            else:
                raise ValueError(f"Expected 4-D array for --roi, got shape {data.shape}")

        if data.shape[0] < 2:
            raise ValueError(f"Expected ≥2 channels, got {data.shape[0]}")

        nuc_stack  = data[0].copy()   # Ch0 = nuclei
        cond_stack = data[1].copy()   # Ch1 = condensate
    else:
        cond_stack = tiff.imread(cond_path)
        nuc_stack  = tiff.imread(nuc_path)

    print(f"Condensate stack : {cond_stack.shape}  dtype={cond_stack.dtype}")
    print(f"Nuclei stack     : {nuc_stack.shape}  dtype={nuc_stack.dtype}")
    return cond_stack, nuc_stack


# ── Step 2: Denoise ───────────────────────────────────────────────────────────

def denoise_stack(stack: np.ndarray, dn_model, label: str) -> np.ndarray:
    """
    Restore a Z-stack with the Cellpose 3 DenoiseModel.

    DenoiseModel.eval expects a list of 2D arrays; returns a list of 2D arrays.
    Removes shot noise before segmentation, which tightens condensate boundaries.
    """
    print(f"  Denoising {label}...")
    restored = dn_model.eval(
        [stack[z] for z in range(stack.shape[0])],
        diameter=None,
        channels=[0, 0],
    )
    return np.stack(restored)


# ── Step 3: Segment ───────────────────────────────────────────────────────────

def segment_condensates(stack: np.ndarray, seg_model, diameter, cellprob_threshold: float = 0.0) -> np.ndarray:
    """Segment condensates with Cellpose 3 in native 3D mode."""
    print(f"  Segmenting condensates (do_3D=True, cellprob_threshold={cellprob_threshold})...")
    masks_3d, _, _ = seg_model.eval(
        stack,
        do_3D=True,
        diameter=diameter,
        cellprob_threshold=cellprob_threshold,
        channels=[0, 0],
    )
    print(f"    condensates: {masks_3d.max()} objects found")
    return masks_3d.astype(np.int32)


# A blob counts as intra-nuclear when at least this fraction of its rendered
# sphere falls inside the nucleus mask. Replaces the older single-centroid-voxel
# test, which rejected condensates sitting on the nuclear boundary (or in voids
# the fill didn't fully close) even when most of their volume was nuclear, and
# in crowded multi-cell fields mis-assigned the target cell. Validated on JABr:
# r 0.926 -> 0.942, LOO-cal MAE 14.6% -> 13.7%, within-20 75% -> 79%.
MIN_NUC_OVERLAP = 0.5


def detect_condensates_blob(cond_stack: np.ndarray, nuc_mask_3d: np.ndarray,
                            threshold: float, min_sigma: float, max_sigma: float,
                            num_sigma: int) -> np.ndarray:
    """Laplacian-of-Gaussian spot detection. One sphere of radius sigma*sqrt(3) per blob.
    Keeps blobs with >= MIN_NUC_OVERLAP of their sphere inside nuc_mask_3d (binary)."""
    print(f"  Segmenting condensates (blob_log, threshold={threshold}, sigma=[{min_sigma},{max_sigma}])...")
    B = float(cond_stack.min())
    norm = np.clip(cond_stack.astype(np.float32) - B, 0, None)
    mx = norm.max()
    if mx <= 0:
        return np.zeros_like(cond_stack, dtype=np.int32)
    norm /= mx

    blobs = blob_log(norm, min_sigma=min_sigma, max_sigma=max_sigma,
                     num_sigma=num_sigma, threshold=threshold, overlap=0.5)
    if len(blobs) == 0:
        print("    blob_log: 0 blobs detected")
        return np.zeros_like(cond_stack, dtype=np.int32)

    Z, Y, X = cond_stack.shape
    mask = np.zeros_like(cond_stack, dtype=np.int32)
    next_id = 1
    for (z, y, x, s) in blobs:
        r = max(1.5, s * np.sqrt(3))
        zi, yi, xi = int(round(z)), int(round(y)), int(round(x))
        rr = int(np.ceil(r))
        z0, z1 = max(0, zi - rr), min(Z, zi + rr + 1)
        y0, y1 = max(0, yi - rr), min(Y, yi + rr + 1)
        x0, x1 = max(0, xi - rr), min(X, xi + rr + 1)
        zz, yy, xx = np.ogrid[z0:z1, y0:y1, x0:x1]
        sphere = (zz - zi) ** 2 + (yy - yi) ** 2 + (xx - xi) ** 2 <= r * r
        n_sphere = int(sphere.sum())
        n_in = int((sphere & (nuc_mask_3d[z0:z1, y0:y1, x0:x1] > 0)).sum())
        if n_sphere == 0 or n_in / n_sphere < MIN_NUC_OVERLAP:
            continue
        sub = mask[z0:z1, y0:y1, x0:x1]
        sub[sphere & (sub == 0)] = next_id
        next_id += 1
    print(f"    blob_log: {mask.max()} condensates (of {len(blobs)} blobs, "
          f">= {MIN_NUC_OVERLAP:.0%} sphere inside nuclei)")
    return mask


def detect_condensates_blob_both(cond_stack: np.ndarray, nuc_mask_3d: np.ndarray,
                                 threshold: float, min_sigma: float, max_sigma: float,
                                 num_sigma: int):
    """blob_log returning two instance masks: (intra-nuclear, cytoplasmic).
    Splits by sphere overlap: nuclear = >= MIN_NUC_OVERLAP of the rendered sphere
    inside nuc_mask_3d, else cytoplasmic. Renders each as a sphere of radius
    sigma*sqrt(3)."""
    print(f"  Segmenting condensates (blob_log split, threshold={threshold})...")
    B = float(cond_stack.min())
    norm = np.clip(cond_stack.astype(np.float32) - B, 0, None)
    mx = norm.max()
    if mx <= 0:
        empty = np.zeros_like(cond_stack, dtype=np.int32)
        return empty, empty.copy()
    norm /= mx
    blobs = blob_log(norm, min_sigma=min_sigma, max_sigma=max_sigma,
                     num_sigma=num_sigma, threshold=threshold, overlap=0.5)
    if len(blobs) == 0:
        empty = np.zeros_like(cond_stack, dtype=np.int32)
        return empty, empty.copy()

    Z, Y, X = cond_stack.shape
    nuc_mask = np.zeros_like(cond_stack, dtype=np.int32)
    cyto_mask = np.zeros_like(cond_stack, dtype=np.int32)
    n_nuc = n_cyto = 0
    for (z, y, x, s) in blobs:
        r = max(1.5, s * np.sqrt(3))
        zi, yi, xi = int(round(z)), int(round(y)), int(round(x))
        rr = int(np.ceil(r))
        z0, z1 = max(0, zi - rr), min(Z, zi + rr + 1)
        y0, y1 = max(0, yi - rr), min(Y, yi + rr + 1)
        x0, x1 = max(0, xi - rr), min(X, xi + rr + 1)
        zz, yy, xx = np.ogrid[z0:z1, y0:y1, x0:x1]
        sphere = (zz - zi) ** 2 + (yy - yi) ** 2 + (xx - xi) ** 2 <= r * r
        n_sphere = int(sphere.sum())
        n_in = int((sphere & (nuc_mask_3d[z0:z1, y0:y1, x0:x1] > 0)).sum())
        is_nuc = n_sphere > 0 and n_in / n_sphere >= MIN_NUC_OVERLAP
        if is_nuc:
            n_nuc += 1
            sub = nuc_mask[z0:z1, y0:y1, x0:x1]
            sub[sphere & (sub == 0)] = n_nuc
        else:
            n_cyto += 1
            sub = cyto_mask[z0:z1, y0:y1, x0:x1]
            sub[sphere & (sub == 0)] = n_cyto
    print(f"    blob_log: {n_nuc} nuclear, {n_cyto} cytoplasmic (of {len(blobs)} total)")
    return nuc_mask, cyto_mask


def segment_nuclei(stack: np.ndarray, seg_model, diameter, cellprob_threshold, close_iters: int = 0) -> np.ndarray:
    """
    Segment nuclei with Cellpose 3 in native 3D mode, then post-process.

    Cellpose over-splits large nuclei into many fragments due to internal
    condensate texture. Post-processing: collapse to a binary mask, re-label
    by 3D connected components, and drop noise fragments < 1000 voxels.
    This gives a clean nucleus count without changing the binary pixel coverage
    (and therefore doesn't affect the PC formula).
    """
    print("  Segmenting nuclei (do_3D=True)...")
    masks_3d, _, _ = seg_model.eval(
        stack,
        do_3D=True,
        diameter=diameter,
        cellprob_threshold=cellprob_threshold,
        channels=[0, 0],
    )
    print(f"    raw Cellpose labels: {masks_3d.max()}")

    # Re-label by connected components and drop noise
    connected = label(masks_3d > 0, connectivity=3)
    sizes     = np.bincount(connected.ravel())
    clean     = np.zeros_like(connected, dtype=np.int32)
    new_lbl   = 0
    for lbl in range(1, connected.max() + 1):
        if sizes[lbl] >= 1000:
            new_lbl += 1
            clean[connected == lbl] = new_lbl
    print(f"    nuclei after connected-component relabeling: {clean.max()}")

    # Fill interior voids per nucleus. Condensates exclude the nuclear stain, so
    # Cellpose carves donut holes around them; those holes sit exactly where
    # nuclear condensates are. Filling per-label (3D + per-slice 2D for z-open
    # tunnels) makes condensates in those voids count as intra-nuclear in the PC.
    filled = np.zeros_like(clean)
    for lbl in range(1, clean.max() + 1):
        m = clean == lbl
        if close_iters > 0:
            # Seal narrow openings ("bays") so they become enclosed holes the
            # fill can close. Per-slice 2D closing avoids bridging across z.
            for z in range(m.shape[0]):
                m[z] = binary_closing(m[z], iterations=close_iters)
        m = binary_fill_holes(m)
        for z in range(m.shape[0]):
            m[z] = binary_fill_holes(m[z])
        filled[m] = lbl
    n_added = int((filled > 0).sum() - (clean > 0).sum())
    label_txt = f" (close={close_iters})" if close_iters else ""
    print(f"    nuclei after void-filling{label_txt}: +{n_added} voxels")
    return filled


# ── Step 4: Per-slice measurements ───────────────────────────────────────────

def extract_slice_measurements(masks_3d: np.ndarray, raw_stack: np.ndarray) -> pd.DataFrame:
    """
    Compute per-slice regionprops from a 3D instance label volume.

    Uses the raw (undenoised) stack for intensity values so measurements
    reflect true signal, not denoising-modified values.

    Returns a DataFrame with columns: label, area, centroid-0, centroid-1,
    mean_intensity, z.
    """
    rows = []
    for z in range(raw_stack.shape[0]):
        props = regionprops_table(
            masks_3d[z],
            intensity_image=raw_stack[z],
            properties=["label", "area", "centroid", "mean_intensity"],
        )
        df = pd.DataFrame(props)
        if not df.empty:
            df["z"] = z
            rows.append(df)
    if rows:
        return pd.concat(rows, ignore_index=True)
    return pd.DataFrame(columns=["label", "area", "centroid-0", "centroid-1", "mean_intensity", "z"])


# ── Step 5: 3D volume estimation ──────────────────────────────────────────────

def compute_volumes(masks_3d: np.ndarray, voxel_xy: float | None, voxel_z: float | None) -> pd.DataFrame:
    """
    Compute per-object 3D volume from a 3D instance label volume.

    Counts voxels per label (label 0 = background, excluded). If voxel_xy and
    voxel_z are provided, also reports physical volume in µm³.

    Returns a DataFrame with columns: label, volume_voxels[, volume_um3].
    """
    props = regionprops_table(masks_3d, properties=["label", "area"])
    df = pd.DataFrame(props).rename(columns={"area": "volume_voxels"})

    if voxel_xy is not None and voxel_z is not None:
        df["volume_um3"] = df["volume_voxels"] * (voxel_xy ** 2) * voxel_z
        print(f"    Voxel size: {voxel_xy} µm (XY) × {voxel_z} µm (Z)")

    return df


# ── Step 6: Partition coefficient ─────────────────────────────────────────────

def compute_partition_coefficient(
    cond_stack: np.ndarray,
    cond_masks_3d: np.ndarray,
    nuc_masks_3d: np.ndarray,
    cond_topx: float = 75.0,
) -> dict:
    """
    Compute the nuclear partition coefficient using the Fabrini et al. method.

    B = minimum voxel intensity across the full FOV (camera background offset).
    Condensed density  = mean clip(pixel - B, 0) over the top-`cond_topx` percent
                         brightest voxels inside (condensate mask AND nucleus
                         mask). The Cellpose 3 do_3D=True mask is more inclusive
                         than a manually-drawn Imaris mask — it adds a dim halo
                         and dark interior pixels around each condensate core.
                         Trimming the bottom 25% of mask pixels (cond_topx=75)
                         removes that fluff and recovers densities that match
                         the manual reference (validated: -16% bias -> +3% bias
                         across 29 JABr cells; see batch_sweep_topx.py).
    Dilute density     = mean of the 50 lowest-intensity valid 10×10×10 patches
                         fully within the nuclear dilute region. Sorting by
                         intensity approximates manual selection of a quiet
                         representative region. Falls back to all dilute voxels
                         if no valid patches are found.
    PC = condensed density / dilute density.

    Returns a dict with keys: pc, background, cond_density, dilute_density.
    """
    B       = float(cond_stack.min())
    cond_3d = cond_masks_3d > 0
    nuc_3d  = nuc_masks_3d  > 0

    # Condensed phase — top-X% brightest voxels in (cond mask & nuc mask).
    nuclear_cond = cond_3d & nuc_3d
    cond_vals    = np.clip(cond_stack[nuclear_cond].astype(np.float64) - B, 0, None)
    cond_vals_sorted = np.sort(cond_vals)
    cutoff       = int(len(cond_vals_sorted) * (1.0 - cond_topx / 100.0))
    cond_density = float(cond_vals_sorted[cutoff:].mean())

    # Dilute phase — mean of the N_PATCHES lowest-intensity valid 10×10×10 patches.
    # Sorting by intensity (ascending) approximates how a researcher manually picks
    # a quiet representative dilute-phase region, giving a stable and reproducible
    # estimate without relying on random sampling.
    dilute_3d  = nuc_3d & ~cond_3d
    PATCH      = 10
    N_PATCHES  = 50
    Z, Y, X    = cond_stack.shape
    candidates = np.argwhere(dilute_3d)
    in_bounds  = candidates[
        (candidates[:, 0] + PATCH <= Z) &
        (candidates[:, 1] + PATCH <= Y) &
        (candidates[:, 2] + PATCH <= X)
    ]

    patch_means = []
    for z0, y0, x0 in in_bounds:
        if dilute_3d[z0:z0+PATCH, y0:y0+PATCH, x0:x0+PATCH].all():
            patch = cond_stack[z0:z0+PATCH, y0:y0+PATCH, x0:x0+PATCH].astype(np.float64) - B
            patch_means.append(np.clip(patch, 0, None).mean())

    if patch_means:
        patch_means.sort()
        dilute_density = float(np.mean(patch_means[:N_PATCHES]))
    else:
        dilute_density = np.clip(
            cond_stack[dilute_3d].astype(np.float64) - B, 0, None
        ).mean()

    pc = cond_density / dilute_density
    return {"pc": pc, "background": B, "cond_density": cond_density, "dilute_density": dilute_density}


def compute_cytoplasmic_partition_coefficient(
    cond_stack: np.ndarray,
    cyto_cond_masks_3d: np.ndarray,
    nuc_masks_3d: np.ndarray,
    cond_topx: float = 75.0,
) -> dict:
    """Cytoplasmic PC, Fabrini methodology.

    Condensate mask: cytoplasmic condensates (centroids outside nuclei).
    Dilute pool: voxels outside any nucleus AND outside any condensate, restricted
    to "signal-bearing" cytoplasm so empty background space doesn't sneak in.
    Signal threshold = 25th percentile of cond_stack (matches Fabrini's
    intensity-above-background definition; not anatomy-based).

    Returns: pc, background, cond_density, dilute_density (None if no signal-
    bearing cytoplasmic region exists, e.g. all blobs were nuclear)."""
    B       = float(cond_stack.min())
    cond_3d = cyto_cond_masks_3d > 0
    nuc_3d  = nuc_masks_3d > 0

    if not cond_3d.any():
        return {"pc": float("nan"), "background": B, "cond_density": 0.0,
                "dilute_density": float("nan")}

    # Condensed phase — top-X% brightest voxels in the cytoplasmic-cond mask.
    cond_vals = np.clip(cond_stack[cond_3d].astype(np.float64) - B, 0, None)
    cond_vals_sorted = np.sort(cond_vals)
    cutoff = int(len(cond_vals_sorted) * (1.0 - cond_topx / 100.0))
    cond_density = float(cond_vals_sorted[cutoff:].mean()) if len(cond_vals_sorted) > cutoff else float(cond_vals_sorted.mean())

    # Cytoplasmic dilute = outside nucleus, outside condensate, above signal floor.
    signal_floor = float(np.percentile(cond_stack, 25))
    dilute_3d = (~nuc_3d) & (~cond_3d) & (cond_stack > signal_floor)
    PATCH, N_PATCHES = 10, 50
    Z, Y, X = cond_stack.shape
    candidates = np.argwhere(dilute_3d)
    in_bounds = candidates[
        (candidates[:, 0] + PATCH <= Z) &
        (candidates[:, 1] + PATCH <= Y) &
        (candidates[:, 2] + PATCH <= X)
    ]
    patch_means = []
    for z0, y0, x0 in in_bounds:
        sub = dilute_3d[z0:z0+PATCH, y0:y0+PATCH, x0:x0+PATCH]
        if sub.all():
            patch = cond_stack[z0:z0+PATCH, y0:y0+PATCH, x0:x0+PATCH].astype(np.float64) - B
            patch_means.append(np.clip(patch, 0, None).mean())

    if patch_means:
        patch_means.sort()
        dilute_density = float(np.mean(patch_means[:N_PATCHES]))
    elif dilute_3d.any():
        dilute_density = float(np.clip(
            cond_stack[dilute_3d].astype(np.float64) - B, 0, None
        ).mean())
    else:
        return {"pc": float("nan"), "background": B,
                "cond_density": cond_density, "dilute_density": float("nan")}

    pc = cond_density / dilute_density if dilute_density > 0 else float("nan")
    return {"pc": pc, "background": B, "cond_density": cond_density,
            "dilute_density": dilute_density, "signal_floor": signal_floor}


# ── Step 7: Save outputs ──────────────────────────────────────────────────────

def save_outputs(
    output_dir: Path,
    cond_restored, nuc_restored,
    cond_masks_3d, nuc_masks_3d,
    cond_df, nuc_df,
    cond_vol_df, nuc_vol_df,
    pc_result: dict,
    voxel_xy, voxel_z,
):
    """Save all masks, tables, and summary CSV to output_dir."""
    tiff.imwrite(output_dir / "cond_restored.tif",    cond_restored)
    tiff.imwrite(output_dir / "nuc_restored.tif",     nuc_restored)
    tiff.imwrite(output_dir / "condensate_masks.tif", cond_masks_3d)
    tiff.imwrite(output_dir / "nuclei_masks.tif",     nuc_masks_3d)

    cond_df.to_csv(output_dir / "condensate_measurements.csv", index=False)
    nuc_df.to_csv(output_dir  / "nuclei_measurements.csv",     index=False)
    cond_vol_df.to_csv(output_dir / "condensate_volumes.csv",  index=False)
    nuc_vol_df.to_csv(output_dir  / "nuclei_volumes.csv",      index=False)

    summary_rows = [
        ("partition_coefficient", pc_result["pc"]),
        ("background",            pc_result["background"]),
        ("condensate_density",    pc_result["cond_density"]),
        ("dilute_density",        pc_result["dilute_density"]),
        ("n_condensates",         int(cond_vol_df["label"].nunique())),
        ("n_nuclei",              int(nuc_vol_df["label"].nunique())),
        ("voxel_xy_um",           voxel_xy if voxel_xy else float("nan")),
        ("voxel_z_um",            voxel_z  if voxel_z  else float("nan")),
    ]
    if "pc_calibrated" in pc_result:
        summary_rows.append(("pc_calibrated", pc_result["pc_calibrated"]))
        summary_rows.append(("construct",     pc_result["construct"]))
    pd.DataFrame(summary_rows, columns=["metric", "value"]).to_csv(
        output_dir / "summary.csv", index=False
    )


# ── Step 8: Visualisation ─────────────────────────────────────────────────────

_COND_COLOR = "#2e8b57"   # green for condensates
_NUC_COLOR  = "#4169e1"   # blue for nuclei
_REFERENCE  = 6.32


def plot_summary(
    output_dir: Path,
    cond_df, nuc_df,
    cond_vol_df, nuc_vol_df,
    pc_result: dict,
):
    """
    Save a 3×3 summary figure to output_dir/results.png.

    Layout:
      Row 0: [PC scorecard] [Objects per Z-slice ── spans 2 cols ──────────]
      Row 1: [Condensate Area] [Condensate Intensity] [Condensate 3D Volume]
      Row 2: [Nuclei Area]     [Nuclei Intensity]     [Nuclei 3D Volume]
    """
    fig = plt.figure(figsize=(15, 10))
    gs  = GridSpec(3, 3, figure=fig, hspace=0.45, wspace=0.35)

    ax_pc  = fig.add_subplot(gs[0, 0])
    ax_obj = fig.add_subplot(gs[0, 1:])
    ax_ca  = fig.add_subplot(gs[1, 0])
    ax_ci  = fig.add_subplot(gs[1, 1])
    ax_cv  = fig.add_subplot(gs[1, 2])
    ax_na  = fig.add_subplot(gs[2, 0])
    ax_ni  = fig.add_subplot(gs[2, 1])
    ax_nv  = fig.add_subplot(gs[2, 2])

    # ── PC scorecard ──────────────────────────────────────────────────────────
    pc_val = pc_result["pc"]
    ax_pc.axis("off")
    ax_pc.text(0.5, 0.96, "Partition Coefficient",
               ha="center", va="top", fontsize=12, fontweight="bold",
               transform=ax_pc.transAxes, color="#333333")
    ax_pc.text(0.5, 0.70, f"{pc_val:.3f}",
               ha="center", va="center", fontsize=42, fontweight="bold",
               color=_COND_COLOR, transform=ax_pc.transAxes)

    bar = ax_pc.inset_axes([0.05, 0.05, 0.90, 0.32])
    bar.barh(1, pc_val,     color=_COND_COLOR, height=0.5)
    bar.barh(0, _REFERENCE, color="#cc3333",   height=0.5)
    bar.set_yticks([0, 1])
    bar.set_yticklabels(
        [f"Reference  {_REFERENCE:.2f}", f"Pipeline  {pc_val:.3f}"],
        fontsize=8.5,
    )
    bar.set_xlim(0, max(pc_val, _REFERENCE) * 1.25)
    bar.xaxis.set_visible(False)
    bar.spines[["top", "right", "left"]].set_visible(False)

    # ── Objects per Z-slice ───────────────────────────────────────────────────
    cond_counts = cond_df.groupby("z")["label"].count()
    nuc_counts  = nuc_df.groupby("z")["label"].count()
    ax_obj.plot(cond_counts.index, cond_counts.values, color=_COND_COLOR, label="Condensates")
    ax_obj.plot(nuc_counts.index,  nuc_counts.values,  color=_NUC_COLOR,  label="Nuclei")
    ax_obj.set_title("Objects per Z-slice")
    ax_obj.set_xlabel("Z-slice")
    ax_obj.set_ylabel("Count")
    ax_obj.legend()

    # ── Condensate row ────────────────────────────────────────────────────────
    ax_ca.hist(cond_df["area"],           bins=40, color=_COND_COLOR)
    ax_ca.set_title("Condensate Area (px²)")
    ax_ca.set_xlabel("Area")
    ax_ca.set_ylabel("Count")

    ax_ci.hist(cond_df["mean_intensity"], bins=40, color=_COND_COLOR)
    ax_ci.set_title("Condensate Intensity")
    ax_ci.set_xlabel("Mean Intensity")

    ax_cv.hist(cond_vol_df["volume_voxels"], bins=40, color=_COND_COLOR)
    ax_cv.set_title("Condensate 3D Volume (voxels)")
    ax_cv.set_xlabel("Volume (voxels)")

    # ── Nuclei row ────────────────────────────────────────────────────────────
    ax_na.hist(nuc_df["area"],           bins=40, color=_NUC_COLOR)
    ax_na.set_title("Nuclei Area (px²)")
    ax_na.set_xlabel("Area")
    ax_na.set_ylabel("Count")

    ax_ni.hist(nuc_df["mean_intensity"], bins=40, color=_NUC_COLOR)
    ax_ni.set_title("Nuclei Intensity")
    ax_ni.set_xlabel("Mean Intensity")

    ax_nv.hist(nuc_vol_df["volume_voxels"], bins=40, color=_NUC_COLOR)
    ax_nv.set_title("Nuclei 3D Volume (voxels)")
    ax_nv.set_xlabel("Volume (voxels)")

    plt.savefig(output_dir / "results.png", dpi=150, bbox_inches="tight")
    plt.close()


# ── Step 9: Interactive Z-stack viewer ────────────────────────────────────────

def _stretch(slice2d: np.ndarray, lo: float, hi: float) -> np.ndarray:
    """Contrast-stretch a 2D slice to [0, 1] using precomputed stack percentiles."""
    if hi <= lo:
        return np.zeros_like(slice2d, dtype=float)
    return np.clip((slice2d.astype(float) - lo) / (hi - lo), 0.0, 1.0)


def _slice_panels_png(
    z: int,
    cond_stack, nuc_stack,
    cond_masks_3d, nuc_masks_3d,
    cond_lo, cond_hi, nuc_lo, nuc_hi,
) -> str:
    """Render the 4-panel row for one Z-slice and return it as a base64 PNG string."""
    import base64, io
    from skimage.segmentation import find_boundaries
    from scipy.ndimage import binary_dilation

    nuc_n = _stretch(nuc_stack[z],  nuc_lo,  nuc_hi)
    cond_n = _stretch(cond_stack[z], cond_lo, cond_hi)

    # Raw composite: nuclei -> blue, condensate -> green
    composite = np.zeros((*nuc_n.shape, 3), dtype=float)
    composite[..., 1] = cond_n   # green
    composite[..., 2] = nuc_n    # blue

    cond_m = cond_masks_3d[z] > 0
    nuc_m  = nuc_masks_3d[z]  > 0

    # Overlay: dim the raw so masks pop, lay translucent fills under thick,
    # high-contrast outlines (magenta = nuclei, yellow = condensate).
    NUC_RGB  = np.array([1.0, 0.2, 0.9])   # magenta
    COND_RGB = np.array([1.0, 0.95, 0.1])  # yellow
    overlay = composite * 0.45             # dim background to ~45%

    def _blend(img, mask, rgb, alpha):
        img[mask] = (1 - alpha) * img[mask] + alpha * rgb

    _blend(overlay, nuc_m,  NUC_RGB,  0.18)
    _blend(overlay, cond_m, COND_RGB, 0.30)

    # Thicken boundaries to ~2 px so they read at a glance.
    nuc_b  = binary_dilation(find_boundaries(nuc_m,  mode="outer"))
    cond_b = binary_dilation(find_boundaries(cond_m, mode="outer"))
    overlay[nuc_b]  = NUC_RGB
    overlay[cond_b] = COND_RGB
    overlay = np.clip(overlay, 0, 1)

    cond_panel = np.zeros((*cond_m.shape, 3), dtype=float)
    cond_panel[..., 1] = cond_m.astype(float)
    nuc_panel = np.zeros((*nuc_m.shape, 3), dtype=float)
    nuc_panel[..., 2] = nuc_m.astype(float)

    panels = [
        (composite,   f"Raw  (nuc=blue, cond=green)"),
        (overlay,     f"Masks on raw  (nuc=magenta, cond=yellow)"),
        (cond_panel,  "Condensate mask"),
        (nuc_panel,   "Nuclei mask"),
    ]

    fig, axes = plt.subplots(1, 4, figsize=(16, 4.2))
    for ax, (img, title) in zip(axes, panels):
        ax.imshow(img)
        ax.set_title(title, fontsize=10)
        ax.axis("off")
    fig.tight_layout()

    buf = io.BytesIO()
    # JPEG keeps the file small: the raw/overlay panels are photographic and
    # compress ~5x better than PNG; the flat mask panels survive at high quality.
    fig.savefig(buf, format="jpg", dpi=90, bbox_inches="tight", pil_kwargs={"quality": 85})
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode("ascii")


def _fmt(v, nd=3):
    """Format a number for the viewer, or '—' for None/NaN."""
    if v is None:
        return "—"
    try:
        if isinstance(v, float) and np.isnan(v):
            return "—"
    except TypeError:
        pass
    if isinstance(v, float):
        return f"{v:.{nd}f}"
    return str(v)


def write_zviewer(
    output_dir: Path,
    cond_stack, nuc_stack,
    cond_masks_3d, nuc_masks_3d,
    pc_result: dict | None = None,
    voxel_xy=None, voxel_z=None,
):
    """
    Write a self-contained zviewer.html into output_dir: a Z-slider over four
    panels (raw composite, masks-on-raw, condensate mask, nuclei mask), a
    stack-level summary (PC, densities, object counts, dimensions), and a live
    per-slice readout. All slice images are embedded as base64 PNGs so the file
    is fully portable.
    """
    n_z, ny, nx = cond_stack.shape
    cond_lo, cond_hi = np.percentile(cond_stack, [2, 99.5])
    nuc_lo,  nuc_hi  = np.percentile(nuc_stack,  [2, 99.5])
    pc = pc_result or {}

    # ── Shared intensity-histogram bins (in-mask voxels, both channels) ─────────
    cond_all = cond_stack[cond_masks_3d > 0]
    nuc_all  = nuc_stack[nuc_masks_3d  > 0]
    pooled   = np.concatenate([cond_all, nuc_all]) if (cond_all.size + nuc_all.size) else np.array([1.0])
    hist_hi  = float(np.percentile(pooled, 99.5)) or 1.0
    NBINS    = 40
    edges    = np.linspace(0, hist_hi, NBINS + 1)
    centers  = ((edges[:-1] + edges[1:]) / 2).round(1).tolist()

    def _hist(vals):
        return np.histogram(vals, bins=edges)[0].astype(int).tolist() if vals.size else [0] * NBINS

    # ── Per-slice stats (computed from masks + raw stacks) ──────────────────────
    print(f"    Rendering {n_z} Z-slices for viewer...")
    slices = []
    for z in range(n_z):
        cm = cond_masks_3d[z] > 0
        nm = nuc_masks_3d[z]  > 0
        cond_ct = int(len(np.unique(cond_masks_3d[z])) - 1)  # minus background
        nuc_ct  = int(len(np.unique(nuc_masks_3d[z]))  - 1)
        png = _slice_panels_png(
            z, cond_stack, nuc_stack, cond_masks_3d, nuc_masks_3d,
            cond_lo, cond_hi, nuc_lo, nuc_hi,
        )
        slices.append({
            "z": z,
            "cond": cond_ct,
            "nuc": nuc_ct,
            "cond_area": int(cm.sum()),
            "nuc_area":  int(nm.sum()),
            "cond_int": round(float(cond_stack[z][cm].mean()), 1) if cm.any() else None,
            "nuc_int":  round(float(nuc_stack[z][nm].mean()),  1) if nm.any() else None,
            "hcond": _hist(cond_stack[z][cm]),
            "hnuc":  _hist(nuc_stack[z][nm]),
            "png": png,
        })

    # ── Stack-level summary ─────────────────────────────────────────────────────
    n_cond_obj = int(len(np.unique(cond_masks_3d)) - 1)
    n_nuc_obj  = int(len(np.unique(nuc_masks_3d))  - 1)
    summary = {
        "pc":            pc.get("pc"),
        "pc_calibrated": pc.get("pc_calibrated"),
        "pc_cytoplasmic": pc.get("pc_cytoplasmic"),
        "construct":     pc.get("construct"),
        "background":    pc.get("background"),
        "cond_density":  pc.get("cond_density"),
        "dilute_density": pc.get("dilute_density"),
        "n_cond_obj":    n_cond_obj,
        "n_nuc_obj":     n_nuc_obj,
        "dims":          f"{n_z} × {ny} × {nx}",
        "voxel_xy":      voxel_xy,
        "voxel_z":       voxel_z,
        "hist_centers":  centers,
        "hist_cond":     _hist(cond_all),
        "hist_nuc":      _hist(nuc_all),
    }
    # Pre-format display strings so the template stays logic-light.
    summary["pc_main"] = _fmt(summary["pc_calibrated"] if summary["pc_calibrated"] is not None else summary["pc"])
    summary["pc_sub"]  = (f"raw {_fmt(summary['pc'])}"
                          + (f" · {summary['construct']} calibrated" if summary["pc_calibrated"] is not None else " (uncalibrated)"))

    html = (_ZVIEWER_TEMPLATE
            .replace("__N_Z__", str(n_z))
            .replace("__SUMMARY__", json.dumps(summary))
            .replace("__SLICES__", json.dumps(slices)))
    out = output_dir / "zviewer.html"
    out.write_text(html, encoding="utf-8")
    print(f"    Z-viewer written to: {out}")


_ZVIEWER_TEMPLATE = """<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<title>Condensate pipeline — Z-stack viewer</title>
<style>
  :root { --cond:#f5f21a; --nuc:#ff33e6; --accent:#3ad07a; }
  * { box-sizing: border-box; }
  body { font-family: -apple-system, Segoe UI, Roboto, sans-serif; margin: 0;
         background: #0e0e12; color: #eee; }
  header { padding: 14px 20px; background: #16161d; border-bottom: 1px solid #2a2a33; }
  h1 { font-size: 16px; margin: 0; font-weight: 600; }
  .wrap { display: flex; gap: 16px; padding: 16px; align-items: flex-start;
          flex-wrap: wrap; }
  .side { flex: 0 0 260px; display: flex; flex-direction: column; gap: 14px; }
  .main { flex: 1 1 640px; min-width: 360px; }
  .card { background: #16161d; border: 1px solid #2a2a33; border-radius: 10px;
          padding: 14px 16px; }
  .card h2 { margin: 0 0 10px; font-size: 12px; letter-spacing: .06em;
             text-transform: uppercase; color: #8a8a99; font-weight: 600; }
  .pc-big { font-size: 46px; font-weight: 700; color: var(--accent); line-height: 1; }
  .pc-sub { font-size: 12px; color: #9aa; margin-top: 6px; }
  table { width: 100%; border-collapse: collapse; font-size: 13px; }
  td { padding: 4px 0; vertical-align: top; }
  td.k { color: #9aa; }
  td.v { text-align: right; font-variant-numeric: tabular-nums; font-weight: 600; }
  .dot { display:inline-block; width:9px; height:9px; border-radius:2px;
         margin-right:6px; vertical-align: baseline; }
  .frame img { width: 100%; border-radius: 8px; background:#000; display:block; }
  .zlabel { text-align:center; font-size:14px; color:#cde; margin: 10px 0 6px;
            font-variant-numeric: tabular-nums; }
  input[type=range] { width: 100%; }
  .hint { font-size: 12px; color: #777; margin-top: 6px; text-align:center; }
  .hist-row { display:flex; gap:16px; margin-top:16px; flex-wrap:wrap; }
  .hist-row .card { flex:1 1 300px; }
  canvas { width:100%; height:150px; display:block; }
  .legend { font-size:11px; color:#9aa; margin-top:6px; }
</style></head>
<body>
<header>
  <h1>Condensate Partition Pipeline — Z-stack viewer</h1>
</header>
<div class="wrap">
  <div class="side">
    <div class="card">
      <h2>Partition Coefficient</h2>
      <div class="pc-big" id="pcBig">—</div>
      <div class="pc-sub" id="pcSub"></div>
    </div>
    <div class="card">
      <h2>Stack summary</h2>
      <table id="stackTbl"></table>
    </div>
    <div class="card">
      <h2>This slice</h2>
      <table id="sliceTbl"></table>
    </div>
  </div>
  <div class="main">
    <div class="card frame">
      <img id="img" alt="z-slice">
      <div class="zlabel" id="zlabel"></div>
      <input type="range" id="slider" min="0" max="__N_Z__" value="0" step="1">
      <div class="hint">Drag the slider or use &larr; / &rarr; arrow keys</div>
    </div>
    <div class="hist-row">
      <div class="card">
        <h2 id="sliceHistTitle">This slice — in-mask intensity</h2>
        <canvas id="sliceHist" width="600" height="150"></canvas>
        <div class="legend">
          <span class="dot" style="background:var(--cond)"></span>Condensate ch
          &nbsp;&nbsp;<span class="dot" style="background:var(--nuc)"></span>Nuclei ch
          &nbsp;— x: intensity, y: voxel count
        </div>
      </div>
      <div class="card">
        <h2>Whole stack — in-mask intensity</h2>
        <canvas id="overallHist" width="600" height="150"></canvas>
        <div class="legend">
          <span class="dot" style="background:var(--cond)"></span>Condensate ch
          &nbsp;&nbsp;<span class="dot" style="background:var(--nuc)"></span>Nuclei ch
          &nbsp;— pooled over all slices
        </div>
      </div>
    </div>
  </div>
</div>
<script>
  const SLICES = __SLICES__;
  const SUM = __SUMMARY__;
  const img = document.getElementById('img');
  const slider = document.getElementById('slider');
  slider.max = SLICES.length - 1;

  const COND = '<span class="dot" style="background:var(--cond)"></span>';
  const NUC  = '<span class="dot" style="background:var(--nuc)"></span>';
  const row = (k, v) => `<tr><td class="k">${k}</td><td class="v">${v}</td></tr>`;

  // Static PC + stack summary
  document.getElementById('pcBig').textContent = SUM.pc_main;
  document.getElementById('pcSub').textContent = SUM.pc_sub;
  const f = (x, d=2) => (x === null || x === undefined) ? '—'
                        : (typeof x === 'number' ? x.toFixed(d) : x);
  document.getElementById('stackTbl').innerHTML =
      row('Cytoplasmic PC', f(SUM.pc_cytoplasmic, 3))
    + row('Cond. density',  f(SUM.cond_density, 1))
    + row('Dilute density', f(SUM.dilute_density, 1))
    + row('Background',      f(SUM.background, 1))
    + row(COND + 'Condensates', SUM.n_cond_obj)
    + row(NUC  + 'Nuclei',      SUM.n_nuc_obj)
    + row('Dimensions (Z×Y×X)', SUM.dims)
    + (SUM.voxel_xy ? row('Voxel XY (µm)', f(SUM.voxel_xy, 3)) : '')
    + (SUM.voxel_z  ? row('Voxel Z (µm)',  f(SUM.voxel_z, 3))  : '');

  const sliceTbl = document.getElementById('sliceTbl');
  const zlabel = document.getElementById('zlabel');

  // ── Histogram drawing (overlaid condensate + nuclei) ───────────────────────
  const CSS = getComputedStyle(document.documentElement);
  const COND_C = CSS.getPropertyValue('--cond').trim();
  const NUC_C  = CSS.getPropertyValue('--nuc').trim();

  function drawHist(canvas, hc, hn, ymax) {
    const ctx = canvas.getContext('2d');
    const W = canvas.width, H = canvas.height, pad = 4;
    ctx.clearRect(0, 0, W, H);
    const n = hc.length, bw = (W - 2 * pad) / n;
    const top = ymax || Math.max(1, ...hc, ...hn);
    const bar = (arr, color) => {
      ctx.fillStyle = color; ctx.globalAlpha = 0.55;
      for (let i = 0; i < n; i++) {
        const h = (arr[i] / top) * (H - 2 * pad);
        ctx.fillRect(pad + i * bw, H - pad - h, Math.max(1, bw - 0.5), h);
      }
      ctx.globalAlpha = 1;
    };
    bar(hn, NUC_C);
    bar(hc, COND_C);
  }

  // Overall histogram: fixed y-scale, drawn once.
  const overallMax = Math.max(1, ...SUM.hist_cond, ...SUM.hist_nuc);
  drawHist(document.getElementById('overallHist'), SUM.hist_cond, SUM.hist_nuc, overallMax);
  const sliceCanvas = document.getElementById('sliceHist');
  const sliceHistTitle = document.getElementById('sliceHistTitle');

  function show(i) {
    i = Math.max(0, Math.min(SLICES.length - 1, i));
    const s = SLICES[i];
    img.src = 'data:image/jpeg;base64,' + s.png;
    zlabel.textContent = `Z ${s.z + 1} / ${SLICES.length}`;
    sliceTbl.innerHTML =
        row(COND + 'Condensates', s.cond)
      + row(NUC  + 'Nuclei',      s.nuc)
      + row(COND + 'Cond. area (px)', s.cond_area)
      + row(NUC  + 'Nuc. area (px)',  s.nuc_area)
      + row(COND + 'Cond. mean int.', f(s.cond_int, 1))
      + row(NUC  + 'Nuc. mean int.',  f(s.nuc_int, 1));
    drawHist(sliceCanvas, s.hcond, s.hnuc);
    sliceHistTitle.textContent = `This slice (Z ${s.z + 1}) — in-mask intensity`;
    slider.value = i;
  }
  slider.addEventListener('input', e => show(+e.target.value));
  document.addEventListener('keydown', e => {
    if (e.key === 'ArrowLeft')  show(+slider.value - 1);
    if (e.key === 'ArrowRight') show(+slider.value + 1);
  });
  show(0);
</script>
</body></html>
"""


# ── Main ──────────────────────────────────────────────────────────────────────

def run_single_channel(args, use_gpu, output_dir):
    """Segment a single-channel image (nuclei OR condensate) on its own.

    No partition coefficient — that requires both channels. Outputs the 3D
    instance mask, object count, per-object volumes, and per-slice measurements.
    Nuclei use Cellpose cyto3; condensate uses the same detector the full
    pipeline would (construct/--detector), with no intra-nuclear gate.
    """
    print(f"\nSingle-channel mode: {args.channel}")
    stack = tiff.imread(args.single)
    print(f"  stack: {stack.shape}  dtype={stack.dtype}")

    print("\n[1/3] Denoising...")
    dn_model = denoise.DenoiseModel(model_type="denoise_cyto3", gpu=use_gpu)
    restored = denoise_stack(stack, dn_model, args.channel)

    print("\n[2/3] Segmenting...")
    if args.channel == "nuclei":
        seg = models.CellposeModel(gpu=use_gpu, model_type="cyto3")
        masks = segment_nuclei(restored, seg, args.nuc_diameter, args.nuc_cellprob)
    else:
        detector = args.detector
        if detector == "auto":
            detector = "blob_log" if args.construct == "JABr" else "cellpose"
        print(f"    condensate detector: {detector}")
        if detector == "blob_log":
            # No nucleus channel here, so the intra-nuclear gate is disabled:
            # an all-true mask keeps every detected blob.
            all_nuc = np.ones(stack.shape, dtype=bool)
            masks = detect_condensates_blob(stack, all_nuc, args.blob_threshold,
                                            args.min_sigma, args.max_sigma, args.num_sigma)
        else:
            if args.cond_model is not None:
                seg = models.CellposeModel(gpu=use_gpu, pretrained_model=str(args.cond_model))
            else:
                seg = models.CellposeModel(gpu=use_gpu, model_type="cyto3")
            masks = segment_condensates(restored, seg, args.diameter, cellprob_threshold=args.cond_cellprob)

    print("\n[3/3] Measuring...")
    meas = extract_slice_measurements(masks, stack)
    vol = compute_volumes(masks, args.voxel_xy, args.voxel_z)
    tiff.imwrite(output_dir / f"{args.channel}_masks.tif", masks)
    meas.to_csv(output_dir / f"{args.channel}_measurements.csv", index=False)
    vol.to_csv(output_dir / f"{args.channel}_volumes.csv", index=False)
    n = int(masks.max())
    pd.DataFrame([("mode", "single_channel"), ("channel", args.channel), ("n_objects", n)],
                 columns=["metric", "value"]).to_csv(output_dir / "summary.csv", index=False)
    print(f"\n  {args.channel}: {n} objects")
    print(f"All outputs saved to: {output_dir}")


def run_full_pipeline(args, use_gpu, output_dir, log=print):
    """Run the full two-channel PC pipeline and write all outputs.

    Shared by the CLI (`main`) and the web app. `args` is any object exposing the
    same attributes as the argparse namespace; `log` is a line-callback so the
    web app can stream progress (defaults to ``print`` for the CLI). Returns the
    ``pc_result`` dict.
    """
    # Load
    log("\n[1/6] Loading stacks...")
    cond_stack, nuc_stack = load_stacks(args.cond, args.nuc, args.roi)

    # Denoise
    log("\n[2/6] Denoising with Cellpose 3 DenoiseModel...")
    dn_model     = denoise.DenoiseModel(model_type="denoise_cyto3", gpu=use_gpu)
    cond_restored = denoise_stack(cond_stack, dn_model, "condensates")
    nuc_restored  = denoise_stack(nuc_stack,  dn_model, "nuclei")

    # Resolve detector: 'auto' routes JABr to blob_log, everything else to Cellpose
    detector = args.detector
    if detector == "auto":
        detector = "blob_log" if args.construct == "JABr" else "cellpose"
    log(f"    Condensate detector: {detector}")

    # Segment nuclei (always Cellpose cyto3)
    log("\n[3/6] Segmenting with Cellpose 3 (do_3D=True)...")
    nuc_seg_model = models.CellposeModel(gpu=use_gpu, model_type="cyto3")
    nuc_masks_3d  = segment_nuclei(nuc_restored,  nuc_seg_model,  args.nuc_diameter, args.nuc_cellprob, close_iters=args.nuc_close)

    # Segment condensates. For blob_log we split into nuclear vs cytoplasmic
    # by centroid; for cellpose the same mask is used and the split happens at
    # PC computation time via intersection with the nucleus mask.
    cyto_cond_masks_3d = None
    if detector == "blob_log":
        cond_masks_3d, cyto_cond_masks_3d = detect_condensates_blob_both(
            cond_stack, nuc_masks_3d > 0,
            threshold=args.blob_threshold, min_sigma=args.min_sigma,
            max_sigma=args.max_sigma, num_sigma=args.num_sigma,
        )
    else:
        if args.cond_model is not None:
            log(f"    Condensate model: {args.cond_model}")
            cond_seg_model = models.CellposeModel(gpu=use_gpu, pretrained_model=str(args.cond_model))
        else:
            cond_seg_model = nuc_seg_model
        all_cond_masks_3d = segment_condensates(cond_restored, cond_seg_model, args.diameter, cellprob_threshold=args.cond_cellprob)
        # Split: nuclear = AND nucleus, cytoplasmic = AND NOT nucleus
        nuc_3d = nuc_masks_3d > 0
        cond_masks_3d = (all_cond_masks_3d * nuc_3d).astype(np.int32)
        cyto_cond_masks_3d = (all_cond_masks_3d * (~nuc_3d)).astype(np.int32)

    # Per-slice measurements
    log("\n[4/6] Extracting per-slice measurements...")
    cond_df = extract_slice_measurements(cond_masks_3d, cond_stack)
    nuc_df  = extract_slice_measurements(nuc_masks_3d,  nuc_stack)
    log(f"    Condensate rows: {len(cond_df)}  |  Nuclei rows: {len(nuc_df)}")

    # 3D volumes
    log("\n[5/6] Computing 3D volumes...")
    cond_vol_df = compute_volumes(cond_masks_3d, args.voxel_xy, args.voxel_z)
    nuc_vol_df  = compute_volumes(nuc_masks_3d,  args.voxel_xy, args.voxel_z)
    log(f"    Condensate objects: {len(cond_vol_df)}  |  Nuclei objects: {len(nuc_vol_df)}")

    # Partition coefficient
    log("\n[6/6] Computing partition coefficient...")
    pc_result = compute_partition_coefficient(cond_stack, cond_masks_3d, nuc_masks_3d, cond_topx=args.cond_topx)
    log(f"    [nuclear]    PC               : {pc_result['pc']:.3f}")
    log(f"                 Background (B)   : {pc_result['background']:.2f}")
    log(f"                 Cond density     : {pc_result['cond_density']:.2f}")
    log(f"                 Dilute density   : {pc_result['dilute_density']:.2f}")

    if cyto_cond_masks_3d is not None:
        pc_cyto = compute_cytoplasmic_partition_coefficient(
            cond_stack, cyto_cond_masks_3d, nuc_masks_3d, cond_topx=args.cond_topx)
        if not np.isnan(pc_cyto["pc"]):
            log(f"    [cytoplasmic] PC              : {pc_cyto['pc']:.3f}")
            log(f"                 Cond density     : {pc_cyto['cond_density']:.2f}")
            log(f"                 Dilute density   : {pc_cyto['dilute_density']:.2f}")
        else:
            log(f"    [cytoplasmic] no cytoplasmic condensates detected")
        pc_result["pc_cytoplasmic"] = pc_cyto["pc"]
        pc_result["cond_density_cytoplasmic"] = pc_cyto["cond_density"]
        pc_result["dilute_density_cytoplasmic"] = pc_cyto["dilute_density"]

    pc_cal, was_cal = apply_calibration(pc_result["pc"], args.construct)
    if was_cal:
        entry = CALIBRATION[args.construct]
        if entry.get("kind") == "linear":
            m, b = entry["slope"], entry["intercept"]
            log(f"    PC calibrated ({args.construct}, {m:.3f}*pc+{b:.3f}): {pc_result['pc']:.3f} -> {pc_cal:.3f}")
        else:
            log(f"    PC calibrated ({args.construct}, {entry.get('kind', 'unknown')}): {pc_result['pc']:.3f} -> {pc_cal:.3f}")
        pc_result["pc_calibrated"] = pc_cal
        pc_result["construct"] = args.construct
    elif args.construct is not None:
        log(f"    [warn] No calibration for '{args.construct}'; returning raw pipeline_pc")

    # Save
    save_outputs(
        output_dir,
        cond_restored, nuc_restored,
        cond_masks_3d, nuc_masks_3d,
        cond_df, nuc_df,
        cond_vol_df, nuc_vol_df,
        pc_result,
        args.voxel_xy, args.voxel_z,
    )
    plot_summary(output_dir, cond_df, nuc_df, cond_vol_df, nuc_vol_df, pc_result)

    if args.view:
        log("\n[viewer] Building interactive Z-stack viewer...")
        write_zviewer(output_dir, cond_stack, nuc_stack, cond_masks_3d, nuc_masks_3d,
                      pc_result=pc_result, voxel_xy=args.voxel_xy, voxel_z=args.voxel_z)

    log(f"\nAll outputs saved to: {output_dir}")
    return pc_result


def max_overlap_nucleus(nuc_masks, cond_masks, log=print):
    """Keep only the nucleus with the most condensate-mask overlap.

    The "target cell" is the one Imaris's manual analysis was done on — almost
    always the cell containing the visible condensates. Picking the nucleus with
    the largest cond ∩ nuc voxel count selects that cell directly, regardless of
    where it sits in the cropped FOV.
    """
    if nuc_masks.max() == 0:
        return nuc_masks
    cond_3d = cond_masks > 0
    labels = np.unique(nuc_masks)
    labels = labels[labels > 0]
    best_label, best_overlap = None, -1
    for lbl in labels:
        overlap = int(((nuc_masks == lbl) & cond_3d).sum())
        if overlap > best_overlap:
            best_overlap, best_label = overlap, int(lbl)
    out = np.zeros_like(nuc_masks)
    if best_label is not None and best_overlap > 0:
        out[nuc_masks == best_label] = 1
        log(f"    target nucleus: label {best_label}  ({best_overlap} cond voxels overlap)")
    else:
        log("    no nucleus has condensate overlap — keeping all nuclei")
        out = (nuc_masks > 0).astype(nuc_masks.dtype)
    return out


def run_batch(folder, output_dir, args, use_gpu, ref_csv=None, log=print):
    """Process every TIF in a folder; write per-cell outputs + comparison.csv.

    Shared by the GUI batch tab and the web app. Each cell gets its own
    subfolder (and its own zviewer.html when ``args.view``). Returns a structured
    result dict (per-cell rows + correlation metrics) for the caller to display.
    """
    folder = Path(folder)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    tif_files = sorted(folder.glob("*.tif")) + sorted(folder.glob("*.tiff"))
    if not tif_files:
        log("No .tif files found in the selected folder.")
        return {"cells": [], "n": 0, "scatter": False}

    detector = args.detector
    if detector == "auto":
        detector = "blob_log" if args.construct == "JABr" else "cellpose"
    log(f"Found {len(tif_files)} TIF files in {folder}")
    log(f"Construct: {args.construct or '(none)'}    Detector: {detector}")
    log(f"GPU: {'enabled' if use_gpu else 'disabled'}\n")

    dn_model      = denoise.DenoiseModel(model_type="denoise_cyto3", gpu=use_gpu)
    nuc_seg_model = models.CellposeModel(gpu=use_gpu, model_type="cyto3")
    cond_seg_model = None
    if detector == "cellpose":
        cond_seg_model = (models.CellposeModel(gpu=use_gpu, pretrained_model=str(args.cond_model))
                          if args.cond_model else nuc_seg_model)

    ref_df = pc_col = None
    if ref_csv:
        ref_df = pd.read_csv(ref_csv)
        ref_df["_stem"] = ref_df[ref_df.columns[0]].apply(lambda x: Path(x).stem)
        pc_col = [c for c in ref_df.columns if "partition" in c.lower()][0]

    rows, cells = [], []
    for i, tif_path in enumerate(tif_files, 1):
        log(f"--- [{i}/{len(tif_files)}] {tif_path.name}")
        cell_dir = output_dir / tif_path.stem
        cell_dir.mkdir(exist_ok=True)
        cell = {"name": tif_path.name, "stem": tif_path.stem,
                "ok": False, "error": None, "has_viewer": False}
        try:
            cond_stack, nuc_stack = load_stacks(None, None, tif_path)
            cond_restored = denoise_stack(cond_stack, dn_model, "condensates")
            nuc_restored  = denoise_stack(nuc_stack,  dn_model, "nuclei")
            nuc_masks_3d  = segment_nuclei(nuc_restored, nuc_seg_model, args.nuc_diameter,
                                           args.nuc_cellprob, close_iters=args.nuc_close)

            if detector == "blob_log":
                cond_masks_3d, cyto_cond_masks_3d = detect_condensates_blob_both(
                    cond_stack, nuc_masks_3d > 0, threshold=args.blob_threshold,
                    min_sigma=args.min_sigma, max_sigma=args.max_sigma, num_sigma=args.num_sigma)
            else:
                all_cond = segment_condensates(cond_restored, cond_seg_model, args.diameter)
                nb = nuc_masks_3d > 0
                cond_masks_3d = (all_cond * nb).astype(np.int32)
                cyto_cond_masks_3d = (all_cond * (~nb)).astype(np.int32)

            nuc_single = max_overlap_nucleus(nuc_masks_3d, cond_masks_3d, log=log)
            pc = compute_partition_coefficient(cond_stack, cond_masks_3d, nuc_single, cond_topx=args.cond_topx)
            pc_cal, was_cal = apply_calibration(pc["pc"], args.construct)
            pc_cyto = compute_cytoplasmic_partition_coefficient(
                cond_stack, cyto_cond_masks_3d, nuc_masks_3d, cond_topx=args.cond_topx)

            row = {"file": tif_path.name, "pipeline_pc_nuclear": round(pc["pc"], 4),
                   "pipeline_pc_cytoplasmic": (round(pc_cyto["pc"], 4)
                                               if not np.isnan(pc_cyto["pc"]) else None)}
            if was_cal:
                row["pipeline_pc_nuclear_calibrated"] = round(pc_cal, 4)
            if ref_df is not None:
                match = ref_df[ref_df["_stem"] == tif_path.stem]
                if not match.empty:
                    ref_val = float(match.iloc[0][pc_col])
                    row["reference_pc"] = round(ref_val, 4)
                    row["error_pct_raw"] = round((pc["pc"] - ref_val) / ref_val * 100, 1)
                    if was_cal:
                        row["error_pct_calibrated"] = round((pc_cal - ref_val) / ref_val * 100, 1)
            rows.append(row)

            msg = f"    nuc PC raw = {pc['pc']:.3f}"
            if was_cal:               msg += f"  cal = {pc_cal:.3f}"
            if not np.isnan(pc_cyto["pc"]): msg += f"   cyto PC = {pc_cyto['pc']:.3f}"
            if ref_df is not None and "reference_pc" in row:
                msg += f"  ref = {row['reference_pc']}"
            log(msg)

            if args.view:
                if was_cal:
                    pc["pc_calibrated"] = pc_cal
                    pc["construct"] = args.construct
                pc["pc_cytoplasmic"] = pc_cyto["pc"] if not np.isnan(pc_cyto["pc"]) else None
                write_zviewer(cell_dir, cond_stack, nuc_stack, cond_masks_3d, nuc_masks_3d, pc_result=pc)
                cell["has_viewer"] = True

            cell.update(ok=True, pc=round(pc["pc"], 4),
                        pc_calibrated=(round(pc_cal, 4) if was_cal else None),
                        pc_cytoplasmic=(round(pc_cyto["pc"], 4) if not np.isnan(pc_cyto["pc"]) else None),
                        reference_pc=row.get("reference_pc"))
        except Exception as e:
            log(f"    ERROR: {e}")
            rows.append({"file": tif_path.name, "pipeline_pc_nuclear": float("nan"), "error": str(e)})
            cell["error"] = str(e)
        cells.append(cell)

    results_df = pd.DataFrame(rows)
    results_df.to_csv(output_dir / "comparison.csv", index=False)
    log(f"\nSaved comparison.csv -> {output_dir}")

    out = {"cells": cells, "n": len(tif_files), "scatter": False,
           "comparison_csv": "comparison.csv"}

    if ref_df is not None and "reference_pc" in results_df.columns:
        y_col = ("pipeline_pc_nuclear_calibrated"
                 if "pipeline_pc_nuclear_calibrated" in results_df.columns
                 else "pipeline_pc_nuclear")
        valid = results_df.dropna(subset=["reference_pc", y_col])
        if len(valid) > 1:
            r = float(np.corrcoef(valid["reference_pc"], valid[y_col])[0, 1])
            rmse = float(np.sqrt(((valid[y_col] - valid["reference_pc"]) ** 2).mean()))
            mae_pct = float((abs(valid[y_col] - valid["reference_pc"]) / valid["reference_pc"]).mean() * 100)
            fig, ax = plt.subplots(figsize=(5, 5))
            ax.scatter(valid["reference_pc"], valid[y_col], color="#1A73E8",
                       alpha=0.8, edgecolors="white", s=60)
            lim = max(valid["reference_pc"].max(), valid[y_col].max()) * 1.1
            ax.plot([0, lim], [0, lim], "k--", lw=0.8, alpha=0.5)
            ax.set_xlabel("Reference PC")
            ax.set_ylabel("Pipeline PC (calibrated)" if "calibrated" in y_col else "Pipeline PC (raw)")
            ax.set_title(f"r = {r:.3f}  |  RMSE = {rmse:.3f}  |  MAE = {mae_pct:.1f}%")
            ax.set_xlim(0, lim); ax.set_ylim(0, lim)
            plt.tight_layout(); plt.savefig(output_dir / "scatter.png", dpi=150); plt.close()
            log(f"Saved scatter.png  (r={r:.3f}, RMSE={rmse:.3f}, MAE={mae_pct:.1f}%)")
            out.update(scatter=True, r=round(r, 3), rmse=round(rmse, 3), mae=round(mae_pct, 1))

    log(f"\nAll outputs saved to: {output_dir}")
    return out


def main():
    args = parse_args()

    output_dir = args.output or (Path(__file__).parent / "outputs")
    output_dir.mkdir(parents=True, exist_ok=True)

    use_gpu = core.use_gpu() and not args.no_gpu
    print(f"GPU: {'enabled — ' + torch.cuda.get_device_name(0) if use_gpu else 'disabled'}")

    # Single-channel segmentation mode (no PC).
    if args.single is not None:
        run_single_channel(args, use_gpu, output_dir)
        return

    if args.roi is None and (args.cond is None or args.nuc is None):
        print("Error: provide either --roi, both --cond and --nuc, or --single")
        raise SystemExit(1)

    run_full_pipeline(args, use_gpu, output_dir)


if __name__ == "__main__":
    main()
