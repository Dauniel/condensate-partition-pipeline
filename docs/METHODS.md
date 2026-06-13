# Methods — How the pipeline measures the partition coefficient

**Franco Lab · UCLA · Spring 2026**

This document describes what the pipeline computes, why each step is designed the
way it is, and — importantly — the boundaries of where it is valid.

---

## 1. The quantity: nuclear partition coefficient

The **partition coefficient (PC)** quantifies how strongly a biomolecule
concentrates into condensates relative to the surrounding dilute phase:

```
PC = condensed-phase density / dilute-phase density
```

We follow the **Fabrini et al.** intensity-based definition. All intensities are
**background-subtracted**, where the background `B` is the minimum voxel intensity
across the field of view (the camera offset).

- **Condensed-phase density** — the mean of `clip(intensity − B, 0)` over the
  **brightest 75 %** of voxels that lie inside *both* a condensate mask *and* a
  nucleus. Automated masks (Cellpose `do_3D`, or a rendered LoG sphere) include a
  dim halo and dark interior pixels that a manual Imaris tracing would exclude;
  trimming the bottom 25 % of mask voxels removes that fluff and recovers
  densities consistent with the manual reference (validated by a top-X% sweep:
  −16 % bias → +3 % bias across 29 JABr cells).
- **Dilute-phase density** — the mean of the **50 lowest-intensity** valid
  10×10×10-voxel patches that fit entirely within the nuclear dilute region
  (nucleus minus condensates). Sorting by intensity and averaging many patches
  approximates how a researcher manually picks a quiet representative region, and
  is far more stable than a single random patch.

Both densities are background-subtracted; the PC is their ratio.

---

## 2. The pipeline, step by step

| Step | Method | Why |
|---|---|---|
| **1. Denoise** | Cellpose 3 `denoise_cyto3` on each Z-slice | Removes shot noise, which tightens segmentation boundaries. Intensities for the PC are read from the **raw** stack, not the denoised one. |
| **2. Segment nuclei** | Cellpose 3 `cyto3`, `do_3D=True` | Native 3D nuclear segmentation. |
| **2b. Clean nuclei** | Connected-component relabel; drop fragments < 1000 voxels | Cellpose over-splits large nuclei around internal condensate texture; relabeling restores whole nuclei without changing pixel coverage. |
| **2c. Fill voids** | Per-nucleus `binary_fill_holes` (3D + per-slice 2D) | Condensates exclude the nuclear stain, so Cellpose carves "donut holes" exactly where nuclear condensates sit. Filling them makes those condensates count as intra-nuclear. |
| **3. Detect condensates** | `blob_log` (Laplacian-of-Gaussian), σ = 1.5–6.0, threshold 0.03 | LoG is a classical, parameter-light blob detector well-matched to roughly-spherical condensates. Each detection is rendered as a sphere of radius σ·√3. |
| **3b. Intra-nuclear gate** | Keep a blob if **≥ 50 %** of its sphere lies inside a nucleus | Replaces an older single-centroid test that rejected boundary condensates and misassigned cells in crowded fields. Validated on JABr: r 0.926 → **0.942**, MAE 14.6 % → **13.7 %**. |
| **4–5. Measure** | `regionprops`, voxel counting | Per-object area, intensity, and 3D volume. |
| **6. Partition coefficient** | Section 1 formula | Nuclear PC (and cytoplasmic PC for blobs outside the nucleus). |
| **7. Calibrate** | Per-construct map (section 3) | Standardizes the automated PC onto the manual Imaris scale. |

---

## 3. Why there is a calibration

The automated pipeline and the manual Imaris workflow **define the condensed and
dilute regions slightly differently**:

- The automated condensed region is defined by the brightest mask voxels, which
  runs **brighter** than a hand-drawn Imaris condensate.
- The automated dilute region is the *quietest* sampled patches, which runs
  **dimmer** than a manually chosen region.

Both effects push the ratio the same direction, so the **raw automated PC is
systematically ~3× the manual reference** — but *consistently* so. Because the
bias is consistent (high correlation, r = 0.94), a single linear map converts
the automated value back onto the manual scale:

```
calibrated_PC = 0.3642 · raw_PC + 3.0405      (JABr)
```

This is **cross-method standardization**, not a fudge factor: it aligns two valid
but differently-defined measurements. Calibration **cannot improve correlation**
(a linear map is rank-preserving) — it only removes the scale/offset bias. The
correlation `r` is therefore the true accuracy ceiling, and calibration is only
meaningful where `r` is already high.

Calibrations for each construct live in
[`outputs/calibration_table.json`](../outputs/calibration_table.json) (linear for
JABr / JABr_4arm / Tornado; isotonic/piecewise for GABr / AABr).

---

## 4. Validation (JABr)

On the 28 JABr cells with a manual Imaris reference:

| Metric | Value |
|---|---|
| Pearson r (raw or calibrated) | **0.942** |
| Mean absolute error (calibrated) | **12.9 %** |
| Cells within ±20 % | **79 %** |

![JABr validation](figures/jabr_validation.png)

---

## 5. Scope and generalization

**The pipeline is validated and production-ready for JABr.** It is
**construct-specific**: the `blob_log` parameters (σ range, threshold) encode
JABr's condensate size and brightness, and do not transfer unchanged to other
constructs.

### Cross-construct performance (raw correlation vs Imaris)

| Construct | Pearson r | Status |
|---|---|---|
| **JABr** | **0.942** | Production |
| GABr | 0.678 | Mediocre; isotonic calibration available, use with caution |
| AABr | 0.122 | Not usable as-is |
| AwtBr | −0.116 | Not usable as-is |
| GwtBr | −0.341 | Not usable as-is |

The root cause is **scale**: the wild-type (wt) constructs have condensates up to
~2× larger than `blob_log`'s σ ceiling can represent, so a fixed detector either
misses them or fragments them.

### Can one model avoid per-construct tuning? (leave-one-construct-out test)

We tested whether a single learned detector (Cellpose fine-tuned on the Imaris
labels) could generalize **zero-shot** to a construct it never saw, via
leave-one-construct-out cross-validation:

| Held-out construct | learned-model r | cells with 0 detections |
|---|---|---|
| JABr | 0.875 | 19 / 30 |
| GABr | −0.207 | 25 / 30 |
| AABr | 0.084 | 17 / 25 |
| AwtBr | 0.369 | 23 / 33 |
| GwtBr | −0.221 | 0 / 31 |

**Conclusion:** a single model trained without examples of the target construct
**fails by under-detection** — it stays silent on most cells of a construct it
has never seen. The one fold with zero empty detections (GwtBr) is the one whose
morphological sibling (AwtBr) was in training, which points to the fix:

- **Per-construct calibration** where correlation is already acceptable (GABr).
- **Few-shot retraining** — adding even a handful of labeled ROIs of a new
  construct should restore detection — far cheaper than a full parameter sweep,
  and the recommended path to extend the tool beyond JABr.

This is a deliberate, evidence-based scope statement: the tool does one thing
(JABr) well, and the route to broadening it is mapped out rather than assumed.

---

## 6. Reference

- Fabrini et al. — intensity-based partition-coefficient methodology
  (10×10×10-voxel patch dilute-density estimate).
- Stringer & Pachitariu — *Cellpose 3* (segmentation + denoising).
- Lindeberg — Laplacian-of-Gaussian scale-space blob detection (`skimage.feature.blob_log`).
