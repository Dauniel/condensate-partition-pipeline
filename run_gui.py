"""
run_gui.py — Simple GUI for the condensate pipeline.

Usage:
    python run_gui.py
"""

import sys
import threading
import tkinter as tk
from tkinter import ttk, filedialog, scrolledtext
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))


# ── Helpers ───────────────────────────────────────────────────────────────────

def browse_file(entry, title="Select file", parent=None):
    path = filedialog.askopenfilename(
        parent=parent,
        title=title,
        filetypes=[("TIF files", "*.tif *.tiff"), ("All files", "*.*")],
    )
    if path:
        entry.delete(0, tk.END)
        entry.insert(0, path)


def browse_folder(entry, title="Select folder", parent=None):
    path = filedialog.askdirectory(parent=parent, title=title)
    if path:
        entry.delete(0, tk.END)
        entry.insert(0, path)


# ── Main window ───────────────────────────────────────────────────────────────

class PipelineGUI:
    def __init__(self, root: tk.Tk):
        self.root = root
        root.title("Condensate Pipeline")
        root.resizable(True, True)
        root.minsize(780, 600)
        root.columnconfigure(0, weight=1)
        root.rowconfigure(0, weight=1)

        ACCENT  = "#1A73E8"
        BG      = "#F8F9FA"
        CARD_BG = "#FFFFFF"

        root.configure(bg=BG)

        style = ttk.Style()
        style.theme_use("clam")
        style.configure("TFrame",        background=BG)
        style.configure("TNotebook",     background=BG, tabposition="n")
        style.configure("TNotebook.Tab", font=("Helvetica", 11), padding=[14, 6])
        style.configure("TLabel",        background=CARD_BG, font=("Helvetica", 11))
        style.configure("BG.TLabel",     background=BG,      font=("Helvetica", 11))
        style.configure("Head.TLabel",   background=BG,      font=("Helvetica", 13, "bold"), foreground=ACCENT)
        style.configure("Sub.TLabel",    background=BG,      font=("Helvetica", 10), foreground="#555555")
        style.configure("TEntry",        font=("Helvetica", 11), fieldbackground=CARD_BG)
        style.configure("TCheckbutton",  background=BG,      font=("Helvetica", 11))
        style.configure("Run.TButton",   font=("Helvetica", 12, "bold"),
                        foreground="white", background=ACCENT, padding=8)
        style.map("Run.TButton",
                  background=[("active", "#1558B0"), ("disabled", "#AAAAAA")])

        outer = ttk.Frame(root, padding=20)
        outer.grid(row=0, column=0, sticky="nsew")
        outer.columnconfigure(0, weight=1)
        outer.rowconfigure(5, weight=1)  # log row expands

        # ── Title ─────────────────────────────────────────────────────────────
        ttk.Label(outer, text="Condensate Partition Coefficient Pipeline",
                  style="Head.TLabel").grid(row=0, column=0, sticky="w")
        ttk.Label(outer, text="Franco Lab  ·  Spring 2026",
                  style="Sub.TLabel").grid(row=1, column=0, sticky="w", pady=(0, 10))

        # ── Notebook tabs ─────────────────────────────────────────────────────
        nb = ttk.Notebook(outer)
        nb.grid(row=2, column=0, sticky="ew")

        self.single_tab = ttk.Frame(nb, padding=14)
        self.batch_tab  = ttk.Frame(nb, padding=14)
        nb.add(self.single_tab, text="  Single File  ")
        nb.add(self.batch_tab,  text="  Batch  ")

        self._build_single_tab()
        self._build_batch_tab()

        # ── Shared settings ───────────────────────────────────────────────────
        card_set = self._card(outer, row=3, title="Settings")

        ttk.Label(card_set, text="Top-X% brightest voxels for cond. density:").grid(
            row=0, column=0, sticky="w", pady=3)
        self.topx_var = tk.StringVar(value="75")
        ttk.Entry(card_set, textvariable=self.topx_var, width=6).grid(
            row=0, column=1, sticky="w", padx=(8, 0))
        ttk.Label(card_set, text="(default 75 — recommended)", foreground="#777").grid(
            row=0, column=2, sticky="w", padx=6)

        ttk.Label(card_set, text="Nuclei cell probability threshold:").grid(
            row=1, column=0, sticky="w", pady=3)
        self.cellprob_var = tk.StringVar(value="-2.0")
        ttk.Entry(card_set, textvariable=self.cellprob_var, width=6).grid(
            row=1, column=1, sticky="w", padx=(8, 0))
        ttk.Label(card_set, text="(default -2.0)", foreground="#777").grid(
            row=1, column=2, sticky="w", padx=6)

        ttk.Label(card_set, text="Construct:").grid(row=2, column=0, sticky="w", pady=3)
        self.construct_var = tk.StringVar(value="(none)")
        ttk.OptionMenu(card_set, self.construct_var, "(none)",
                       "(none)", "JABr", "GABr", "AABr", "JABr_4arm", "Tornado"
                       ).grid(row=2, column=1, sticky="w", padx=(8, 0))
        ttk.Label(card_set, text="(picks calibration + auto-routes detector)",
                  foreground="#777").grid(row=2, column=2, sticky="w", padx=6)

        ttk.Label(card_set, text="Condensate detector:").grid(row=3, column=0, sticky="w", pady=3)
        self.detector_var = tk.StringVar(value="auto")
        ttk.OptionMenu(card_set, self.detector_var, "auto",
                       "auto", "cellpose", "blob_log"
                       ).grid(row=3, column=1, sticky="w", padx=(8, 0))
        ttk.Label(card_set, text="(auto: JABr→blob_log, else Cellpose)",
                  foreground="#777").grid(row=3, column=2, sticky="w", padx=6)

        ttk.Label(card_set, text="Cellpose cond. model path:").grid(row=4, column=0, sticky="w", pady=3)
        self.cond_model_var = tk.StringVar(value="")
        ttk.Entry(card_set, textvariable=self.cond_model_var, width=42).grid(
            row=4, column=1, columnspan=2, sticky="w", padx=(8, 0))
        ttk.Label(card_set, text="(blank = built-in cyto3; recommended: V3 epoch 35)",
                  foreground="#777").grid(row=5, column=1, columnspan=2, sticky="w", padx=(8, 0))

        ttk.Label(card_set, text="blob_log threshold:").grid(row=6, column=0, sticky="w", pady=3)
        self.blob_thresh_var = tk.StringVar(value="0.03")
        ttk.Entry(card_set, textvariable=self.blob_thresh_var, width=6).grid(
            row=6, column=1, sticky="w", padx=(8, 0))
        ttk.Label(card_set, text="(default 0.03 — only used when detector=blob_log)",
                  foreground="#777").grid(row=6, column=2, sticky="w", padx=6)

        self.gpu_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(card_set, text="Disable GPU  (use on laptop / no CUDA)",
                        variable=self.gpu_var).grid(row=7, column=0, columnspan=3,
                                                    sticky="w", pady=(6, 0))

        self.view_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(card_set, text="Write interactive Z-viewer (zviewer.html)",
                        variable=self.view_var).grid(row=8, column=0, columnspan=3,
                                                     sticky="w", pady=(2, 0))

        # ── Run button + status ───────────────────────────────────────────────
        btn_frame = ttk.Frame(outer, style="TFrame")
        btn_frame.grid(row=4, column=0, sticky="ew", pady=(14, 4))
        self.run_btn = ttk.Button(btn_frame, text="▶  Run Pipeline",
                                  style="Run.TButton", command=self._run)
        self.run_btn.pack(side="left")
        self.status_lbl = ttk.Label(btn_frame, text="", style="BG.TLabel",
                                    font=("Helvetica", 11), foreground="#555")
        self.status_lbl.pack(side="left", padx=14)

        # ── Shared log ────────────────────────────────────────────────────────
        log_card = self._card(outer, row=5, title="Output Log")
        log_card.master.columnconfigure(0, weight=1)
        log_card.master.rowconfigure(0, weight=1)
        log_card.columnconfigure(0, weight=1)
        log_card.rowconfigure(1, weight=1)
        self.log = scrolledtext.ScrolledText(
            log_card, width=72, height=18,
            font=("Courier", 10), bg="#1A1A2E", fg="#E0E0E0",
            insertbackground="white", state="disabled",
        )
        self.log.grid(row=1, column=0, sticky="nsew")

        self._toggle_mode()

    # ── Single tab ────────────────────────────────────────────────────────────

    def _build_single_tab(self):
        p = self.single_tab
        BG = "#F8F9FA"

        mode_card = self._card(p, row=0, title="Input Mode")
        self.mode = tk.StringVar(value="roi")
        ttk.Radiobutton(mode_card, text="Single multi-channel TIF  (OME / Cut ROI)",
                        variable=self.mode, value="roi",
                        command=self._toggle_mode).grid(row=0, column=0, sticky="w", pady=2)
        ttk.Radiobutton(mode_card, text="Separate condensate + nuclei files",
                        variable=self.mode, value="split",
                        command=self._toggle_mode).grid(row=1, column=0, sticky="w", pady=2)
        ttk.Radiobutton(mode_card, text="Single channel — segment one image, nuclei or condensate (no PC)",
                        variable=self.mode, value="single",
                        command=self._toggle_mode).grid(row=2, column=0, sticky="w", pady=2)

        self.roi_frame   = self._card(p, row=1, title="Multi-Channel TIF")
        self.roi_entry   = self._file_row(self.roi_frame, 0, "TIF file:", "Select TIF")

        self.split_frame = self._card(p, row=2, title="Channel Files")
        self.cond_entry  = self._file_row(self.split_frame, 0, "Condensate (Ch2):", "Select condensate TIF")
        self.nuc_entry   = self._file_row(self.split_frame, 1, "Nuclei (Ch1):",     "Select nuclei TIF")

        self.single_frame  = self._card(p, row=3, title="Single Channel")
        self.single_ch_entry = self._file_row(self.single_frame, 0, "Image file:", "Select single-channel TIF")
        ttk.Label(self.single_frame, text="This image is:").grid(row=2, column=0, sticky="w", pady=3)
        self.single_ch_type = tk.StringVar(value="condensate")
        ttk.OptionMenu(self.single_frame, self.single_ch_type, "condensate",
                       "condensate", "nuclei").grid(row=2, column=1, sticky="w", padx=(8, 0))
        ttk.Label(self.single_frame, text="(segments this one channel; outputs masks + measurements, no PC)",
                  foreground="#777").grid(row=2, column=2, sticky="w", padx=6)

        out_card = self._card(p, row=4, title="Output")
        self.single_out_entry = self._folder_row(out_card, 0, "Output folder:")

    # ── Batch tab ─────────────────────────────────────────────────────────────

    def _build_batch_tab(self):
        p = self.batch_tab

        dir_card = self._card(p, row=0, title="Input")
        self.batch_dir_entry = self._folder_row(
            dir_card, 0, "Folder of TIF files:",
            title="Select folder containing TIF files")
        ttk.Label(dir_card, text="All .tif files in this folder will be processed.",
                  foreground="#777", font=("Helvetica", 10)).grid(
            row=2, column=0, columnspan=3, sticky="w", pady=(4, 0))

        ref_card = self._card(p, row=1, title="Reference CSV  (optional)")
        self.ref_csv_entry = self._file_row(
            ref_card, 0, "Nuclear PC CSV:", "Select reference CSV",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")])
        ttk.Label(ref_card,
                  text="If provided, outputs comparison.csv and scatter.png vs manual reference.",
                  foreground="#777", font=("Helvetica", 10)).grid(
            row=2, column=0, columnspan=3, sticky="w", pady=(4, 0))

        out_card = self._card(p, row=2, title="Output")
        self.batch_out_entry = self._folder_row(out_card, 0, "Output folder:")

    # ── Layout helpers ────────────────────────────────────────────────────────

    def _card(self, parent, row, title):
        CARD_BG = "#FFFFFF"
        wrapper = tk.Frame(parent, bg="#DADCE0", padx=1, pady=1)
        wrapper.grid(row=row, column=0, sticky="ew", pady=5)
        inner = tk.Frame(wrapper, bg=CARD_BG, padx=14, pady=10)
        inner.pack(fill="both", expand=True)
        tk.Label(inner, text=title, bg=CARD_BG,
                 font=("Helvetica", 10, "bold"), fg="#1A73E8").grid(
            row=0, column=0, columnspan=3, sticky="w", pady=(0, 6))
        return inner

    def _file_row(self, parent, row, label, dialog_title="Select file", filetypes=None):
        r = row + 1
        ttk.Label(parent, text=label).grid(row=r, column=0, sticky="w", pady=2)
        entry = ttk.Entry(parent, width=46)
        entry.grid(row=r, column=1, padx=8)
        ft = filetypes or [("TIF files", "*.tif *.tiff"), ("All files", "*.*")]
        ttk.Button(parent, text="Browse…",
                   command=lambda e=entry, t=dialog_title, f=ft: self._browse_file(e, t, f)
                   ).grid(row=r, column=2)
        return entry

    def _folder_row(self, parent, row, label, title="Select folder"):
        r = row + 1
        ttk.Label(parent, text=label).grid(row=r, column=0, sticky="w", pady=2)
        entry = ttk.Entry(parent, width=46)
        entry.grid(row=r, column=1, padx=8)
        ttk.Button(parent, text="Browse…",
                   command=lambda e=entry, t=title: browse_folder(e, t, self.root)
                   ).grid(row=r, column=2)
        return entry

    def _browse_file(self, entry, title, filetypes):
        path = filedialog.askopenfilename(
            parent=self.root, title=title, filetypes=filetypes)
        if path:
            entry.delete(0, tk.END)
            entry.insert(0, path)

    def _toggle_mode(self):
        m = self.mode.get()
        self.roi_frame.master.grid_remove()
        self.split_frame.master.grid_remove()
        self.single_frame.master.grid_remove()
        if m == "roi":
            self.roi_frame.master.grid()
        elif m == "split":
            self.split_frame.master.grid()
        else:
            self.single_frame.master.grid()

    # ── Logging ───────────────────────────────────────────────────────────────

    def _log(self, msg: str):
        self.log.configure(state="normal")
        self.log.insert(tk.END, msg + "\n")
        self.log.see(tk.END)
        self.log.configure(state="disabled")

    def _set_status(self, msg, color="#555555"):
        self.status_lbl.configure(text=msg, foreground=color)

    # ── Dispatch run to correct tab ───────────────────────────────────────────

    def _run(self):
        try:
            topx        = float(self.topx_var.get())
            cellprob    = float(self.cellprob_var.get())
            blob_thresh = float(self.blob_thresh_var.get())
        except ValueError:
            self._set_status("Invalid setting — check numeric fields.", "#C62828")
            return

        no_gpu     = self.gpu_var.get()
        construct  = self.construct_var.get()
        construct  = None if construct in ("(none)", "") else construct
        detector   = self.detector_var.get()
        if detector == "auto":
            detector = "blob_log" if construct == "JABr" else "cellpose"
        cond_model = self.cond_model_var.get().strip() or None

        cfg = dict(topx=topx, cellprob=cellprob, no_gpu=no_gpu,
                   construct=construct, detector=detector,
                   cond_model=cond_model, blob_thresh=blob_thresh,
                   view=self.view_var.get())

        # Detect active tab
        active = self.root.nametowidget(
            self.root.nametowidget(
                str(self.single_tab.winfo_parent())
            ).select()
        )
        if active is self.single_tab:
            self._run_single(cfg)
        else:
            self._run_batch(cfg)

    # ── Single-cell run ───────────────────────────────────────────────────────

    def _run_single(self, cfg):
        mode = self.mode.get()
        if mode == "single":
            path = self.single_ch_entry.get().strip()
            if not path:
                self._set_status("Please select an image file.", "#C62828"); return
            out = self.single_out_entry.get().strip() or None
            self._start_worker(
                lambda: self._worker_single_channel(path, self.single_ch_type.get(), out, cfg))
            return
        if mode == "roi":
            roi = self.roi_entry.get().strip()
            if not roi:
                self._set_status("Please select a TIF file.", "#C62828"); return
            cond = nuc = None
        else:
            cond = self.cond_entry.get().strip()
            nuc  = self.nuc_entry.get().strip()
            if not cond or not nuc:
                self._set_status("Please select both channel files.", "#C62828"); return
            roi = None

        out = self.single_out_entry.get().strip() or None
        self._start_worker(lambda: self._worker_single(roi, cond, nuc, out, cfg))

    # ── Single-channel run (one image, nuclei or condensate; no PC) ─────────────

    def _worker_single_channel(self, path, channel, out, cfg):
        import io, contextlib

        class LogWriter(io.TextIOBase):
            def __init__(self, gui): self.gui = gui
            def write(self, s):
                if s.strip(): self.gui._log(s.rstrip())
                return len(s)

        import matplotlib
        matplotlib.use("Agg")
        with contextlib.redirect_stdout(LogWriter(self)):
            from pipeline import (denoise_stack, segment_nuclei, segment_condensates,
                                  detect_condensates_blob, extract_slice_measurements,
                                  compute_volumes)
            import numpy as np
            import pandas as pd
            import tifffile as tiff
            from cellpose import models, core, denoise as cp_denoise

            output_dir = Path(out) if out else Path(__file__).parent / "outputs"
            output_dir.mkdir(parents=True, exist_ok=True)

            use_gpu = core.use_gpu() and not cfg["no_gpu"]
            self._log(f"GPU: {'enabled' if use_gpu else 'disabled'}")
            self._log(f"Single-channel mode: {channel}  (no partition coefficient)")

            self._log("\n[1/3] Loading + denoising...")
            stack = tiff.imread(path)
            self._log(f"  stack {stack.shape}  dtype={stack.dtype}")
            dn_model = cp_denoise.DenoiseModel(model_type="denoise_cyto3", gpu=use_gpu)
            restored = denoise_stack(stack, dn_model, channel)

            self._log("\n[2/3] Segmenting...")
            if channel == "nuclei":
                seg = models.CellposeModel(gpu=use_gpu, model_type="cyto3")
                masks = segment_nuclei(restored, seg, None, cfg["cellprob"])
            else:
                detector = cfg["detector"]   # already resolved (auto -> blob_log/cellpose)
                self._log(f"    condensate detector: {detector}")
                if detector == "blob_log":
                    # No nucleus channel -> intra-nuclear gate disabled (keep all blobs).
                    all_nuc = np.ones(stack.shape, dtype=bool)
                    masks = detect_condensates_blob(stack, all_nuc, cfg["blob_thresh"],
                                                    1.5, 6.0, 8)
                else:
                    if cfg["cond_model"]:
                        seg = models.CellposeModel(gpu=use_gpu, pretrained_model=cfg["cond_model"])
                    else:
                        seg = models.CellposeModel(gpu=use_gpu, model_type="cyto3")
                    masks = segment_condensates(restored, seg, None, cellprob_threshold=0.0)

            self._log("\n[3/3] Measuring...")
            meas = extract_slice_measurements(masks, stack)
            vol = compute_volumes(masks, None, None)
            tiff.imwrite(output_dir / f"{channel}_masks.tif", masks)
            meas.to_csv(output_dir / f"{channel}_measurements.csv", index=False)
            vol.to_csv(output_dir / f"{channel}_volumes.csv", index=False)
            n = int(masks.max())
            pd.DataFrame([("mode", "single_channel"), ("channel", channel), ("n_objects", n)],
                         columns=["metric", "value"]).to_csv(output_dir / "summary.csv", index=False)
            self._log(f"\n  {channel}: {n} objects detected")
            self._log(f"All outputs saved to: {output_dir}")

    def _worker_single(self, roi, cond, nuc, out, cfg):
        import io, contextlib

        class LogWriter(io.TextIOBase):
            def __init__(self, gui): self.gui = gui
            def write(self, s):
                if s.strip(): self.gui._log(s.rstrip())
                return len(s)

        import matplotlib
        matplotlib.use("Agg")
        with contextlib.redirect_stdout(LogWriter(self)):
            from pipeline import (
                load_stacks, denoise_stack, segment_condensates,
                segment_nuclei, detect_condensates_blob, detect_condensates_blob_both,
                extract_slice_measurements,
                compute_volumes, compute_partition_coefficient,
                compute_cytoplasmic_partition_coefficient,
                save_outputs, plot_summary, apply_calibration,
            )
            import numpy as np
            from cellpose import models, core, denoise as cp_denoise

            output_dir = Path(out) if out else Path(__file__).parent / "outputs"
            output_dir.mkdir(parents=True, exist_ok=True)

            use_gpu = core.use_gpu() and not cfg["no_gpu"]
            self._log(f"GPU: {'enabled' if use_gpu else 'disabled'}")
            self._log(f"Construct: {cfg['construct'] or '(none)'}    Detector: {cfg['detector']}")

            self._log("\n[1/6] Loading stacks...")
            cond_stack, nuc_stack = load_stacks(
                Path(cond) if cond else None,
                Path(nuc)  if nuc  else None,
                Path(roi)  if roi  else None,
            )

            self._log("\n[2/6] Denoising...")
            dn_model      = cp_denoise.DenoiseModel(model_type="denoise_cyto3", gpu=use_gpu)
            cond_restored = denoise_stack(cond_stack, dn_model, "condensates")
            nuc_restored  = denoise_stack(nuc_stack,  dn_model, "nuclei")

            self._log("\n[3/6] Segmenting nuclei (Cellpose cyto3)...")
            nuc_seg_model = models.CellposeModel(gpu=use_gpu, model_type="cyto3")
            nuc_masks_3d  = segment_nuclei(nuc_restored, nuc_seg_model, None, cfg["cellprob"])

            cyto_cond_masks_3d = None
            if cfg["detector"] == "blob_log":
                self._log(f"\n[3b/6] Detecting condensates (blob_log split, threshold={cfg['blob_thresh']})...")
                cond_masks_3d, cyto_cond_masks_3d = detect_condensates_blob_both(
                    cond_stack, nuc_masks_3d > 0,
                    threshold=cfg["blob_thresh"], min_sigma=1.5, max_sigma=6.0, num_sigma=8,
                )
            else:
                self._log(f"\n[3b/6] Segmenting condensates (Cellpose, model={cfg['cond_model'] or 'cyto3'})...")
                if cfg["cond_model"]:
                    cond_seg_model = models.CellposeModel(gpu=use_gpu, pretrained_model=cfg["cond_model"])
                else:
                    cond_seg_model = nuc_seg_model
                all_cond_masks_3d = segment_condensates(cond_restored, cond_seg_model, None)
                nuc_3d = nuc_masks_3d > 0
                cond_masks_3d = (all_cond_masks_3d * nuc_3d).astype(np.int32)
                cyto_cond_masks_3d = (all_cond_masks_3d * (~nuc_3d)).astype(np.int32)

            self._log("\n[4/6] Measuring...")
            cond_df = extract_slice_measurements(cond_masks_3d, cond_stack)
            nuc_df  = extract_slice_measurements(nuc_masks_3d,  nuc_stack)

            self._log("\n[5/6] Computing 3D volumes...")
            cond_vol_df = compute_volumes(cond_masks_3d, None, None)
            nuc_vol_df  = compute_volumes(nuc_masks_3d,  None, None)

            self._log("\n[6/6] Computing partition coefficient...")
            pc = compute_partition_coefficient(cond_stack, cond_masks_3d, nuc_masks_3d, cond_topx=cfg["topx"])
            self._log(f"\n  [nuclear]      PC raw         : {pc['pc']:.3f}")
            self._log(f"                 Background (B) : {pc['background']:.2f}")
            self._log(f"                 Cond density   : {pc['cond_density']:.2f}")
            self._log(f"                 Dilute density : {pc['dilute_density']:.2f}")
            pc_cal, was_cal = apply_calibration(pc["pc"], cfg["construct"])
            if was_cal:
                self._log(f"                 PC calibrated  : {pc_cal:.3f}   ({cfg['construct']})")
                pc["pc_calibrated"] = pc_cal
                pc["construct"] = cfg["construct"]
            elif cfg["construct"]:
                self._log(f"                 [warn] no calibration entry for '{cfg['construct']}'")

            if cyto_cond_masks_3d is not None:
                pc_cyto = compute_cytoplasmic_partition_coefficient(
                    cond_stack, cyto_cond_masks_3d, nuc_masks_3d, cond_topx=cfg["topx"])
                if not np.isnan(pc_cyto["pc"]):
                    self._log(f"\n  [cytoplasmic]  PC raw         : {pc_cyto['pc']:.3f}")
                    self._log(f"                 Cond density   : {pc_cyto['cond_density']:.2f}")
                    self._log(f"                 Dilute density : {pc_cyto['dilute_density']:.2f}")
                else:
                    self._log("\n  [cytoplasmic]  no cytoplasmic condensates detected")
                pc["pc_cytoplasmic"] = pc_cyto["pc"]
                pc["cond_density_cytoplasmic"] = pc_cyto["cond_density"]
                pc["dilute_density_cytoplasmic"] = pc_cyto["dilute_density"]

            save_outputs(output_dir, cond_restored, nuc_restored,
                         cond_masks_3d, nuc_masks_3d,
                         cond_df, nuc_df, cond_vol_df, nuc_vol_df,
                         pc, None, None)
            plot_summary(output_dir, cond_df, nuc_df, cond_vol_df, nuc_vol_df, pc)
            if cfg.get("view"):
                from pipeline import write_zviewer
                self._log("\nBuilding interactive Z-viewer (zviewer.html)...")
                write_zviewer(output_dir, cond_stack, nuc_stack, cond_masks_3d, nuc_masks_3d,
                              pc_result=pc)
                self._log(f"  Open {output_dir / 'zviewer.html'} in any browser.")
            self._log(f"\nAll outputs saved to: {output_dir}")

    # ── Batch run ─────────────────────────────────────────────────────────────

    def _run_batch(self, cfg):
        folder = self.batch_dir_entry.get().strip()
        if not folder:
            self._set_status("Please select a folder.", "#C62828"); return

        ref_csv = self.ref_csv_entry.get().strip() or None
        out     = self.batch_out_entry.get().strip() or None
        self._start_worker(lambda: self._worker_batch(folder, ref_csv, out, cfg))

    def _worker_batch(self, folder, ref_csv, out, cfg):
        import io, contextlib
        from types import SimpleNamespace

        class LogWriter(io.TextIOBase):
            def __init__(self, gui): self.gui = gui
            def write(self, s):
                if s.strip(): self.gui._log(s.rstrip())
                return len(s)

        import matplotlib
        matplotlib.use("Agg")
        with contextlib.redirect_stdout(LogWriter(self)):
            from pipeline import run_batch
            from cellpose import core

            output_dir = Path(out) if out else Path(__file__).parent / "outputs" / "batch"
            use_gpu = core.use_gpu() and not cfg["no_gpu"]

            args = SimpleNamespace(
                construct=cfg["construct"] or None, detector=cfg["detector"],
                cond_topx=cfg["topx"], nuc_cellprob=cfg["cellprob"],
                nuc_diameter=None, nuc_close=0,
                blob_threshold=cfg["blob_thresh"], min_sigma=1.5, max_sigma=6.0, num_sigma=8,
                diameter=None, cond_model=cfg["cond_model"], cond_cellprob=0.0,
                voxel_xy=None, voxel_z=None, view=cfg.get("view", True),
            )
            run_batch(folder, output_dir, args, use_gpu, ref_csv=ref_csv)

    # ── Worker harness ────────────────────────────────────────────────────────

    def _start_worker(self, fn):
        self.run_btn.configure(state="disabled")
        self._set_status("Running…", "#1A73E8")
        self.log.configure(state="normal"); self.log.delete("1.0", tk.END)
        self.log.configure(state="disabled")

        def worker():
            try:
                fn()
                self.root.after(0, lambda: self._set_status("Done ✓", "#18965C"))
            except Exception as e:
                import traceback
                msg = str(e)
                self._log(f"\nERROR: {msg}\n{traceback.format_exc()}")
                self.root.after(0, lambda m=msg: self._set_status(f"Error: {m}", "#C62828"))
            finally:
                self.root.after(0, lambda: self.run_btn.configure(state="normal"))

        threading.Thread(target=worker, daemon=True).start()


# ── Entry point ───────────────────────────────────────────────────────────────

if __name__ == "__main__":
    root = tk.Tk()
    PipelineGUI(root)
    root.mainloop()
