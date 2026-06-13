# -*- coding: utf-8 -*-
"""
capture_gui.py — Build the pipeline GUI, populate it realistically, and grab
screenshots to docs/figures/ for the manual.

Produces:
  docs/figures/gui_main.png   configured Single-File run (sample loaded, Construct=JABr)
  docs/figures/gui_done.png   the same after a run, with the output log + Done status

Run on a machine with a desktop session:  python docs/capture_gui.py
"""
import sys, time
from pathlib import Path
import tkinter as tk
from PIL import ImageGrab

ROOT = Path(__file__).resolve().parent.parent          # repo root
FIG = ROOT / "docs" / "figures"
sys.path.insert(0, str(ROOT))
from run_gui import PipelineGUI

SAMPLE = ROOT / "sample_data" / "JABr_Sample2_5_3.tif"

DONE_LOG = """GPU: enabled
Construct: JABr    Detector: blob_log

[1/6] Loading stacks...
Condensate stack : (47, 161, 176)  dtype=uint16
Nuclei stack     : (47, 161, 176)  dtype=uint16

[2/6] Denoising...
[3/6] Segmenting nuclei (Cellpose cyto3)...
    nuclei after connected-component relabeling: 6
    nuclei after void-filling: +4120 voxels
[3b/6] Detecting condensates (blob_log split, threshold=0.03)...
    blob_log: 10 nuclear, 7 cytoplasmic (of 17 total)
[4/6] Measuring...
[5/6] Computing 3D volumes...
[6/6] Computing partition coefficient...

  [nuclear]  PC raw         : 4.867
             Background (B) : 79.00
             Cond density   : 319.46
             Dilute density : 65.64
             PC calibrated  : 4.813   (JABr)
All outputs saved to: my_output
"""


def grab(root, path):
    root.update_idletasks()
    root.lift()
    root.attributes("-topmost", True)
    for _ in range(6):
        root.update()
        time.sleep(0.15)
    x, y = root.winfo_rootx(), root.winfo_rooty()
    w, h = root.winfo_width(), root.winfo_height()
    img = ImageGrab.grab(bbox=(x, y, x + w, y + h))
    img.save(path)
    root.attributes("-topmost", False)
    print(f"saved {path}  ({img.size[0]}x{img.size[1]})")


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    root = tk.Tk()
    root.geometry("840x928+60+6")
    gui = PipelineGUI(root)
    gui.log.configure(height=6)   # shrink log so the whole window fits a 949px screen

    # Configure a realistic single-file run
    gui.mode.set("roi")
    gui._toggle_mode()
    gui.roi_entry.delete(0, tk.END)
    gui.roi_entry.insert(0, str(SAMPLE))
    gui.single_out_entry.delete(0, tk.END)
    gui.single_out_entry.insert(0, "my_output")
    gui.construct_var.set("JABr")
    gui.detector_var.set("auto")
    grab(root, FIG / "gui_main.png")

    # Simulate a completed run: fill the log + Done status
    gui.log.configure(state="normal")
    gui.log.delete("1.0", tk.END)
    gui.log.insert(tk.END, DONE_LOG)
    gui.log.see(tk.END)
    gui.log.configure(state="disabled")
    gui._set_status("Done ✓", "#18965C")
    grab(root, FIG / "gui_done.png")

    # Single-channel mode: segment one image (here a condensate-only file)
    gui._set_status("", "#555555")
    gui.log.configure(state="normal"); gui.log.delete("1.0", tk.END); gui.log.configure(state="disabled")
    gui.mode.set("single")
    gui._toggle_mode()
    gui.single_ch_entry.delete(0, tk.END)
    gui.single_ch_entry.insert(0, str(SAMPLE.parent / "condensate_only.tif"))
    gui.single_ch_type.set("condensate")
    grab(root, FIG / "gui_single.png")

    root.destroy()


if __name__ == "__main__":
    main()
