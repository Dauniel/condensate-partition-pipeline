# -*- coding: utf-8 -*-
"""
webapp.py — Single-page web app that unifies the pipeline GUI and the Z-viewer.

One local browser app: configure a run, watch live progress, and scroll the
resulting Z-stack viewer in the same page (the viewer is the existing
`zviewer.html`, embedded in an iframe — no duplicated viewer logic).

Run:    python webapp.py
Then open http://127.0.0.1:5000 (it auto-opens a browser tab).

The heavy pipeline import (torch/cellpose) is deferred until the first run so the
server starts instantly.
"""
import sys, subprocess, threading, traceback, webbrowser
from pathlib import Path
from types import SimpleNamespace

from flask import (Flask, render_template, request, jsonify,
                   send_from_directory, abort)

ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT / "outputs" / "web_run"

app = Flask(__name__, template_folder=str(ROOT / "templates"))

# Single-user local tool → one shared run-state object guarded by a lock.
RUN = {
    "running": False, "done": False, "error": None,
    "log": [], "output_dir": None, "viewer": False, "summary": None,
    "batch": None,
}
LOCK = threading.Lock()


def _log(line):
    # run_full_pipeline emits one call per logical line (may contain \n).
    for part in str(line).split("\n"):
        RUN["log"].append(part)


def _gpu_available(disable):
    if disable:
        return False
    try:
        from cellpose import core
        return bool(core.use_gpu())
    except Exception:
        return False


def _worker(cfg):
    """Background thread: run the pipeline, stream progress into RUN['log']."""
    try:
        import pipeline as P  # deferred heavy import (torch/cellpose)

        out = Path(cfg["output"]).expanduser()
        if not out.is_absolute():
            out = ROOT / out
        out.mkdir(parents=True, exist_ok=True)
        RUN["output_dir"] = str(out)

        use_gpu = _gpu_available(cfg["no_gpu"])
        _log(f"GPU: {'enabled' if use_gpu else 'disabled'}")

        mode = cfg["mode"]
        if mode == "batch":
            args = _make_args(cfg)
            res = P.run_batch(cfg["folder"], out, args, use_gpu,
                              ref_csv=(cfg.get("ref_csv") or None), log=_log)
            RUN["batch"] = res
            RUN["viewer"] = False
        elif mode == "single":
            args = _make_args(cfg, single=cfg["path"])
            P.run_single_channel(args, use_gpu, out)
            RUN["viewer"] = False
        else:
            if mode == "roi":
                args = _make_args(cfg, roi=cfg["path"])
            else:  # separate channels
                args = _make_args(cfg, cond=cfg["cond"], nuc=cfg["nuc"])
            pc = P.run_full_pipeline(args, use_gpu, out, log=_log)
            RUN["viewer"] = bool(cfg["view"])
            RUN["summary"] = {
                "pc": pc.get("pc"),
                "pc_calibrated": pc.get("pc_calibrated"),
                "construct": pc.get("construct"),
            }
        RUN["done"] = True
    except Exception as e:
        RUN["error"] = f"{type(e).__name__}: {e}"
        _log("\nERROR: " + RUN["error"])
        _log(traceback.format_exc())
    finally:
        RUN["running"] = False


def _make_args(cfg, roi=None, cond=None, nuc=None, single=None):
    """Build an argparse-namespace-shaped object for run_full_pipeline."""
    return SimpleNamespace(
        roi=Path(roi) if roi else None,
        cond=Path(cond) if cond else None,
        nuc=Path(nuc) if nuc else None,
        single=Path(single) if single else None,
        channel=cfg.get("channel", "condensate"),
        construct=cfg["construct"] or None,
        detector=cfg["detector"],
        cond_topx=float(cfg["topx"]),
        nuc_cellprob=float(cfg["cellprob"]),
        nuc_diameter=None, nuc_close=0,
        blob_threshold=float(cfg["blob_thresh"]),
        min_sigma=1.5, max_sigma=6.0, num_sigma=8,
        diameter=None, cond_model=None, cond_cellprob=0.0,
        voxel_xy=None, voxel_z=None,
        view=bool(cfg["view"]),
    )


# ── Routes ────────────────────────────────────────────────────────────────────

@app.route("/")
def index():
    return render_template("index.html")


@app.route("/api/run", methods=["POST"])
def api_run():
    if RUN["running"]:
        return jsonify(error="A run is already in progress."), 409
    cfg = request.get_json(force=True)

    # Validate paths up front for a friendly error.
    mode = cfg.get("mode")
    if mode == "roi" and not cfg.get("path"):
        return jsonify(error="Select a multi-channel TIF."), 400
    if mode == "single" and not cfg.get("path"):
        return jsonify(error="Select an image file."), 400
    if mode == "separate" and not (cfg.get("cond") and cfg.get("nuc")):
        return jsonify(error="Select both condensate and nuclei files."), 400
    if mode == "batch" and not cfg.get("folder"):
        return jsonify(error="Select a folder of TIF files."), 400

    # Reset state.
    RUN.update(running=True, done=False, error=None, log=[],
               output_dir=None, viewer=False, summary=None, batch=None)
    threading.Thread(target=_worker, args=(cfg,), daemon=True).start()
    return jsonify(ok=True)


@app.route("/api/status")
def api_status():
    viewer_url = None
    if RUN["done"] and RUN["viewer"]:
        viewer_url = "/files/zviewer.html"

    batch = None
    if RUN["done"] and RUN["batch"]:
        b = RUN["batch"]
        cells = [{
            **c,
            "viewer_url": (f"/files/{c['stem']}/zviewer.html" if c.get("has_viewer") else None),
        } for c in b["cells"]]
        batch = {**b, "cells": cells,
                 "scatter_url": ("/files/scatter.png" if b.get("scatter") else None)}

    return jsonify(
        running=RUN["running"], done=RUN["done"], error=RUN["error"],
        log=RUN["log"], viewer_url=viewer_url, summary=RUN["summary"], batch=batch,
    )


@app.route("/files/<path:relpath>")
def files(relpath):
    """Serve files from the current run's output directory (no traversal)."""
    if not RUN["output_dir"]:
        abort(404)
    base = Path(RUN["output_dir"]).resolve()
    target = (base / relpath).resolve()
    if base not in target.parents and target != base:
        abort(403)
    if not target.exists():
        abort(404)
    return send_from_directory(str(base), relpath)


# A tiny script run as a subprocess so the OS file dialog lives in its own
# process — avoids tkinter's "must run on the main thread" issues under Flask.
_DIALOG_SRC = r"""
import sys, tkinter as tk
from tkinter import filedialog
kind = sys.argv[1]
root = tk.Tk(); root.withdraw()
root.attributes("-topmost", True); root.update()
if kind == "folder":
    p = filedialog.askdirectory(title="Select a folder")
elif kind == "csv":
    p = filedialog.askopenfilename(title="Select a CSV file",
        filetypes=[("CSV files", "*.csv"), ("All files", "*.*")])
else:
    p = filedialog.askopenfilename(title="Select a TIFF image",
        filetypes=[("TIFF images", "*.tif *.tiff"), ("All files", "*.*")])
root.destroy()
sys.stdout.write(p or "")
"""


@app.route("/api/native-pick")
def api_native_pick():
    """Open the operating system's native file/folder dialog and return the path."""
    kind = request.args.get("kind", "file")
    try:
        out = subprocess.run([sys.executable, "-c", _DIALOG_SRC, kind],
                             capture_output=True, text=True, timeout=600)
        return jsonify(path=out.stdout.strip())
    except Exception as e:
        return jsonify(path="", error=str(e))


def main():
    url = "http://127.0.0.1:5000"
    print(f"Condensate Pipeline web app -> {url}")
    if "--no-browser" not in sys.argv:
        threading.Timer(1.0, lambda: webbrowser.open(url)).start()
    # threaded=True so the pipeline run thread doesn't block status polling.
    app.run(host="127.0.0.1", port=5000, debug=False, threaded=True)


if __name__ == "__main__":
    main()
