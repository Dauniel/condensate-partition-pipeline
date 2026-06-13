# -*- coding: utf-8 -*-
"""
build_pdf.py — Render a polished PDF instruction manual with reportlab.

Usage:  python docs/build_pdf.py
Output: docs/Manual.pdf

This is a convenience renderer so a PDF exists without a LaTeX install. The
canonical sources are docs/MANUAL.md (GitHub) and docs/manual_latex/manual.tex
(Overleaf). Keep this in sync with those if the content changes.
"""
from pathlib import Path
from reportlab.lib.pagesizes import letter
from reportlab.lib.units import inch
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_CENTER
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Image,
                                Table, TableStyle, PageBreak)

HERE = Path(__file__).parent
FIG = HERE / "figures"
ACCENT = colors.HexColor("#1A73E8")
INK = colors.HexColor("#222222")
MUTED = colors.HexColor("#555555")

ss = getSampleStyleSheet()
H1 = ParagraphStyle("H1", parent=ss["Heading1"], textColor=ACCENT, fontSize=16,
                    spaceBefore=14, spaceAfter=6)
H2 = ParagraphStyle("H2", parent=ss["Heading2"], textColor=INK, fontSize=12.5,
                    spaceBefore=8, spaceAfter=4)
BODY = ParagraphStyle("Body", parent=ss["BodyText"], fontSize=10.3, leading=14.5,
                      spaceAfter=6)
CODE = ParagraphStyle("Code", parent=ss["BodyText"], fontName="Courier", fontSize=9,
                      backColor=colors.HexColor("#F2F4F7"), borderPadding=5,
                      leading=12, spaceAfter=6)
CAP = ParagraphStyle("Cap", parent=BODY, fontSize=8.5, textColor=MUTED,
                     alignment=TA_CENTER)


def tbl(data, widths):
    t = Table(data, colWidths=widths, repeatRows=1)
    t.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), ACCENT),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTSIZE", (0, 0), (-1, -1), 8.8),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("GRID", (0, 0), (-1, -1), 0.5, colors.HexColor("#D0D5DD")),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#F8F9FB")]),
        ("TOPPADDING", (0, 0), (-1, -1), 5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
        ("LEFTPADDING", (0, 0), (-1, -1), 6),
        ("RIGHTPADDING", (0, 0), (-1, -1), 6),
    ]))
    return t


def P(s, style=BODY):
    return Paragraph(s, style)


def fig(name, width=5.6 * inch, caption=None):
    p = FIG / name
    out = []
    if p.exists():
        out.append(Image(str(p), width=width, height=width * _ratio(p)))
        if caption:
            out.append(Spacer(1, 3))
            out.append(P(caption, CAP))
        out.append(Spacer(1, 8))
    return out


def _ratio(path):
    from PIL import Image as PImage
    w, h = PImage.open(path).size
    return h / w


def build():
    story = []

    # ── Title page ────────────────────────────────────────────────────────────
    story += [Spacer(1, 1.6 * inch)]
    story.append(P("Condensate Partition Coefficient Pipeline",
                   ParagraphStyle("T", parent=H1, fontSize=26, alignment=TA_CENTER,
                                  textColor=ACCENT, leading=30)))
    story.append(Spacer(1, 0.2 * inch))
    story.append(P("Instruction Manual",
                   ParagraphStyle("ST", parent=BODY, fontSize=15, alignment=TA_CENTER)))
    story.append(Spacer(1, 0.5 * inch))
    story.append(P("Automated 3D segmentation of nuclear biomolecular condensates "
                   "and measurement of the partition coefficient from confocal Z-stacks",
                   ParagraphStyle("sub", parent=BODY, alignment=TA_CENTER, fontSize=11,
                                  textColor=MUTED)))
    story.append(Spacer(1, 1.4 * inch))
    story.append(P("<b>Franco Lab &nbsp;&bull;&nbsp; University of California, Los Angeles</b>",
                   ParagraphStyle("lab", parent=BODY, alignment=TA_CENTER, fontSize=12)))
    story.append(P("Author: Daniel Chang", ParagraphStyle("a", parent=BODY,
                   alignment=TA_CENTER, fontSize=11)))
    story.append(P("Principal Investigator: Elisa Franco",
                   ParagraphStyle("pi", parent=BODY, alignment=TA_CENTER, fontSize=11,
                                  textColor=MUTED)))
    story.append(Spacer(1, 0.3 * inch))
    story.append(P("Spring 2026", ParagraphStyle("d", parent=BODY,
                   alignment=TA_CENTER, fontSize=11, textColor=MUTED)))
    story.append(PageBreak())

    # ── Overview ────────────────────────────────────────────────────────────────
    story.append(P("Overview", H1))
    story.append(P("This pipeline replaces the manual Imaris workflow for measuring the "
                   "<b>partition coefficient (PC)</b> of nuclear biomolecular condensates. "
                   "Given a two-channel confocal Z-stack, it segments nuclei and condensates "
                   "in 3D, computes the PC, and reports a calibrated value on the manual "
                   "reference scale."))
    story.append(tbl([
        ["Item", "Description"],
        ["Inputs", "Two-channel Z-stack: channel 0 = nuclei, channel 1 = condensate."],
        ["Outputs", "Nuclear & cytoplasmic PC, 3D masks, per-object volumes, summary figure."],
        ["Validated on", "JABr: Pearson r = 0.942, MAE 12.9%, 79% of cells within ±20% (n = 28)."],
        ["Interfaces", "Graphical (run_gui.py) or command line (pipeline.py)."],
    ], [1.4 * inch, 4.7 * inch]))
    story.append(Spacer(1, 10))
    story += fig("jabr_validation.png", width=3.9 * inch,
                 caption="JABr validation: calibrated pipeline PC vs manual Imaris reference.")

    # ── Installation ──────────────────────────────────────────────────────────
    story.append(P("Installation (one time)", H1))
    story.append(P("1. Install Python 3.10 or 3.11.<br/>"
                   "2. Install PyTorch for your machine from pytorch.org (CPU or your CUDA "
                   "version), e.g. CUDA 12.1:"))
    story.append(P("pip install torch --index-url https://download.pytorch.org/whl/cu121", CODE))
    story.append(P("3. Install the remaining dependencies:"))
    story.append(P("pip install -r requirements.txt", CODE))
    story.append(P("On first run, Cellpose downloads its model weights automatically "
                   "(one-time ~100 MB). A GPU is recommended (~10–15 s per cell); CPU works "
                   "but is slower."))

    # ── Input prep ──────────────────────────────────────────────────────────────
    story.append(P("Preparing your image files", H1))
    story.append(P("The pipeline accepts either one multi-channel TIF (both channels in one "
                   "file) or two separate single-channel TIFs. <b>Channel order matters:</b> "
                   "channel 0 must be nuclei, channel 1 the condensate. Most axis orderings "
                   "(ZCYX, CZYX, ...) are auto-detected. If results look wrong, swapped "
                   "channels are the most common cause."))

    # ── GUI ─────────────────────────────────────────────────────────────────────
    story.append(P("Running with the graphical interface", H1))
    story.append(P("Launch with <font face='Courier'>python run_gui.py</font>. A window with "
                   "<b>Single File</b> and <b>Batch</b> tabs opens."))
    story.append(P("Step by step (single cell)", H2))
    story.append(P("1. On the <b>Single File</b> tab, choose the input mode and Browse to your "
                   "data.<br/>"
                   "2. Choose an Output folder (optional).<br/>"
                   "3. In <b>Settings</b>, set <b>Construct = JABr</b>; leave all other settings "
                   "at their defaults.<br/>"
                   "4. Click <b>▶ Run Pipeline</b>. Progress streams in the Output Log; the "
                   "status reads <b>Done</b> when finished.<br/>"
                   "5. Open your output folder for the masks, tables, and results.png."))
    story.append(P("Live GUI screenshots should be captured on the machine where it runs "
                   "(Windows: Win+Shift+S) and saved into docs/figures/. The actual numerical "
                   "output of a real run is shown under “Reading the outputs.”",
                   ParagraphStyle("note", parent=BODY, textColor=MUTED, fontSize=9.3)))

    # ── Settings ──────────────────────────────────────────────────────────────
    story.append(P("Understanding the settings", H1))
    story.append(P("For routine JABr analysis you only need to set <b>Construct</b>; the rest "
                   "have sensible defaults."))
    story.append(tbl([
        ["Setting", "Default", "What it does"],
        ["Construct", "(none)", "Pick JABr for the validated workflow. Auto-selects the blob_log "
         "detector and the JABr calibration so the reported PC is on the manual scale."],
        ["Condensate detector", "auto", "auto routes JABr→blob_log, else Cellpose. Normally left on auto."],
        ["Top-X% brightest", "75", "Defines condensed-phase density; trims the dim mask halo to "
         "match the manual reference. Leave at 75."],
        ["Nuclei cell-prob.", "-2.0", "How aggressively nuclei are detected. -2.0 suits the lab's stain."],
        ["blob_log threshold", "0.03", "Condensate spot sensitivity (used when detector = blob_log). "
         "Leave at 0.03 for JABr."],
        ["Cellpose cond. model", "(blank)", "Advanced; a fine-tuned condensate model. Leave blank for JABr."],
        ["Disable GPU", "off", "Tick on a laptop with no NVIDIA GPU."],
    ], [1.5 * inch, 0.75 * inch, 3.85 * inch]))
    story.append(Spacer(1, 6))
    story.append(P("<b>The one rule for routine use:</b> set Construct = JABr and run. Everything "
                   "else is for experimentation or other constructs."))

    # ── Outputs ───────────────────────────────────────────────────────────────
    story.append(P("Reading the outputs", H1))
    story.append(P("Every run writes summary.csv (headline numbers incl. the calibrated PC), "
                   "results.png (summary figure), condensate_masks.tif and nuclei_masks.tif "
                   "(3D labeled masks), and per-object/per-slice CSVs."))
    story.append(P("For the bundled sample JABr_Sample2_5_3.tif, a real run produces "
                   "<b>nuclear PC raw 4.867 → calibrated 4.813</b> (manual reference 4.558, +5.6%)."))
    story.append(P("<b>Report the calibrated nuclear PC.</b> The raw value is systematically "
                   "~3× the manual scale by construction; calibration standardizes it. "
                   "<b>Always glance at condensate_masks.tif</b> over the raw channel in "
                   "Fiji/ImageJ — a 10-second check catches the rare misfire."))
    story += fig("example_results.png", width=6.0 * inch,
                 caption="results.png from the bundled-sample run.")

    # ── Batch + CLI ─────────────────────────────────────────────────────────────
    story.append(P("Batch mode and command line", H1))
    story.append(P("On the <b>Batch</b> tab, select a folder of .tif files; every file is "
                   "processed. Optionally supply the manual Imaris nuclear-PC CSV as a Reference "
                   "CSV — the pipeline then writes comparison.csv and scatter.png with the "
                   "correlation, RMSE, and MAE."))
    story.append(P("Command line:"))
    story.append(P("python pipeline.py --roi sample_data/JABr_Sample2_5_3.tif --construct JABr --output my_output", CODE))
    story.append(P("Separate channels: --nuc nuclei.tif --cond condensate.tif. Physical volumes "
                   "in µm³: add --voxel-xy 0.065 --voxel-z 0.3. Force CPU: --no-gpu."))

    # ── Troubleshooting ───────────────────────────────────────────────────────
    story.append(P("Troubleshooting", H1))
    story.append(tbl([
        ["Symptom", "Likely cause / fix"],
        ["PC wildly off / nuclei empty", "Channels swapped: confirm ch0 = nuclei, ch1 = condensate."],
        ["“0 condensates detected”", "Not JABr, or very dim signal. Lower blob_log threshold "
         "(e.g. 0.02); confirm the raw channel has spots."],
        ["Very slow (minutes/cell)", "Running on CPU. Install a CUDA PyTorch build; ensure Disable GPU is unchecked."],
        ["CUDA out of memory", "Tick Disable GPU, or close other GPU programs."],
        ["Cellpose download fails", "First run needs internet for model weights."],
        ["One cell's calibrated PC off", "Calibration is population-level; cells vary (±20% for "
         "~80% of cells). Verify against the mask."],
    ], [2.0 * inch, 4.1 * inch]))

    # ── Method summary ──────────────────────────────────────────────────────────
    story.append(P("Method summary and scope", H1))
    story.append(P("The PC is the ratio of condensed-phase to dilute-phase density (Fabrini "
                   "et al.), both camera-background-subtracted. Condensed density uses the "
                   "brightest 75% of voxels inside (condensate ∩ nucleus); dilute density is "
                   "the mean of the 50 lowest-intensity 10×10×10 patches in the nuclear "
                   "dilute region. Nuclei are segmented with Cellpose 3 (cyto3, 3D), cleaned, "
                   "and void-filled; condensates are detected with Laplacian-of-Gaussian blob_log "
                   "and kept when ≥50% of the rendered sphere lies inside a nucleus. A "
                   "per-construct linear (JABr) or isotonic calibration standardizes the "
                   "automated PC onto the manual scale."))
    story.append(P("<b>Scope.</b> Validated and production-ready for JABr. The detector is "
                   "construct-specific; a leave-one-construct-out test confirmed a single "
                   "detector does not generalize zero-shot (it under-detects on unseen "
                   "constructs). Extending to other constructs needs per-construct calibration "
                   "or few-shot retraining. Full rationale and the cross-construct evidence are "
                   "in docs/METHODS.md."))

    doc = SimpleDocTemplate(str(HERE / "Manual.pdf"), pagesize=letter,
                            topMargin=0.9 * inch, bottomMargin=0.9 * inch,
                            leftMargin=0.95 * inch, rightMargin=0.95 * inch,
                            title="Condensate Pipeline — Instruction Manual",
                            author="Daniel Chang")
    doc.build(story)
    print("wrote", HERE / "Manual.pdf")


if __name__ == "__main__":
    build()
