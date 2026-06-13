# -*- coding: utf-8 -*-
"""
build_combined_pdf.py — Render the combined report (User Manual + Methods +
Research Timeline + JABr Reference Data) as a single black-and-white PDF.

Usage:  python docs/combined/build_combined_pdf.py
Output: docs/combined/Combined_Report.pdf

Black ink only (no color). Reference-data tables are parsed directly from the
Obsidian note so the numbers stay exact; a curated panel gallery is embedded
from docs/combined/figures/.
"""
from pathlib import Path
from reportlab.lib.pagesizes import letter
from reportlab.lib.units import inch
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_CENTER, TA_JUSTIFY
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Image,
                                Table, TableStyle, PageBreak)

HERE = Path(__file__).parent
FIG = HERE / "figures"
NOTE = Path(r"C:\storage\code\research\Research\reference_data\JABr Reference Data.md")

BLACK = colors.black
GREY_HDR = colors.HexColor("#E6E6E6")   # grayscale only — no color
GREY_ALT = colors.HexColor("#F4F4F4")
GREY_RULE = colors.HexColor("#999999")

ss = getSampleStyleSheet()
TITLE = ParagraphStyle("Title", parent=ss["Title"], textColor=BLACK,
                       fontName="Times-Bold", fontSize=26, leading=30, alignment=TA_CENTER)
SUBTITLE = ParagraphStyle("Sub", parent=ss["Normal"], textColor=BLACK,
                          fontName="Times-Roman", fontSize=14, alignment=TA_CENTER, leading=18)
PART = ParagraphStyle("Part", parent=ss["Heading1"], textColor=BLACK,
                      fontName="Times-Bold", fontSize=18, spaceBefore=18, spaceAfter=10)
H1 = ParagraphStyle("H1", parent=ss["Heading2"], textColor=BLACK,
                    fontName="Times-Bold", fontSize=13.5, spaceBefore=12, spaceAfter=5)
H2 = ParagraphStyle("H2", parent=ss["Heading3"], textColor=BLACK,
                    fontName="Times-Bold", fontSize=11.5, spaceBefore=8, spaceAfter=3)
BODY = ParagraphStyle("Body", parent=ss["BodyText"], textColor=BLACK,
                      fontName="Times-Roman", fontSize=10.5, leading=14.5,
                      spaceAfter=6, alignment=TA_JUSTIFY)
BULLET = ParagraphStyle("Bullet", parent=BODY, leftIndent=14, spaceAfter=3)
CODE = ParagraphStyle("Code", parent=ss["BodyText"], fontName="Courier", fontSize=8.8,
                      textColor=BLACK, backColor=colors.HexColor("#F0F0F0"),
                      borderPadding=5, leading=11.5, spaceAfter=6)
CAP = ParagraphStyle("Cap", parent=BODY, fontSize=8.6, alignment=TA_CENTER, spaceBefore=2)


def P(s, style=BODY):
    return Paragraph(s, style)


def _ratio(path):
    from PIL import Image as PImage
    w, h = PImage.open(path).size
    return h / w


def fig(name, width, caption=None):
    p = FIG / name
    out = []
    if p.exists():
        out += [Spacer(1, 4), Image(str(p), width=width, height=width * _ratio(p))]
        if caption:
            out.append(P(caption, CAP))
        out.append(Spacer(1, 8))
    return out


def table(rows, widths, font=8.2, header=True):
    t = Table(rows, colWidths=widths, repeatRows=1 if header else 0)
    style = [
        ("FONTNAME", (0, 0), (-1, -1), "Times-Roman"),
        ("FONTSIZE", (0, 0), (-1, -1), font),
        ("TEXTCOLOR", (0, 0), (-1, -1), BLACK),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("GRID", (0, 0), (-1, -1), 0.4, GREY_RULE),
        ("TOPPADDING", (0, 0), (-1, -1), 2.5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 2.5),
        ("LEFTPADDING", (0, 0), (-1, -1), 4),
        ("RIGHTPADDING", (0, 0), (-1, -1), 4),
    ]
    if header:
        style += [
            ("FONTNAME", (0, 0), (-1, 0), "Times-Bold"),
            ("BACKGROUND", (0, 0), (-1, 0), GREY_HDR),
            ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, GREY_ALT]),
        ]
    t.setStyle(TableStyle(style))
    return t


# ── Parse markdown tables from the note ────────────────────────────────────────

def parse_note_tables():
    """Return {section_title: [ [cells...], ... ]} for every markdown table."""
    text = NOTE.read_text(encoding="utf-8")
    sections, cur, rows = {}, None, []
    for line in text.splitlines():
        if line.startswith("## "):
            if cur and rows:
                sections.setdefault(cur, []).append(rows)
            cur, rows = line[3:].strip(), []
            continue
        s = line.strip()
        if s.startswith("|"):
            # protect escaped pipes
            s2 = s.replace(r"\|", "\x00")
            cells = [c.strip().replace("\x00", "|") for c in s2.strip("|").split("|")]
            if set("".join(cells)) <= set("-: "):   # separator row
                continue
            rows.append(cells)
        else:
            if cur and rows:
                sections.setdefault(cur, []).append(rows)
                rows = []
    if cur and rows:
        sections.setdefault(cur, []).append(rows)
    return sections


def build():
    sec = parse_note_tables()
    story = []

    # ── Title page ──────────────────────────────────────────────────────────────
    story += [Spacer(1, 1.7 * inch)]
    story.append(P("Condensate Partition Coefficient Pipeline", TITLE))
    story.append(Spacer(1, 0.25 * inch))
    story.append(P("Combined Report", SUBTITLE))
    story.append(P("User Manual &mdash; Methods &mdash; Research Timeline &mdash; JABr Reference Data", SUBTITLE))
    story.append(Spacer(1, 1.6 * inch))
    story.append(P("<b>Franco Lab &mdash; University of California, Los Angeles</b>",
                   ParagraphStyle("lab", parent=SUBTITLE, fontSize=13)))
    story.append(Spacer(1, 0.15 * inch))
    story.append(P("Author: Daniel Chang", ParagraphStyle("a", parent=SUBTITLE, fontSize=12)))
    story.append(P("Principal Investigator: ______________________",
                   ParagraphStyle("pi", parent=SUBTITLE, fontSize=12)))
    story.append(Spacer(1, 0.3 * inch))
    story.append(P("Spring 2026", ParagraphStyle("d", parent=SUBTITLE, fontSize=11)))
    story.append(PageBreak())

    # ── Contents ──────────────────────────────────────────────────────────────
    story.append(P("Contents", PART))
    for i, t in enumerate(["Part I &mdash; User Manual",
                           "Part II &mdash; Methods",
                           "Part III &mdash; Research Timeline",
                           "Part IV &mdash; JABr Reference Data"], 1):
        story.append(P(t, ParagraphStyle("toc", parent=BODY, fontSize=12, spaceAfter=8)))
    story.append(PageBreak())

    # ════════════════════════════ PART I — MANUAL ════════════════════════════
    story.append(P("Part I &mdash; User Manual", PART))
    story.append(P("This pipeline replaces the manual Imaris workflow for measuring the "
                   "partition coefficient (PC) of nuclear biomolecular condensates. Given a "
                   "two-channel confocal Z-stack (channel 0 = nuclei, channel 1 = condensate), "
                   "it segments nuclei and condensates in 3D, computes the PC, and reports a "
                   "calibrated value on the manual reference scale. It runs from a point-and-click "
                   "graphical interface or the command line."))
    story += fig("jabr_validation.png", 3.7 * inch,
                 "Validation on JABr: calibrated pipeline PC vs the manual Imaris reference.")

    story.append(P("1. Installation (one time)", H1))
    story.append(P("Install Python 3.10 or 3.11, then install PyTorch for your machine "
                   "(CPU or your CUDA version) from pytorch.org. Example for CUDA 12.1:"))
    story.append(P("pip install torch --index-url https://download.pytorch.org/whl/cu121", CODE))
    story.append(P("Then install the remaining dependencies:"))
    story.append(P("pip install -r requirements.txt", CODE))
    story.append(P("On first run, Cellpose downloads its model weights automatically (one-time "
                   "~100 MB). A GPU is recommended (~10&ndash;15 s per cell); CPU works but is slower."))

    story.append(P("2. Preparing image files", H1))
    story.append(P("The pipeline accepts either one multi-channel TIF (both channels in one file) "
                   "or two separate single-channel TIFs. Channel order matters: channel 0 must be "
                   "nuclei, channel 1 the condensate. Most axis orderings (ZCYX, CZYX, &hellip;) are "
                   "auto-detected; if results look wrong, swapped channels are the most common cause."))

    story.append(P("3. Running with the graphical interface", H1))
    story.append(P("Launch with <font face='Courier'>python run_gui.py</font>. A window with "
                   "<b>Single File</b> and <b>Batch</b> tabs opens.", BODY))
    for s in ["1. On the <b>Single File</b> tab, choose the input mode and browse to your data.",
              "2. Choose an output folder (optional).",
              "3. In <b>Settings</b>, set <b>Construct = JABr</b>; leave all other settings at their defaults.",
              "4. Click <b>Run Pipeline</b>. Progress streams in the Output Log; the status reads <b>Done</b> when finished.",
              "5. Open your output folder for the masks, tables, and results.png."]:
        story.append(P(s, BULLET))

    story.append(P("4. Settings", H1))
    story.append(P("For routine JABr analysis you only need to set <b>Construct</b>; the rest have "
                   "sensible defaults."))
    story.append(table([
        ["Setting", "Default", "What it does"],
        ["Construct", "(none)", "Pick JABr for the validated workflow. Auto-selects the blob_log detector and the JABr calibration so the reported PC is on the manual scale."],
        ["Condensate detector", "auto", "auto routes JABr to blob_log, else Cellpose. Normally left on auto."],
        ["Top-X% brightest", "75", "Defines condensed-phase density; trims the dim mask halo to match the manual reference. Leave at 75."],
        ["Nuclei cell-prob.", "-2.0", "How aggressively nuclei are detected. -2.0 suits the lab's nuclear stain."],
        ["blob_log threshold", "0.03", "Condensate spot sensitivity (used when detector = blob_log). Leave at 0.03 for JABr."],
        ["Disable GPU", "off", "Tick on a laptop with no NVIDIA GPU."],
    ], [1.4 * inch, 0.7 * inch, 4.0 * inch], font=9))
    story.append(Spacer(1, 5))
    story.append(P("<b>The one rule for routine use:</b> set Construct = JABr and run.", BODY))

    story.append(P("5. Reading the outputs", H1))
    story.append(P("Every run writes summary.csv (headline numbers incl. the calibrated PC), "
                   "results.png, condensate_masks.tif and nuclei_masks.tif (3D labeled masks), and "
                   "per-object/per-slice CSVs. For the bundled sample JABr_Sample2_5_3.tif, a real "
                   "run produces nuclear PC raw 4.867 calibrated to 4.813 (manual reference 4.558, "
                   "+5.6%). Report the calibrated nuclear PC, and always glance at "
                   "condensate_masks.tif over the raw channel in Fiji/ImageJ to catch the rare misfire."))
    story += fig("example_results.png", 6.0 * inch, "results.png from the bundled-sample run.")

    story.append(P("6. Command line", H1))
    story.append(P("python pipeline.py --roi sample_data/JABr_Sample2_5_3.tif --construct JABr --output my_output", CODE))
    story.append(P("Separate channels: --nuc nuclei.tif --cond condensate.tif. Physical volumes in "
                   "cubic microns: add --voxel-xy 0.065 --voxel-z 0.3. Force CPU: --no-gpu. For many "
                   "files (and optional comparison vs a manual reference CSV), use the GUI Batch tab."))

    story.append(P("7. Troubleshooting", H1))
    story.append(table([
        ["Symptom", "Likely cause / fix"],
        ["PC wildly off / nuclei empty", "Channels swapped: confirm ch0 = nuclei, ch1 = condensate."],
        ["“0 condensates detected”", "Not JABr, or very dim signal. Lower blob_log threshold (e.g. 0.02)."],
        ["Very slow (minutes/cell)", "Running on CPU. Install a CUDA PyTorch build; ensure Disable GPU is unchecked."],
        ["CUDA out of memory", "Tick Disable GPU, or close other GPU programs."],
        ["One cell's calibrated PC off", "Calibration is population-level; cells vary (+/-20% for ~80% of cells)."],
    ], [2.0 * inch, 4.1 * inch], font=9))
    story.append(PageBreak())

    # ════════════════════════════ PART II — METHODS ════════════════════════════
    story.append(P("Part II &mdash; Methods", PART))
    story.append(P("1. The quantity", H1))
    story.append(P("The partition coefficient is PC = condensed-phase density / dilute-phase "
                   "density (Fabrini et al.). All intensities are background-subtracted, where the "
                   "background B is the minimum voxel intensity across the field of view (camera offset)."))
    for s in ["<b>Condensed-phase density</b> &mdash; mean of clip(intensity &minus; B, 0) over the brightest 75% of voxels inside both a condensate mask and a nucleus. Automated masks include a dim halo a manual tracing would exclude; trimming the bottom 25% recovers densities consistent with the manual reference.",
              "<b>Dilute-phase density</b> &mdash; mean of the 50 lowest-intensity valid 10&times;10&times;10-voxel patches inside the nuclear dilute region. Averaging many low patches approximates a manual quiet-region selection and is far more stable than a single random patch."]:
        story.append(P(s, BULLET))

    story.append(P("2. Pipeline steps", H1))
    story.append(table([
        ["Step", "Method", "Why"],
        ["Denoise", "Cellpose 3 denoise_cyto3 per slice", "Removes shot noise; PC intensities are read from the raw stack."],
        ["Segment nuclei", "Cellpose 3 cyto3, do_3D", "Native 3D nuclear segmentation."],
        ["Clean nuclei", "Connected-component relabel; drop < 1000-voxel fragments", "Cellpose over-splits large nuclei; relabeling restores whole nuclei without changing pixel coverage."],
        ["Fill voids", "Per-nucleus binary_fill_holes (3D + 2D)", "Condensates carve donut holes in the nuclear stain exactly where nuclear condensates sit; filling makes them count as intra-nuclear."],
        ["Detect condensates", "blob_log (LoG), sigma 1.5-6.0, threshold 0.03", "Classical, parameter-light blob detector matched to roughly-spherical condensates; each detection rendered as a sphere of radius sigma*sqrt(3)."],
        ["Intra-nuclear gate", "Keep a blob if >= 50% of its sphere lies inside a nucleus", "Replaces an older single-centroid test that mis-handled boundary blobs and crowded fields."],
        ["Calibrate", "Per-construct linear (JABr) / isotonic map", "Standardizes the automated PC onto the manual Imaris scale."],
    ], [1.1 * inch, 2.1 * inch, 2.9 * inch], font=8.4))

    story.append(P("3. Why there is a calibration", H1))
    story.append(P("The automated and manual workflows define the condensed and dilute regions "
                   "slightly differently: the automated condensed region runs brighter (brightest "
                   "mask voxels) and the automated dilute region runs dimmer (quietest sampled "
                   "patches) than a hand-drawn Imaris region. Both push the ratio the same way, so "
                   "the raw automated PC is systematically about 3x the manual reference &mdash; but "
                   "consistently so (r = 0.94). A single linear map converts it back:"))
    story.append(P("calibrated_PC = 0.3642 * raw_PC + 3.0405      (JABr)", CODE))
    story.append(P("This is cross-method standardization, not a fudge factor. A linear map is "
                   "rank-preserving, so calibration cannot improve correlation &mdash; it only removes "
                   "the scale/offset bias. The correlation r is therefore the true accuracy ceiling, "
                   "and calibration is meaningful only where r is already high."))

    story.append(P("4. Scope and generalization", H1))
    story.append(P("The pipeline is validated and production-ready for JABr. It is "
                   "construct-specific: the blob_log parameters encode JABr's condensate size and "
                   "brightness and do not transfer unchanged. Cross-construct correlation vs Imaris:"))
    story.append(table([
        ["Construct", "Pearson r", "Status"],
        ["JABr", "0.942", "Production"],
        ["GABr", "0.678", "Mediocre; isotonic calibration, use with caution"],
        ["AABr", "0.122", "Not usable as-is"],
        ["AwtBr", "-0.116", "Not usable as-is"],
        ["GwtBr", "-0.341", "Not usable as-is"],
    ], [1.3 * inch, 1.2 * inch, 3.6 * inch], font=9))
    story.append(Spacer(1, 5))
    story.append(P("A leave-one-construct-out test (a single learned Cellpose detector, evaluated "
                   "on a construct it never trained on) confirmed that one model does not generalize "
                   "zero-shot: it fails by under-detection, staying silent on most cells of an unseen "
                   "construct. The route to extend the tool is per-construct calibration where "
                   "correlation is already acceptable, or few-shot retraining with a handful of "
                   "labeled ROIs of the new construct &mdash; far cheaper than a manual parameter sweep."))
    story.append(PageBreak())

    # ════════════════════════════ PART III — TIMELINE ════════════════════════════
    story.append(P("Part III &mdash; Research Timeline", PART))
    story.append(P("A chronological record of how the pipeline was built. Dates are from the "
                   "project's version-control history."))
    timeline = [
        ("Winter 2026 — Foundation", [
            ("2026-03-17", "Initial clean pipeline established."),
            ("2026-03-24", "Winter 2026 report finalized; baseline segmentation + PC concept."),
        ]),
        ("April 2026 — Segmentation survey & the PC formula", [
            ("2026-04-22", "Spring segmentation-model survey; 4 candidate models compared."),
            ("2026-04-25", "Background-subtracted PC formula adopted across all scripts."),
            ("2026-04-30", "PC now matches reference: connected-component nuclei fix + lowest-patch dilute density; 50-patch dilute stability fix."),
        ]),
        ("May 2026 — Calibration, GUI, generalization", [
            ("2026-05-06", "Batch cross-reference pipeline; top-75% condensate-density calibration."),
            ("2026-05-07", "GUI built; laptop (CPU) validation; PI meeting."),
            ("2026-05-14", "Crowded-field debugging; max-overlap-nucleus target-cell heuristic."),
            ("2026-05-16", "First fine-tuned Cellpose model; 5-construct accuracy sweep; per-construct calibration."),
            ("2026-05-18", "M187 research poster (PC = 6.297 vs 6.32 reference)."),
            ("2026-05-28", "blob_log detector wired into the pipeline; JABr calibration refit; GUI gains construct/detector controls; cytoplasmic PC added."),
        ]),
        ("June 2026 — Production hardening & honest scope", [
            ("2026-06-10", "Discussion section; nuclei void-filling so condensates in donut holes count toward the PC."),
            ("2026-06-11", "Sphere-overlap intra-nuclear gate promoted to production: JABr r 0.926 to 0.942, MAE 14.6% to 13.7%."),
            ("2026-06-12", "Leave-one-construct-out experiment: a single learned detector does not generalize zero-shot. Few-shot path mapped."),
            ("2026-06-13", "Clean public release: repository, instruction manual, methods writeup, and this combined report."),
        ]),
    ]
    for title, items in timeline:
        story.append(P(title, H1))
        rows = [["Date", "Milestone"]] + [[d, m] for d, m in items]
        story.append(table(rows, [1.1 * inch, 5.0 * inch], font=9.2))
        story.append(Spacer(1, 4))
    story.append(P("Where it stands", H1))
    story.append(P("Validated, production-ready for JABr (r = 0.942, MAE 12.9%, 79% of cells within "
                   "+/-20%, n = 28). Delivered: a CLI, a GUI for non-coding lab use, per-construct "
                   "calibration, and full documentation. Known boundary, with a plan: the detector "
                   "is JABr-specific; extending it needs per-construct calibration or few-shot retraining."))
    story.append(PageBreak())

    # ════════════════════════════ PART IV — JABr REFERENCE DATA ════════════════
    story.append(P("Part IV &mdash; JABr Reference Data", PART))
    story.append(P("Source: Box &rsaquo; Condensate Volume Quantification &rsaquo; JABr. The manual "
                   "Imaris measurements (cytoplasmic and nuclear) are the ground truth against which "
                   "the pipeline is validated; the pipeline results table that follows is the "
                   "production blob_log V4 output (sphere-overlap gate, linear calibration)."))

    # Reference: Nuclear (the validated region) + Cytoplasmic
    def ref_table(title, key):
        rows_groups = sec.get(key, [])
        if not rows_groups:
            return
        rows = rows_groups[0]
        story.append(P(title, H1))
        story.append(table(rows, [1.5 * inch, 1.7 * inch, 1.5 * inch, 1.4 * inch], font=8.2))
        story.append(Spacer(1, 6))

    ref_table("Manual reference &mdash; Nuclear (n = 30)", "Nuclear")
    story.append(PageBreak())
    ref_table("Manual reference &mdash; Cytoplasmic (n = 30)", "Cytoplasmic")

    # Pipeline production results table
    pipe_key = "Pipeline — Nuclear (blob_log, current production)"
    groups = sec.get(pipe_key, [])
    if groups:
        story.append(P("Pipeline results &mdash; Nuclear (blob_log V4, sphere-overlap gate)", H1))
        story.append(P("Detector blob_log (threshold 0.03, sigma 1.5-6.0); nuclei via Cellpose cyto3 "
                       "with void-filling; intra-nuclear gate >= 50% sphere overlap; calibration "
                       "0.3642 * raw + 3.0405. n = 28 (Sample3_3_10 and Sample3_3_15 find 0 "
                       "intra-nuclear condensates and are excluded)."))
        rows = groups[0]
        widths = [1.0 * inch, 0.85 * inch, 0.85 * inch, 0.6 * inch, 0.6 * inch, 0.7 * inch, 0.6 * inch, 0.55 * inch]
        story.append(table(rows, widths, font=7.4))
        story.append(Spacer(1, 5))
        story.append(P("<b>Summary:</b> n = 28, mean |err| = 13.7%, 22/28 within +/-20%, Pearson r = 0.942.", BODY))
        story.append(P("Notable failures (> 25% error)", H2))
        for s in ["<b>Sample1_2_3 (+63%)</b>: smallest reference PC (2.55); the linear calibration "
                  "floor (intercept ~3.0) dominates, so the percentage error is large though the "
                  "absolute error is only ~1.6 PC units.",
                  "<b>Sample1_1_1 (+43%)</b>: a poor nucleus mask drifting into a dim region inflates "
                  "the cond/dilute ratio; both blob_log and Cellpose fail here.",
                  "<b>Sample3_3_3 (+29%)</b>: an over-bright blob set keeps cond density high relative "
                  "to the reference PC; not fixed by the gate."]:
            story.append(P(s, BULLET))

    # Curated panel gallery
    story.append(PageBreak())
    story.append(P("Visual panels (representative)", H1))
    story.append(P("Per-ROI max-intensity Z-projection panels: merged reference (nuclei/condensate) "
                   "| nuclei channel | condensate channel | pipeline nuclei mask | pipeline "
                   "condensate mask | reference (Imaris) condensate mask | classification overlay. "
                   "A representative subset is shown; the full set of 28 panels is in the project "
                   "repository (docs / reference_data)."))
    gallery = [
        ("Sample1_4_2", "best case: ref PC 5.05, calibrated 5.06, err 0.3%"),
        ("Sample2_5_3", "bundled sample: ref PC 4.56, calibrated 4.67, err 2.4%"),
        ("Sample1_1_3", "high PC: ref 19.77, calibrated 19.53, err 1.2%"),
        ("Sample2_5_5", "rescued by the sphere-overlap gate: was +58%, now ref 6.73 / calibrated 5.72, err 15%"),
        ("Sample1_2_3", "failure: ref 2.55, calibrated 4.15, err 63% (calibration-floor effect)"),
        ("Sample3_3_3", "failure: ref 8.82, calibrated 11.42, err 29% (over-bright blobs)"),
    ]
    for name, cap in gallery:
        story += fig(f"{name}.png", 6.4 * inch, f"{name} &mdash; {cap}")

    # ── Page numbers ──────────────────────────────────────────────────────────
    def footer(canvas, doc):
        canvas.saveState()
        canvas.setFont("Times-Roman", 9)
        canvas.setFillColor(BLACK)
        canvas.drawCentredString(letter[0] / 2, 0.55 * inch, str(doc.page))
        canvas.restoreState()

    doc = SimpleDocTemplate(str(HERE / "Combined_Report.pdf"), pagesize=letter,
                            topMargin=0.85 * inch, bottomMargin=0.85 * inch,
                            leftMargin=0.95 * inch, rightMargin=0.95 * inch,
                            title="Condensate Pipeline — Combined Report", author="Daniel Chang")
    doc.build(story, onFirstPage=footer, onLaterPages=footer)
    print("wrote", HERE / "Combined_Report.pdf")


if __name__ == "__main__":
    build()
