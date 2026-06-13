# -*- coding: utf-8 -*-
"""
gen_combined_tex.py — Generate combined_report.tex (Overleaf-ready, black ink
only) for the combined report: Manual + Methods + Timeline + JABr Reference Data.

Table rows are parsed from the Obsidian note so the numbers stay exact.
Run:  python docs/combined/gen_combined_tex.py
Out:  docs/combined/combined_report.tex   (compile with pdfLaTeX)
"""
from pathlib import Path

HERE = Path(__file__).parent
OUT = HERE / "combined_report.tex"
NOTE = Path(r"C:\storage\code\research\Research\reference_data\JABr Reference Data.md")


def parse_note_tables():
    """Return {section_title: [ [cells...], ... ]} for every markdown table in the note."""
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
            s2 = s.replace(r"\|", "\x00")                      # protect escaped pipes
            cells = [c.strip().replace("\x00", "|") for c in s2.strip("|").split("|")]
            if set("".join(cells)) <= set("-: "):              # separator row
                continue
            rows.append(cells)
        else:
            if cur and rows:
                sections.setdefault(cur, []).append(rows)
                rows = []
    if cur and rows:
        sections.setdefault(cur, []).append(rows)
    return sections


def esc(s):
    """Escape LaTeX special characters in a table cell."""
    return (s.replace("\\", r"\textbackslash{}")
             .replace("&", r"\&").replace("%", r"\%").replace("$", r"\$")
             .replace("#", r"\#").replace("_", r"\_").replace("{", r"\{")
             .replace("}", r"\}").replace("~", r"\textasciitilde{}")
             .replace("^", r"\textasciicircum{}"))


def longtable(header_cells, data_rows, colspec, caption=None):
    out = []
    out.append(r"\begin{center}")
    out.append(r"\begin{longtable}{%s}" % colspec)
    out.append(r"\toprule")
    out.append(" & ".join(r"\textbf{%s}" % h for h in header_cells) + r" \\")
    out.append(r"\midrule")
    out.append(r"\endfirsthead")
    out.append(r"\toprule")
    out.append(" & ".join(r"\textbf{%s}" % h for h in header_cells) + r" \\")
    out.append(r"\midrule")
    out.append(r"\endhead")
    out.append(r"\bottomrule")
    out.append(r"\endlastfoot")
    for row in data_rows:
        out.append(" & ".join(esc(c) for c in row) + r" \\")
    out.append(r"\end{longtable}")
    out.append(r"\end{center}")
    return "\n".join(out)


def main():
    sec = parse_note_tables()
    nuclear = sec["Nuclear"][0][1:]            # drop parsed header
    cyto = sec["Cytoplasmic"][0][1:]
    pipe_key = "Pipeline — Nuclear (blob_log, current production)"
    pipe = sec[pipe_key][0][1:]

    L = []
    A = L.append

    # ── Preamble ────────────────────────────────────────────────────────────────
    A(r"""% Combined Report — Condensate Partition Coefficient Pipeline
% Black ink only (no color). Compile with pdfLaTeX (e.g. upload this folder to
% Overleaf and press Recompile). Figures are expected in ./figures.
\documentclass[11pt]{article}
\usepackage[margin=1in]{geometry}
\usepackage{amsmath}
\usepackage{mathptmx}            % Times serif — clean, no color
\usepackage[T1]{fontenc}
\usepackage{graphicx}
\usepackage{booktabs}
\usepackage{longtable}
\usepackage{array}
\usepackage{ragged2e}
\usepackage{titlesec}
\usepackage{fancyhdr}
\usepackage[hidelinks]{hyperref}
\graphicspath{{figures/}}

\titleformat{\section}{\Large\bfseries}{\thesection}{0.6em}{}
\titleformat{\subsection}{\large\bfseries}{\thesubsection}{0.6em}{}
\titleformat{\subsubsection}{\normalsize\bfseries}{}{0em}{}
\setlength{\parskip}{0.5em}
\setlength{\parindent}{0pt}
\setlength{\emergencystretch}{3em}   % absorb long unbreakable code tokens in justified text
\newcommand{\code}[1]{\texttt{\small #1}}

\pagestyle{fancy}\fancyhf{}
\lhead{\small Condensate Partition Coefficient Pipeline}
\rhead{\small Franco Lab \textbullet\ UCLA}
\cfoot{\small\thepage}
\renewcommand{\headrulewidth}{0.4pt}

\begin{document}

% ── Title page ────────────────────────────────────────────────────────────────
\begin{titlepage}
\centering
\vspace*{3cm}
{\Huge\bfseries Condensate Partition\\[0.3em] Coefficient Pipeline\par}
\vspace{1cm}
{\LARGE Combined Report\par}
\vspace{0.5cm}
{\large User Manual \quad\textbullet\quad Methods \quad\textbullet\quad Research Timeline \quad\textbullet\quad JABr Reference Data\par}
\vfill
{\large\bfseries Franco Lab \textbullet\ University of California, Los Angeles\par}
\vspace{0.4cm}
{\large Author: Daniel Chang\par}
{\large Principal Investigator: Elisa Franco\par}
\vspace{0.8cm}
{Spring 2026\par}
\end{titlepage}

\tableofcontents
\newpage
""")

    # ── Part I — Manual ───────────────────────────────────────────────────────
    A(r"""\section{Part I --- User Manual}

This pipeline replaces the manual Imaris workflow for measuring the partition
coefficient (PC) of nuclear biomolecular condensates. Given a two-channel confocal
Z-stack (channel 0 = nuclei, channel 1 = condensate), it segments nuclei and
condensates in 3D, computes the PC, and reports a calibrated value on the manual
reference scale. It runs from a point-and-click graphical interface or the command line.

\begin{center}
\includegraphics[width=0.6\textwidth]{jabr_validation.png}\\
{\small Validation on JABr: calibrated pipeline PC vs the manual Imaris reference.}
\end{center}

\subsection{Installation (one time)}
Install Python 3.10 or 3.11, then install PyTorch for your machine (CPU or your CUDA
version) from \href{https://pytorch.org/get-started/locally/}{pytorch.org}. Example for CUDA 12.1:
\begin{quote}\ttfamily\small pip install torch -{}-index-url\\ \hspace*{1.5em}https://download.pytorch.org/whl/cu121\end{quote}
Then install the remaining dependencies:
\begin{quote}\code{pip install -r requirements.txt}\end{quote}
On first run, Cellpose downloads its model weights automatically (one-time
$\sim$100\,MB). A GPU is recommended ($\sim$10--15\,s per cell); CPU works but is slower.

\subsection{Preparing image files}
The pipeline accepts either one multi-channel TIF (both channels in one file) or two
separate single-channel TIFs. \textbf{Channel order matters:} channel 0 must be nuclei,
channel 1 the condensate. Most axis orderings (ZCYX, CZYX, \dots) are auto-detected; if
results look wrong, swapped channels are the most common cause.

\subsection{Running with the graphical interface}
Launch with \code{python run\_gui.py}. A window with \textbf{Single File} and \textbf{Batch}
tabs opens.
\begin{enumerate}\setlength{\itemsep}{1pt}
\item On the \textbf{Single File} tab, choose the input mode and browse to your data.
\item Choose an output folder (optional).
\item In \textbf{Settings}, set \textbf{Construct = JABr}; leave all other settings at their defaults.
\item Click \textbf{Run Pipeline}. Progress streams in the Output Log; the status reads \textbf{Done} when finished.
\item Open your output folder for the masks, tables, and \code{results.png}.
\end{enumerate}

\subsection{Settings}
For routine JABr analysis you only need to set \textbf{Construct}; the rest have sensible defaults.

\begin{center}
\renewcommand{\arraystretch}{1.2}
\begin{tabular}{>{\RaggedRight}p{3.2cm} >{\RaggedRight}p{1.5cm} >{\RaggedRight}p{8.3cm}}
\toprule
\textbf{Setting} & \textbf{Default} & \textbf{What it does}\\
\midrule
Construct & (none) & Pick JABr for the validated workflow. Auto-selects the \code{blob\_log} detector and the JABr calibration so the reported PC is on the manual scale.\\
Condensate detector & auto & \code{auto} routes JABr to \code{blob\_log}, else Cellpose. Normally left on auto.\\
Top-X\% brightest & 75 & Defines condensed-phase density; trims the dim mask halo to match the manual reference. Leave at 75.\\
Nuclei cell-prob. & $-2.0$ & How aggressively nuclei are detected. $-2.0$ suits the lab's nuclear stain.\\
blob\_log threshold & 0.03 & Condensate spot sensitivity (used when detector = \code{blob\_log}). Leave at 0.03 for JABr.\\
Disable GPU & off & Tick on a laptop with no NVIDIA GPU.\\
\bottomrule
\end{tabular}
\end{center}
\textbf{The one rule for routine use:} set Construct = JABr and run.

\subsection{Reading the outputs}
Every run writes \code{summary.csv} (headline numbers including the calibrated PC),
\code{results.png}, \code{condensate\_masks.tif} and \code{nuclei\_masks.tif} (3D labeled
masks), and per-object/per-slice CSVs. For the bundled sample \code{JABr\_Sample2\_5\_3.tif},
a real run produces nuclear PC raw 4.867 calibrated to 4.813 (manual reference 4.558,
$+5.6\%$). Report the calibrated nuclear PC, and always glance at \code{condensate\_masks.tif}
over the raw channel in Fiji/ImageJ to catch the rare misfire.

\begin{center}
\includegraphics[width=0.95\textwidth]{example_results.png}\\
{\small \code{results.png} from the bundled-sample run.}
\end{center}

\subsection{Command line}
\begin{quote}\ttfamily\small python pipeline.py -{}-roi sample\_data/JABr\_Sample2\_5\_3.tif\\ \hspace*{1.5em}-{}-construct JABr -{}-output my\_output\end{quote}
Separate channels: \code{-{}-nuc nuclei.tif -{}-cond condensate.tif}. Physical volumes in
$\mu m^3$: add \code{-{}-voxel-xy 0.065 -{}-voxel-z 0.3}. Force CPU: \code{-{}-no-gpu}. For
many files (and optional comparison vs a manual reference CSV), use the GUI \textbf{Batch} tab.

\subsection{Troubleshooting}
\begin{center}
\begin{tabular}{>{\RaggedRight}p{5cm} >{\RaggedRight}p{8cm}}
\toprule
\textbf{Symptom} & \textbf{Likely cause / fix}\\
\midrule
PC wildly off / nuclei empty & Channels swapped: confirm ch0 = nuclei, ch1 = condensate.\\
``0 condensates detected'' & Not JABr, or very dim signal. Lower \code{blob\_log} threshold (e.g.\ 0.02).\\
Very slow (minutes/cell) & Running on CPU. Install a CUDA PyTorch build; ensure Disable GPU is unchecked.\\
\code{CUDA out of memory} & Tick Disable GPU, or close other GPU programs.\\
One cell's calibrated PC off & Calibration is population-level; cells vary ($\pm 20\%$ for $\sim$80\% of cells).\\
\bottomrule
\end{tabular}
\end{center}
\newpage
""")

    # ── Part II — Methods ───────────────────────────────────────────────────────
    A(r"""\section{Part II --- Methods}

\subsection{The quantity}
The partition coefficient is PC = condensed-phase density / dilute-phase density
(Fabrini et al.). All intensities are background-subtracted, where the background $B$
is the minimum voxel intensity across the field of view (camera offset).
\begin{itemize}\setlength{\itemsep}{1pt}
\item \textbf{Condensed-phase density} --- mean of $\mathrm{clip}(\text{intensity}-B,0)$
over the brightest 75\% of voxels inside both a condensate mask and a nucleus. Automated
masks include a dim halo a manual tracing would exclude; trimming the bottom 25\% recovers
densities consistent with the manual reference.
\item \textbf{Dilute-phase density} --- mean of the 50 lowest-intensity valid
$10\times10\times10$-voxel patches inside the nuclear dilute region. Averaging many low
patches approximates a manual quiet-region selection and is far more stable than a single random patch.
\end{itemize}

\subsection{Pipeline steps}
\begin{center}
\begin{tabular}{>{\RaggedRight}p{2.4cm} >{\RaggedRight}p{4.6cm} >{\RaggedRight}p{6cm}}
\toprule
\textbf{Step} & \textbf{Method} & \textbf{Why}\\
\midrule
Denoise & Cellpose 3 \code{denoise\_cyto3} per slice & Removes shot noise; PC intensities are read from the raw stack.\\
Segment nuclei & Cellpose 3 \code{cyto3}, \code{do\_3D} & Native 3D nuclear segmentation.\\
Clean nuclei & Connected-component relabel; drop $<1000$-voxel fragments & Cellpose over-splits large nuclei; relabeling restores whole nuclei.\\
Fill voids & Per-nucleus \code{binary\_fill\_holes} (3D + 2D) & Condensates carve donut holes in the nuclear stain exactly where nuclear condensates sit; filling makes them count as intra-nuclear.\\
Detect condensates & \code{blob\_log} (LoG), $\sigma$ 1.5--6.0, threshold 0.03 & Classical, parameter-light detector matched to roughly-spherical condensates; each rendered as a sphere of radius $\sigma\sqrt{3}$.\\
Intra-nuclear gate & Keep a blob if $\geq 50\%$ of its sphere lies inside a nucleus & Replaces an older single-centroid test that mis-handled boundary blobs and crowded fields.\\
Calibrate & Per-construct linear (JABr) / isotonic map & Standardizes the automated PC onto the manual Imaris scale.\\
\bottomrule
\end{tabular}
\end{center}

\subsection{Why there is a calibration}
The automated and manual workflows define the condensed and dilute regions slightly
differently: the automated condensed region runs brighter (brightest mask voxels) and the
automated dilute region runs dimmer (quietest sampled patches) than a hand-drawn Imaris
region. Both push the ratio the same way, so the raw automated PC is systematically about
$3\times$ the manual reference --- but consistently so ($r = 0.94$). A single linear map converts it back:
\begin{quote}\code{calibrated\_PC = 0.3642 * raw\_PC + 3.0405}\quad (JABr)\end{quote}
This is cross-method standardization, not a fudge factor. A linear map is rank-preserving,
so calibration cannot improve correlation --- it only removes the scale/offset bias. The
correlation $r$ is therefore the true accuracy ceiling, and calibration is meaningful only where $r$ is already high.

\subsection{Scope and generalization}
The pipeline is validated and production-ready for JABr. It is construct-specific: the
\code{blob\_log} parameters encode JABr's condensate size and brightness and do not transfer
unchanged. Cross-construct correlation vs Imaris:
\begin{center}
\begin{tabular}{l c >{\RaggedRight}p{7cm}}
\toprule
\textbf{Construct} & \textbf{Pearson r} & \textbf{Status}\\
\midrule
JABr & 0.942 & Production\\
GABr & 0.678 & Mediocre; isotonic calibration, use with caution\\
AABr & 0.122 & Not usable as-is\\
AwtBr & $-0.116$ & Not usable as-is\\
GwtBr & $-0.341$ & Not usable as-is\\
\bottomrule
\end{tabular}
\end{center}
A leave-one-construct-out test (a single learned Cellpose detector, evaluated on a construct
it never trained on) confirmed that one model does not generalize zero-shot: it fails by
under-detection, staying silent on most cells of an unseen construct. The route to extend the
tool is per-construct calibration where correlation is already acceptable, or few-shot
retraining with a handful of labeled ROIs of the new construct --- far cheaper than a manual parameter sweep.
\newpage
""")

    # ── Part III — Timeline ───────────────────────────────────────────────────
    timeline = [
        ("Winter 2026 --- Foundation", [
            ("2026-03-17", "Initial clean pipeline established."),
            ("2026-03-24", "Winter 2026 report finalized; baseline segmentation + PC concept."),
        ]),
        ("April 2026 --- Segmentation survey \\& the PC formula", [
            ("2026-04-22", "Spring segmentation-model survey; 4 candidate models compared."),
            ("2026-04-25", "Background-subtracted PC formula adopted across all scripts."),
            ("2026-04-30", "PC now matches reference: connected-component nuclei fix + lowest-patch dilute density; 50-patch stability fix."),
        ]),
        ("May 2026 --- Calibration, GUI, generalization", [
            ("2026-05-06", "Batch cross-reference pipeline; top-75\\% condensate-density calibration."),
            ("2026-05-07", "GUI built; laptop (CPU) validation; PI meeting."),
            ("2026-05-14", "Crowded-field debugging; max-overlap-nucleus target-cell heuristic."),
            ("2026-05-16", "First fine-tuned Cellpose model; 5-construct accuracy sweep; per-construct calibration."),
            ("2026-05-18", "M187 research poster (PC = 6.297 vs 6.32 reference)."),
            ("2026-05-28", "blob\\_log detector wired into the pipeline; JABr calibration refit; GUI gains construct/detector controls; cytoplasmic PC added."),
        ]),
        ("June 2026 --- Production hardening \\& honest scope", [
            ("2026-06-10", "Discussion section; nuclei void-filling so condensates in donut holes count toward the PC."),
            ("2026-06-11", "Sphere-overlap intra-nuclear gate promoted to production: JABr r 0.926 to 0.942, MAE 14.6\\% to 13.7\\%."),
            ("2026-06-12", "Leave-one-construct-out experiment: a single learned detector does not generalize zero-shot. Few-shot path mapped."),
            ("2026-06-13", "Clean public release: repository, instruction manual, methods writeup, and this combined report."),
        ]),
    ]
    A(r"\section{Part III --- Research Timeline}")
    A("A chronological record of how the pipeline was built. Dates are from the project's version-control history.\n")
    for title, items in timeline:
        A(r"\subsection*{%s}" % title)
        A(r"\begin{center}")
        A(r"\begin{tabular}{>{\RaggedRight}p{2cm} >{\RaggedRight}p{11.5cm}}")
        A(r"\toprule")
        A(r"\textbf{Date} & \textbf{Milestone}\\")
        A(r"\midrule")
        for d, m in items:
            A("%s & %s \\\\" % (d, m))
        A(r"\bottomrule")
        A(r"\end{tabular}")
        A(r"\end{center}")
    A(r"""\subsection*{Where it stands}
Validated, production-ready for JABr ($r = 0.942$, MAE 12.9\%, 79\% of cells within
$\pm 20\%$, n = 28). Delivered: a CLI, a GUI for non-coding lab use, per-construct
calibration, and full documentation. Known boundary, with a plan: the detector is
JABr-specific; extending it needs per-construct calibration or few-shot retraining.
\newpage
""")

    # ── Part IV — JABr Reference Data ───────────────────────────────────────────
    A(r"""\section{Part IV --- JABr Reference Data}
Source: Box \textrightarrow\ Condensate Volume Quantification \textrightarrow\ JABr. The
manual Imaris measurements (nuclear and cytoplasmic) are the ground truth against which the
pipeline is validated; the pipeline results table is the production \code{blob\_log} V4 output
(sphere-overlap gate, linear calibration).
""")
    A(r"\subsection*{Manual reference --- Nuclear (n = 30)}")
    A(longtable(["File", "Condensate Density", "Dilute Density", "Partition Coefficient"],
                nuclear, "l r r r"))
    A(r"\subsection*{Manual reference --- Cytoplasmic (n = 30)}")
    A(longtable(["File", "Condensate Density", "Dilute Density", "Partition Coefficient"],
                cyto, "l r r r"))

    A(r"\subsection*{Pipeline results --- Nuclear (blob\_log V4, sphere-overlap gate)}")
    A(r"""Detector \code{blob\_log} (threshold 0.03, $\sigma$ 1.5--6.0); nuclei via Cellpose
\code{cyto3} with void-filling; intra-nuclear gate $\geq 50\%$ sphere overlap; calibration
$0.3642 \cdot \mathrm{raw} + 3.0405$. n = 28 (Sample3\_3\_10 and Sample3\_3\_15 find 0
intra-nuclear condensates and are excluded).
""")
    A(longtable(["File", "Cond Dens.", "Dilute Dens.", "Bg", "Raw PC", "Cal PC", "Ref PC",
                 r"\textbar err\textbar\,\%"], pipe, "l r r r r r r r"))
    A(r"""\vspace{0.5em}
\textbf{Summary:} n = 28, mean $|\text{err}|$ = 13.7\%, 22/28 within $\pm 20\%$, Pearson $r = 0.942$.

\subsubsection*{Notable failures ($>25\%$ error)}
\begin{itemize}\setlength{\itemsep}{1pt}
\item \textbf{Sample1\_2\_3 ($+63\%$):} smallest reference PC (2.55); the linear calibration
floor (intercept $\sim$3.0) dominates, so the percentage error is large though the absolute
error is only $\sim$1.6 PC units.
\item \textbf{Sample1\_1\_1 ($+43\%$):} a poor nucleus mask drifting into a dim region inflates
the cond/dilute ratio; both \code{blob\_log} and Cellpose fail here.
\item \textbf{Sample3\_3\_3 ($+29\%$):} an over-bright blob set keeps cond density high relative
to the reference PC; not fixed by the gate.
\end{itemize}
\newpage
""")

    A(r"""\subsection*{Visual panels (representative)}
Per-ROI max-intensity Z-projection panels: merged reference (nuclei/condensate) $|$ nuclei
channel $|$ condensate channel $|$ pipeline nuclei mask $|$ pipeline condensate mask $|$
reference (Imaris) condensate mask $|$ classification overlay. A representative subset is
shown; the full set of 28 panels is in the project repository.
""")
    gallery = [
        ("Sample1_4_2", "best case: ref PC 5.05, calibrated 5.06, err 0.3\\%"),
        ("Sample2_5_3", "bundled sample: ref PC 4.56, calibrated 4.67, err 2.4\\%"),
        ("Sample1_1_3", "high PC: ref 19.77, calibrated 19.53, err 1.2\\%"),
        ("Sample2_5_5", "rescued by the sphere-overlap gate: was $+58\\%$, now ref 6.73 / calibrated 5.72, err 15\\%"),
        ("Sample1_2_3", "failure: ref 2.55, calibrated 4.15, err 63\\% (calibration-floor effect)"),
        ("Sample3_3_3", "failure: ref 8.82, calibrated 11.42, err 29\\% (over-bright blobs)"),
    ]
    for name, cap in gallery:
        tex_name = name.replace("_", r"\_")
        A(r"\begin{center}")
        A(r"\includegraphics[width=\textwidth]{%s.png}\\" % name)
        A(r"{\small %s --- %s}" % (tex_name, cap))
        A(r"\end{center}")
        A(r"\vspace{0.4em}")

    A(r"\end{document}")

    OUT.write_text("\n".join(L), encoding="utf-8")
    print("wrote", OUT, "(", len("\n".join(L)), "chars )")


if __name__ == "__main__":
    main()
