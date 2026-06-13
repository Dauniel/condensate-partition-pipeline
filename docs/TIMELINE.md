# Research Timeline — Condensate Partition Coefficient Pipeline

**Franco Lab · UCLA · Daniel Chang · Winter–Spring 2026**

A chronological record of how the pipeline was built, the key decisions, and the
results at each stage. Dates are from the project's version-control history.

---

## Winter 2026 — Foundation

| Date | Milestone |
|---|---|
| 2026-03-17 | Initial clean pipeline established. |
| 2026-03-24 | Winter 2026 report finalized — baseline segmentation + PC concept. |

The Winter work established the core problem (automating the manual Imaris
partition-coefficient measurement) and a first Cellpose-based segmentation.

---

## Spring 2026 — Method development

### April — Segmentation survey & the PC formula

| Date | Milestone |
|---|---|
| 2026-04-22 | Spring segmentation-model survey begun; 4 candidate models compared. |
| 2026-04-23 | Reference-PC investigation; old pipeline archived into `winter_implementation/`. |
| 2026-04-25 | **Background-subtracted PC formula** adopted across all scripts. |
| 2026-04-28 | `spring_implementation` pipeline created; nuclei over-segmentation fix. |
| 2026-04-30 | **PC now matches reference** — connected-component nuclei fix + lowest-patch dilute density. |
| 2026-04-30 | Dilute-density stability fix: average **50** patches instead of one random patch. |

This month produced the measurement that works: the Fabrini-style PC with
background subtraction, a 50-patch dilute estimate, and clean connected-component
nuclei. First single-cell agreement with the manual reference.

### May — Calibration, GUI, and the generalization question

| Date | Milestone |
|---|---|
| 2026-05-05 | Pipeline diagram + design notes. |
| 2026-05-06 | Batch cross-reference pipeline; **top-75 % condensate-density** calibration. |
| 2026-05-07 | **GUI built**; laptop (CPU) validation; PI meeting prep. |
| 2026-05-14 | Crowded-field debugging; **max-overlap-nucleus** target-cell heuristic. |
| 2026-05-16 | First fine-tuned Cellpose model; **5-construct accuracy sweep**; per-construct isotonic calibration scaffolding. |
| 2026-05-18 | v3 construct-balanced training; M187 research **poster** (PC = 6.297 vs 6.32 reference). |
| 2026-05-19 | Cross-construct head-to-head (blob_log vs Cellpose V3). |
| 2026-05-21 | Spring 2026 research report draft + PDF. |
| 2026-05-22 | Methods revisions: dataset rationale, Fabrini citation, limitations. |
| 2026-05-28 | **`blob_log` detector wired into the pipeline**; JABr/JABr_4arm/Tornado calibrations refit; GUI gains construct/detector controls; cytoplasmic PC added. |

May moved from "works on one cell" to "works across the JABr batch," added the
point-and-click GUI for the lab, and surfaced the central scientific finding that
the method is **construct-specific**.

### June — Production hardening & honest scope

| Date | Milestone |
|---|---|
| 2026-06-10 | Discussion section + references; **nuclei void-filling** so condensates in donut holes count toward the PC. |
| 2026-06-11 | **Sphere-overlap (≥50 %) intra-nuclear gate** promoted to production: JABr r 0.926 → **0.942**, MAE 14.6 % → **13.7 %**. |
| 2026-06-11 | V4 cross-construct test confirms construct-specificity (JABr 0.94, GABr 0.68, AABr/AwtBr/GwtBr fail). |
| 2026-06-12 | **Leave-one-construct-out experiment**: a single learned detector does not generalize zero-shot (fails by under-detection). Mapped the few-shot path to extend beyond JABr. |
| 2026-06-13 | Clean public release: this repository, instruction manual, and methods writeup. |

---

## Where it stands

- **Validated, production-ready for JABr**: r = 0.942, MAE 12.9 %, 79 % of cells
  within ±20 % of the manual Imaris reference (n = 28).
- **Delivered**: a CLI (`pipeline.py`), a GUI (`run_gui.py`) for non-coding lab
  use, per-construct calibration, and full documentation.
- **Known boundary, with a plan**: the detector is JABr-specific; extending it to
  other constructs needs per-construct calibration or few-shot retraining, not a
  manual parameter sweep. See [METHODS.md](METHODS.md) § 5.
