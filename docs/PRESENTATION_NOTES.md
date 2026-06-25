# Presentation Notes — Condensate Partition Coefficient Pipeline

# June 19

**Daniel Chang · Franco Lab · UCLA · PI: Elisa Franco**

---

## 0. Start

> "I built an automated pipeline that measures how strongly a protein concentrates
> into nuclear condensates — the *partition coefficient* — directly from confocal
> Z-stacks. It replaces a slow, manual Imaris workflow. On our main construct,
> JABr, it matches the human-traced reference with a correlation of **0.94** and a
> mean error of about **13%**, in seconds per cell instead of minutes of manual
> tracing."

Three takeaways:
1. **What** — automated partition-coefficient measurement from microscopy images.
2. **How good** — r = 0.94, MAE ~13%, 79% of cells within ±20% (n = 28, JABr).
3. **Scope** — production-ready for JABr; construct-specific by design, with a mapped-out path to extend it.

---

## 1. The problem (why this matters)

- Biomolecular **condensates** are membraneless droplets where biomolecules
  concentrate. How strongly a molecule partitions into them is a key biophysical
  readout.
- The **partition coefficient (PC)** = how concentrated the molecule is *inside*
  the condensate vs. *outside* (the dilute phase).
- **Today this is measured by hand in Imaris** — a researcher traces condensates
  and picks dilute regions slice by slice. It's slow, subjective, and doesn't scale.
- **Goal:** reproduce that measurement automatically, objectively, and at batch scale.

---

## 2. What the PC actually is (the one equation)

```
PC = condensed-phase density / dilute-phase density
```

Following the **Fabrini et al.** intensity-based definition. Key details:

- All intensities are **background-subtracted** (background = the camera's minimum
  voxel value, i.e. the offset).
- **Condensed-phase density** = mean intensity of the **brightest 75%** of voxels
  that are inside *both* a condensate *and* a nucleus.
  - *Why top 75%?* Automated masks include a dim halo + dark interior pixels a
    human would exclude. Trimming the bottom 25% recovers densities that match the
    manual tracing. (Validated by a sweep: bias went from −16% to +3%.)
- **Dilute-phase density** = mean of the **50 quietest** 10×10×10-voxel patches
  inside the nucleus but outside condensates.
  - *Why 50 patches?* Mimics how a researcher picks a quiet representative region;
    averaging many patches is far more stable than one random patch.

On "is picking the brightest/quietest cherry-picking?" — it mirrors what the manual
Imaris operator does by eye, and both choices are applied *consistently*, which is
what makes the calibration valid (§5).

---

## 3. The pipeline, step by step

1. **Denoise** — Cellpose 3 `denoise_cyto3` on each slice. The PC
   intensities are read from the **raw** stack, not the denoised one — denoising
   only helps draw cleaner boundaries.
2. **Segment nuclei in 3D** — Cellpose 3 `cyto3`, native 3D.
3. **Clean nuclei** — relabel connected components, drop fragments < 1000 voxels
   (Cellpose over-splits big nuclei around internal texture).
4. **Fill voids** — condensates push out the nuclear stain, leaving "donut holes."
   Filling them lets nuclear condensates count as intra-nuclear.
5. **Detect condensates** — `blob_log` (Laplacian-of-Gaussian blob detector),
   σ = 1.5–6.0. A classical, parameter-light detector well-matched to roughly
   spherical condensates. Each blob is rendered as a sphere.
6. **Intra-nuclear gate** — keep a blob only if **≥ 50% of its volume** is inside a
   nucleus. This replaced an older single-center-point test and bumped r from
   0.926 → **0.942** while cutting error.
7. **Measure** — per-object volumes, intensities (regionprops + voxel counting).
8. **Compute PC**, then **calibrate** onto the manual Imaris scale.

---

## 4. Why there's a calibration

- The automated method and manual Imaris **define the regions slightly
  differently**: automated condensed region runs brighter, automated dilute region
  runs dimmer. Both push the ratio the same way.
- Result: the **raw automated PC is ~3× the manual value — but consistently so**
  (r = 0.94).
- Because the bias is consistent, **one linear map** converts it back:

```
calibrated_PC = 0.3642 · raw_PC + 3.0405        (JABr)
```

This is an ordinary least-squares fit of the manual Imaris PC against the raw
automated PC over the 28 JABr cells. Calibration is cross-method *standardization*,
not a fudge factor: a linear map is rank-preserving, so it **cannot improve
correlation** — it only removes the scale/offset. The correlation r is the true
accuracy ceiling, and calibration only makes sense where r is already high.

---

## 5. Validation results

On the **28 JABr cells** with a manual Imaris reference:

| Metric | Value |
|---|---|
| Pearson r | **0.942** |
| Mean absolute error (calibrated) | **12.9%** |
| Cells within ±20% of manual | **79%** |

Bundled sample sanity check: `JABr_Sample2_5_3.tif` → raw PC 4.867 → calibrated
**4.813** vs. manual reference 4.558 (**+5.6%**).

---

## 6. Scope & limitations

The tool is validated and production-ready for JABr. It is **construct-specific by
design** — the detector parameters encode JABr's condensate size/brightness.

Cross-construct raw correlation vs. Imaris:

| Construct | r | Verdict |
|---|---|---|
| **JABr** | **0.94** | Production |
| GABr | 0.68 | Mediocre — use with caution |
| AABr | 0.12 | Not usable as-is |
| AwtBr / GwtBr | ≈ 0 / negative | Not usable as-is |

- **Root cause:** *scale.* Wild-type (wt) constructs have condensates up to ~2×
  larger than the detector's σ ceiling, so a fixed detector misses or fragments them.
- A leave-one-construct-out test (detector fine-tuned on Imaris labels, held out one
  construct) showed it **fails by under-detection** — staying silent on most cells
  of a construct it never saw. So a single detector does not generalize zero-shot.
- **Path forward:**
  - Per-construct **calibration** where correlation is already acceptable.
  - **Few-shot retraining** — adding a handful of labeled ROIs of a new construct
    should restore detection; far cheaper than a manual parameter sweep.

On "so it only works on one construct?" — it's validated on one, and the limits are
*measured*, not guessed. Scoped to do JABr well, with exactly what it would take to
extend it documented. That's more useful than a tool that claims to work everywhere
and silently fails.

---

## 7. Deliverables

Three interfaces, one shared code path (CLI, GUI, and web app all call the same
`run_full_pipeline` / `run_batch`):

1. **Web app** (`webapp.py`, Flask) — configure, run, and view in one browser page.
   Single cells *and* whole batches. Live progress streaming; batch results as
   clickable chips with r / RMSE / MAE; native OS file dialogs.
2. **Desktop GUI** (`run_gui.py`, Tkinter) — point-and-click for non-coders in the lab.
3. **CLI** (`pipeline.py`) — scriptable, one bundled sample to reproduce.

**The interactive Z-viewer (`zviewer.html`)** — the demo highlight:
- Self-contained HTML; scroll through every Z-slice.
- Four panels: **raw image | mask overlay | condensate masks | nuclei masks**.
- Live **per-slice stats** + a **stack summary** (the PC, densities, object counts).
- **Intensity histograms** that update per slice + an overall one.
- Portable and shareable — just a file you open in any browser.

---

## 8. Anticipated questions

- **How long per cell?** → ~10–15 s on a GPU; CPU works but slower.
- **Why blob_log instead of a learned detector for condensates?** → It's classical,
  parameter-light, and matches spherical condensates well; and the leave-one-out
  test showed a learned detector doesn't generalize zero-shot anyway. Fewer moving
  parts, easier to trust.
- **Why Cellpose for nuclei?** → Strong general cell/nucleus segmentation, works in
  native 3D, with a denoising model in the same toolkit. Note: it's Cellpose `cyto3`
  plus post-processing (connected-component relabel + void-filling), not a raw call.
- **What's the human-in-the-loop part?** → Glance at the saved masks / zviewer to
  confirm segmentation is sane. It's an assistant, not a black box.
- **Could the calibration drift over time / microscopes?** → Possibly — it's tied
  to acquisition. The right move is a small re-validation set per setup; the
  framework supports re-fitting trivially.
- **What's next?** → Unify calibration across condensate-forming constructs
  (excluding non-condensate ones like Tornado), and few-shot retraining to extend
  detection beyond JABr.

---

## 9. Closing

> "So: a manual, subjective, hours-long Imaris measurement is now an automated,
> reproducible, seconds-per-cell pipeline that matches the human reference at
> r = 0.94 on JABr — with the limits measured and a clear path to extend it. Happy
> to demo the live viewer."

---

### Important numbers
- **r = 0.942**, **MAE 12.9%**, **79% within ±20%**, **n = 28** (JABr).
- Raw PC ≈ **3×** manual; linear calibration `0.3642·x + 3.0405`.
- Detector: `blob_log`, σ 1.5–6.0; intra-nuclear gate **≥ 50%**.
- Condensed = brightest **75%**; dilute = **50** quietest 10³ patches.
- Gate improvement: r **0.926 → 0.942**.
</content>
