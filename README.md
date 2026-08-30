# TreeTimeAlign

**A stem-based workflow for intertemporal co-registration of TLS forest point clouds**

## Abstract

TreeTimeAlign computes a rigid transform with all 6 degrees of freedom that brings a LiDAR scan of a forest plot into alignment with a second scan of the same plot taken at a different time. It is **target-less** and **ground-independent**: the alignment is derived entirely from the geometry of individual tree stems, which are modeled by fitting ellipses to the stem sections within horizontal and axis-perpendicular slices of each instance segmented tree. The workflow begins by automatically identifying corresponding trees between epochs for subsequent misalignment analysis. The workflow is written as six Quarto documents that hand off to one another through checkpoint files, plus a batch driver that runs the whole sequence unattended over any number of scan pairs. This README describes what the program does and how it is organised.

## Table of Contents

- [What this repository contains](#what-this-repository-contains)
- [Structure of the workflow](#structure-of-the-workflow)
- [The batch driver](#the-batch-driver)
- [Checkpoints](#checkpoints)
- [Stage 1: Initial alignment and tree identification](#stage-1-initial-alignment-and-tree-identification)
- [Stage 2: Full 3D rotation](#stage-2-full-3d-rotation)
- [Stage 3: Vertical translation](#stage-3-vertical-translation)
- [Stage 4: Yaw from the inter-stem line network](#stage-4-yaw-from-the-inter-stem-line-network)
- [Stage 5: Residual horizontal translation](#stage-5-residual-horizontal-translation)
- [Stage 6: Composition and export](#stage-6-composition-and-export)
- [Truth values and testing](#truth-values-and-testing)
- [Running the workflow](#running-the-workflow)
- [Requirements](#requirements)

## What this repository contains

| File | Role |
|------|------|
| `RunR.R` | Batch driver. Runs all six stages over a list of scan pairs without supervision, with logging, checkpointing and a report. |
| `Stage1_TemporalPairs.qmd` | Rough vertical alignment; ellipse fitting; temporal pair identification; the rotation pivot. |
| `Stage2_XYZRotation.qmd` | Full 3D rotation from stem axis vectors. |
| `Stage3_ZTranslation.qmd` | Axis-perpendicular ellipse refit; vertical offset from stem-property curve matching. |
| `Stage4_ZRotation.qmd` | Yaw from the inter-stem line network. |
| `Stage5_XYTranslation.qmd` | Residual horizontal translation. |
| `Stage6_Export.qmd` | Composes every solved transform and applies it once, writing the aligned cloud. |
| `RANSACellipse.R` | Shared ellipse-fitting and point-matching library, sourced by Stages 1 and 3. |
| `TestSampleGen.qmd` | Generates synthetic test samples with a known transform, for evaluating the workflow. |

## Structure of the workflow

Each stage solves one component of the transform and records it. Each subsequent stage applies the components solved before it, so that it measures the residual misalignment rather than the original: Stage 3 applies the Stage 2 rotation before refitting ellipses, Stage 4 applies the Stage 3 vertical offset before solving the yaw, and Stage 5 applies the Stage 4 yaw before measuring the horizontal displacement. What each stage applies is applied to the fitted geometry it works on, not written back to the cloud.

Stage 6 then applies the whole composed transform to the original point cloud in a single pass. No stage writes a transformed cloud to disk, so the exported result is one matrix applied to the file as it was read, not an accumulation of five separate rewrites.

### Two key design consequences:

1. **Efficient computation**: Only Stages 1, 3 and 6 load LAS data. Stages 2, 4 and 5 work on fitted ellipse geometry, so re-solving a component costs seconds rather than gigabytes of input.

2. **Auditability**: Because no stage consumes its inputs in place, the rotation matrices, ellipse tables, per-pair residuals and filter verdicts all survive to the end of the run and can be inspected.

### Single rotation pivot

A single rotation pivot, `pivot1`, is fixed once in Stage 1 and never recomputed. By default it is read from a marker point carried in the input file; where none is present it is computed as the mid-point of the extent of the first cloud in all three axes. Every rotation in the workflow turns about that point, so no rotation introduces a translation of its own.

## The batch driver

`RunR.R` runs the six stages over a list of scan pairs without supervision. It reads each stage's Quarto document, extracts its executable chunks, and evaluates them in a fresh environment per scan pair. It is not a wrapper around `quarto render`: it parses the chunks directly, which is what allows it to skip figure and diagnostic chunks, inject parameters, and continue after a failure on one pair.

### What the driver provides:

- **Preflight validation**: Every stage file, the ellipse library, the input clouds, the output directories and the parameter map are checked *before* any computation starts. A typo in a path or a chunk label stops the run in seconds rather than after an hour of fitting.

- **Parameter injection**: Each stage's parameters live in dedicated `params-*` chunks placed immediately before the chunk that uses them. The driver overrides named parameters after the defaults are assigned and before they are read. The defaults remain in the `.qmd` files as the single definition; the driver names only what a batch changes.

- **Chunk classification**: Figure and diagnostic chunks are identified by label and skipped in batch mode.

- **Isolation and continuation**: Each pair runs in its own environment. A failure is recorded, the report is rewritten, and the next pair starts.

- **Progress and logging**: Per-pair console progress, and optionally a per-pair log file capturing every stage printout.

- **Reporting**: A CSV row per pair, rewritten after every pair so an interrupted batch still leaves a valid report.

## Checkpoints

Each stage writes a checkpoint, a self-contained `.rds` file containing the cloud, ellipses, solved components, and anything else later stages might need. Checkpoints are named `*_S1`, `*_S3` and `*_S6` because only Stages 1, 3 and 6 load point cloud data; Stages 2, 4 and 5 take a checkpoint from the previous stage as input and update it in-place in memory.

The checkpoint from Stage 1 includes the instance-segmented cloud, the rough alignment matrix, the pooled ellipse table, and a tree-pair summary with geospatial and structural scores. Stage 3 inherits it and adds refit ellipses from the axis-perpendicular slices. Every stage appends its own output to the checkpoint, so no information is lost and nothing is read twice.

If you want to re-solve just one component (e.g. Stage 2) after tuning its parameters, you load the Stage 1 checkpoint, apply everything before Stage 2, and run Stage 2 forward—without touching Stages 1, 3–5, or the cloud.

## Stage 1: Initial alignment and tree identification

### Stage 1A: Rough vertical alignment

Stage 1 begins by detecting whether the two clouds are grossly misaligned (rotated sideways or upside-down) and applying a coarse vertical alignment. This is done by:

1. Fitting the two clouds' stem density distributions along the vertical and along the forest axis (the axes of greatest and second-greatest variance in the horizontal plane).
2. Aligning them by rotation about the vertical only.

### Ellipse fitting

Ellipses are fitted independently to the stem sections in each cloud. For each tree and each cloud:

1. Every stem point is sliced into horizontal sections every 50 cm, plus one axis-perpendicular section (i.e. perpendicular to the forest axis).
2. Each section is fitted with a RANSAC algorithm (Fischler & Bolles, 1981) that scores points by how far they deviate from the nearest point on an ellipse and selects the ellipse that maximizes the count of inliers.
3. Ellipses are retained if they pass outlier and aspect-ratio checks.

The result is a table of ellipses with fitted centre, axis lengths, rotation, and slice metadata, one row per ellipse, one to ten ellipses per tree depending on its height and stem straightness.

### Tree identification

Trees present in both epochs are identified by matching their ellipse sets:

1. For each tree in epoch 1, every candidate tree in epoch 2 is scored on the spatial proximity and structural similarity of their ellipse sets.
2. Candidates are ranked, and the top-ranked candidate (if present) is returned as the match.
3. Isolated trees or trees misaligned by more than ~3 m are flagged as singletons and excluded from subsequent fitting.

The matching yields a pairwise table of scores in multiple dimensions: distance, structural similarity (Euclidean distances between ellipses), and others. These scores are carried forward to Stage 3, where they feed a random forest trained on samples with known truth, to weight the stem contributions to the vertical offset.

## Stage 2: Full 3D rotation

Stage 2 applies the rotation pivot, then computes the mean stem axis vector for each matched tree in each epoch. The stem axis is the principal axis of the cloud of all ellipse centres (in 3D) for that tree—a regression line through a weighted point cloud.

The axis vector for each tree is a direction in 3D; for each matched pair, the rotation that aligns the epoch-1 axis to the epoch-2 axis is computed. Axes with low confidence (because the ellipse set is small or the points are poorly aligned to the line) are downweighted, and the final rotation is a weighted mean of the per-tree rotations, iteratively refined by MAD (median absolute deviation) outlier removal.

The rotation matrix, `R_xyz`, is applied to ellipse centres and stem geometry only; the point cloud itself is not transformed at this stage.

## Stage 3: Vertical translation

Stage 3 applies the rotation solved in Stage 2 to the cloud, then refits ellipses in slices that are perpendicular to the rotated forest axis. This axis-perpendicular refit improves ellipse stability after rotation.

Once refitted, the epoch-1 and epoch-2 ellipse sets for each matched tree are thinned to a common height range. For each epoch, a representative stem property (mean major-axis length, minor-axis length, or area) is fitted as a function of height using robust regression with Tukey biweights (Beaton & Tukey, 1974). The vertical offset that aligns these two curves is solved by minimizing the sum of squared residuals.

This produces a single vertical offset, `t_z`, that applies to all matched trees.

### Random forest weighting

Where ground truth is available (e.g. from synthetic samples or measurements), a random forest is trained offline on a set of samples to predict the error in each tree's implied vertical offset. The forest is trained on the twelve scores from the tree-pair summary (spatial and structural similarity), with response variable the log-transformed error:

$$y_i = \log\!\big(\max(|s^* - s_i|,\; \varepsilon)\big)$$

The log scale is used because errors span millimetres to tens of metres, and a forest fitted on the raw scale would spend its capacity separating very wrong candidates from extremely wrong ones.

At run time, the model supplies a predicted error for each tree. For a pool of candidates with predicted errors $\hat{e}_i$, each is weighted:

$$w_i = \exp\!\left(-k \left(\frac{\hat{e}_i - \min_j \hat{e}_j}{\delta}\right)^2\right)$$

The reduction is the weighted mean of the pool's implied shifts.

**Important**: The model reads no truth at run time. It is trained offline on samples that have one and thereafter predicts from the twelve scores alone, so it applies unchanged to a pair with no known truth.

## Stage 4: Yaw from the inter-stem line network

Stage 4 applies the vertical offset, then averages each stem's ellipse centres to one representative horizontal point per stem per epoch. Every pair of stems defines an edge, and the angular difference between the epoch-1 and epoch-2 versions of that edge is measured. A pure yaw rotates every edge by the same angle, so each edge is an independent estimate of the same quantity. Stems whose edges disagree with the consensus are removed iteratively by MAD, and the rotation is the mean of the surviving edge-wise differences.

## Stage 5: Residual horizontal translation

The Stage 4 yaw is applied to the averaged epoch-1 centres, the horizontal displacement to each epoch-2 partner is measured, stems disagreeing with the consensus are removed by MAD on each component, and the surviving displacements are averaged into `t_trans_xy`.

## Stage 6: Composition and export

Stage 6 rebuilds every solved component and applies the composed transform to the cloud in one pass:

$$p' = T_z\big(R_{\text{align+xyz}}\,p + t_{\text{align+xyz}}\big) + t$$

where $R_{\text{align+xyz}}$ is the Stage 1A rough vertical alignment composed with the Stage 2 rotation, both about `pivot1`; $T_z$ is the Stage 4 yaw about the same pivot; and $t$ carries the Stage 5 horizontal translation and the Stage 4 vertical offset. The aligned cloud is written to disk, and the composed matrix is recorded alongside it.

## Truth values and testing

### How a known truth is constructed

`TestSampleGen.qmd` builds test samples whose correct answer is known by construction. It takes a reference cloud, applies a chosen rigid transform (a translation, a tilt of chosen magnitude and direction, and a yaw) and writes the result as a new sample. Because the workflow's task is to undo that transform, its inverse is the answer the workflow should produce.

The generator also writes eight **marker points** into the sample, at the vertices of the reference cloud's bounding box, twice over: once at their true positions, and once carried through the same forward transform. Each marker carries an identifier so that the two copies pair unambiguously. These markers give a second, independent way to score a run: apply the transform the workflow solved to the displaced markers and measure how far each lands from its true twin. This is a displacement in metres measured on geometry, rather than a comparison of matrices, so it is unaffected by any difference in pivot convention or factorisation between the generator and the workflow.

The generator can sweep a list of reference clouds against a list of tilt magnitudes, producing one sample per combination together with a single consolidated table holding every sample's known transform and marker coordinates. Marker points carry no `TreeID` and are removed by Stage 1 before any stem-wise computation, so their presence does not influence the alignment.

### What happens when there is no truth

**On a real scan pair the workflow has no truth and does not invent one.** This is worth stating plainly, because the same code paths run in both cases.

- Truth values are read from an external reference table written by the generator, or from marker points present in the sample. Neither exists for an ordinary scan pair.

- When neither is available, every truth-derived quantity is reported as `NA` and every error column is empty. The workflow does not substitute a default, an estimate, or a value carried over from another sample.

- **No solved transform depends on a truth value.** Every component, meaning the rough alignment, the rotation, the vertical offset, the yaw and the horizontal translation, is derived from the two clouds alone. Truth is used for *scoring* a run, never for producing one. A pair with no truth yields exactly the same transform it would have yielded with one.

- Where the source of a truth value could be ambiguous, the run records which source actually supplied it, distinct from which source was requested, so a silent fallback cannot be mistaken for a successful lookup.

## Running the workflow

### Batch mode

Set the scan pairs and the output locations at the top of `RunR.R` and source it:

```r
source("RunR.R")
```

The driver validates every input, then runs each pair through all six stages, writing checkpoints, an aligned cloud, a per-pair log and a report row.

### Interactively, one stage at a time

Open the stage you want, point its load path at the previous stage's checkpoint, and run the document. Because checkpoints are self-describing, this works in a fresh R session with nothing else loaded. This is the path to use when diagnosing a stage or tuning its parameters: the documents render every figure and diagnostic table that the batch skips.

## Requirements

R ≥ 4.3, with the following packages:

- `lidR`
- `data.table`
- `dplyr`
- `tidyr`
- `tibble`
- `ggplot2`
- `plotly`
- `dbscan`
- `cli`
- `htmltools`
- `ranger` (for the optional random forest weighting)

Quarto is needed only to render the documents; the batch driver parses them directly and does not require it.

## References

Beaton, A. E., & Tukey, J. W. (1974). The fitting of power series, meaning polynomials, illustrated on band-spectroscopic data. *Technometrics*, 16(2), 147–185.

Breiman, L. (2001). Random Forests. *Machine Learning*, 45(1), 5–32. https://doi.org/10.1023/A:1010933404324

Fischler, M. A., & Bolles, R. C. (1981). Random Sample Consensus: A Paradigm for Model Fitting with Applications to Image Analysis and Automated Cartography. *Communications of the ACM*, 24(6), 381–395. https://doi.org/10.1145/358669.358692

Sen, P. K. (1968). Estimates of the Regression Coefficient Based on Kendall's Tau. *Journal of the American Statistical Association*, 63(324), 1379–1389. https://doi.org/10.1080/01621459.1968.10480934

Theil, H. (1950). A Rank-Invariant Method of Linear and Polynomial Regression Analysis. *Indagationes Mathematicae*, 12(85), 173.

Wright, M. N., & Ziegler, A. (2017). ranger: A Fast Implementation of Random Forests for High Dimensional Data in C++ and R. *Journal of Statistical Software*, 77(1), 1–17. https://doi.org/10.18637/jss.v077.i01
