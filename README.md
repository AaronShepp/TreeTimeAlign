# TreeTimeAlign

### A stem-based workflow for intertemporal co-registration of TLS forest point clouds

**Aaron Sheppard**

## Abstract

TreeTimeAlign computes a rigid transform with all 6 degrees of freedom that brings a LiDAR scan of a forest plot into alignment with a second scan of the same plot taken at a different time. It is *target-less* and *ground-independent*: the alignment is derived entirely from the geometry of individual tree stems, which are modeled by fitting ellipses to the stem sections within horizontal and axis-perpendicular slices of each instance segmented tree. The workflow begins by automatically identifying corresponding trees between epochs for subsequent misalignment analysis. The workflow is written as six Quarto documents that hand off to one another through checkpoint files, plus a batch driver that runs the whole sequence unattended over any number of scan pairs. This READme document describes what the program does and how it is organised.

## Table of Contents

- [What this repository contains](#what-this-repository-contains)
- [Structure of the workflow](#structure-of-the-workflow)
- [The batch driver](#the-batch-driver)
  - [Checkpoints](#checkpoints)
- [Method](#method)
  - [Inputs and prerequisites](#inputs-and-prerequisites)
  - [Outlier rejection](#outlier-rejection)
  - [Ellipse fitting](#ellipse-fitting)
  - [Stage 1: rough vertical alignment, ellipse fitting, temporal pairing](#stage-1-rough-vertical-alignment-ellipse-fitting-temporal-pairing)
  - [Stage 2: three-dimensional rotation from stem axes](#stage-2-three-dimensional-rotation-from-stem-axes)
  - [Stage 3: vertical offset from stem-property curves](#stage-3-vertical-offset-from-stem-property-curves)
  - [Learned error prediction](#learned-error-prediction)
  - [Stage 4: yaw from the inter-stem line network](#stage-4-yaw-from-the-inter-stem-line-network)
  - [Stage 5: residual horizontal translation](#stage-5-residual-horizontal-translation)
  - [Stage 6: composition and export](#stage-6-composition-and-export)
- [Truth values and testing](#truth-values-and-testing)
- [Running the workflow](#running-the-workflow)
- [Requirements](#requirements)
- [References](#references)

---

## What this repository contains

| File | Role |
|---|---|
| `RunR.R` | Batch driver. Runs all six stages over a list of scan pairs without supervision, with logging, checkpointing and a report. |
| `Stage1_TemporalPairs.qmd` | Rough vertical alignment; ellipse fitting; temporal pair identification; the rotation pivot. |
| `Stage2_XYZRotation.qmd` | Full 3D rotation from stem axis vectors. |
| `Stage3_ZTranslation.qmd` | Axis-perpendicular ellipse refit; vertical offset from stem-property curve matching. |
| `Stage4_ZRotation.qmd` | Yaw from the inter-stem line network. |
| `Stage5_XYTranslation.qmd` | Residual horizontal translation. |
| `Stage6_Export.qmd` | Composes every solved transform and applies it once, writing the aligned cloud. |
| `RANSACellipse.R` | Shared ellipse-fitting and point-matching library, sourced by Stages 1 and 3. |
| `TestSampleGen.qmd` | Generates synthetic test samples with a known transform, for evaluating the workflow ([Truth values and testing](#truth-values-and-testing)). |

---

## Structure of the workflow

Each stage solves one component of the transform and records it. Each subsequent stage applies the components solved before it, so that it measures the residual misalignment rather than the original: Stage 3 applies the Stage 2 rotation before refitting ellipses, Stage 4 applies the Stage 3 vertical offset before solving the yaw, and Stage 5 applies the Stage 4 yaw before measuring the horizontal displacement. What each stage applies is applied to the fitted geometry it works on, not written back to the cloud.

Stage 6 then applies the whole composed transform to the original point cloud in a single pass. No stage writes a transformed cloud to disk, so the exported result is one matrix applied to the file as it was read, not an accumulation of five separate rewrites.

Two consequences follow:

- Only Stages 1, 3 and 6 load LAS data. Stages 2, 4 and 5 work on fitted ellipse geometry, so re-solving a component costs seconds rather than gigabytes of input.
- Because no stage consumes its inputs in place, the rotation matrices, ellipse tables, per-pair residuals and filter verdicts all survive to the end of the run and can be inspected.

A single rotation pivot, `pivot1`, is fixed once in Stage 1 and never recomputed. By default it is read from a marker point carried in the input file; where none is present it is computed as the mid-point of the extent of the first cloud in all three axes. Every rotation in the workflow turns about that point, so no rotation introduces a translation of its own.

---

## The batch driver

`RunR.R` runs the six stages over a list of scan pairs without supervision. It reads each stage's Quarto document, extracts its executable chunks, and evaluates them in a fresh environment per scan pair. It is not a wrapper around `quarto render`: it parses the chunks directly, which is what allows it to skip figure and diagnostic chunks, inject parameters, and continue after a failure on one pair.

What the driver provides over running the documents by hand:

**Preflight validation**
Every stage file, the ellipse library, the input clouds, the output directories and the parameter map are checked *before* any computation starts. A typo in a path or a chunk label stops the run in seconds rather than after an hour of fitting.

**Parameter injection**
Each stage's parameters live in dedicated `params-*` chunks placed immediately before the chunk that uses them. The driver overrides named parameters after the defaults are assigned and before they are read. The defaults remain in the `.qmd` files as the single definition; the driver names only what a batch changes.

**Chunk classification**
Figure and diagnostic chunks are identified by label and skipped in batch mode.

**Isolation and continuation**
Each pair runs in its own environment. A failure is recorded, the report is rewritten, and the next pair starts.

**Progress and logging**
Per-pair console progress, and optionally a per-pair log file capturing every stage printout.

**Reporting**
A CSV row per pair, rewritten after every pair so an interrupted batch still leaves a valid report.

### Checkpoints

Each stage writes an `.rds` checkpoint containing every object the downstream stages need. The checkpoints are what allow the six Quarto documents to be run individually.

- **Any stage can be opened and run on its own.** Load the previous stage's checkpoint and the document reproduces that stage exactly, including every figure, diagnostic table and intermediate object, without repeating the work that came before it.
- **Checkpoints are self-describing.** The manifest of what a checkpoint contains is stored *inside* the file, so a checkpoint can be loaded in a fresh R session with nothing else defined. Resuming at an arbitrary point does not require first reconstructing state from the stages that preceded it.
- **Stage 3 writes two checkpoints.** The first is taken immediately after the axis-perpendicular ellipse fit, as ellipse fitting is the most expensive process in the workflow. So the whole matching and filtering half of Stage 3 can be re-tuned without ellipse refitting.

The driver runs the stages over many pairs; the .qmd files are where a single stage is inspected, diagnosed and tuned. Both paths execute the same code, so a change in a .qmd is automatically included by the driver.

---

## Method

### Inputs and prerequisites

The workflow takes two `.las`/`.laz` files of the same plot, each carrying a `TreeID` scalar field that isolates individual stems. Only *instance* segmentation is required: stems must be separated from one another, but individual returns need not be labelled as stem rather than branch or foliage. No ground classification is used at any point.

Point cloud input and output use `lidR` (Roussel et al., 2020); neighbourhood searches use `dbscan` (Hahsler et al., 2019); data manipulation uses `dplyr` (Wickham et al., 2026); figures use `ggplot2` (Wickham, 2016) and `plotly` (Sievert, 2020).

### Outlier filtering

Several steps reduce a set of scalar estimates to one value, and each first removes outliers using the median absolute deviation. For a set $\{x_i\}$ with median $\tilde{x}$ and $\mathrm{MAD} = \mathrm{median}(|x_i - \tilde{x}|)$ (Rousseeuw & Croux, 1993), each value is scored

$$
z_i = \frac{0.6745\,(x_i - \tilde{x})}{\mathrm{MAD}},
$$

and flagged when $|z_i|$ exceeds a per-step multiplier (Iglewicz & Hoaglin, 1993). The median and MAD are themselves resistant to extreme values, so the criterion is not distorted by the outliers it is meant to find. The multiplier is set separately at each filtering step, since the tolerable dispersion differs between, for example, angles measured across a stem network and vertical offsets measured across signals.

### Ellipse fitting

Stem cross-sections are recovered by fitting ellipses to thin slices of each segmented stem. Fitting uses RANSAC (Frey et al., 2026; Fischler & Bolles, 1981): a minimal subset of slice points is drawn, an ellipse is fitted to it by the direct least squares method of Fitzgibbon et al. (1999), the remaining points are scored against it, and the model with the largest consensus set is retained.

#### Inliers and the distance threshold

A point is an inlier when its distance to the candidate ellipse is below `distance_threshold`. The distance used is a genuine metre distance to the ellipse boundary rather than the algebraic residual of the conic, so the threshold means the same thing at every point on the curve. An algebraic residual varies with position around an eccentric ellipse, which would make one threshold behave as several.

The threshold is therefore the half-width of an inlier band around the candidate boundary, and its size is a statement about the stem surface rather than about the fit. It only helps while it stays small relative to the stem radius: for stems with semi-minor axes of roughly 5 to 15 cm, a band of ±15 mm is already 10 to 30% of the radius, and beyond that the fit stops being constrained by the stem surface and begins accepting whatever the slice contains.

#### Iteration budget and early exit

The number of hypotheses drawn is adaptive rather than fixed. When `iter_per_point` is set, the budget is

$$
\operatorname{clamp}\big(\texttt{iter\_per\_point} \times n_{\text{points}},\; \texttt{iter\_min},\; \texttt{iter\_max}\big),
$$

so a sparse slice is not given the same effort as a dense one, and no slice can consume unbounded time. A fixed `n_iterations` is available as an alternative.

Iteration also stops early, by default, once a hypothesis covers a sufficient share of the cross-section: the boundary is divided into 36 arc segments of 10° each, and a candidate whose inliers occupy at least `early_exit_segs` of them is accepted without exhausting the budget. The criterion is angular coverage rather than inlier count, because a cross-section is only well determined when its inliers are distributed around the boundary; a hypothesis supported by many points along one arc constrains the opposite side of the ellipse hardly at all. Early exit can be disabled, in which case every slice runs its full budget.

#### Slice thickness and the adaptive density ramp

Slice thickness and the distance threshold are the two settings that decide whether a slice contains a resolvable cross-section at all, and the right values depend on how densely the stems were sampled. A TLS epoch and a ULS epoch of the same plot can differ by an order of magnitude in points per slice.

Stage 1 therefore runs a **density probe** before any production fitting. It fits a small set of stems in both epochs and reduces the result to a single scarcity value $t$, computed from points per *kept* slice rather than per slice, since an epoch can average plenty of points across all slices and still retain almost nothing after cleaning. The reference points are $t = 0$ at or above 1000 points per kept slice and $t = 1$ at or below 150, interpolated linearly between.

Both epochs then take the *scarcer* epoch's $t$. They must be fitted identically or their ellipses are not comparable across time, so the thinner epoch decides for both. Every affected parameter interpolates from its base value at $t = 0$ toward its cap at $t = 1$:

| Parameter | Base ($t=0$) | Cap ($t=1$) |
|---|---|---|
| Slice thickness | 5 cm | 20 cm |
| Distance threshold | 5 mm | 15 mm |
| Consensus fraction, centre | base | 1.00 |
| Consensus fraction, semi-minor | base | 0.60 |
| Consensus fraction, semi-major | base | 0.60 |

The base values are the TLS settings. The distance threshold cap is deliberately held at three times the base rather than ramped further, for the reason given above: widening the inlier band past a certain fraction of the stem radius cannot recover a slice that does not contain a resolvable cross-section, and an over-wide band shows up downstream as inflated semi-axes.

The probe is optional and can be switched off, in which case the base values are used unchanged. When it runs, the recommended thickness and threshold are carried into Stage 3 so that the axis-perpendicular refit uses the same geometry as the Stage 1 fit.

### Stage 1: rough vertical alignment, ellipse fitting, temporal pairing

#### Rough vertical alignment

A scan may sit in its file in an arbitrary orientation. Stage 1 first tests the per-stem bounding boxes for a cloud lying on its side and, if it finds one, tips the first cloud upright before anything is fitted. It then fits a coarse ellipse set on every epoch-1 stem, diagnoses upright versus inverted from the trend of cross-sectional area against height, builds a forest-average stem axis, and rotates the first cloud so that axis points along global $+Z$.

Only the first cloud is touched, only in memory, and the composed rotation is recorded. Because the file on disk is never rewritten, this rough rotation is carried through to Stage 6 and composed into the final transform.

#### Temporal pairing by triangle congruence

Stems must be matched between epochs before any of their geometry can be compared. The default is automatic matching by **triangle congruence**, described below. A manually supplied mapping table of corresponding `TreeID`s is available as an alternative, and restricts the run to the pairs it lists.

Automatic matching works by triangle congruence, which is the established basis for target-less registration of forest scans. This registration methods does not require a stem in one epoch to be near its partner in the other for temporal pairing. Instead, only the relative arrangement of stems must be comparable between epochs. Triangles formed by triplets of stem positions are invariant under the rigid motion separating the two epochs, so a triangle in one epoch whose side lengths match a triangle in the other, within tolerance, is evidence of a correspondence between their vertices. Triangle congruence is the established basis for marker-free registration of forest terrestrial scans, and the method used here follows a specific line of that literature.

Kelbe et al. (2016) recovered the transformation between two scans from stem maps by matching triplets of stem positions, eliminating dissimilar triplets using DBH and geometric descriptors derived by principal component analysis, and attaching a confidence metric to each correspondence. Their approach becomes expensive as the number of stems grows, and Tremblay & Béland (2018) produced a parallelised variant that compares triangle side lengths and DBH instead of the PCA descriptors.

**Dai et al. (2020)** removed the dependence on DBH by sorting each triangle's vertices by the length of the edge opposite each one, so that matching needs only the 3D positions. They discard triangles that are isosceles, equilateral, highly collinear, or too small to be informative, then discard triangle pairs whose corresponding edges differ by more than a tolerance. For each surviving pair they solve a transform by singular value decomposition, apply it, count how many key points fall within a distance tolerance of a partner, and keep the transform with the most correspondences, breaking ties on RMSE.

The scoring used here follows **Wang et al. (2023)** instead, whose *inlier-grouping* framework establishes all correspondences once rather than iterating trial-and-error over candidate transforms. Their method builds triangles only between each stem and its $K$ nearest neighbours, matches triangles locally on edge length, and then *grows a consensus set*: two locally matched triangle pairs join the same set when they are mutually consistent, and the largest consensus set yields the correspondences.

The implementation here is of that family. All congruent candidates are retained rather than only the best per triplet, and each grows outward into an *island* of mutually consistent correspondences. The consistency test differs from GlobalMatch's: two correspondences are neighbours when they share two identical stem pairs — an edge — and a candidate that would claim a stem pair already assigned differently by the island is refused outright, so contradictory mappings cannot merge. The largest island wins; mean side residual breaks ties only. Each surviving pair records how many triangles support it.

A reporting residual is computed afterwards by fitting the winning island's correspondences with the closed-form rigid superposition of Kabsch (1976), using the singular value decomposition and determinant correction of Umeyama (1991) — the same least-squares rigid fit both Dai et al. (2020) and Wang et al. (2023) use to recover their transforms. Here it is diagnostic only, and does not enter the pairing decision.

The implementation here differs from a nearest-match scheme in how a candidate correspondence is *scored*. Rather than accepting the single best-matching triangle per triplet, all congruent candidates are retained, and each candidate is scored by how much other evidence agrees with it: two correspondences are neighbours when they share two identical stem pairs (an edge), and a correspondence grows outward into an *island* of mutually consistent correspondences over the same stems in the same arrangement. A candidate that brings a stem pair contradicting one already claimed by the island is refused, so incompatible mappings cannot merge. The largest island wins; mean side residual breaks ties only. Each surviving pair carries the number of triangles supporting it.

A reporting residual is computed afterwards by fitting the winning island's stem correspondences with the closed-form rigid superposition of Kabsch (1976), using the singular value decomposition and determinant correction described by Umeyama (1991). This residual is diagnostic; it does not enter the pairing decision.

### Stage 2: three-dimensional rotation from stem axes

#### The straightness gate

Pairs are collected into two groups using a gold length, 8 m by default,

**Group A**
Within a pair, both stems reach the gold length. Each will be measured over its bottom 8 m only.

**Group B**
Either stem falls short. Each will be measured over its full kept range, capped at the gold length, against a threshold scaled to that range. Thresholds scale as $t(L) = t_{\text{gold}}(L/\text{gold})^{p}$, with $p = 2$ for lateral deviation, which grows as $L^2$ at constant curvature, and $p = 1$ for direction change, which grows as $L$. Both hold curvature constant across window lengths, so a short window cannot pass simply by having had less stem over which to bend.

Stems are given straightness scores on three metrics; gate thresholds are set for each metric and adjusted by that stem's own measured noise floor. The three metrics are: lateral deviation from a fitted line, sag, and total direction change along the stem. A stem is called straight when at least two of the three metrics pass their threshold. So no single metric decides alone.

Every pair is given a straightness score based on the lower of the three metrics. Passing pairs are those where both stems are straight.

When too few stems reach the gold length, a **silver-length** pass calculates the same metrics over a shorter window, 4 m by default. When silver-lengths are implemented, gates will not be applied and all stems will be used downstream.

#### When gold-lengths are implemented, which pairs are used downstream?

The number of passing pairs determines which set is handed to the Stage 2 rotation and the Stage 3 vertical offset solve. Three regimes arise, and the workflow reports which one was taken:

| Regime | Condition | What each solve receives |
|---|---|---|
| `gate` | Enough pairs pass the vote | Both the rotation and the Stage 3 vertical offset use the gate survivors. Nothing is added. |
| `topup` | Some pass, but fewer than the target | The rotation still uses only the gate survivors. The vertical offset is topped up to the target with the best-ranked pairs that failed the vote. |
| `irls_all` | Too few pass to gate meaningfully | The gate filters nothing and every usable pair enters the rotation solve. The vertical offset uses this group as well. |

#### Stem lengths, and why they differ

Three different lengths of stem are in use at once, and they are independent by design:

- **The kept range** is whatever survived ellipse cleaning for that stem, and differs from stem to stem. A stem must span at least 4 m to be usable at all.
- **The straightness window** is the gold length, or the stem's capped kept range in group B, or the silver length under the fallback. So the two stems of one pair can be judged over different lengths from each other when one falls short.
- **The axis window** is the bottom 4 m of each stem's kept centres.

#### Axis fitting

A stem axis is fitted to each stem's ellipse centres over the axis window by Theil–Sen regression of $x$ on $z$ and $y$ on $z$ (Theil, 1950; Sen, 1968). A stem can therefore be judged straight over 8 m while its axis is fitted over the bottom 4 m.

The Theil–Sen regression is a median-of-slopes estimator that resists the leverage a single bad ellipse centre would exert on ordinary least squares. Ordinary least squares is available as an alternative but is not the default, since one bad centre with a long lever arm rotates the whole axis.

#### The rotation solve

The rotation carrying the epoch-1 axis vectors onto their epoch-2 partners is solved as a weighted orthogonal Procrustes problem by the Kabsch algorithm (Kabsch, 1976) with the determinant correction of Umeyama (1991), which guarantees a proper rotation rather than a reflection. The solve is iterated with reweighting (Holland & Welsch, 1977) under Tukey's biweight (Beaton & Tukey, 1974), so that stems whose axes disagree with the emerging consensus lose influence progressively rather than being cut by a hard threshold.

The rotation is written to the checkpoint and is not applied here.

### Stage 3: vertical offset from stem-property curves

Stage 3 applies the Stage 2 rotation to the first cloud about `pivot1`, rotates the retained stem axis vectors through the same transform, and refits ellipses on planes *perpendicular to each pair's average stem axis* rather than on horizontal planes. Fitting perpendicular to the axis means a leaning stem yields a true cross-section rather than an oblique one.

The vertical offset is then recovered by matching how stem properties vary along each stem. For every pair, a set of *signals* is computed as a function of distance along the stem axis, defined in [The signals](#the-signals). A smoothed LOESS curve is fitted for each stem-signal (Cleveland, 1979; Cleveland & Devlin, 1988). Features (peaks and troughs) are detected on each smoothed curve, and every plausible epoch-1 to epoch-2 feature match becomes a *candidate* carrying the vertical shift it implies.

#### Fitting, cleaning, refitting, recleaning and trimming

A single RANSAC pass over every slice produces both good cross-sections and confident nonsense, and the two are not separable by looking at a slice alone: an ellipse fitted to a branch stub is a perfectly valid ellipse. What distinguishes them is context. A stem is a continuous object, so a slice that describes it should resemble the slices immediately above and below it. The sequence below uses that fact, first to remove slices, then to give the removed ones a second chance under better conditions, then to judge the result again.

**Pass 1: fit.** Every slice along every stem is fitted independently by RANSAC (see [Ellipse fitting](#ellipse-fitting)). Slices with too few points, or where no hypothesis reaches consensus, produce no ellipse and are marked unfitted.

**Pass 1: clean.** Cleaning runs in three stages, in order.

**Stage A, neighbour consensus**
Each fitted slice is compared with the nearest fitted slice above and below it along the axis. Unfitted slices are skipped over, so a gap in the stack does not break the chain. Two slices agree only when *all three* of the following hold: the perpendicular distance between their centres is within a fraction of their mean semi-minor axis; their semi-minor axes differ by within a fraction of the smaller; and their semi-major axes likewise. A slice with no agreeing neighbour is removed. The tolerances are fractions rather than absolute distances so that the same setting applies to a thick stem and a thin one. Axes are compared largest-to-largest rather than by label, since which axis is named major can swap between adjacent slices on a near-circular section.

**Stage B, MAD on the semi-minor axis**
Across the surviving slices of one stem, semi-minor axes are screened by the median absolute deviation (see [Outlier rejection](#outlier-rejection)).

**Stage C, MAD on the semi-major axis**
The same, on the semi-major axis.

Consensus runs first because it is local: it removes slices that disagree with their immediate neighbours, which is where most bad fits are. The MAD stages are global to the stem and catch what survives that, such as a run of several adjacent slices that agree with each other and with nothing else.

**Pass 2: masked refit.** A slice removed by consensus has not necessarily failed. It may have found the wrong points: a branch, a neighbouring stem, or foliage inside the slice window can dominate a fit that had no reason to prefer the stem. Those slices are therefore refitted with the surrounding distraction masked away.

Each target slice takes an *anchor*, the nearest available slice below it, looking above only when none exists below. Available means kept by pass-1 cleaning, or successfully refit earlier in the same loop, so targets are processed from the bottom upward and a refit cascades through a run of consecutive failures. The anchor's centre is projected up the local $Z$ axis to the target slice, and a circular mask of the anchor's semi-major axis scaled by a factor, 1.3 by default, is applied. Points outside the mask are invisible to the refit.

Because Stage 3 fits in the axis frame, the local $Z$ axis *is* the stem axis, so projecting the anchor's centre up to the target slice introduces no lean error at all. Stage 1 has to absorb that error in its mask factor because it has no axis yet. This is one of the things fitting in the axis frame buys.

The refit inherits pass 1's distance threshold and iteration rule by default, but the iteration budget is computed on the *masked* point count, which is smaller, so a refit slice receives a smaller budget than its first pass did. By default only consensus removals are retried; extending this to MAD removals is available.

**Pass 2: reclean.** The refit slices join the kept set and the whole stack is put through consensus and both MAD stages again, with the pass-1 tolerances inherited unchanged by default.

Running the same tolerances is the conservative choice here, not a lenient one. A refit slice was fitted against a masked point set anchored on its own neighbour, so it is biased *toward* agreeing with that neighbour. Loosening the second pass would admit slices whose only evidence is an agreement the masking helped produce.

**Trimming the ends.** Cleaning judges each slice against its neighbours, which is the right question in the body of a stem and the wrong one at its ends. Near the base and up in the crown the survivors thin into a scatter: long runs of nothing, broken by the occasional slice that found enough points to agree with something. Those stragglers pass every filter so far, since no slice is removed merely for being isolated, but a handful of ellipses separated from the stem by a metre of gap describe nothing continuous, and a curve fitted through them fits noise exactly where it has least support.

Trimming finds long gaps, defined as runs of consecutive slices carrying no kept ellipse, and cuts off any side of a gap holding only a few kept ellipses. Why each slice is empty does not matter, whether never fitted, dropped by consensus, dropped by either MAD stage, or refit and dropped again. What is being detected is a break in continuity, and a break is a break whatever produced it.

**What survives.** The signals of [The signals](#the-signals) are computed on the trimmed, twice-cleaned stack. Each stage can be switched off independently, in which case the preceding output passes through unchanged, and the counts removed at each stage are reported so that a stem losing most of its slices is visible rather than silent.

#### The signals

Every signal is a per-slice scalar, computed in each stem's own local frame where the stem axis is the $z$ axis, and indexed by distance along that axis.

**Five ellipse signals.** These describe the fitted ellipse itself, and follow from the semi-major axis $a$, the semi-minor axis $b$, and the fitted centre.

| Signal | Definition |
|---|---|
| `area` | $\pi a b$, the area of the fitted ellipse |
| `sa` | $a$, the semi-major axis |
| `sb` | $b$, the semi-minor axis |
| `ecc` | $\sqrt{1 - (b/a)^2}$, eccentricity, zero for a circular section and approaching one as the section flattens |
| `center_offset` | The perpendicular distance of the fitted ellipse centre from the stem axis, $\sqrt{c_x^2 + c_y^2}$ in the local frame. It measures how far the section sits off the axis rather than the shape of the section, so it responds to a stem that wanders where the shape signals do not |

**Eight protrusion signals.** The five signals above all describe the fitted ellipse, and none of them says anything about the points the ellipse did *not* capture. That is where a branch stub, a bark plate or a buttress lives. The protrusion signals fill that gap: the axis-orthogonal plane is cut into eight sectors, and each records the fraction of that sector's points lying outside the fitted ellipse.

There are eight rather than one because the position of a protrusion around the stem is the informative part. A single fraction-outside would average a one-sided bulge away to nothing.

**The sector frame**
The sectors lie in the plane perpendicular to the pair axis, so they are *not* horizontal and their labels are *not* compass directions; N, NE, E and so on are used only because eight labelled directions read more easily than indices. North is defined as the in-plane direction with the greatest global $Z$ component, that is, the projection of the global $Z$ axis onto the orthogonal plane. This depends only on the pair axis, so both stems of a pair receive the same sector boundaries and their signals can be compared slice for slice. Anchoring on anything derived from the points themselves, such as the widest radius or the first principal direction, would let the two epochs rotate independently and destroy that comparability. The remaining seven directions sit at 45° intervals turning right-handed about the axis, and the sectors are the wedges between adjacent directions, so no point can fall in two.

**The mask**
Counting every point in the slice window would count leaves: a slice cut through a crown catches foliage metres from the stem, all of it outside the fitted ellipse, so a leafy stem would read as protruding in every direction at once and the signal would measure canopy density rather than stem shape. Each slice is therefore masked by its own fitted ellipse scaled by a factor, 1.3 by default. The mask is elliptical rather than circular because the ellipse is what the points are being scored against; a circular mask on a flattened section would admit more of the plane along the minor axis than the major one and bias the sectors accordingly. The same mask is applied to every slice, refit or not, so that two slices' protrusion values are always computed under comparable masks.

**The apex**
The wedges radiate from the fitted ellipse centre by default, so a bulge lands in the wedge it bulges into, measured from the shape it is bulging out of. Radiating from the stem axis instead is available, and is the better choice where the fits themselves are unstable, since the axis does not move when a fit wobbles.

**Reading the values**
A point counts as outside when its radial excess exceeds a tolerance, zero by default. Zero is the literal definition, and it has a consequence worth knowing: the fit places the ellipse *through* the point cloud, so roughly half the inlier shell sits marginally outside it and every sector starts from a baseline near 0.5. These values should be read as departures from about 0.5 rather than from 0. The baseline is shared by both epochs and cancels in the comparison. Setting the tolerance to the fit's own distance threshold instead measures the fraction outside the fitted *shell*, which starts near 0 and responds only to genuine protrusions.

**Minimum points**
A sector with fewer than three points returns `NA` rather than a fraction, since a quotient over one or two points swings between 0 and 1 on noise and would look like a measurement.

Every candidate is scored on twelve factors, each clamped independently to $[0,1]$. They are defined in [The twelve candidate factors](#the-twelve-candidate-factors).

Candidates are aggregated to produce a single shift per pair-signal. The default aggregation mode takes the implied shift of the single highest-scoring candidate in the pool, scored according to a predicted-error derived from the 12 scoring factors.

The eight protrusion signals are aggregated into a single protrusion signal by applying a MAD filter (see [Outlier rejection](#outlier-rejection)) to the 8 sectors then taking the median. This is done because each sector alone constitutes 8 unreliable signals.

Each pair-signal passes through two MAD stages: across signals within a pair, then across pairs within a signal, before being averaged into `t_trans_z`. The shift is recorded and applied in Stage 4.

#### The twelve candidate factors

**Three underlying measurements.** Six of the twelve factors are built from three underlying measurements, so here they are defined.

**Majesty**
Walk from the feature outward until a *larger feature of the same type* is met on both sides; the max elevation change along the way is recorded, then the lesser of these two height differences is the majesty.

**Major prominence**
The *lesser* of the two differences between a feature and its nearest opposite-type kept feature on each side, falling back to the curve's opposite global extreme where neither side has one.

**Minor prominence**
The drop on the *other* side, the one major prominence did not take. Where that side carries no opposite-type feature, the drop to the curve's own extreme over that side.

All are computed *after* weak features have been screened out, and re-measured against the survivors alone, so that a feature whose larger neighbour was filtered away is re-judged against the structure that actually remains.

**Six relative factors.** These come in three `_tot` / `_sim` pairs, one pair per underlying measurement. Each quantity is first normalised by that curve's own maximum *of the same kind and the same type*. The normalisation is not cosmetic: the signals live on incomparable scales (square metres of cross-section, metres of centre offset, a dimensionless protrusion fraction) and a curve with large features must not be judged against the same absolute bar as a flat one.

| Factor | What it asks |
|---|---|
| `rel_maj_tot` | Mean of the two features' relative majesty. *Are these two strong features?* |
| `rel_maj_sim` | Gaussian similarity of those two relative majesties (≤0.03 difference → 1, ≥0.50 → ≈0). *Are they strong to the same degree?* |
| `rel_majorprom_tot` | The first question, on major prominence |
| `rel_majorprom_sim` | The second question, on major prominence |
| `rel_minorprom_tot` | The first question, on minor prominence |
| `rel_minorprom_sim` | The second question, on minor prominence |

**Two shape factors.** Both compare the curves themselves over a window around the two features.

**`pearson`**
Correlation of the two curves over that window, mapped $\exp(-\lambda(1-r))$ so that $r = 1 \rightarrow 1$, $r = 0 \rightarrow 0.1$, and negative correlations decay to a floor.

**`slope_agree`**
Both curves are interpolated onto a common grid, subtracted, and the residual differentiated; the factor is the standard deviation of that derivative, normalised by the curves' own mean slope-SD to stay scale-free, and mapped as $(1-r)$. A *flat* residual means identical shapes at *any* vertical offset, which is exactly what a pure vertical shift produces.

**Two pool factors.** These are the only two that depend on the *other* candidates rather than on the candidate alone, and that distinction matters in [Learned error prediction](#learned-error-prediction).

**`median_prox`**
Gaussian proximity of this candidate's implied shift to the median of its pool, with the width set by the pool's own interquartile range. A pool of fewer than two candidates receives a neutral value.

**`consensus`**
How strongly the other candidates corroborate this one's implied shift, with near-agreement counting fully and disagreement beyond roughly 20 cm counting for nothing, summed and then min-max normalised within the pool.

**Two spacing factors.** Both test global vertical distance to the *opposite* feature type as opposed to signal geometric relationship with other features.

**`left_opp`**
Similarity of the distance to the nearest opposite-type feature to the *left* of each feature: for a peak candidate the nearest trough, for a trough candidate the nearest peak. Where one side has such a feature and the other does not, the factor takes a low value; where neither does, a neutral one.

**`right_opp`**
The same, to the right.

Two features that are genuinely the same feature should sit at similar distances from the structure surrounding them, which is what these two measure.

### Learned error prediction

A random forest (Breiman, 2001) model predicts, from a candidate's twelve scores, how far that candidate's implied shift falls from the truth. Candidates are then weighted by their predicted error. The forest model was fitted with `ranger` (Wright & Ziegler, 2017).

#### Training

Training used the candidate-level tables written by previous runs on samples whose true shift is known (see [Truth values and testing](#truth-values-and-testing)). For a candidate with implied shift $s_i$ and known true shift $s^{\ast}$, the response is

$$
y_i = \log\!\big(\max(|s^{\ast} - s_i|,\; \varepsilon)\big),
$$

on the log scale because these errors span millimetres to tens of metres, and a forest fitted on the raw scale would spend its capacity separating very wrong candidates from extremely wrong ones while the weighting only distinguishes the first few centimetres. The floor $\varepsilon$ prevents a candidate that lands exactly right from producing $-\infty$. The predictors are the twelve raw factor scores of [The twelve candidate factors](#the-twelve-candidate-factors).

One forest is fitted per signal, except that the eight protrusion sectors share a single model, and a fallback forest is fitted on all signals pooled for any signal a model file does not cover. A signal with too few training rows or too few distinct stem pairs is not given its own model and is routed to the fallback.

#### Application

At run time the model file supplies, for each signal, the forest to use, the predictor names in the order the forest expects, and the scale of its response. Predictions are exponentiated back to metres. For a pool of candidates with predicted errors $\hat{e}_i$, each is weighted

$$
w_i = \exp\!\left(-k \left(\frac{\hat{e}_i - \min_j \hat{e}_j}{\delta}\right)^{\!2}\right),
$$

so that a candidate whose predicted error exceeds the pool's best by $\delta$ carries weight $e^{-k}$. The reduction is the weighted mean of the pool's implied shifts.

The anchor is the pool's own minimum rather than zero. A pool in which every candidate is predicted to be poor still yields an answer weighted toward its own best, instead of collapsing to zero weights and dropping out; and pools of differing overall difficulty remain comparable. The kernel is applied within each pool separately, never across pools, since that is the set the minimum is taken over.

**The model reads no truth at run time.** It is trained offline on samples that have one and thereafter predicts from the twelve scores alone, so it applies unchanged to a pair with no known truth (see [What happens when there is no truth](#what-happens-when-there-is-no-truth)).

### Stage 4: yaw from the inter-stem line network

Stage 4 applies the vertical offset, then averages each stem's ellipse centres to one representative horizontal point per stem per epoch. Every pair of stems defines an edge, and the angular difference between the epoch-1 and epoch-2 versions of that edge is measured. A pure yaw rotates every edge by the same angle, so each edge is an independent estimate of the same quantity. Stems whose edges disagree with the consensus are removed iteratively by MAD, and the rotation is the mean of the surviving edge-wise differences.

### Stage 5: residual horizontal translation

The Stage 4 yaw is applied to the averaged epoch-1 centres, the horizontal displacement to each epoch-2 partner is measured, stems disagreeing with the consensus are removed by MAD on each component, and the surviving displacements are averaged into `t_trans_xy`.

### Stage 6: composition and export

Stage 6 rebuilds every solved component and applies the composed transform to the cloud in one pass:

$$
p' = T_z\big(R_{\text{align+xyz}}\,p + t_{\text{align+xyz}}\big) + t,
$$

where $R_{\text{align+xyz}}$ is the Stage 1A rough vertical alignment composed with the Stage 2 rotation, both about `pivot1`; $T_z$ is the Stage 4 yaw about the same pivot; and $t$ carries the Stage 5 horizontal translation and the Stage 4 vertical offset. The aligned cloud is written to disk, and the composed matrix is recorded alongside it.

---

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

---

## Running the workflow

### Batch

Set the scan pairs and the output locations at the top of `RunR.R` and source it:

```r
source("RunR.R")
```

The driver validates every input, then runs each pair through all six stages, writing checkpoints, an aligned cloud, a per-pair log and a report row.

### Interactively, one stage at a time

Open the stage you want, point its load path at the previous stage's checkpoint, and run the document. Because checkpoints are self-describing, this works in a fresh R session with nothing else loaded. This is the path to use when diagnosing a stage or tuning its parameters: the documents render every figure and diagnostic table that the batch skips.

---

## Requirements

R ≥ 4.3, with `lidR`, `data.table`, `dplyr`, `tidyr`, `tibble`, `ggplot2`, `plotly`, `dbscan`, `cli`, `htmltools`, and `ranger` for the optional random forest weighting. Quarto is needed only to render the documents; the batch driver parses them directly and does not require it.

---

## References

Beaton, A. E., & Tukey, J. W. (1974). The fitting of power series, meaning polynomials, illustrated on band-spectroscopic data. *Technometrics, 16*(2), 147–185. https://doi.org/10.1080/00401706.1974.10489171

Breiman, L. (2001). Random forests. *Machine Learning, 45*(1), 5–32. https://doi.org/10.1023/A:1010933404324

Cleveland, W. S. (1979). Robust locally weighted regression and smoothing scatterplots. *Journal of the American Statistical Association, 74*(368), 829–836. https://doi.org/10.1080/01621459.1979.10481038

Cleveland, W. S., & Devlin, S. J. (1988). Locally weighted regression: An approach to regression analysis by local fitting. *Journal of the American Statistical Association, 83*(403), 596–610. https://doi.org/10.1080/01621459.1988.10478639

Dai, W., Yang, B., Liang, X., Dong, Z., Huang, R., Wang, Y., Pyörälä, J., & Kukko, A. (2020).** Fast registration of forest terrestrial laser scans using key points detected from crowns and stems. *International Journal of Digital Earth, 13*(12), 1585–1603. https://doi.org/10.1080/17538947.2020.1764118

Fischler, M. A., & Bolles, R. C. (1981). Random sample consensus: A paradigm for model fitting with applications to image analysis and automated cartography. *Communications of the ACM, 24*(6), 381–395. https://doi.org/10.1145/358669.358692

Fitzgibbon, A., Pilu, M., & Fisher, R. B. (1999). Direct least square fitting of ellipses. *IEEE Transactions on Pattern Analysis and Machine Intelligence, 21*(5), 476–480. https://doi.org/10.1109/34.765658

Frey, J., Schindler, Z., & Kröner, K. (2026). Cspstandsegmentation: Stand segmenta￾tion for individual tree detection from airborne laser scanning data [R package
version 0.2.0]. https://CRAN.R-project.org/package=CspStandSegmentation

Hahsler, M., Piekenbrock, M., & Doran, D. (2019). dbscan: Fast density-based clustering with R. *Journal of Statistical Software, 91*(1). https://doi.org/10.18637/jss.v091.i01

Holland, P. W., & Welsch, R. E. (1977). Robust regression using iteratively reweighted least-squares. *Communications in Statistics - Theory and Methods, 6*(9), 813–827. https://doi.org/10.1080/03610927708827533

Iglewicz, B., & Hoaglin, D. C. (1993). *How to detect and handle outliers* (The ASQC Basic References in Quality Control: Statistical Techniques, Vol. 16). ASQC Quality Press.

Kabsch, W. (1976). A solution for the best rotation to relate two sets of vectors. *Acta Crystallographica Section A, 32*(5), 922–923. https://doi.org/10.1107/S0567739476001873

Kelbe, D., van Aardt, J., Romanczyk, P., van Leeuwen, M., & Cawse-Nicholson, K. (2016). Marker-free registration of forest terrestrial laser scanner data pairs with embedded confidence metrics. *IEEE Transactions on Geoscience and Remote Sensing, 54*(7), 4314–4330. https://doi.org/10.1109/TGRS.2016.2539219

Roussel, J.-R., Auty, D., Coops, N. C., Tompalski, P., Goodbody, T. R. H., Sánchez Meador, A., Bourdon, J.-F., de Boissieu, F., & Achim, A. (2020). lidR: An R package for analysis of Airborne Laser Scanning (ALS) data. *Remote Sensing of Environment, 251*, 112061. https://doi.org/10.1016/j.rse.2020.112061

Rousseeuw, P. J., & Croux, C. (1993). Alternatives to the median absolute deviation. *Journal of the American Statistical Association, 88*(424), 1273–1283. https://doi.org/10.1080/01621459.1993.10476408

Sen, P. K. (1968). Estimates of the regression coefficient based on Kendall's tau. *Journal of the American Statistical Association, 63*(324), 1379–1389. https://doi.org/10.1080/01621459.1968.10480934

Sievert, C. (2020). *Interactive web-based data visualization with R, plotly, and shiny*. Chapman and Hall/CRC. https://plotly-r.com

Theil, H. (1950). A rank-invariant method of linear and polynomial regression analysis. *Indagationes Mathematicae, 12*(85), 173.

Tremblay, J.-F., & Béland, M. (2018). Towards operational marker-free registration of terrestrial lidar data in forests. *ISPRS Journal of Photogrammetry and Remote Sensing, 146*, 430–435. https://doi.org/10.1016/j.isprsjprs.2018.10.011

Umeyama, S. (1991). Least-squares estimation of transformation parameters between two point patterns. *IEEE Transactions on Pattern Analysis and Machine Intelligence, 13*(4), 376–380. https://doi.org/10.1109/34.88573

Wang, X., Yang, Z., Cheng, X., Stoter, J., Xu, W., Wu, Z., & Nan, L. (2023).** GlobalMatch: Registration of forest terrestrial point clouds by global matching of relative stem positions. *ISPRS Journal of Photogrammetry and Remote Sensing, 197*, 71–86. https://doi.org/10.1016/j.isprsjprs.2023.01.013

Wickham, H. (2016). *ggplot2: Elegant graphics for data analysis*. Springer-Verlag New York. https://ggplot2.tidyverse.org

Wickham, H., François, R., Henry, L., Müller, K., & Vaughan, D. (2026). *dplyr: A grammar of data manipulation* (R package version 1.2.1). https://CRAN.R-project.org/package=dplyr

Wright, M. N., & Ziegler, A. (2017). ranger: A fast implementation of random forests for high dimensional data in C++ and R. *Journal of Statistical Software, 77*(1), 1–17. https://doi.org/10.18637/jss.v077.i01
