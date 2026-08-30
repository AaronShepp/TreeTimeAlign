# =============================================================================
# RunAlignmentBatch - run the full six-stage alignment over a list of LAS pairs
# =============================================================================
#
# WHAT THIS IS
#
# A driver, not a copy. It contains none of the alignment logic: it reads the
# six Stage .qmd files, pulls their R chunks out, and evaluates those chunks in
# order. Every calculation still lives in exactly one place - the .qmd - and an
# edit made there takes effect here on the next run with nothing to keep in
# sync. The alternative, pasting the stage code into a script, would have given
# two copies of a 21,000-line pipeline that drift apart within a week.
#
# HOW THE PARAMETERS GET IN
#
# The .qmd files hardcode their paths, because they are also meant to be opened
# and rendered by hand. So this driver does not try to pre-set variables and
# hope the chunk leaves them alone - the chunk would overwrite them. Instead it
# runs the chunk, lets it assign whatever it likes, and then OVERWRITES the
# specific variables afterwards, before the next chunk sees them. Overrides are
# therefore keyed by chunk label: "after this chunk runs, set these."
#
# That is why a change to a .qmd can break this file: if a parameter is renamed
# or moves to a different chunk, the override silently stops landing where it
# is needed. Every override target is checked for existence after it is applied
# and the run stops if one is missing, so a rename fails loudly on the first
# pair rather than quietly producing six wrong results.
#
# WHAT IT PRODUCES
#
#   per pair:  the aligned .las               (Stage 6)
#              the Z-translation optimisation .csv (Stage 3.16)
#              the temporal pairs .csv        (Stage 1)
#              four .rds checkpoints          (Stages 1, 2, 3, 4, 5 - see below)
#   once:      one report .csv covering every pair, with the nine error terms,
#              every rotation matrix, and the pair counts from each filter
#
# The per-pair result rows come from Stage 6 Step 6.6, which is the only place
# that knows how to compute the errors. This driver collects them; it does not
# recompute anything.
#
# THE PIVOT
#
# Nothing has to be told which pivot to use. Each generated sample carries its
# own, as the point with PivotPoint = 1, so Stage 1 is put in pivot1_mode =
# "las" and reads it out of the cloud it was handed. A batch of twenty samples
# needs twenty different pivots and no bookkeeping at all.
#
# =============================================================================


# =============================================================================
# 1. INPUT PARAMETERS - the only part meant to be edited
# =============================================================================

## -- 1.1 The .qmd files ------------------------------------------------------
## Order matters and is not inferred from the names: the driver runs them in
## the order given here.
qmd_paths <- c(
  stage1 = file.path("C:/Users/ascsh/Documents/1_Artemis/11_PaqLab/115_TreeTimeAlign/Stage1_TemporalPairing", "Stage1_2026-08-29_3.qmd"),
  stage2 = file.path("C:/Users/ascsh/Documents/1_Artemis/11_PaqLab/115_TreeTimeAlign/Stage2_XYZRot", "Stage2_2026-08-29_3.qmd"),
  stage3 = file.path("C:/Users/ascsh/Documents/1_Artemis/11_PaqLab/115_TreeTimeAlign/Stage3_ZShift", "Stage3_2026-08-29_7.qmd"),
  stage4 = file.path("C:/Users/ascsh/Documents/1_Artemis/11_PaqLab/115_TreeTimeAlign/Stage4_ZRot", "Stage4_2026-08-29_3.qmd"),
  stage5 = file.path("C:/Users/ascsh/Documents/1_Artemis/11_PaqLab/115_TreeTimeAlign/Stage5_XYShift", "Stage5_2026-08-29_3.qmd"),
  stage6 = file.path("C:/Users/ascsh/Documents/1_Artemis/11_PaqLab/115_TreeTimeAlign/Stage6_Export", "Stage6_2026-08-29_3.qmd")
)

## -- 1.2 The RANSAC ellipse library ------------------------------------------
## Sourced by Stages 1 and 3. Given here once and pushed into both.
ransac_library_path <-
  "C:/Users/ascsh/Documents/1_Artemis/11_PaqLab/115_TreeTimeAlign/RANSAC/RANSACellipse_2026-08-19_1_PointMatch.R"

## -- 1.3 The consolidated truth table ----------------------------------------
## One row per LAS pair, in the column layout the test-sample generator writes
## (the per-sample CSVs from TestSampleGen, rbind-ed together). Matched to each
## pair on `sample_name`.
truth_csv_path <-
  "C:/Users/ascsh/Documents/1_Artemis/11_PaqLab/115_TreeTimeAlign/ReferenceCSV/SNBG_F1_F1006_F1125_a_2020_2024_TLS_TLS.csv"

## -- 1.4 The pairs to run ----------------------------------------------------
## Just the names. Every path is resolved from the consolidated truth table in
## Section 2.5, so a sample is named ONCE here instead of three times - and the
## paths the driver uses are the ones the generator actually wrote, not a
## transcription of them.
##
## Each name must match a sample_name row in truth_csv_path.
sample_names <- c(
  "I_SNBG_F1_a_0p00_2020_2024_TLS_TLS",
  "I_SNBG_F1_a_29p71_2020_2024_TLS_TLS",
  "I_SNBG_F1_a_45p59_2020_2024_TLS_TLS",
  "I_SNBG_F1_a_90p00_2020_2024_TLS_TLS",
  "I_SNBG_F1_a_134p77_2020_2024_TLS_TLS",
  "I_SNBG_F1_a_180p00_2020_2024_TLS_TLS",
  "I_SNBG_F1006_a_0p00_2020_2024_TLS_TLS",
  "I_SNBG_F1006_a_29p71_2020_2024_TLS_TLS",
  "I_SNBG_F1006_a_45p59_2020_2024_TLS_TLS",
  "I_SNBG_F1006_a_90p00_2020_2024_TLS_TLS",
  "I_SNBG_F1006_a_134p77_2020_2024_TLS_TLS",
  "I_SNBG_F1006_a_180p00_2020_2024_TLS_TLS",
  "I_SNBG_F1125_a_0p00_2020_2024_TLS_TLS",
  "I_SNBG_F1125_a_29p71_2020_2024_TLS_TLS",
  "I_SNBG_F1125_a_45p59_2020_2024_TLS_TLS",
  "I_SNBG_F1125_a_90p00_2020_2024_TLS_TLS",
  "I_SNBG_F1125_a_134p77_2020_2024_TLS_TLS",
  "I_SNBG_F1125_a_180p00_2020_2024_TLS_TLS"
)

## las_1 - the generated sample. The truth table records the path the generator
## wrote it to, in `sample_las`, so this is a lookup and not a guess.
##
## las_2 - the SECOND-EPOCH TARGET, and a different file from the reference the
## sample was made from. `reference_las_1` in the truth table is the R_ file the
## damage was applied to; pairing against that would be scoring the pipeline on
## a cloud it has already seen. So las_2 is built from a template over the truth
## row's own identifier fields:
##
##   T_<Sample_location_ID>_<Cluster_ID>_<Sample_Type_ID>_<Epoch_2>_<Scanner_2>.las
##
## Anything in {braces} is replaced by that column of the sample's truth row.
## Change the template, not twenty-eight paths, when the convention changes.
las_2_dir      <- "C:/Users/ascsh/Documents/1_Artemis/11_PaqLab/115_TreeTimeAlign/Samples"
las_2_template <- "T_{Sample_location_ID}_{Cluster_ID}_{Sample_Type_ID}_{Epoch_2}_{Scanner_2}.las"

## Overrides, for a pair that does not follow the convention. Any sample not
## named here is resolved as above; anything named here wins. Leave empty to
## resolve everything automatically.
##   las_1_overrides <- c("I_SNBG_F1_a_0p00_2020_2024_TLS_TLS" = "D:/somewhere/odd.las")
las_1_overrides <- character(0)
las_2_overrides <- character(0)

## -- 1.5 Output directories --------------------------------------------------
out_root <- "C:/Users/ascsh/Documents/1_Artemis/11_PaqLab/115_TreeTimeAlign/RunR/exports"

out_dirs <- list(
  aligned      = file.path(out_root, "AlignedLAS"),        # per pair: the .las
  optimisation = file.path(out_root, "Optimization"),   # per pair: <sample>_fmatch_optim.csv
  pairs        = file.path(out_root, "TemporalPairs"),  # per pair: pairs csv
  checkpoints  = file.path(out_root, "SaveLoad"),       # per pair: the .rds set
  logs         = file.path(out_root, "Logs"),           # per pair: console log
  report       = out_root                               # once: the batch report
)

report_csv_name <- sprintf("Report_%s.csv", format(Sys.Date(), "%Y-%m-%d"))

## -- 1.6 Run behaviour -------------------------------------------------------

# Figures are the single largest cost in a batch and nothing downstream reads
# them. FALSE sets params$make_plots and skips every chunk whose label marks it
# as a figure.
make_plots <- FALSE

# Diagnostics - MAD multiplier sweeps, convergence plots, IQR summaries and the
# like - are read by a person deciding how to set a parameter. A batch is not
# that; it runs with the parameters already chosen, and no diagnostic product
# appears in any stage's save manifest. So they are off, and the cost is not
# only their runtime: a diagnostic that stops on a missing precondition fails
# the whole pair, which is how a completed 148-minute run was recorded as a
# failure because no RANSAC convergence history had been stored.
#
# TRUE runs them anyway, but even then they can no longer fail a pair - see
# .run_qmd, which evaluates every non-core chunk under its own error handler.
run_diagnostics <- FALSE

# One bad pair should not cost the other nineteen. On FALSE the driver stops at
# the first failure instead, which is what you want while setting a batch up.
continue_on_error <- TRUE

# Mirror each pair's console output to its own log file. Worth keeping on: the
# stage printouts are the only record of what the filters did, and in a batch
# they otherwise scroll past unread.
write_logs <- TRUE

# Keep the .rds checkpoints per pair. They are large. FALSE deletes them once
# the pair finishes, which loses the ability to resume a single stage.
keep_checkpoints <- TRUE

# On a failure, hand the pair's environment back instead of dropping it. The
# stages write their checkpoint at the END of a stage, so a failure part-way
# through one loses everything that stage computed - which for Stage 1 is an
# hour of RANSAC per cloud. With this on, the environment is left in
# .last_failed_env, and the fitted objects can be inspected or the run resumed
# from the chunk that died (see .run_qmd's start_after, and the resume recipe
# printed after a failure).
#
# It costs memory: one pair's full object set, clouds included, stays resident
# for the rest of the batch. On a long batch of large clouds, or if a pair is
# expected to fail for a reason you do not intend to recover from, set FALSE.
keep_failed_env <- FALSE

## -- 1.7 Resuming a failed run -----------------------------------------------
#
# Both NULL for a normal batch. Set them together to carry on from where a
# failure stopped, without repeating the stages and chunks that already ran:
#
#   resume_env  <- .last_failed_env
#   resume_from <- list(stage = "stage1", after = "diag-ransac-convergence")
#
# `after` is a chunk LABEL, and the resumed stage starts at the chunk following
# it. The failure message prints both values ready to paste, so this is normally
# a copy rather than something to work out.
#
# Applies to the FIRST pair only - the environment being resumed belongs to one
# sample. Later pairs in the same batch run normally in fresh environments, so
# put the failed sample first in las_pairs if there are several.
#
# A resume trusts the environment it is given. If it is not the one that failed,
# the run will stop at the first object it cannot find, which is the intended
# behaviour: half a result is worse than none.
resume_env  <- NULL  #"NULL": restart, ".last_failed_env": from checkpoint
resume_from <- NULL # list(stage = "stage3", after = "helpers")

# ── Stage 3's two checkpoints ────────────────────────────────────────────────
# Stage 3 saves twice: 3.1 straight after the axis-perpendicular RANSAC fit, and
# 3.2 at the end of the stage. Only one is worth keeping, and it is 3.1.
#
# 3.1 sits immediately after the only expensive thing Stage 3 does. Everything
# after it - peak finding, candidate scoring, the random-forest weighting, the
# sector pooling, both MAD passes - runs in seconds from that checkpoint. So 3.1
# is the file that lets the whole match-and-reduce end of the stage be re-tuned
# without refitting a single ellipse, which is exactly what tuning the forest
# weights requires.
#
# 3.2 is a strict superset of 3.1 plus a few seconds of arithmetic. Keeping it
# buys the ability to start at Stage 4 without re-running 3.12-3.17, which costs
# less than storing the file.
#
# This is the reverse of the earlier default, which kept 3.2 and skipped 3.1.
# That made sense when nothing downstream was being re-tuned. It does not now.
#
# Note both chunks still RUN either way: 3.1's manifest defines
# .stage3_1_required, which 3.2's save builds on, so the chunk cannot be skipped
# outright - only its write suppressed.
keep_stage3_1_checkpoint <- TRUE     # post-fit: re-tune matching without refitting
keep_stage3_2_checkpoint <- FALSE    # end of stage: cheap to rebuild from 3.1

# ── Stage 4 and Stage 5 checkpoints ──────────────────────────────────────────
# Neither stage fits an ellipse. Stage 4 applies t_trans_z to geometry Stage 3
# already produced; Stage 5 filters and reports on it. Their checkpoints are
# therefore near-copies of the Stage 3.2 one - cumulative manifests, so each is
# a strict superset of the last - carrying no new expensive computation. On a
# batch of twenty pairs that is forty large files written and never read.
#
# What is actually lost by switching them off: the ability to resume a run at
# Stage 5 or Stage 6 without re-running Stage 4. Since Stage 4 and Stage 5 take
# seconds once Stage 3.2 exists, resuming from the Stage 3.2 checkpoint costs
# less than storing the other two.
#
# TRUE restores them, per stage, if a batch is being used to produce
# checkpoints for later interactive work rather than to produce results.
keep_stage4_checkpoint <- FALSE
keep_stage5_checkpoint <- FALSE

# Stage 3's post-rotation and Stage 4's post-Z-shift .las snapshots. Diagnostic
# exports meant for opening in CloudCompare one at a time; both default to TRUE
# in the .qmds because that is the right default for someone running a stage by
# hand. In a batch they are a full point cloud written per pair per stage, for
# nobody, and Stage 4's has to reload both clouds from disk first. TRUE if you
# want the intermediate frames of every pair on disk to inspect later; the final
# aligned cloud is exported either way.
export_stage_clouds <- FALSE

# Score each pair against the generator's truth CSV (Stage 6 §6.6) and write its
# result row. This is normally the point of running a batch. FALSE only makes
# sense for pairs with no known truth - and note that with it off, Stage 6 never
# builds `alignment_result`, so every pair is reported as "no_result" and the
# batch report carries nothing but sample names and runtimes.
score_truth <- TRUE

# Where Stage 3's DIAGNOSTIC true Z-shift comes from - the yardstick §3.13.3
# scores the shift search against. Separate from score_truth above, which is
# Stage 6 scoring the final transform.
#
#   "updatedBox" re-derive it from the marker box AFTER this run's rotation:
#                the eight box-2 corners take the same two legs las_1 took, and
#                the shift is the mean of box1_z - rotated_box2_z. Asks what the
#                stage can actually answer - given the rotation I ended up with,
#                what Z-shift is correct - instead of charging it for rotation
#                error it did not cause.
#   "csv"        the generator's target_shift_z_m, from the per-pair truth row.
#   "manual"     the constant in the .qmd (s3_true_shift_manual).
#   "none"       no truth; §3.13.3 is skipped and its error columns come out NA.
#
# Set here rather than left to the .qmd because it is the one input to that
# chunk the driver could not previously see: an interactive edit left at "none"
# would otherwise follow every pair of the next batch in silence, which is the
# same failure that moved the truth CSV path onto the driver to begin with.
# Each source falls through to the next if it cannot produce a number, and none
# of them can stop a pair - a missing marker box costs the diagnostic only.
truth_shift_source <- "updatedBox"


# =============================================================================
# 2. MACHINERY - reading chunks out of a .qmd and running them
# =============================================================================

suppressPackageStartupMessages({
  library(lidR)
  library(data.table)
  library(dplyr)
})

#' Split a .qmd into its executable R chunks.
#'
#' Deliberately a parser rather than knitr::purl(): purl flattens everything
#' into one script and throws away the chunk labels, and the labels are the only
#' handle this driver has for placing an override or skipping a figure.
#'
#' Returns a list of chunks, each with: label, code, and the `eval` option if
#' the chunk declares one.
.read_qmd_chunks <- function(path) {

  if (!file.exists(path)) stop("qmd not found:\n  ", path)
  lines <- readLines(path, warn = FALSE)

  # Opening fences for R chunks only. ```{python}, ```{=html} and plain ```
  # blocks inside prose are all left alone.
  open_i  <- grep("^```\\{r[ ,}]", lines)
  close_i <- grep("^```\\s*$", lines)

  chunks <- vector("list", length(open_i))

  for (k in seq_along(open_i)) {
    start <- open_i[k]
    ends  <- close_i[close_i > start]
    if (!length(ends)) stop("Unclosed chunk at line ", start, " of ", basename(path))
    end <- ends[1]

    body <- if (end > start + 1L) lines[(start + 1L):(end - 1L)] else character(0)

    # #| options form an unbroken run at the TOP of the body and are not code.
    # Only that leading run is stripped: a later "#|" would be an ordinary
    # comment inside the code and removing it would change nothing, but
    # removing it by mistake from inside a string literal would.
    n_opt <- 0L
    while (n_opt < length(body) && grepl("^#\\|", body[n_opt + 1L])) n_opt <- n_opt + 1L
    opt_lines <- if (n_opt) body[seq_len(n_opt)] else character(0)
    code      <- if (n_opt) body[-seq_len(n_opt)] else body

    .opt <- function(nm) {
      hit <- grep(sprintf("^#\\|\\s*%s\\s*:", nm), opt_lines, value = TRUE)
      if (!length(hit)) return(NA_character_)
      trimws(sub(sprintf("^#\\|\\s*%s\\s*:", nm), "", hit[1]))
    }

    lab <- .opt("label")
    chunks[[k]] <- list(
      label = if (is.na(lab)) sprintf("unlabelled_%03d", k) else lab,
      eval  = .opt("eval"),
      code  = code,
      line  = start
    )
  }
  chunks
}

# ── The driver's own progress reporting ───────────────────────────────────────
#
# Self-contained on purpose. The stages draw cli progress bars inside their long
# loops, and those work when a .qmd is run on its own - but cli draws on stderr,
# and this driver holds a message sink (which, unlike the stdout sink, CANNOT be
# split), so nothing cli writes reaches the terminal during a batch. Rather than
# fight that, the driver reports progress at its own level, in plain cat() lines
# on stdout, which the split sink carries to both the console and the log.
#
# Two levels, and neither reads anything out of the stage code:
#
#   pair-level   how many of the batch's pairs are done, and roughly how long
#                the rest will take
#   chunk-level  which stage and which chunk is running right now, as a count
#                out of that file's total
#
# The chunk count comes from the parser, which already knows how many chunks a
# .qmd holds before it evaluates any of them, so this needs no cooperation from
# the stages at all. A long chunk still prints nothing while it runs - the line
# for the chunk appears when it STARTS, so a stalled run shows the name of the
# chunk it is stuck in, which is the thing worth knowing.

.fmt_dur <- function(sec) {
  if (!is.finite(sec)) return("--")
  if (sec < 90) sprintf("%.0fs", sec) else sprintf("%.1fm", sec / 60)
}

.bar <- function(frac, width = 20L) {
  n <- as.integer(round(max(0, min(1, frac)) * width))
  paste0("[", strrep("#", n), strrep("-", width - n), "]")
}

.say <- function(...) { cat(..., sep = ""); utils::flush.console() }

# ── Chunks that are always skipped, figures on or off ─────────────────────────
#
# The override mechanism sets a variable AFTER a chunk runs, which works for a
# parameter chunk because the parameter is consumed later. It does NOT work for
# a chunk that assigns a path and then immediately uses it - the chunk has
# already read the hardcoded value by the time the override lands.
#
# Stage 2's fig-3d-stem-with-vectors is the one chunk in the pipeline that does
# this: it re-declares las_1_path / las_2_path for the on-demand cloud loader
# and reads them in the same breath. It is a figure and nothing downstream
# consumes its output, so it is skipped unconditionally rather than left to
# read some other sample's clouds. Turning make_plots on does not bring it back.
always_skip <- c(
  stage2 = "fig-3d-stem-with-vectors"
)

#' Classify a chunk from its label: "figure", "diagnostic" or "core".
#'
#' Matching is on hyphen/underscore-delimited TOKENS, not substrings, so
#' "align-vis-slice-status" is caught by its middle token while a core chunk
#' that merely contains the letters is not. Labels are the only classification
#' the .qmd files offer, and they are a convention rather than a contract - so
#' anything unrecognised is treated as CORE. Skipping a core chunk produces a
#' failure a long way from its cause; running a diagnostic costs seconds.
#'
#' Not classified as diagnostic: summary-*. Several of those derive values that
#' later chunks read, and they are cheap. If a summary chunk is ever made to
#' stop on a missing input, it belongs in always_skip, not here.
.chunk_class <- function(label) {
  # A params-<X> chunk is whatever <X> is. It holds only <X>'s inputs, so if <X>
  # is a figure that this batch skips, its parameters must be skipped too -
  # otherwise the assignments run with the functional chunk absent. That is not
  # hypothetical: params-fig-3d-stem-with-vectors reassigns las_1_path and
  # las_2_path, and letting it run would point the rest of the stage at whatever
  # collection that figure happens to name.
  label <- sub("^params-", "", label)
  tok <- unlist(strsplit(tolower(label), "[-_.]+"))
  if (any(tok %in% c("fig", "figure", "figures", "plot", "plots",
                     "vis", "viz", "visual", "visualise", "visualize",
                     "visualisation", "visualization",
                     "inspect", "gallery", "preview")))            return("figure")
  if (any(tok %in% c("diag", "diagnostic", "diagnostics",
                     "debug", "sweep")))                           return("diagnostic")
  "core"
}

#' Run one .qmd's chunks in an environment, applying overrides as they come due.
#'
#' @param overrides named list: names are chunk labels, values are named lists
#'   of variable = value applied AFTER that chunk runs.
#' @param force_eval chunk labels to run even though they carry `eval: false`.
#' @param skip       chunk labels to skip entirely.
#' @param skip_figures     skip chunks classed "figure" (see .chunk_class).
#' @param skip_diagnostics skip chunks classed "diagnostic".
#' @param start_after resume point: skip every chunk up to AND INCLUDING this
#'   label, then run the rest. For restarting a stage that failed part-way in an
#'   environment that already holds everything the earlier chunks produced -
#'   which is what .last_failed_env is. NA/NULL runs the whole file.
#'
#'   Two things it does NOT do, both deliberate. It does not apply the overrides
#'   belonging to the chunks it skipped: those chunks already ran, and their
#'   overrides already landed, in the run being resumed. And it does not check
#'   that the environment actually contains what the skipped chunks would have
#'   made - resuming into a fresh environment will fail somewhere downstream
#'   with an object-not-found error, and that is the intended outcome rather
#'   than a silent half-run.
.run_qmd <- function(path, env,
                     overrides  = list(),
                     force_eval = character(0),
                     skip       = character(0),
                     skip_figures = TRUE,
                     skip_diagnostics = TRUE,
                     start_after = NULL,
                     verbose    = TRUE) {

  chunks <- .read_qmd_chunks(path)
  nm <- basename(path)

  # A start_after that matches nothing would silently run the entire file and
  # repeat the expensive work the resume was meant to avoid, so it is an error.
  resuming <- !is.null(start_after) && !is.na(start_after) && nzchar(start_after)
  if (resuming) {
    .labs <- vapply(chunks, `[[`, character(1), "label")
    .at   <- match(start_after, .labs)
    if (is.na(.at))
      stop(sprintf("start_after = '%s' is not a chunk label in %s. Labels are:\n  %s",
                   start_after, nm, paste(.labs, collapse = "\n  ")))
    if (verbose)
      .say(sprintf("    resuming %s after chunk '%s' (%d of %d chunk(s) held back)\n",
                   nm, start_after, .at, length(chunks)))
  }

  # Counted up front so the running commentary can say "12/47" rather than just
  # naming a chunk with no sense of how much is left.
  n_total <- length(chunks)
  k       <- 0L
  t_stage <- proc.time()[["elapsed"]]
  n_skip_fig <- 0L
  n_skip_dia <- 0L
  .prev_label <- NA_character_

  for (ch in chunks) {
    k <- k + 1L

    # ── Should this chunk run? ────────────────────────────────────────────────
    cls    <- .chunk_class(ch$label)
    is_fig <- skip_figures     && identical(cls, "figure")
    is_dia <- skip_diagnostics && identical(cls, "diagnostic")

    if (resuming && k <= .at)             next
    # Same rule as .chunk_class: skipping <X> by name skips params-<X> with it.
    if (ch$label %in% skip ||
        sub("^params-", "", ch$label) %in% skip) next
    if (is_fig)               { n_skip_fig <- n_skip_fig + 1L; next }
    if (is_dia)               { n_skip_dia <- n_skip_dia + 1L; next }
    if (identical(ch$eval, "false") && !(ch$label %in% force_eval)) next
    if (!length(ch$code) || !any(nzchar(trimws(ch$code)))) next

    if (verbose) {
      .el <- proc.time()[["elapsed"]] - t_stage
      .say(sprintf("    %s %2d/%2d %3.0f%%  %-34s  %s\n",
                   .bar(k / n_total, 12L), k, n_total, 100 * k / n_total,
                   ch$label, .fmt_dur(.el)))
    }

    # Where the run currently is, kept outside the pair environment so it
    # survives a failure. The error handler in section 4 turns this into the
    # exact start_after value needed to resume - which is the PREVIOUS chunk,
    # not this one: a chunk that died did not finish its work, so a resume has
    # to run it again (after whatever caused the failure has been fixed).
    assign(".last_chunk",
           list(stage = nm, label = ch$label, prev = .prev_label),
           envir = globalenv())

    # Chunks are evaluated as ONE expression block, not line by line, so a
    # multi-line function definition or an if/else spanning lines parses.
    expr <- tryCatch(parse(text = paste(ch$code, collapse = "\n")),
                     error = function(e)
                       stop(sprintf("Parse error in %s chunk '%s' (line %d):\n  %s",
                                    nm, ch$label, ch$line, conditionMessage(e))))

    # A figure or diagnostic that is running anyway (because it was switched
    # back on) still must not be able to fail the pair: nothing downstream and
    # no save manifest reads what it produces, so its error is reported and the
    # stage carries on. Core chunks are NOT caught - a real failure there has to
    # stop the pair, or the run would save half-computed results.
    if (identical(cls, "core")) {
      eval(expr, envir = env)
    } else {
      tryCatch(eval(expr, envir = env),
               error = function(e)
                 .say(sprintf("           !! %s chunk '%s' failed and was skipped: %s\n",
                              cls, ch$label, conditionMessage(e))))
    }

    # This chunk is done, so it becomes the resume point for whatever fails next.
    .prev_label <- ch$label

    # ── Overrides, applied after the chunk has had its say ───────────────────
    if (!is.null(overrides[[ch$label]])) {
      ov <- overrides[[ch$label]]
      for (v in names(ov)) {
        if (!exists(v, envir = env, inherits = FALSE))
          stop(sprintf(paste0("Override target '%s' does not exist after chunk '%s' of %s. ",
                              "The .qmd has been edited and the variable renamed or moved; ",
                              "fix the override map in section 3 of this driver."),
                       v, ch$label, nm))
        assign(v, ov[[v]], envir = env)
      }
      if (verbose) .say(sprintf("           -> set: %s\n", paste(names(ov), collapse = ", ")))
    }
  }

  # Printed so that if a later stage ever complains about a missing object, the
  # log says how many chunks were held back and of which kind.
  if (verbose && (n_skip_fig + n_skip_dia) > 0L)
    .say(sprintf("    (skipped %d figure + %d diagnostic chunk(s) in %s)\n",
                 n_skip_fig, n_skip_dia, nm))

  invisible(env)
}


# =============================================================================
# 3. THE OVERRIDE MAP - where each parameter lives
# =============================================================================
#
# One function so it can be rebuilt per pair with that pair's paths. Keyed by
# chunk label; see .run_qmd above for how they are applied.
#
# The .rds paths chain the stages together: each stage writes one and the next
# reads it, so both ends of every link are set here from the same expression
# and cannot drift.

.build_params <- function(p) list(

  # ── Stage 1 ────────────────────────────────────────────────────────────────
  "params-stage1" = list(
    las_1_path            = p$las_1,
    las_2_path            = p$las_2,
    ransac_ellipse_script = ransac_library_path
  ),
  # Each generated sample carries its own pivot as the PivotPoint = 1 marker, so
  # this needs no per-pair value - which is exactly why it is the mode used here.
  "params-pivot1" = list(
    pivot1_mode = "las"
  ),
  "params-export-pairs-csv" = list(
    export_pairs = TRUE
  ),

  # ── Stage 2 ────────────────────────────────────────────────────────────────
  # No entries. Stage 2 reads everything from the session, and its params chunks
  # run with their .qmd defaults - which for a batch are inert, since the
  # checkpoint path in params-load_Stage2 is never read (.stage_load is a no-op
  # under .batch_mode).

  # ── Stage 3 ────────────────────────────────────────────────────────────────
  # Belt and braces. params-stage3 now defaults to s3_paths = "inherit", so it
  # already carries what Stage 1 measured and these two are the same values
  # arriving by a second route. Kept because they cost nothing and they make the
  # driver's intent explicit: THIS pair's clouds, whatever a stage may think it
  # remembers. Set s3_paths <- "custom" in the .qmd and this entry still wins,
  # which is the correct precedence for a batch.
  "params-stage3" = list(
    las_1_path            = p$las_1,
    las_2_path            = p$las_2
  ),
  # The diagnostic truth source for §3.13.3. The truth CSV PATH still arrives by
  # .batch_opt(".batch_truth_csv") inside the chunk, because that chunk opens the
  # file it names; only the source selector is a parameter, and it lands here so
  # a batch is never at the mercy of whatever the .qmd was last left set to.
  "params-peak-score-residual-glm" = list(
    s3_truth_source = truth_shift_source
  ),
  # ransac_ellipse_script is deliberately absent: Stage 3 never assigns it, it
  # arrives from Stage 1 in the same session (and from the checkpoint by hand),
  # and it is already this driver's value by then.

  # ── Stages 4 and 5 ─────────────────────────────────────────────────────────
  # No entries. Everything a batch changes in these two - whether the diagnostic
  # .las snapshots are written, and where - is read inside the chunk through
  # .batch_opt(), because those chunks also decide a filename and suppress
  # .unique_path, which are not parameters. One mechanism per job: parameters
  # here, export destinations there.

  # ── Stage 6 ────────────────────────────────────────────────────────────────
  # Re-asserted rather than needed: both paths arrive from Stage 1 in this same
  # session, and Stage 6 assigns neither itself. Keyed to params-load-stage6,
  # the chunk that actually READS las_1, so the value lands immediately before
  # the read rather than beside an unrelated checkpoint path. .par_apply allows
  # this because it checks that the name exists in the environment, not that the
  # chunk authored it. 6.6's scoring switch and its two paths stay on
  # .batch_opt(), as above.
  "params-load-stage6" = list(
    las_1_path = p$las_1,
    las_2_path = p$las_2
  )
)


# ── The legacy override map - now empty, and kept only as a mechanism ────────
#
# Overrides here are applied AFTER the named chunk. Every parameter all six
# stages need has moved to .build_params() above, which lands BEFORE the
# functional chunk and so can reach a value that is read in the same breath as
# it is set - the thing an after-override structurally cannot do.
#
# The mechanism stays because it is still the right tool for a value that a
# LATER chunk reads, and because it existence-checks its targets, which makes it
# a usable escape hatch for a one-off. If you find yourself adding an entry
# here, check first whether the chunk in question has a params-* half: if it
# does, the entry belongs in .build_params() instead.
.build_overrides <- function(p) list(

  # Stages 1-3 are empty: their parameters live in params-* chunks now and are
  # set through .build_params() above, before the chunk that reads them.
  stage1 = list(),

  stage2 = list(),

  stage3 = list(),

  stage4 = list(),

  stage5 = list(),

  # Stage 6's las paths moved to .build_params("params-load_Stage6") when its
  # parameters were split out, so this map is now empty for every stage.
  stage6 = list()
)

# ── Why there are no save/load or EXPORT overrides above ─────────────────────
#
# The same trap caught two families of override, and the second was found only
# after a batch had already run on it.
#
# FIRST, the checkpoints. There used to be one override per stage, setting
# stage<N>_rds_path, stage<N>_save and
# stage<N>_load. Not one of them could ever take effect: overrides are applied
# AFTER a chunk runs, and every one of those chunks assigns the path and the
# toggles and then acts on them, all in the same chunk. The override landed on
# variables nothing would read again.
#
# What replaced them, in section 4.2's per-pair setup:
#
#   .batch_rds_dir / .batch_sample   the per-sample checkpoint name, derived
#                                    inside .stage_save() where the writing is
#   .batch_skip_saves                switching an individual checkpoint off
#   .batch_mode                      loads become no-ops, so the load paths
#                                    that these overrides used to set are not
#                                    read by anything
#
# The .qmd files keep their hardcoded interactive paths, untouched, and a stage
# run by hand still behaves exactly as it always did.
#
# SECOND, the exports - removed for the same reason, once a batch failed on one
# of them. Gone from above: "export-pairs-csv" (.pairs_out_dir,
# .pairs_basename), "pm-optim-export" (.pm_out_dir, .pm_basename),
# "export-aligned-cloud" (output_path) and "alignment-error-vs-truth"
# (score_against_truth, truth_csv_path, result_csv_dir). Every one of those
# chunks writes its file in the same chunk that names it.
#
# They failed in three different ways, which is the argument for not leaving a
# dead override in place looking like a handled path:
#
#   * the pairs CSV pointed at another machine's home directory and stopped the
#     pair outright, two hours in
#   * the aligned cloud silently went to one filename for every sample
#   * 6.6's score_against_truth stayed FALSE, so the chunk ran, wrote nothing,
#     and every pair came back "no_result" with no error anywhere
#
# What replaced them: .batch_pairs_dir, .batch_aligned_las, .batch_truth_csv,
# .batch_results_dir, .batch_score_truth, .batch_export_stage_clouds and
# .batch_export_dir, all assigned in 4.2 before the first chunk and read inside
# the chunks with .batch_opt(name, <the chunk's own default>).
#
# pm-optim-export got no replacement: that chunk no longer writes a CSV at all,
# so its override had been dead twice over.

# =============================================================================
# 4. RUN
# =============================================================================

# Output directories are created and write-tested by the preflight below.

## -- 4.1 Preflight - everything checkable, checked before anything runs ------
##
## A batch that dies forty minutes in because of a typo in a path has wasted
## forty minutes, and the typo was visible from the first second.
##
## Every problem is COLLECTED and reported together rather than stopping at the
## first one. Fixing six paths one failed run at a time is six runs; this way it
## is one. Checked here: the six stage .qmd files, the RANSAC ellipse library,
## the truth CSV (and its contents), both input clouds of every pair, and every
## output directory the run will write a .las or a .csv into - each one created
## if it is not there and then actually written to, because a directory that
## exists and cannot be written to fails exactly as late as one that does not
## exist at all.

.problems <- character(0)
.note <- function(...) .problems <<- c(.problems, paste0(...))

# ── The six stage documents ──────────────────────────────────────────────────
for (.st in names(qmd_paths)) {
  .p <- qmd_paths[[.st]]
  if (!nzchar(.p))            .note(sprintf("%s: no .qmd path given", .st))
  else if (!file.exists(.p))  .note(sprintf("%s .qmd not found:\n      %s", .st, .p))
  else if (dir.exists(.p))    .note(sprintf("%s .qmd is a directory, not a file:\n      %s", .st, .p))
}
if (!all(names(qmd_paths) %in% c("stage1", "stage2", "stage3", "stage4", "stage5", "stage6")) ||
    length(qmd_paths) != 6L)
  .note("qmd_paths must name all six stages (stage1 ... stage6); it currently has: ",
        paste(names(qmd_paths), collapse = ", "))

# ── The RANSAC ellipse library ───────────────────────────────────────────────
# Sourced by Stages 1 and 3 and by this driver. Missing, the fit cannot run at
# all, so it is worth its own line rather than surfacing as an obscure
# "could not find function" an hour in.
if (!nzchar(ransac_library_path)) {
  .note("ransac_library_path is empty")
} else if (!file.exists(ransac_library_path)) {
  .note("RANSAC ellipse library not found:\n      ", ransac_library_path)
}

# ── The truth CSV, and what is in it ─────────────────────────────────────────
truth_all <- NULL
if (!nzchar(truth_csv_path)) {
  .note("truth_csv_path is empty")
} else if (!file.exists(truth_csv_path)) {
  .note("Truth CSV not found:\n      ", truth_csv_path)
} else {
  truth_all <- tryCatch(utils::read.csv(truth_csv_path, stringsAsFactors = FALSE),
                        error = function(e) { .note("Truth CSV could not be read: ",
                                                    conditionMessage(e)); NULL })
  if (!is.null(truth_all)) {
    if (!"sample_name" %in% names(truth_all)) {
      .note("The truth CSV has no sample_name column - it is not a TestSampleGen table:\n      ",
            truth_csv_path)
    } else {
      .dupes <- unique(truth_all$sample_name[duplicated(truth_all$sample_name)])
      if (length(.dupes))
        .note("The truth CSV has more than one row for: ", paste(.dupes, collapse = ", "),
              "\n      Each sample must appear once, or there is no way to know which row is the truth.")
      .no_truth <- setdiff(sample_names, truth_all$sample_name)
      if (length(.no_truth))
        .note("No truth row for: ", paste(.no_truth, collapse = ", "),
              "\n      Add them to the consolidated CSV, or remove them from las_pairs.")
    }
  }
}

# ── Resolve las_pairs from the truth table ───────────────────────────────────
# Done here, after the truth CSV has been read and checked, because every path
# comes out of it. A name with no truth row was already reported above; it is
# dropped here rather than carried forward as a row of NA paths.
.fill_template <- function(tmpl, row) {
  out <- tmpl
  for (f in regmatches(tmpl, gregexpr("\\{[^}]+\\}", tmpl))[[1]]) {
    col <- substr(f, 2L, nchar(f) - 1L)
    val <- if (col %in% names(row)) as.character(row[[col]][1]) else NA_character_
    if (is.na(val) || !nzchar(val)) {
      .note(sprintf("las_2 template field {%s} is missing from the truth row for %s.",
                    col, row$sample_name[1]))
      return(NA_character_)
    }
    out <- sub(f, val, out, fixed = TRUE)
  }
  out
}

las_pairs <- data.frame(sample_name = character(0), las_1_path = character(0),
                        las_2_path = character(0), stringsAsFactors = FALSE)

if (!length(sample_names)) {
  .note("sample_names is empty - there is nothing to run.")
} else if (is.null(truth_all) || !"sample_name" %in% names(truth_all)) {
  .note("Without a readable truth table there is no way to resolve las_1 and las_2 paths.")
} else {
  if (anyDuplicated(sample_names))
    .note("sample_names has repeated entries: ",
          paste(unique(sample_names[duplicated(sample_names)]), collapse = ", "))

  .keep <- intersect(unique(sample_names), truth_all$sample_name)
  .rows <- lapply(.keep, function(nm) {
    tr <- truth_all[truth_all$sample_name == nm, , drop = FALSE]

    p1 <- if (nm %in% names(las_1_overrides)) unname(las_1_overrides[[nm]])
          else if ("sample_las" %in% names(tr)) as.character(tr$sample_las[1])
          else NA_character_
    if (is.na(p1) || !nzchar(p1))
      .note(sprintf("%s: no sample_las column in the truth table and no las_1 override.", nm))

    p2 <- if (nm %in% names(las_2_overrides)) unname(las_2_overrides[[nm]])
          else {
            f <- .fill_template(las_2_template, tr)
            if (is.na(f)) NA_character_ else file.path(las_2_dir, f)
          }

    data.frame(sample_name = nm, las_1_path = p1, las_2_path = p2,
               stringsAsFactors = FALSE)
  })
  if (length(.rows)) las_pairs <- do.call(rbind, .rows)

  cat(sprintf("\n  Resolved %d of %d requested sample(s) from the truth table.\n",
              nrow(las_pairs), length(unique(sample_names))))
  cat(sprintf("    las_1: the generator's own sample_las column\n"))
  cat(sprintf("    las_2: %s  in  %s\n", las_2_template, las_2_dir))
  if (length(las_1_overrides) || length(las_2_overrides))
    cat(sprintf("    overrides in force: %d las_1, %d las_2\n",
                length(las_1_overrides), length(las_2_overrides)))
  if (nrow(las_pairs))
    cat(sprintf("    e.g. %s\n           las_1 %s\n           las_2 %s\n",
                las_pairs$sample_name[1], las_pairs$las_1_path[1], las_pairs$las_2_path[1]))
}

# ── The input clouds, per pair ───────────────────────────────────────────────
if (!nrow(las_pairs)) .note("las_pairs is empty - there is nothing to run.")
for (.i in seq_len(nrow(las_pairs))) {
  for (.col in c("las_1_path", "las_2_path")) {
    .p <- las_pairs[[.col]][.i]
    if (is.na(.p) || !nzchar(.p))
      .note(sprintf("%s: %s is blank", las_pairs$sample_name[.i], .col))
    else if (!file.exists(.p))
      .note(sprintf("%s: %s not found:\n      %s", las_pairs$sample_name[.i], .col, .p))
  }
}
if (anyDuplicated(las_pairs$sample_name))
  .note("las_pairs has repeated sample_name(s): ",
        paste(unique(las_pairs$sample_name[duplicated(las_pairs$sample_name)]), collapse = ", "),
        "\n      Outputs are named from sample_name, so the second run would overwrite the first.")

# ── The output directories - the .las and .csv exports ───────────────────────
# Created here rather than merely tested, since creating one is the fix for it
# being absent. What is actually being checked is that the path can be made and
# then written to.
.out_targets <- c(
  `aligned .las export`      = out_dirs$aligned,
  `match-score optimisation .csv export` = out_dirs$optimisation,
  `temporal pairs .csv export` = out_dirs$pairs,
  `checkpoint .rds`          = out_dirs$checkpoints,
  `per-pair log`             = out_dirs$logs,
  `batch report .csv`        = out_dirs$report,
  `per-pair result .csv`     = file.path(out_dirs$report, "PerPairResults")
)
for (.lab in names(.out_targets)) {
  .d <- .out_targets[[.lab]]
  if (!dir.exists(.d))
    dir.create(.d, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(.d)) {
    .note(sprintf("%s directory could not be created:\n      %s", .lab, .d))
  } else {
    .probe <- file.path(.d, ".write_test")
    .ok <- tryCatch({ cat("", file = .probe); unlink(.probe); TRUE },
                    error = function(e) FALSE, warning = function(w) FALSE)
    if (!.ok) .note(sprintf("%s directory is not writable:\n      %s", .lab, .d))
  }
}

# ── The parameter map points at chunks that exist ────────────────────────────
# .par_apply() catches a bad VARIABLE name at run time, but a bad CHUNK label is
# silent by construction: no chunk answers to it, so nothing is applied and the
# batch runs on the .qmd default. Cheap to check here, invisible otherwise.
.pmap <- .build_params(list(name = "x", las_1 = "x", las_2 = "x",
                            truth_row_csv = "x", aligned_las = "x"))
.pmap <- .pmap[!vapply(.pmap, is.null, logical(1))]
if (length(.pmap)) {
  .all_labels <- unlist(lapply(qmd_paths, function(q)
    if (file.exists(q)) vapply(.read_qmd_chunks(q), `[[`, character(1), "label")
    else character(0)), use.names = FALSE)
  .no_chunk <- setdiff(names(.pmap), .all_labels)
  if (length(.no_chunk))
    .note("These .build_params() entries name no chunk in any stage .qmd:\n      ",
          paste(.no_chunk, collapse = "\n      "),
          "\n      They would be silently ignored. Check the label, or the chunk's name in the .qmd.")
}

if (length(.problems))
  stop(sprintf("Preflight found %d problem(s) - nothing has been run:\n\n  - %s\n",
               length(.problems), paste(.problems, collapse = "\n\n  - ")),
       call. = FALSE)

cat(sprintf("Preflight OK: 6 stage .qmd, RANSAC library, truth table (%d row(s)), %d cloud(s), %d output dir(s).\n",
            nrow(truth_all), 2L * nrow(las_pairs), length(.out_targets)))

## -- 4.1b Per-tree progress from inside the RANSAC library -------------------
## The ellipse fit is the longest step in the pipeline and its per-tree loop is
## inside the library, so nothing in this driver or in the .qmds can wrap a bar
## around it. The library emits progress events instead, and whatever handler is
## registered here receives them.
##
## Registering ransac_progress_cat_handler() switches the library from its
## built-in cli bar (stderr, invisible under this driver's message sink) to
## plain lines on stdout, which the split sink carries to both the console and
## the per-pair log. The .qmds are untouched: run one on its own with no handler
## registered and it still draws the cli bar it always did.
##
## The library is sourced by the STAGES, in their own per-pair environments, so
## it is sourced here as well purely to make the registration functions
## available. The handler is stored as an OPTION, not in an environment defined
## by the library, so it survives every one of those re-sourcings - an
## environment would be recreated empty by each source() call and the handler
## silently shadowed. One registration therefore covers the whole batch.
source(ransac_library_path)

if (exists("ransac_set_progress_handler")) {
  ransac_set_progress_handler(
    ransac_progress_cat_handler(step_frac = 0.05, min_seconds = 15))
  .say("Per-tree RANSAC progress: routed to stdout via the library handler.\n")
} else {
  .say("Per-tree RANSAC progress: NOT available - this copy of the RANSAC library\n",
       "  predates ransac_set_progress_handler(). The fit will run silently.\n")
}

cat(sprintf("Batch: %d pair(s) | figures %s | diagnostics %s | truth table %d row(s)\n",
            nrow(las_pairs), if (make_plots) "ON" else "OFF",
            if (run_diagnostics) "ON" else "OFF", nrow(truth_all)))
cat(sprintf("       truth scoring %s | stage cloud snapshots %s\n",
            if (score_truth) "ON" else "OFF (every pair will report no_result)",
            if (export_stage_clouds) "ON" else "OFF"))

# ── Resume settings, checked before anything expensive runs ──────────────────
if (xor(is.null(resume_env), is.null(resume_from)))
  stop("resume_env and resume_from must be set together, or both left NULL.")

.resume_skip_stage <- function(st) FALSE   # replaced below when resuming

if (!is.null(resume_from)) {
  if (!all(c("stage", "after") %in% names(resume_from)))
    stop("resume_from must be list(stage = ..., after = ...).")
  if (!resume_from$stage %in% names(qmd_paths))
    stop(sprintf("resume_from$stage '%s' is not one of: %s",
                 resume_from$stage, paste(names(qmd_paths), collapse = ", ")))
  if (!is.environment(resume_env))
    stop("resume_env is not an environment. Pass .last_failed_env, not its name.")

  # Stages BEFORE the resume stage are not run at all: their products are
  # already in the environment being carried over.
  .resume_at <- match(resume_from$stage, names(qmd_paths))
  .resume_skip_stage <- function(st) match(st, names(qmd_paths)) < .resume_at

  cat(sprintf("  RESUME: %s onward, starting after chunk '%s'%s\n",
              resume_from$stage, resume_from$after,
              if (nrow(las_pairs) > 1L)
                sprintf(" (first pair only - the other %d run normally)",
                        nrow(las_pairs) - 1L) else ""))
}
# Named rather than counted. A bare number invites the question "which four?",
# and the answer is the thing worth printing.
.kept_cp <- c("Stage 1", "Stage 2",
              if (keep_stage3_1_checkpoint) "Stage 3.1",
              if (keep_stage3_2_checkpoint) "Stage 3.2",
              if (keep_stage4_checkpoint)   "Stage 4",
              if (keep_stage5_checkpoint)   "Stage 5")
cat(sprintf("  checkpoints per pair: %d  (%s)\n",
            length(.kept_cp), paste(.kept_cp, collapse = ", ")))
if (!keep_stage3_2_checkpoint)
  cat("    Stage 3.2 is not written, so a resume must start at Stage 3.1 and re-run\n",
      "    3.12-3.17 - seconds of arithmetic, no refitting.\n", sep = "")
rm(.kept_cp)
if (make_plots)
  cat("  Figures are ON. They add minutes per pair and nothing downstream reads them;\n",
      "  the chunks in always_skip stay skipped regardless.\n", sep = "")
cat("\n")

## -- 4.2 The report file, and how it is kept ---------------------------------
##
## The report is rewritten after EVERY pair rather than once at the end. A batch
## of twenty pairs is most of a day, and a crash on pair nineteen used to take
## the other eighteen results with it - they existed only as a list in memory.
## Now each pair's row is on disk the moment that pair finishes.
##
## The name is resolved ONCE, here, and the same file is overwritten from then
## on, each time carrying every pair completed so far. So the file grows a row
## at a time and there is never more than one report per batch to reconcile. An
## existing report from an EARLIER batch is still not clobbered - that is what
## the suffix loop below is for; it just runs before the loop instead of after
## it, so the name cannot change midway.

report_path <- file.path(out_dirs$report, report_csv_name)
if (file.exists(report_path)) {
  .k <- 1L
  repeat {
    cand <- sub("\\.csv$", sprintf("_%03d.csv", .k), report_path)
    if (!file.exists(cand)) { report_path <- cand; break }
    .k <- .k + 1L
  }
  cat(sprintf("  a report of today's date already exists - writing to %s\n",
              basename(report_path)))
}

#' Build the report from the rows collected so far and write it.
#'
#' bind_rows rather than rbind: a failed pair contributes a five-column row and
#' a successful one contributes two hundred, and rbind cannot reconcile those.
#' The failures come out with NA in every error column, which reads correctly.
#'
#' Pairs that have not run yet are NULL in `results` and are dropped, so this is
#' the same function whether it is called after pair one or after the last.
.write_report <- function(rows, path) {
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(invisible(NULL))
  rep <- dplyr::bind_rows(rows)

  # sample_name and status first - the two columns anyone opening this looks at.
  # The marker-box residuals lead the shift columns because they are the primary
  # shift measurement (Stage 6 §6.6): they are metres measured between two sets
  # of points, where the err_shift_* terms are read off the solved matrix. Both
  # are kept - a run can be right about where the box went and wrong about how it
  # is turned, and only the rotation terms say so.
  .lead <- intersect(c("sample_name", "batch_status", "batch_error", "runtime_min",
                       "box_res_xyz_m", "box_res_xy_m", "box_res_x_m",
                       "box_res_y_m", "box_res_z_m",
                       "box_percorner_mean_m", "box_percorner_max_m",
                       "err_rot_total_deg", "err_rot_xy_mag_deg", "err_rot_xy_dir_deg",
                       "err_rot_z_deg", "err_shift_x_m", "err_shift_y_m",
                       "err_shift_z_m", "err_shift_xy_m", "err_shift_xyz_m"),
                     names(rep))
  rep <- rep[, c(.lead, setdiff(names(rep), .lead)), drop = FALSE]

  # Written to a temporary file and renamed, so an interrupt during the write
  # cannot leave a half-row of CSV where the previous complete report was.
  tmp <- paste0(path, ".tmp")
  ok <- tryCatch({
    utils::write.csv(rep, tmp, row.names = FALSE)
    if (file.exists(path)) unlink(path)
    file.rename(tmp, path)
    TRUE
  }, error = function(e) {
    cat(sprintf("  report NOT written (%s) - the rows are still in `results`\n",
                conditionMessage(e)))
    unlink(tmp)
    FALSE
  })
  invisible(if (ok) rep else NULL)
}

## -- 4.3 The loop ------------------------------------------------------------

results   <- vector("list", nrow(las_pairs))
run_log   <- vector("list", nrow(las_pairs))
t_batch   <- Sys.time()

for (i in seq_len(nrow(las_pairs))) {

  nm <- las_pairs$sample_name[i]
  t0 <- Sys.time()

  # Pair-level line. The ETA is the mean time per completed pair times the
  # number left - crude, but the pairs are broadly comparable in size, so it is
  # good enough to answer "finished before lunch or not". Absent on the first
  # pair, since there is nothing yet to average.
  .done_n   <- i - 1L
  .el_batch <- as.numeric(difftime(Sys.time(), t_batch, units = "secs"))
  .eta_batch <- if (.done_n > 0L) .el_batch / .done_n * (nrow(las_pairs) - .done_n) else NA_real_
  .say(sprintf("\n%s PAIR %d/%d  %s\n", .bar(.done_n / nrow(las_pairs)),
               i, nrow(las_pairs), nm))
  .say(sprintf("  batch elapsed %s | ~%s left\n",
               .fmt_dur(.el_batch), .fmt_dur(.eta_batch)))

  # The per-pair truth row, written out as its own one-row CSV. Stage 6 §6.6
  # and Stage 1's "csv" pivot mode both require exactly one row and refuse a
  # table, which is the correct behaviour for them - a stage should not have to
  # guess which sample it is looking at. Splitting the table here is what lets
  # the consolidated file be the thing the user maintains.
  truth_row_csv <- file.path(out_dirs$report, "PerPairResults",
                             paste0(nm, "_truth.csv"))
  utils::write.csv(truth_all[truth_all$sample_name == nm, , drop = FALSE],
                   truth_row_csv, row.names = FALSE)

  p <- list(
    name          = nm,
    las_1         = las_pairs$las_1_path[i],
    las_2         = las_pairs$las_2_path[i],
    truth_row_csv = truth_row_csv,
    aligned_las   = file.path(out_dirs$aligned, paste0(nm, "_aligned.las")),
    rds1          = file.path(out_dirs$checkpoints, paste0(nm, "_Stage1.rds")),
    rds2          = file.path(out_dirs$checkpoints, paste0(nm, "_Stage2.rds")),
    rds3a         = file.path(out_dirs$checkpoints, paste0(nm, "_Stage3_1.rds")),
    rds3b         = file.path(out_dirs$checkpoints, paste0(nm, "_Stage3_2.rds")),
    rds4          = file.path(out_dirs$checkpoints, paste0(nm, "_Stage4.rds")),
    rds5          = file.path(out_dirs$checkpoints, paste0(nm, "_Stage5.rds"))
  )

  ov <- .build_overrides(p)

  # ── A fresh environment per pair ────────────────────────────────────────────
  # Not the global environment, and not reused between pairs. The stages leave
  # hundreds of objects behind and many are guarded with exists() - a leftover
  # pairs_stage2 or run_counts from the previous sample would be picked up as
  # though it belonged to this one, and the result would be wrong in a way that
  # nothing would flag. Parented to globalenv() so the packages and this
  # driver's own helpers are still visible.
  #
  # The one exception is a resume: there the whole point is to keep the objects
  # the failed run had already computed, so its environment is adopted as-is.
  .resume_on <- i == 1L && !is.null(resume_env) && !is.null(resume_from)
  if (.resume_on) {
    env <- resume_env
    .say(sprintf("  RESUMING in the supplied environment (%d object(s)) at %s, after '%s'\n",
                 length(ls(envir = env, all.names = TRUE)),
                 resume_from$stage, resume_from$after))
  } else {
    env <- new.env(parent = globalenv())
  }
  assign("params", list(make_plots = make_plots), envir = env)

  # ── This pair's parameter overrides ────────────────────────────────────────
  # Read by .par_apply() at the end of every params-* chunk, keyed by that
  # chunk's label. This is the mechanism that replaces the after-the-fact
  # overrides for Stages 1-3: the params chunk assigns the .qmd's defaults, this
  # list overwrites the ones the batch is changing, and only then does the
  # functional chunk run. A parameter can no longer be read in the same chunk
  # that sets it, so there is no longer a window for the override to miss.
  #
  # The defaults stay in the .qmds. Naming a parameter here is how the BATCH
  # differs from a by-hand run, not a second copy of the parameter set - copying
  # all ~670 would just move the drift problem, which is what produced a Stage 3
  # warning that its hardcoded collection did not match the renamed sample.
  #
  # .par_apply() stops the pair if a name here is not one the chunk assigns, so
  # a rename in a .qmd surfaces immediately rather than as a batch that quietly
  # ran on the default.
  assign(".batch_params", .build_params(p), envir = env)

  # ── The one flag the stages read from this driver ──────────────────────────
  # .batch_mode tells the stage helpers that the whole pipeline is running in a
  # single session, and it changes exactly one thing: the checkpoints stop being
  # load-bearing.
  #
  #   * every save manifest becomes advisory - a missing object is reported and
  #     left out instead of stopping the pair
  #   * a failed write is reported, not thrown
  #   * .stage_load() does nothing at all
  #
  # Which is correct here and nowhere else: Stage N+1 already holds in memory
  # everything Stage N computed, so there is nothing for it to read back, and
  # the .rds files are a record of the run rather than the thing carrying it.
  # Stage 1's save chunk asking for thirty objects it could not reach must never
  # again be able to discard an hour of RANSAC.
  #
  # Run any stage .qmd by hand and this is absent, so the checkpoints are strict
  # again - which is right, because then the checkpoint IS the handover.
  assign(".batch_mode", TRUE, envir = env)

  # ── Where this pair's checkpoints go ───────────────────────────────────────
  # Not an override, and it cannot be one. Every save chunk assigns its .rds
  # path and then USES it a few lines later in the same chunk, so an override
  # keyed to that chunk arrives after the file has been written - which for a
  # batch meant every sample writing over the same hardcoded name (or, with
  # .unique_path, piling up _001, _002 next to it). These two values are in
  # place before the first chunk runs, and .stage_save() derives the per-sample
  # name from them plus the label the chunk already passes it:
  #
  #     <checkpoints>/<sample>_Stage1.rds, ..._Stage3_2.rds, ...
  #
  # which is exactly the naming this driver documents, now produced by the one
  # piece of code that actually does the writing.
  assign(".batch_rds_dir", out_dirs$checkpoints, envir = env)
  assign(".batch_sample",  nm,                   envir = env)

  # Checkpoints to switch off BY LABEL. Every chunk named here still RUNS -
  # later manifests inherit from earlier ones, so a chunk cannot be skipped
  # outright - only its write is suppressed. See keep_stage3_1_checkpoint and
  # keep_stage4_checkpoint in section 1 for which and why.
  assign(".batch_skip_saves",
         c(if (keep_stage3_1_checkpoint) character(0) else "Stage 3.1",
           if (keep_stage3_2_checkpoint) character(0) else "Stage 3.2",
           if (keep_stage4_checkpoint)   character(0) else "Stage 4",
           if (keep_stage5_checkpoint)   character(0) else "Stage 5"),
         envir = env)

  # ── Where this pair's EXPORTS go ───────────────────────────────────────────
  # Same problem as the checkpoints, same shape of fix. Five chunks across four
  # stages write a file to a path they set a few lines earlier in the same
  # chunk, so the overrides that used to be keyed to them in section 3 never
  # landed. What that cost, before these were read inside the chunks instead:
  #
  #   Stage 1 export-pairs-csv   wrote to a hardcoded path on ANOTHER machine
  #                              and failed the pair outright
  #   Stage 6 export-aligned-cloud   every sample queued behind one filename
  #   Stage 6 6.6 scoring        score_against_truth stayed FALSE, so the chunk
  #                              ran, did nothing, and every pair came back
  #                              "no_result" with no error to show for it
  #
  # The stage side reads these with .batch_opt(name, <its own default>), so an
  # older .qmd that has not been updated keeps its hardcoded default and this
  # driver still runs - it just does not redirect that one export.
  assign(".batch_pairs_dir",   out_dirs$pairs,        envir = env)
  assign(".batch_aligned_las", p$aligned_las,         envir = env)
  assign(".batch_truth_csv",   p$truth_row_csv,       envir = env)
  assign(".batch_results_dir", file.path(out_dirs$report, "PerPairResults"),
         envir = env)
  assign(".batch_score_truth", score_truth,           envir = env)

  # Stage 3's post-rotation and Stage 4's post-Z-shift .las snapshots. Both are
  # eyeballing aids for CloudCompare, both default to TRUE in the .qmds, and a
  # batch has nobody eyeballing - so off unless asked. They are not cheap:
  # Stage 4's reloads both clouds from disk to rebuild the transform before it
  # can write. Stage 6's aligned cloud, which is the actual product, is
  # unaffected by this switch.
  assign(".batch_export_stage_clouds", export_stage_clouds, envir = env)
  assign(".batch_export_dir",          out_dirs$aligned,    envir = env)

  # Stage 3's candidate-grain match-scoring CSV - one row per peak/trough match
  # candidate, carrying every raw factor score, the weights applied to it, the
  # aggregate, and the shift it fed. This is the table the match weighting is
  # tuned on, so it wants to exist for EVERY pair, named per pair.
  assign(".batch_optim_dir", out_dirs$optimisation, envir = env)

  # ── Logging, and why the two streams are treated differently ───────────────
  # sink(split = TRUE) on stdout sends output to the log AND the terminal.
  # sink(type = "message") has no split argument - it can only redirect - so
  # anything written to stderr vanishes from the terminal for the whole run.
  #
  # That is why the driver's own progress lines use cat() rather than message(),
  # and why the stages' .pb_* progress helpers cat() to stdout rather than using
  # cli (which draws on stderr and would be swallowed here). Both therefore
  # reach the terminal and the log. Warnings still go to the log only, which is
  # the right trade - they are for reading afterwards, the progress is for
  # watching now.
  log_path <- file.path(out_dirs$logs, paste0(nm, "_run.log"))
  if (write_logs) {
    con <- file(log_path, open = "wt")
    sink(con, split = TRUE)
    sink(con, type = "message")
  }

  ok <- TRUE
  err_msg <- NA_character_

  tryCatch({

    for (st in names(qmd_paths)) {
      # resume_from lets a re-run pick up where a failure left off: stages
      # before the named one are not run at all, and the named one starts after
      # the chunk that died. It only makes sense with an environment that
      # already holds the earlier work, so it is ignored unless one was handed
      # in via resume_env (see section 1.6).
      if (.resume_on && .resume_skip_stage(st)) {
        .say(sprintf("  -- %s (%d/%d)  skipped: resuming at %s\n",
                     st, which(names(qmd_paths) == st), length(qmd_paths),
                     resume_from$stage))
        next
      }
      .say(sprintf("  -- %s (%d/%d)\n", st, which(names(qmd_paths) == st), length(qmd_paths)))
      .run_qmd(qmd_paths[[st]], env,
               overrides    = ov[[st]],
               force_eval   = if (st == "stage6") "alignment-error-vs-truth" else character(0),
               skip         = unname(always_skip[names(always_skip) == st]),
               skip_figures = !make_plots,
               skip_diagnostics = !run_diagnostics,
               start_after  = if (.resume_on && identical(st, resume_from$stage))
                                resume_from$after else NULL,
               verbose      = TRUE)
    }

  }, error = function(e) {
    ok      <<- FALSE
    err_msg <<- conditionMessage(e)
    cat("  FAILED: ", err_msg, "\n", sep = "")

    # Where it died, and the two lines that pick up from there. Printed rather
    # than only stored, because the log is what gets read the next morning.
    .lc <- if (exists(".last_chunk", envir = globalenv())) get(".last_chunk", envir = globalenv()) else NULL
    if (!is.null(.lc)) {
      .st <- names(qmd_paths)[match(.lc$stage, basename(unlist(qmd_paths)))]
      cat(sprintf("  died in %s, chunk '%s'\n", .lc$stage, .lc$label))
      if (!is.na(.st) && keep_failed_env && !is.na(.lc$prev))
        cat(sprintf(paste0("  to resume without repeating that stage's earlier work:\n",
                           "    resume_env  <- .last_failed_env\n",
                           "    resume_from <- list(stage = \"%s\", after = \"%s\")\n",
                           "  ('%s' is the last chunk that COMPLETED, so '%s' runs again.)\n"),
                    .st, .lc$prev, .lc$prev, .lc$label))
    }
  })

  if (write_logs) { sink(type = "message"); sink(); close(con) }

  # ── Collect the result row ──────────────────────────────────────────────────
  # Stage 6 §6.6 builds `alignment_result`. Read it out of the environment
  # rather than off disk: the CSV it wrote went through .unique_path and may
  # carry a _001 suffix, and guessing the name back is a needless failure mode.
  if (ok && exists("alignment_result", envir = env, inherits = FALSE)) {
    row <- get("alignment_result", envir = env)
    row$batch_status   <- "ok"
    row$batch_error    <- NA_character_
    row$runtime_min    <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    row$log_file       <- if (write_logs) log_path else NA_character_
    results[[i]] <- row
  } else {
    # A failed pair still gets a row. A batch report with a silently missing
    # sample is worse than one that says which sample died and why.
    results[[i]] <- data.frame(
      sample_name  = nm,
      batch_status = if (ok) "no_result" else "failed",
      batch_error  = if (ok) "Stage 6 finished but wrote no alignment_result - check score_against_truth."
                     else err_msg,
      runtime_min  = as.numeric(difftime(Sys.time(), t0, units = "mins")),
      log_file     = if (write_logs) log_path else NA_character_,
      stringsAsFactors = FALSE
    )
  }

  if (!keep_checkpoints)
    unlink(unlist(p[c("rds1", "rds2", "rds3a", "rds3b", "rds4", "rds5")]))

  # A successful pair's environment is dropped: its results are on disk and in
  # `results`, and holding it would keep a full point cloud resident for the
  # rest of the batch. A FAILED pair's is kept when keep_failed_env is on -
  # everything the stages had computed up to the failure lives only here, and
  # for Stage 1 that is an hour of RANSAC per cloud.
  if (!ok && keep_failed_env) {
    assign(".last_failed_env", env, envir = globalenv())
    assign(".last_failed_pair", nm, envir = globalenv())
    .say(sprintf("  environment kept in .last_failed_env (%d object(s), sample '%s')\n",
                 length(ls(envir = env, all.names = TRUE)), nm))
  }
  rm(env); gc(verbose = FALSE)

  # ── The report, rewritten now that this pair has a row ─────────────────────
  # Every pair completed so far, over the top of the same file. If the next pair
  # brings the session down, this is what survives.
  .rep_now <- .write_report(results, report_path)
  if (!is.null(.rep_now))
    .say(sprintf("  report updated: %d pair(s), %d column(s) -> %s\n",
                 nrow(.rep_now), ncol(.rep_now), basename(report_path)))

  .say(sprintf("  %s %s in %s\n",
               if (ok) "OK  " else "FAIL",
               nm,
               .fmt_dur(as.numeric(difftime(Sys.time(), t0, units = "secs")))))

  if (!ok && !continue_on_error)
    stop("Stopping: continue_on_error is FALSE and '", nm, "' failed.")
}

## -- 4.4 Close out -----------------------------------------------------------
## The report is already on disk - it was rewritten after each pair. This is the
## same call once more, so the summary below reports the file as it finally
## stands rather than as it stood one pair ago.

report <- .write_report(results, report_path)
if (is.null(report)) report <- dplyr::bind_rows(Filter(Negate(is.null), results))

cat("═══ Batch complete ═══════════════════════════════════════════════════════\n")
cat(sprintf("  pairs run   : %d (%d ok, %d failed)\n",
            nrow(report),
            sum(report$batch_status == "ok"),
            sum(report$batch_status != "ok")))
cat(sprintf("  elapsed     : %.1f min\n",
            as.numeric(difftime(Sys.time(), t_batch, units = "mins"))))
cat(sprintf("  report      : %s (%d columns)\n", report_path, ncol(report)))
cat(sprintf("  aligned las : %s\n", out_dirs$aligned))
cat(sprintf("  optim csv   : %s  (<sample>_fmatch_optim.csv per pair)\n",
            out_dirs$optimisation))
cat(sprintf("  pairs csv   : %s\n", out_dirs$pairs))

if (any(report$batch_status != "ok")) {
  cat("\n  Failures:\n")
  .f <- report[report$batch_status != "ok", c("sample_name", "batch_error")]
  for (r in seq_len(nrow(.f)))
    cat(sprintf("    %s: %s\n", .f$sample_name[r], .f$batch_error[r]))
}

if ("err_shift_xyz_m" %in% names(report)) {
  .o <- report[report$batch_status == "ok", ]
  if (nrow(.o)) {
    cat("\n  Across successful pairs:\n")
    # The marker-box figures first, when the samples carried a box. na.rm
    # throughout: a batch can legitimately mix generated samples with real ones,
    # and one real pair should not turn the whole summary into NA.
    if ("box_res_xyz_m" %in% names(.o) && any(!is.na(.o$box_res_xyz_m))) {
      cat(sprintf("    box XYZ residual: median %.4f m | max %.4f m  (%d of %d pair(s) carried a box)\n",
                  median(.o$box_res_xyz_m, na.rm = TRUE),
                  max(.o$box_res_xyz_m, na.rm = TRUE),
                  sum(!is.na(.o$box_res_xyz_m)), nrow(.o)))
      cat(sprintf("    box per-corner  : median %.4f m | max %.4f m  (does NOT cancel - a rotation\n",
                  median(.o$box_percorner_mean_m, na.rm = TRUE),
                  max(.o$box_percorner_max_m, na.rm = TRUE)))
      cat("                      error shows here while the residual above averages to zero)\n")
    }
    cat(sprintf("    XYZ shift error : median %.4f m | max %.4f m  (from the solved matrix)\n",
                median(.o$err_shift_xyz_m), max(.o$err_shift_xyz_m)))
    cat(sprintf("    full rot error  : median %.4f deg | max %.4f deg\n",
                median(.o$err_rot_total_deg), max(.o$err_rot_total_deg)))
  }
}
cat("══════════════════════════════════════════════════════════════════════════\n")
