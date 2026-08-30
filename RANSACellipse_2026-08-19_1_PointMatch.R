# =============================================================================
# ransac_ellipse.R
#
# Pure function library — no code runs on source().
# Drop-in companion to ransac_circles.R that replaces circle fits with
# ellipse fits. Only the geometric primitive (circle → ellipse) and the
# associated parameter/output fields differ.
#
# ── Public API ────────────────────────────────────────────────────────────────
#
#   ransac_ellipse_fit_params(...)
#     Build / override a fit-parameter list.  All arguments have defaults.
#
#   ransac_ellipse_clean_params(...)
#     Build / override a cleaning-parameter list.  Identical signature to
#     ransac_clean_params(); cleaning operates on semi-major axis (a),
#     semi-minor axis (b), and centre XY instead of radius + centre XY.
#
#   custom_ransac_ellipse_fit(data, n_iterations, distance_threshold)
#     RANSAC ellipse fit with inlier ID tracking.
#     Samples 5 points per iteration, fits a general conic with
#     conicfit::EllipseDirectFit, converts to geometric parameters, and
#     scores hypotheses by angular coverage then inlier count.
#     Returns: list($ellipse, $inliers, $angle_segs, $n_iter)
#       $ellipse = c(cx, cy, a, b, angle_deg)   (a >= b, angle in [0, 180))
#
#   ransac_set_progress_handler(fn) / ransac_clear_progress_handler()
#     Register a function to receive per-tree progress events during a fit,
#     instead of the built-in cli bar. See SECTION 0 for the event contract.
#
#   ransac_progress_cat_handler(...)
#     A ready-made handler that prints throttled plain lines to stdout - the
#     one to register from a batch driver, where cli's stderr bar is invisible.
#
#   ransac_ellipse_fit_one(las, tree_ids, params, label)
#     Fit raw RANSAC ellipses for one LAS object.
#     Output columns per slice: slice_z, center_x, center_y, semi_a, semi_b,
#     angle_deg, inlier_ids, slice_point_ids.
#
#   ransac_ellipse_filter_window(raw_ellipses, bot, top)
#     Restrict raw ellipses to a Z sub-window.
#
#   ransac_ellipse_clean(raw_ellipses, params)
#     Spike-removal on centre XY, semi_a, semi_b independently.
#
#   ransac_ellipse_extract_ids(raw_ellipses, type)
#     Extract per-tree point_id sets ("inliers" | "appeared").
#
#   ransac_ellipse_run(las, tree_ids, fit_params, clean_params,
#                      window_bot, window_top, clean, seed, label, verbose)
#     Full pipeline: fit → filter window → [clean] → extract IDs.
#     Returns: $raw, $windowed_raw, $clean, $valid, $inlier_ids, $appeared_ids.
#
#   ransac_ellipse_run_multi(las_list, fit_params, clean_params,
#                            window_bot, window_top, clean, seed, verbose)
#     Convenience wrapper over ransac_ellipse_run() for named LAS lists.
#
#   ransac_ellipse_print_clean_summary(clean_result, label)
#
#   ellipse_point_distance(x, y, cx, cy, a, b, angle_deg)
#     Distance (m) from points to a fitted ellipse — the fitter's own metric,
#     exposed for diagnostics.
#
#   ellipse_count_inliers(x, y, cx, cy, a, b, angle_deg, distance_threshold)
#     Inlier count for an already-fitted ellipse, using the identical
#     classification path as the fitter (dead zone + confocal distance).
#     Use this for threshold-sensitivity sweeps so they match the fit.
#
# ── Key differences vs circle version ────────────────────────────────────────
#
#  Geometric primitive
#    Circle  : cx, cy, radius  (1 shape parameter)
#    Ellipse : cx, cy, semi_a, semi_b, angle_deg  (3 shape parameters)
#
#  Minimum sample size
#    Circle  : 3 points sufficient for a unique solution
#    Ellipse : 5 points required (general conic has 5 degrees of freedom)
#              → n_pts_sample fixed at 5 (was min(5, n) for circles)
#
#  Fit algorithm
#    Circle  : conicfit::CircleFitByPratt  (algebraic, Pratt constraint)
#    Ellipse : conicfit::EllipseDirectFit  (Bookstein algebraic fit)
#              followed by conicfit::AtoG() to obtain geometric parameters
#              (cx, cy, semi-axes a & b, tilt angle)
#
#  Validity guards
#    Ellipses add one extra rejection rule:
#      Both axes must be positive finite numbers
#
#  Inlier distance metric
#    Circle  : |sqrt((x-cx)²+(y-cy)²) - r| < threshold
#    Ellipse : closed-form confocal-hyperbola distance from the point to the
#              ellipse boundary (Rosin 1998; closed form Maalek & Lichti 2021),
#              screened by a certified dead zone that resolves ~99% of points
#              without evaluating any distance. See SECTION 1b.
#              A point is an inlier when that distance < distance_threshold.
#
#              distance_threshold is a genuine metre distance with the same
#              meaning for every stem, at every eccentricity. The v1 metric (a
#              scaled algebraic residual) was dimensionless and imposed a
#              tolerance of about 0.02·a metres, silently tightening as stems
#              got thinner. The v5 metric (bisection on Eberly's equation) was a
#              true distance but cost 50 vector iterations per hypothesis for
#              precision ~1e-13 m, twelve orders below the noise floor of the
#              point cloud; the confocal distance reproduces it to well under
#              0.05 mm for stem-like aspect ratios in closed form.
#
#  Cleaning
#    Circles : spike removal on cx, cy, radius (3 series)
#    Ellipses: spike removal on cx, cy, semi_a, semi_b (4 series);
#              orientation angle (angle_deg) is NOT cleaned because it wraps
#              at 180° and MAD-based spike detection is unreliable on circular
#              quantities.
#
#  Z frame for slicing
#    bot/top are height ABOVE EACH TREE'S OWN LOWEST POINT, not literal
#    absolute Z — no ground/DTM normalisation of the point cloud itself is
#    performed; each tree's own minimum Z stands in for its base and is added
#    back in per tree. `slice_z` in the output is therefore absolute Z (each
#    tree's own base + the offset within [bot, top]), so it is directly
#    comparable across trees whose bases sit at different absolute
#    elevations, without requiring the input LAS to be pre-normalised.
#
# =============================================================================


# =============================================================================
# SECTION 0 — Progress reporting
#
# THE PROBLEM THIS SOLVES
#
# The per-tree fit loop lives in here, so this file is the only place that
# knows how far along a fit is. Until now it drew that progress itself, with a
# cli bar. cli writes to stderr and redraws one line, which works when a .qmd
# is run interactively and fails silently in every other context that matters:
# under a message sink (a batch driver's log - sink(type = "message") cannot be
# split, unlike the stdout one), with stderr redirected to a file, while
# knitting, or from a worker process.
#
# The fix is not a different bar. It is to stop deciding. This library now
# EMITS progress events and lets the caller decide what to do with them:
#
#   register a handler  -> the library calls it once per tree and draws nothing
#   register nothing    -> the built-in cli bar, exactly as before
#
# So an interactive .qmd keeps the cli bar it has always had with no change at
# all, and a batch driver registers a handler that cat()s to stdout - which a
# split sink does carry - and gets a per-tree bar in its console and its log.
#
# THE HANDLER CONTRACT
#
# A handler is a function of one argument, a list with:
#
#   $event    "start" | "tick" | "done"
#   $label    the fit label ("las_1", "las_1 (align)", ...)
#   $i        trees completed so far (0 at "start")
#   $total    trees to be fitted
#   $tree_id  the TreeID just finished, NA at "start" and "done"
#   $elapsed  seconds since "start"
#
# It is called for its side effect; its return value is ignored. It is called
# inside tryCatch: a handler that errors is reported once and then disabled for
# the rest of the fit. A fit that has been running for six minutes must not die
# because a print statement was malformed.
# =============================================================================

# ── Where the handler is kept, and why it is an option ───────────────────────
#
# NOT in an environment defined by this file. This file gets source()d more than
# once per run and, under the batch driver, into a DIFFERENT environment each
# time - one per pair. An environment created here would be recreated by every
# source() call, and the fit functions (also redefined by that call) would
# resolve to the newest, empty one. A handler registered by the driver would be
# silently shadowed and never seen: exactly the failure this section exists to
# prevent.
#
# options() is process-wide and survives any number of source() calls into any
# number of environments, which is precisely the property needed.

#' Register a progress handler, replacing any previous one.
#'
#' @param fn A function of one argument (the event list described above), or
#'           NULL to remove the handler and restore the automatic default.
#' @return The previous handler, invisibly - so a caller can restore it.
ransac_set_progress_handler <- function(fn) {
  if (!is.null(fn) && !is.function(fn))
    stop("ransac_set_progress_handler(): fn must be a function or NULL.")
  old <- getOption("ransac.progress_handler")
  options(ransac.progress_handler = fn)
  invisible(old)
}

ransac_get_progress_handler <- function() getOption("ransac.progress_handler")

ransac_clear_progress_handler <- function() ransac_set_progress_handler(NULL)

#' Can cli actually draw a progress bar right now?
#'
#' cli redraws one line on the message stream, which requires a terminal it can
#' rewrite. That is false while knitting, with stderr redirected or sunk, and
#' under any driver holding a message sink. cli then draws NOTHING rather than
#' degrading, which is why fits have been running silently.
.ransac_cli_usable <- function() {
  requireNamespace("cli", quietly = TRUE) &&
    is.null(getOption("knitr.in.progress")) &&
    sink.number(type = "message") == 0L &&
    isTRUE(tryCatch(cli::is_dynamic_tty(), error = function(e) FALSE))
}

#' The handler actually used for a fit.
#'
#' Explicit registration wins. Otherwise cli is used when it can draw, and when
#' it cannot the printed-line handler is used instead. There is no longer a
#' configuration in which a fit reports nothing at all - that was the bug.
.ransac_resolve_handler <- function() {
  h <- getOption("ransac.progress_handler")
  if (!is.null(h)) return(h)
  if (.ransac_cli_usable()) return(NULL)          # NULL => use the cli bar
  ransac_progress_cat_handler()                   # cli cannot draw: print lines
}

#' A ready-made handler that prints plain lines to stdout.
#'
#' The one to use from a batch driver. Returns a CLOSURE so each fit gets its
#' own throttling state rather than sharing one counter across every call.
#'
#' Printing is throttled or a 400-tree fit buries the console: a line goes out
#' when either `step_frac` of the trees have been done or `min_seconds` have
#' passed, whichever comes first, plus always the first and last tree.
#'
#' @param step_frac   Print at most this often, as a fraction of the total.
#' @param min_seconds ...and at least this often on a slow fit.
#' @param width       Characters in the bar.
#' @param con         Where to write. stdout() by default, which is what a
#'                    split sink carries to both console and log. Anything on
#'                    stderr() would be swallowed by a message sink, which is
#'                    the whole reason this exists.
ransac_progress_cat_handler <- function(step_frac   = 0.05,
                                        min_seconds = 15,
                                        width       = 20L,
                                        con         = stdout()) {

  st <- new.env(parent = emptyenv())
  st$i_next <- 1L
  st$t_last <- 0
  st$step   <- 1L

  .dur <- function(sec) {
    if (!is.finite(sec)) return("--")
    if (sec < 90) sprintf("%.0fs", sec) else sprintf("%.1fm", sec / 60)
  }

  function(ev) {

    if (identical(ev$event, "start")) {
      st$step   <- max(1L, as.integer(ceiling(ev$total * step_frac)))
      st$i_next <- 1L
      st$t_last <- 0
      cat(sprintf("\n%s: fitting %d tree(s)\n", ev$label, ev$total), file = con)
      utils::flush.console()
      return(invisible(NULL))
    }

    if (identical(ev$event, "done")) {
      cat(sprintf("  %s: %d tree(s) fitted in %s\n",
                  ev$label, ev$i, .dur(ev$elapsed)), file = con)
      utils::flush.console()
      return(invisible(NULL))
    }

    due <- ev$i >= st$i_next ||
           (ev$elapsed - st$t_last) >= min_seconds ||
           ev$i >= ev$total
    if (!due) return(invisible(NULL))

    st$i_next <- ev$i + st$step
    st$t_last <- ev$elapsed

    frac <- if (ev$total > 0L) ev$i / ev$total else 1
    # ETA from the mean rate so far. Deliberately crude: per-tree cost varies by
    # an order of magnitude between a sparse stem and a dense one, so anything
    # more elaborate would be false precision. It answers "minutes or hours".
    eta  <- if (ev$i > 0L) ev$elapsed * (ev$total - ev$i) / ev$i else NA_real_
    n    <- as.integer(round(frac * width))

    cat(sprintf("  %s [%s] %d/%d %3.0f%% | %s elapsed | ~%s left\n",
                ev$label, paste0(strrep("#", n), strrep("-", width - n)),
                ev$i, ev$total, frac * 100, .dur(ev$elapsed), .dur(eta)),
        file = con)
    utils::flush.console()
    invisible(NULL)
  }
}

# ── Internal: the progress object used by the fit loop ────────────────────────
# One object whether the handler or the cli bar is in use, so the loop below has
# a single code path and cannot drift between the two.
#
# NOTE on cli and .envir: cli ties a bar's lifetime to a frame and destroys it
# when that frame exits. The bar is created in .ransac_pb_start(), so it MUST be
# handed to the caller's frame - the fit function - or it would be torn down the
# instant this returns and nothing would ever be drawn.
.ransac_pb_start <- function(label, total, show_progress, .envir = parent.frame()) {
  if (!isTRUE(show_progress) || !is.finite(total) || total < 1L) return(NULL)

  pb <- new.env(parent = emptyenv())
  pb$label   <- label
  pb$total   <- as.integer(total)
  pb$i       <- 0L
  pb$t0      <- proc.time()[["elapsed"]]
  pb$handler <- .ransac_resolve_handler()
  pb$broken  <- FALSE
  pb$cli_id  <- NULL

  if (is.null(pb$handler)) {
    pb$cli_id <- cli::cli_progress_bar(
      name   = label,
      total  = pb$total,
      format = paste0(label,
                      " {cli::pb_bar} {cli::pb_current}/{cli::pb_total}",
                      " trees | {cli::pb_elapsed}"),
      .envir = .envir)
  } else {
    .ransac_pb_emit(pb, "start", NA_integer_)
  }
  pb
}

.ransac_pb_emit <- function(pb, event, tree_id) {
  if (is.null(pb$handler) || isTRUE(pb$broken)) return(invisible(NULL))
  ev <- list(event = event, label = pb$label, i = pb$i, total = pb$total,
             tree_id = tree_id, elapsed = proc.time()[["elapsed"]] - pb$t0)
  tryCatch(pb$handler(ev), error = function(e) {
    pb$broken <- TRUE
    warning(sprintf("RANSAC progress handler failed (%s) - progress reporting disabled for this fit.",
                    conditionMessage(e)), call. = FALSE)
  })
  invisible(NULL)
}

.ransac_pb_tick <- function(pb, tree_id = NA_integer_) {
  if (is.null(pb)) return(invisible(NULL))
  pb$i <- pb$i + 1L
  if (!is.null(pb$cli_id)) cli::cli_progress_update(id = pb$cli_id)
  else .ransac_pb_emit(pb, "tick", tree_id)
  invisible(NULL)
}

.ransac_pb_done <- function(pb) {
  if (is.null(pb)) return(invisible(NULL))
  if (!is.null(pb$cli_id)) cli::cli_progress_done(id = pb$cli_id)
  else .ransac_pb_emit(pb, "done", NA_integer_)
  invisible(NULL)
}


# =============================================================================
# SECTION 1 — Parameter constructors
# =============================================================================

#' Build an ellipse fit-parameter list with documented defaults.
#'
#' @param bot                Lower edge of the first slice, as height (m)
#'                           ABOVE THIS TREE'S OWN LOWEST POINT — not an
#'                           absolute Z value. Applied per tree: each tree's
#'                           own minimum Z is used as its base and `bot` is
#'                           added to it. No ground/DTM normalisation of the
#'                           LAS is required for this to be correct.
#'                           NULL = adaptive: equivalent to bot = 0, i.e. use
#'                           each tree's own lowest point directly.
#' @param top                Upper edge of the last slice, as height (m)
#'                           above the same per-tree base as `bot`.
#'                           NULL = adaptive: use each tree's own highest point
#'                           instead of a shared fixed upper bound. Also
#'                           disables the height pre-filter that otherwise
#'                           drops trees whose own height range doesn't reach
#'                           a fixed `top`, since every tree is then fit up to
#'                           its own maximum height.
#' @param step               Distance between consecutive slice centres (m).
#' @param width              Thickness of each slice (m).
#' @param distance_threshold Max distance (metres) from a point to the ellipse
#'                           boundary for it to count as a RANSAC inlier,
#'                           measured by the confocal-hyperbola distance of
#'                           SECTION 1b. Directly comparable to the radial
#'                           threshold used for circles, for any eccentricity
#'                           and any stem size; for a circular fit the two
#'                           coincide exactly.
#' @param n_iterations       Number of RANSAC iterations per slice.
#' @param early_exit         Logical. If TRUE (default), stop iterating once a
#'                           hypothesis covers early_exit_segs 10-degree arc
#'                           segments. Set FALSE to always run all n_iterations.
#' @param early_exit_segs    Number of 10-degree arc segments (out of 36) that
#'                           trigger the early exit. Default 30 (~83% of arc).
#'                           Ignored when early_exit = FALSE.
ransac_ellipse_fit_params <- function(
    bot                = 0.5,
    top                = 5.5,
    height_min         = NULL,
    step               = 0.10,
    width              = 0.20,
    distance_threshold = 0.02,
    n_iterations       = 700L,
    iter_per_point     = NULL,
    iter_max           = 1000L,
    iter_min           = 1L,
    early_exit         = TRUE,
    early_exit_segs    = 30L,
    track_history      = FALSE
) {
  stopifnot(is.null(bot) || (is.numeric(bot) && length(bot) == 1L))
  stopifnot(is.null(top) || (is.numeric(top) && length(top) == 1L))
  stopifnot(is.null(height_min) || (is.numeric(height_min) && length(height_min) == 1L))
  stopifnot(is.logical(early_exit), length(early_exit) == 1L)
  stopifnot(is.numeric(early_exit_segs),
            early_exit_segs >= 1L, early_exit_segs <= 36L)
  list(
    bot                = bot,
    top                = top,
    # Minimum height above its own base a tree must reach to be fitted at all.
    # Independent of `top`: `top` caps how far UP a tree is sliced, height_min
    # decides WHICH trees are sliced. A tree taller than height_min but shorter
    # than top is fitted up to its own highest point.
    height_min         = height_min,
    step               = step,
    width              = width,
    distance_threshold = distance_threshold,
    n_iterations       = n_iterations,
    # Adaptive iteration budget. When iter_per_point is a number, the fixed
    # n_iterations is ignored and each slice instead gets
    #   clamp(iter_per_point * (points available to that fit), iter_min, iter_max)
    # so sparse slices are not given a budget far larger than their own sample
    # space, and dense ones are not starved. NULL keeps the fixed budget.
    iter_per_point     = iter_per_point,
    iter_max           = as.integer(iter_max),
    iter_min           = as.integer(iter_min),
    early_exit         = early_exit,
    early_exit_segs    = as.integer(early_exit_segs),
    track_history      = isTRUE(track_history)
  )
}

#' Build an ellipse cleaning-parameter list with documented defaults.
#'
#' Cleaning follows a two-stage flagging-first process:
#'   1. XY consensus filter — a slice is flagged spike_xy only when BOTH
#'      center_x AND center_y are independently identified as spikes
#'      (consensus requirement).
#'   2. Semi-axis filters   — spike_a and spike_b are detected independently;
#'      a slice is flagged spike_radii only when BOTH axes are spikes
#'      (XY-flagged rows are excluded from axis fits).
#' The composite flag  spike = spike_xy | spike_radii  is used by
#' ransac_ellipse_run() to produce the $valid list.  Orientation angle is not
#' cleaned because circular-quantity statistics are unreliable with MAD.
#'
#' @param mad_tol_xy    MAD multiplier for XY-centre spike removal.
#' @param local_k       Neighbours each side for the XY local-quadratic fit.
#' @param mad_tol_r     MAD multiplier for semi-axis spike removal (applied to
#'                      semi_a and semi_b independently).
#' @param local_k_r     Neighbours each side for the semi-axis local-quadratic fit.
ransac_ellipse_clean_params <- function(
    mad_tol_xy    = 1.5,
    local_k       = 2L,
    mad_tol_r     = 5.0,
    local_k_r     = 3L
) {
  list(
    mad_tol_xy    = mad_tol_xy,
    local_k       = local_k,
    mad_tol_r     = mad_tol_r,
    local_k_r     = local_k_r
  )
}


# =============================================================================
# SECTION 1b — Point-to-ellipse distance and inlier classification
# =============================================================================
#
# Two components, used together:
#   .ellipse_confocal_distance()  closed-form distance from a point to the
#                                 ellipse boundary (confocal hyperbola method)
#   .ellipse_inlier_mask()        threshold test, with a certified dead-zone
#                                 that answers most points without evaluating
#                                 the distance at all
#
# WHY NOT THE ALGEBRAIC RESIDUAL (the original v1 metric)
# -------------------------------------------------------
# Writing s = rho / r(theta), where rho is the point's distance from the centre
# and r(theta) = ab/sqrt((b cos)^2 + (a sin)^2) is the ellipse's polar radius in
# the point's own direction, the v1 quantity
#
#     |u^2/a^2 + v^2/b^2 - 1| * ab/sqrt((bu)^2 + (av)^2)   ==   |s^2 - 1| / s
#
# is a RATIO OF TWO LENGTHS and therefore dimensionless. Its threshold was a
# pure number, not a distance in metres: near the boundary it evaluates to
# approximately (2/r)*d_true, so a fixed cutoff silently tightened as stems got
# thinner (about 1 mm at a = 5 cm, about 6 mm at a = 30 cm) and drifted further
# with eccentricity. It applied the r/rho angular correction once where a
# further factor of about 2/r was also needed.
#
# WHY NOT THE EXACT (EBERLY) DISTANCE
# -----------------------------------
# Exact point-to-ellipse distance has no numerically stable closed form. The
# quartic of Safaee-Rad (1991) is exact but its analytic roots are unstable; the
# distance equation of Uteshev & Goncharova (2018) is exact but still requires
# root isolation with a delicate root-selection rule; the eigenvalue method of
# Chernov & Wijewickrema (2013) is the stable reference but is iterative and
# roughly 900x slower than closed-form alternatives. Bisection on Eberly's
# equation is a defensible exact method, but at 50 vector iterations per RANSAC
# hypothesis it dominated the run time — and precision far below the noise floor
# of the point cloud buys nothing. Every hypothesis is now scored by the same
# metric, including the winner: no hypothesis is measured more finely than the
# ones it is being compared against.
#
# THE CONFOCAL HYPERBOLA DISTANCE
# -------------------------------
# Confocal ellipses and hyperbolas intersect orthogonally. The confocal
# hyperbola through the point therefore crosses the ellipse very nearly along
# the boundary normal, and its intersection with the ellipse is an excellent
# stand-in for the true foot of the perpendicular. Introduced as an error-of-fit
# measure by Rosin (1998); closed form and validated against cylindrical pipe
# point clouds by Maalek & Lichti (2021).
#
# In the ellipse's own frame with X = |u|, Y = |v|, semi-axes ae >= be and
# f^2 = ae^2 - be^2:
#
#     T  = X^2 + Y^2 + f^2
#     D  = T^2 - 4 X^2 f^2
#     XI = ae X / sqrt((T + sqrt(D))/2)
#     YI = be sqrt(1 - (XI/ae)^2)
#     d  = sqrt((X - XI)^2 + (Y - YI)^2)
#
# Properties that matter here:
#   * Closed form. No loop, no convergence test, fully vectorised.
#   * At f = 0 (a circle) it reduces exactly to |rho - R|, so it agrees with the
#     circle pipeline's radial test to machine precision.
#   * No singularity at the centre or at the foci — unlike the Sampson distance
#     (infinite at the centre) or Harker & O'Leary (infinite at the foci).
#   * Points on an axis need no special case: at Y = 0 the formula returns
#     |X| - ae for |X| >= f and the correct off-axis foot inside the evolute.
#
# Measured against bisection on Eberly's equation, 300k points scattered within
# +/-25 mm of the boundary (99th percentile absolute error):
#
#     a = 0.10 m, b/a = 0.90 :  0.004 mm        a = 0.30 m, b/a = 0.90 : 0.0003 mm
#     a = 0.10 m, b/a = 0.80 :  0.027 mm        a = 0.30 m, b/a = 0.80 : 0.0014 mm
#     a = 0.10 m, b/a = 0.60 :  0.830 mm        a = 0.30 m, b/a = 0.60 : 0.0136 mm
#
# Its worst behaviour is for points deep INSIDE the ellipse on the major axis,
# far outside any plausible inlier threshold; restricted to points within 20 mm
# of the boundary the maximum error is 0.37 mm even at b/a = 0.6. This is three
# to four orders of magnitude below the v1 residual's error and about an order
# of magnitude below the point cloud's own noise.
#
# THE CERTIFIED DEAD ZONE
# -----------------------
# The inlier test only asks whether d < tau, not what d is. Two rigorous bounds
# follow from convexity, both expressed through the single scalar
#
#     s^2 = u^2/a^2 + v^2/b^2 :
#
#   Upper: the radial projection of the point onto the boundary is a point ON
#          the curve, at distance |s - 1| * r(theta) <= |s - 1| * a.
#   Lower: outside, the tangent line at that projection separates the point from
#          the whole (convex) ellipse, giving d >= |s - 1| / ||n|| >= |s - 1| * b;
#          inside, d >= (1 - s) * b by the same norm bound on the gradient.
#
# Hence, with no square root and no distance evaluation whatsoever:
#
#     |s - 1| <  tau/a   ==>   certified INLIER
#     |s - 1| >  tau/b   ==>   certified OUTLIER
#     otherwise          ==>   ambiguous, evaluate the confocal distance
#
# Implemented by comparing s^2 against the four constants (1 +/- tau/a)^2 and
# (1 +/- tau/b)^2, all precomputed once per hypothesis. The ambiguous band has
# relative width tau/b - tau/a, which closes as the stem becomes circular.
#
# Both bounds were checked against bisection on Eberly's equation over 7.2
# million random points spanning a in [0.05, 0.30] m and b/a in [0.3, 1.0]:
# zero violations. On a simulated slice (70% stem points with 5 mm noise, 30%
# clutter, tau = 20 mm) the share of points reaching the confocal evaluation is:
#
#     a = 0.05 m : 0.49% (b/a 0.95) to 2.44% (b/a 0.75)
#     a = 0.10 m : 0.23%            to 1.27%
#     a = 0.30 m : 0.08%            to 0.47%
#
# The composite predicate is therefore "true distance < tau" for the ~99% of
# points the bounds certify, and "confocal distance < tau" for the remainder —
# i.e. exactly the confocal criterion, up to disagreements confined to points
# lying within a few hundredths of a millimetre of the threshold itself.
#
# Net effect versus the previous exact implementation: the 50-iteration
# bisection is gone, and the metric that replaces it runs on well under 5% of
# the points. Measured end to end on 2000-point slices, the inlier test is about
# 50x faster than bisection on all points and about 1.8x faster than evaluating
# the confocal distance on all points.

#' Closed-form distance from points to an axis-aligned ellipse (confocal
#' hyperbola approximation; exact for circles).
#'
#' @param u,v  Numeric vectors: point coordinates in the ellipse's own frame
#'             (centre at the origin, axes along u and v).
#' @param a,b  Positive semi-axes. Either order is accepted.
#'
#' @return Numeric vector of distances in metres, same length as u.
.ellipse_confocal_distance <- function(u, v, a, b) {

  # Orient so ae is the major semi-axis; swap the coordinates to match, and
  # fold into the first quadrant (the nearest boundary point is always in the
  # point's own quadrant, so signs carry no information).
  if (a >= b) { ae <- a; be <- b; X <- abs(u); Y <- abs(v) }
  else        { ae <- b; be <- a; X <- abs(v); Y <- abs(u) }

  f2 <- ae * ae - be * be

  # ── Circle: the confocal construction degenerates, but its limit is simply
  # the radial distance, which is exact. Handled explicitly to avoid 0/0.
  if (f2 <= 0) return(abs(sqrt(X * X + Y * Y) - ae))

  Tt <- X * X + Y * Y + f2
  Dd <- Tt * Tt - 4 * X * X * f2
  Dd[Dd < 0] <- 0                       # guard rounding at Y = 0, |X| = f

  XI <- ae * X / sqrt(0.5 * (Tt + sqrt(Dd)))
  L2 <- (XI / ae)^2
  L2[L2 > 1] <- 1                       # guard rounding on the major axis
  YI <- be * sqrt(1 - L2)

  sqrt((X - XI)^2 + (Y - YI)^2)
}


#' Inlier mask for one ellipse hypothesis, with the certified dead zone.
#'
#' Returns the same answer as thresholding .ellipse_confocal_distance() over
#' every point, but evaluates that distance only for points the cheap bounds
#' cannot resolve (typically well under 5%).
#'
#' @param u,v  Point coordinates in the ellipse's own frame.
#' @param a,b  Positive semi-axes.
#' @param tau  Distance threshold in metres.
#'
#' @return Logical vector, TRUE for inliers.
.ellipse_inlier_mask <- function(u, v, a, b, tau) {

  s2 <- (u / a)^2 + (v / b)^2

  al <- tau / a                          # narrower band  -> certified inliers
  bl <- tau / b                          # wider band     -> certified outliers

  hi_in  <- (1 + al)^2
  lo_in  <- if (al < 1) (1 - al)^2 else -1  # -1: test is vacuously satisfied
  hi_out <- (1 + bl)^2
  lo_out <- if (bl < 1) (1 - bl)^2 else -1  # -1: no point is ever below it

  cert_in  <- s2 < hi_in  & s2 > lo_in
  cert_out <- s2 > hi_out | s2 < lo_out

  # Non-finite coordinates can never be inliers; keep them out of the band so
  # the subsetted distance call never sees them.
  ok       <- is.finite(s2)
  cert_in  <- cert_in  & ok
  cert_out <- cert_out | !ok

  mask <- cert_in
  band <- !(cert_in | cert_out)
  if (any(band))
    mask[band] <- .ellipse_confocal_distance(u[band], v[band], a, b) < tau

  mask
}


#' Distance from points to a fitted ellipse, in the original XY frame.
#'
#' Public wrapper around .ellipse_confocal_distance() for diagnostics and for
#' downstream code (e.g. threshold-sensitivity sweeps) that needs to reproduce
#' the fitter's metric on already-fitted ellipses. Use this rather than
#' re-deriving a distance by hand, so that the sweep and the fit agree.
#'
#' @param x,y        Point coordinates.
#' @param cx,cy      Ellipse centre.
#' @param a,b        Semi-axes (metres).
#' @param angle_deg  Ellipse tilt in degrees.
#'
#' @return Numeric vector of distances in metres.
ellipse_point_distance <- function(x, y, cx, cy, a, b, angle_deg) {
  tht <- angle_deg * pi / 180
  u   <- (x - cx) *  cos(tht) + (y - cy) * sin(tht)
  v   <- (x - cx) * -sin(tht) + (y - cy) * cos(tht)
  .ellipse_confocal_distance(u, v, a, b)
}


#' Count inliers for an already-fitted ellipse at a given threshold.
#'
#' Uses the identical classification path as the fitter (dead zone + confocal),
#' so threshold sweeps are directly comparable to the fit itself.
#'
#' @inheritParams ellipse_point_distance
#' @param distance_threshold Max distance (m) from the boundary, in metres.
#'
#' @return Integer count of inliers.
ellipse_count_inliers <- function(x, y, cx, cy, a, b, angle_deg,
                                  distance_threshold) {
  tht <- angle_deg * pi / 180
  u   <- (x - cx) *  cos(tht) + (y - cy) * sin(tht)
  v   <- (x - cx) * -sin(tht) + (y - cy) * cos(tht)
  sum(.ellipse_inlier_mask(u, v, a, b, distance_threshold), na.rm = TRUE)
}


# =============================================================================
# SECTION 1c — RANSAC ellipse fitter
# =============================================================================

#' Robust RANSAC ellipse fit with inlier ID tracking.
#'
#' Samples exactly 5 points per iteration (the minimum for a unique conic),
#' fits an ellipse with conicfit::EllipseDirectFit + conicfit::AtoG, scores
#' hypotheses by angular coverage (10-degree segments) then inlier count.
#'
#' Inlier criterion uses the closed-form confocal-hyperbola distance from the
#' point to the ellipse boundary (SECTION 1b):
#'   d((x,y), ellipse) < distance_threshold
#' A certified dead zone derived from s^2 = (u/a)^2 + (v/b)^2 classifies most
#' points before any distance is evaluated. Both the dead-zone bounds and the
#' distance are applied identically to every hypothesis, the winner included.
#'
#' @param data               Two-column matrix or data.frame of (X, Y) points.
#' @param n_iterations       Number of RANSAC iterations.
#' @param distance_threshold Max true distance (m) from the boundary to count
#'                           as an inlier.
#' @param early_exit         Logical. TRUE (default) breaks the loop as soon as
#'                           a hypothesis covers early_exit_segs segments.
#'                           FALSE always runs all n_iterations.
#' @param early_exit_segs    Segment threshold for the early exit (default 30L).
#'                           Must be in [1, 36]. Ignored when early_exit = FALSE.
#'
#' @return A list with:
#'   $ellipse        — numeric(5): cx, cy, semi_a, semi_b, angle_deg
#'                     (semi_a >= semi_b, angle_deg in [0,180), all NA if no fit)
#'   $inliers        — integer vector of inlier row indices into data
#'   $angle_segs     — integer: 10-degree arc segments covered by best inliers
#'   $n_iter         — integer: iteration at which best fit was found
#'   $early_exit_hit — logical: TRUE if the loop exited via the early-exit rule
custom_ransac_ellipse_fit <- function(data,
                                      n_iterations       = 500L,
                                      iter_per_point     = NULL,
                                      iter_max           = 1000L,
                                      iter_min           = 1L,
                                      distance_threshold = 0.06,
                                      early_exit         = TRUE,
                                      early_exit_segs    = 30L,
                                      track_history      = FALSE) {
  
  data <- as.matrix(data)
  n    <- nrow(data)
  
  if (n < 5L) {
    na5 <- c(mean(data[, 1]), mean(data[, 2]), NA_real_, NA_real_, NA_real_)
    return(list(ellipse = na5, inliers = integer(0),
                angle_segs = 0L, n_iter = 0L, early_exit_hit = FALSE,
                inlier_history = integer(0), segs_history = integer(0),
                n_points = n))
  }
  
  # ── Adaptive iteration budget ───────────────────────────────────────────────
  # Resolved here rather than by the caller because it depends on n, the number
  # of candidate points THIS fit can draw from — for a masked refit that is the
  # masked subset, not the whole slice. Every downstream artefact (the history
  # vectors, n_iter, early_exit_hit) refers to this resolved budget.
  if (!is.null(iter_per_point)) {
    n_iterations <- as.integer(max(iter_min, min(iter_max,
                                                 ceiling(iter_per_point * n))))
  }

  # ── Centre coordinates for numerical stability of EllipseDirectFit/AtoG ──
  x_off <- mean(data[, 1])
  y_off <- mean(data[, 2])
  data[, 1] <- data[, 1] - x_off
  data[, 2] <- data[, 2] - y_off
  
  x <- data[, 1]
  y <- data[, 2]
  
  best_ellipse    <- NULL
  best_inliers    <- 0L
  best_angle_segs <- 0L
  best_inlier_idx <- integer(0)
  best_n_iter     <- 0L
  early_exit_hit  <- FALSE

  # ── Convergence history (SECTION 1c option) ─────────────────────────────────
  # When track_history = TRUE, the state of the CURRENT BEST hypothesis is
  # recorded after every iteration: its inlier count and its angular coverage.
  # Both are needed downstream, because the acceptance rule ranks coverage above
  # inlier count — so the inlier series is not a running maximum, and the
  # iteration at which the winning hypothesis was adopted can only be identified
  # from the two together.
  hist_inliers <- if (track_history) rep(NA_integer_, n_iterations) else NULL
  hist_segs    <- if (track_history) rep(NA_integer_, n_iterations) else NULL
  iters_run    <- 0L
  
  p_densities <- dbscan::pointdensity(data, eps = 0.05)
  
  for (i in seq_len(n_iterations)) {
    
    idx           <- sample.int(n, 5L, prob = p_densities + 1)
    sample_points <- data[idx, , drop = FALSE]
    
    conic <- tryCatch({
      suppressWarnings(conicfit::EllipseDirectFit(sample_points))
    }, error = function(e) NULL)
    
    if (is.null(conic) || length(conic) < 6L || any(!is.finite(conic))) {
      if (track_history) { hist_inliers[i] <- best_inliers; hist_segs[i] <- best_angle_segs }
      iters_run <- i; next
    }
    
    geom <- tryCatch({
      suppressWarnings(conicfit::AtoG(matrix(conic, nrow = 1)))
    }, error = function(e) NULL)
    
    if (is.null(geom) || is.null(geom$ParG) || geom$exitCode != 1L) {
      if (track_history) { hist_inliers[i] <- best_inliers; hist_segs[i] <- best_angle_segs }
      iters_run <- i; next
    }
    
    # ── Read from $ParG: [cx, cy, a, b, t] ──────────────────────────────────
    pg  <- as.numeric(geom$ParG)
    cx  <- pg[1L]
    cy  <- pg[2L]
    a   <- pg[3L]
    b   <- pg[4L]
    tht <- pg[5L]   # radians
    
    if (any(!is.finite(c(cx, cy, a, b, tht))) || a <= 0 || b <= 0) {
      if (track_history) { hist_inliers[i] <- best_inliers; hist_segs[i] <- best_angle_segs }
      iters_run <- i; next
    }
    
    # Ensure a >= b (swap if AtoG returns them the other way)
    if (b > a) { tmp <- a; a <- b; b <- tmp; tht <- tht + pi / 2 }
    
    u  <- (x - cx) *  cos(tht) + (y - cy) * sin(tht)
    v  <- (x - cx) * -sin(tht) + (y - cy) * cos(tht)

    # ── Inlier test: dead zone first, confocal distance only where needed ────
    # .ellipse_inlier_mask() resolves the great majority of points from the
    # scalar s^2 = (u/a)^2 + (v/b)^2 alone, using the rigorous two-sided bound
    # |s - 1| * b <= d <= |s - 1| * a, and evaluates the confocal distance only
    # inside the residual ambiguous band. See SECTION 1b.
    inlier_mask <- .ellipse_inlier_mask(u, v, a, b, distance_threshold)
    inliers     <- sum(inlier_mask)
    
    ang          <- atan2(v[inlier_mask], u[inlier_mask]) * 180 / pi
    ang[ang < 0] <- ang[ang < 0] + 360
    angle_segs   <- length(unique(floor(ang / 10)))
    
    if (angle_segs > best_angle_segs ||
        (angle_segs == best_angle_segs && inliers >= best_inliers)) {
      angle_deg       <- tht * 180 / pi
      angle_deg       <- angle_deg %% 180
      # Add offsets back so centre is in original coordinate space
      best_ellipse    <- c(cx + x_off, cy + y_off, a, b, angle_deg)
      best_angle_segs <- angle_segs
      best_inliers    <- inliers
      best_inlier_idx <- which(inlier_mask)
      best_n_iter     <- i
    }
    
    if (track_history) {
      hist_inliers[i] <- best_inliers
      hist_segs[i]    <- best_angle_segs
    }
    iters_run <- i

    if (early_exit && best_angle_segs >= early_exit_segs) {
      early_exit_hit <- TRUE
      break
    }
  }
  
  if (is.null(best_ellipse)) {
    na5 <- c(mean(x) + x_off, mean(y) + y_off, NA_real_, NA_real_, NA_real_)
    return(list(ellipse        = na5,
                inliers        = integer(0),
                angle_segs     = 0L,
                n_iter         = 0L,
                early_exit_hit = FALSE,
                inlier_history = if (track_history)
                                   hist_inliers[seq_len(iters_run)] else integer(0),
                segs_history   = if (track_history)
                                   hist_segs[seq_len(iters_run)] else integer(0),
                n_points       = n))
  }
  
  list(ellipse        = best_ellipse,
       inliers        = best_inlier_idx,
       angle_segs     = best_angle_segs,
       n_iter         = best_n_iter,
       early_exit_hit = early_exit_hit,
       # Per-iteration state of the best hypothesis; empty unless track_history.
       inlier_history = if (track_history)
                          hist_inliers[seq_len(iters_run)] else integer(0),
       segs_history   = if (track_history)
                          hist_segs[seq_len(iters_run)] else integer(0),
       # Candidate points the fit drew from — the denominator for a fraction.
       n_points       = n)
}


# =============================================================================
# SECTION 2 — Low-level fitting (single LAS)
# =============================================================================

#' Fit RANSAC ellipses for every tree in one LAS object.
#'
#' Identical pipeline to ransac_fit_one() with the following differences:
#'   • custom_ransac_ellipse_fit() replaces custom_ransac_circle_fit()
#'   • Output data.frame has columns center_x, center_y, semi_a, semi_b,
#'     angle_deg instead of center_x, center_y, radius
#'
#' @param las       lidR LAS object; must have $TreeID and $point_id columns.
#' @param tree_ids  Integer vector of TreeIDs to process. NULL = all non-zero.
#' @param params    List from ransac_ellipse_fit_params().
#' @param label     Character label used in the progress bar.
#' @param seed      Integer RNG seed base. IMPORTANT: this is combined with
#'                  each tree's TreeID and re-applied via set.seed() at the
#'                  START of every tree's fit, so that a given tree's RANSAC
#'                  draws are reproducible regardless of which other trees
#'                  are present in `las` or what order `tree_ids` is
#'                  processed in. A single set.seed() call made once by the
#'                  caller before this function runs is NOT sufficient: it
#'                  only fixes the *first* tree's stream position, and every
#'                  subsequent tree's draws shift depending on how many
#'                  trees (and slices/iterations) were fitted before it.
#'                  That positional dependence is exactly why the same tree
#'                  could fit a different ellipse when processed as part of
#'                  a large "all trees" run vs. a small pre-filtered subset
#'                  (e.g. CSV-restricted trees) even with the same seed.
#'
#' @return Named list (names = as.character(tree_ids)); each element is a
#'         data.table with columns:
#'         slice_z, center_x, center_y, semi_a, semi_b, angle_deg,
#'         inlier_ids (list-col), slice_point_ids (list-col).
ransac_ellipse_fit_one <- function(las,
                                   tree_ids      = NULL,
                                   params        = ransac_ellipse_fit_params(),
                                   label         = "las",
                                   show_progress = TRUE,
                                   seed          = 42L) {

  if (is.null(tree_ids)) {
    tree_ids <- unique(las@data$TreeID)
    tree_ids <- tree_ids[!is.na(tree_ids) & tree_ids != 0L]
  }

  bot                <- params$bot
  top                <- params$top
  step               <- params$step
  width              <- params$width
  distance_threshold <- params$distance_threshold
  n_iterations       <- params$n_iterations
  early_exit         <- if (!is.null(params$early_exit))      params$early_exit      else TRUE
  early_exit_segs    <- if (!is.null(params$early_exit_segs)) params$early_exit_segs else 30L
  track_history      <- isTRUE(params$track_history)
  iter_per_point     <- params$iter_per_point
  iter_max           <- if (!is.null(params$iter_max)) params$iter_max else 1000L
  iter_min           <- if (!is.null(params$iter_min)) params$iter_min else 1L

  # ── Pre-filter: remove trees too short to be worth fitting ─────────────────
  # Keyed on `height_min`, NOT on `top`. The two are independent:
  #   height_min  decides WHICH trees are fitted
  #   top         caps how far UP each fitted tree is sliced
  # A tree taller than height_min but shorter than top is still fitted — it is
  # simply sliced up to its own highest point (see the per-tree bounds below).
  # Both are measured above each tree's OWN lowest point, not as absolute Z, so
  # the comparison is against (max_z - min_z).
  height_min <- params$height_min

  if (is.null(height_min)) {

    kept_ids <- tree_ids
    if (show_progress) cat(sprintf(
      "\n[%s] Height filter: skipped (height_min = NULL, every tree fitted).\n",
      label))

  } else {

    las_dt    <- data.table::as.data.table(las@data)[, .(TreeID, Z)]
    height_ok <- las_dt[
      TreeID %in% tree_ids,
      .(min_z = min(Z), max_z = max(Z)),
      by = TreeID
    ]
    height_ok$height_above_base <- height_ok$max_z - height_ok$min_z
    filtered_ids <- height_ok$TreeID[height_ok$height_above_base < height_min]
    kept_ids     <- height_ok$TreeID[height_ok$height_above_base >= height_min]

    if (show_progress) {
      if (length(filtered_ids) > 0L) {
        cat(sprintf(
          "\n[%s] Height filter (height_min = %.2f m above own base):\n",
          label, height_min))
        cat(sprintf(
          "  Trees REMOVED (height above base < %.2f m): %d  [TreeIDs: %s]\n",
          height_min, length(filtered_ids),
          paste(sort(filtered_ids), collapse = ", ")))
        cat(sprintf(
          "  Trees KEPT   (height above base >= %.2f m): %d\n",
          height_min, length(kept_ids)))
      } else {
        cat(sprintf(
          "\n[%s] Height filter: all %d trees reach %.2f m above their own base.\n",
          label, length(tree_ids), height_min))
      }
    }
  }

  tree_ids <- as.integer(kept_ids)
  n_trees  <- length(tree_ids)

  # Per-tree progress. Draws a cli bar when no handler is registered (the
  # original behaviour, unchanged), and otherwise hands each completed tree to
  # the registered handler and draws nothing itself. See SECTION 0.
  pb <- .ransac_pb_start(label, n_trees, show_progress)
  on.exit(.ransac_pb_done(pb), add = TRUE)

  lapply(
    setNames(as.character(tree_ids), as.character(tree_ids)),
    function(tid) {

      # ── Per-tree deterministic reseed ─────────────────────────────────────
      # Fixes cross-run irreproducibility: a single set.seed() call made once
      # for the whole `las` object (as the caller used to do) makes every
      # tree's RANSAC draws depend on the RNG stream position left behind by
      # every tree fitted before it. That position changes whenever the set
      # or order of tree_ids changes (e.g. fitting all trees in automatic
      # spatial-pairing mode vs. fitting only a CSV-restricted subset in
      # manual mode) even though the seed itself is identical. Reseeding
      # here, keyed on the tree's own TreeID, makes each tree's ellipse fit
      # depend only on (seed, TreeID) and nothing else — identical results
      # no matter what else is in `las` or what order tree_ids is processed.
      set.seed(as.integer(seed) + as.integer(tid))

      # !is.na() first, and a data.table subset rather than las[idx]. A bare
      # `TreeID == tid` yields NA on every point whose TreeID is NA - which the
      # test-sample generator's marker points are, along with any unsegmented
      # return - and a logical index carrying NA is not reliably read as FALSE
      # by every subsetting method it might reach. Filtering the coordinates
      # directly also skips rebuilding a LAS object per tree, which is header,
      # extent and validation work for a table this code immediately discards.
      idx <- !is.na(las@data$TreeID) & las@data$TreeID == as.integer(tid)
      pts <- data.table::as.data.table(las@data)[idx, .(X, Y, Z, point_id)]

      # ── Adaptive bounds ────────────────────────────────────────────────────
      # bot/top are height ABOVE THIS TREE'S OWN LOWEST POINT, not literal
      # absolute Z — no ground/DTM normalisation of the point cloud itself is
      # performed or required; each tree's own minimum Z stands in for its
      # base. This makes bot/top behave the way the calling code's parameter
      # comments describe them ("m above tree base") regardless of whether
      # the LAS Z values are true elevation or already ground-normalised.
      # NULL bot  → 0 m above this tree's own lowest point (same as before).
      # NULL top  → this tree's own highest point, so every tree is fit
      #             across its full available height.
      # A fixed top is a CEILING, not a requirement: a tree shorter than top is
      # sliced up to its own highest point rather than being skipped. Whether a
      # tree is fitted at all is decided by height_min in the filter above.
      z_min_tree <- min(pts$Z)
      bot_i <- z_min_tree + (if (is.null(bot)) 0 else bot)
      top_i <- if (is.null(top)) max(pts$Z) else min(z_min_tree + top, max(pts$Z))

      z_seq <- if (top_i - step >= bot_i) seq(bot_i, top_i - step, by = step)
               else numeric(0)

      # NA row helper — five shape fields instead of one radius field
      .make_na_row <- function(zi) {
        data.frame(
          slice_z         = zi + width / 2,
          center_x        = NA_real_,
          center_y        = NA_real_,
          semi_a          = NA_real_,
          semi_b          = NA_real_,
          angle_deg       = NA_real_,
          inlier_ids      = I(list(integer(0))),
          slice_point_ids = I(list(integer(0))),
          n_cand_pts      = NA_integer_,
          inlier_history  = I(list(integer(0))),
          segs_history    = I(list(integer(0)))
        )
      }

      out <- lapply(z_seq, function(zi) {
        pts_i <- pts[Z >= zi & Z < zi + width]
        # Ellipse needs >= 5 points (not >= 3 as for circles)
        if (nrow(pts_i) < 5L) return(.make_na_row(zi))

        full_slice_pids <- if ("point_id" %in% names(pts_i)) pts_i$point_id else integer(0)

        # ── RANSAC ellipse fit ───────────────────────────────────────────────
        ransac_result <- tryCatch(
          custom_ransac_ellipse_fit(
            pts_i[, .(X, Y)],
            n_iterations       = n_iterations,
            iter_per_point     = iter_per_point,
            iter_max           = iter_max,
            iter_min           = iter_min,
            distance_threshold = distance_threshold,
            early_exit         = early_exit,
            early_exit_segs    = early_exit_segs,
            track_history      = track_history
          ),
          error = function(e) list(ellipse        = rep(NA_real_, 5),
                                   inliers        = integer(0),
                                   early_exit_hit = FALSE,
                                   inlier_history = integer(0),
                                   segs_history   = integer(0),
                                   n_points       = nrow(pts_i))
        )
        ell <- ransac_result$ellipse

        inlier_ids <- if (length(ransac_result$inliers) > 0 &&
                          "point_id" %in% names(pts_i))
          pts_i$point_id[ransac_result$inliers] else integer(0)

        slice_pids <- if (!is.na(ell[1])) full_slice_pids else integer(0)

        data.frame(
          slice_z         = zi + width / 2,
          center_x        = ell[1],
          center_y        = ell[2],
          semi_a          = ell[3],
          semi_b          = ell[4],
          angle_deg       = ell[5],
          early_exit_hit  = isTRUE(ransac_result$early_exit_hit),
          inlier_ids      = I(list(inlier_ids)),
          slice_point_ids = I(list(slice_pids)),
          # Convergence record — empty vectors unless params$track_history.
          # n_cand_pts is the denominator for the inlier fraction: the points
          # this fit could draw from, which for a masked refit is the masked
          # subset rather than the whole slice.
          n_cand_pts      = if (is.null(ransac_result$n_points))
                              nrow(pts_i) else as.integer(ransac_result$n_points),
          inlier_history  = I(list(if (is.null(ransac_result$inlier_history))
                                     integer(0) else ransac_result$inlier_history)),
          segs_history    = I(list(if (is.null(ransac_result$segs_history))
                                     integer(0) else ransac_result$segs_history))
        )
      })

      .ransac_pb_tick(pb, tree_id = as.integer(tid))
      data.table::rbindlist(out, fill = TRUE)
    }
  )
}


# =============================================================================
# SECTION 3 — Window filtering
# =============================================================================

#' Restrict a raw-ellipses list to slices whose slice_z falls in [bot, top].
#'
#' @param raw_ellipses Named list as returned by ransac_ellipse_fit_one().
#' @param bot          Lower Z bound, inclusive (same absolute Z units as
#'                     `slice_z`, e.g. height above ground).
#' @param top          Upper Z bound, inclusive (same absolute Z units as
#'                     `slice_z`).
#'
#' @return Filtered named list (same structure, subset of rows per tree).
ransac_ellipse_filter_window <- function(raw_ellipses, bot, top) {
  lapply(raw_ellipses, function(df) {
    if (is.null(df) || nrow(df) == 0L) return(df)
    df[!is.na(df$slice_z) & df$slice_z >= bot & df$slice_z <= top, ,
       drop = FALSE]
  })
}


# =============================================================================
# SECTION 4 — Spike removal and cleaning
# =============================================================================
#
# Design — flagging-first, filter-second
# ──────────────────────────────────────
# .ellipse_flag_one_tree() annotates each row of a single tree's ellipse
# data.frame with boolean flag columns:
#
#   spike_xy     TRUE when BOTH center_x AND center_y are local spikes
#                (consensus requirement, same as doc reference).
#
#   spike_a      TRUE when semi_a is a local spike (independent of spike_b).
#   spike_b      TRUE when semi_b is a local spike (independent of spike_a).
#   spike_radii  TRUE when BOTH spike_a AND spike_b are TRUE, AND the row
#                was not already flagged by spike_xy.
#                XY-flagged rows are excluded from axis fits.
#
#   spike        Composite: spike_xy | spike_radii.
#                Used by ransac_ellipse_run() to produce the $valid list.
#
# Counter columns (scalar, same value on every row — store once, read easily):
#   n_removed_xy, n_removed_radii
#
# This design mirrors the reference chunk (clean_ellipse_df_ax_flagged) with
# the following adaptations for the library context:
#   • Position axis is slice_z (not stem_len; same meaning, different name).
#   • The function returns the fully-annotated data.frame, not a filtered one;
#     filtering happens one level up in ransac_ellipse_clean().
#   • Error handling is in ransac_ellipse_clean(), not inside the flag function.
# =============================================================================

# Internal: two-pass local-quadratic spike detector.
# Matches compute_local_dev / two_pass from the reference chunk exactly,
# with the only change that pos is centred per-observation (p <- p - p0)
# to improve numerical stability of the quadratic fit.
.ellipse_two_pass_spikes <- function(val, pos, mad_tol, local_k) {
  n <- length(val)
  if (n < 5L) return(rep(FALSE, n))

  compute_dev <- function(exclude) {
    vapply(seq_len(n), function(i) {
      if (is.na(val[i])) return(NA_real_)
      left_idx  <- tail(which(seq_len(n) <  i & !exclude), local_k)
      right_idx <- head(which(seq_len(n) >  i & !exclude), local_k)
      nbr_idx   <- c(left_idx, right_idx)
      if (length(nbr_idx) < 2L) return(NA_real_)
      p  <- pos[nbr_idx];  v <- val[nbr_idx]
      ok <- is.finite(p) & is.finite(v)
      p  <- p[ok];          v <- v[ok]
      if (length(p) < 2L) return(NA_real_)
      p0       <- pos[i]
      p        <- p - p0          # centre for numerical stability
      use_quad <- length(p) >= 3L && length(unique(p)) >= 3L
      X        <- if (use_quad) cbind(1, p, p^2) else cbind(1, p)
      if (any(!is.finite(X))) return(NA_real_)
      coefs <- lm.fit(X, v)$coefficients
      if (any(!is.finite(coefs))) return(NA_real_)
      # Prediction at centred position 0 (= the current observation's position)
      pred_i <- coefs[1]
      val[i] - pred_i
    }, numeric(1))
  }

  dev1  <- compute_dev(exclude = rep(FALSE, n))
  mad1  <- mad(dev1, na.rm = TRUE, constant = 1)
  gross <- !is.na(dev1) & !is.na(mad1) & mad1 > 0 & abs(dev1) > mad_tol * 3 * mad1
  dev2  <- compute_dev(exclude = gross)
  mad2  <- mad(dev2, na.rm = TRUE, constant = 1)
  if (!is.na(mad2) && mad2 > 0)
    !is.na(dev2) & abs(dev2) > mad_tol * mad2
  else
    rep(FALSE, n)
}

# Internal: annotate one tree's ellipse data.frame with spike flags.
# Returns the data.frame with flag columns added; does NOT filter rows.
.ellipse_flag_one_tree <- function(df, params) {

  df <- df[!is.na(df$center_x) & !is.na(df$semi_a) & !is.na(df$semi_b), ]
  df <- df[order(df$slice_z), ]

  # Initialise all flag columns to FALSE
  df$spike_xy     <- FALSE
  df$spike_a      <- FALSE
  df$spike_b      <- FALSE
  df$spike_radii  <- FALSE
  df$spike        <- FALSE
  df$n_removed_xy     <- 0L
  df$n_removed_radii  <- 0L
  df$spike_cx     <- FALSE
  df$spike_cy     <- FALSE

  # Early return: too few rows for any meaningful spike detection
  if (nrow(df) < 5L) {
    df$spike <- FALSE
    return(df)
  }

  pos <- df$slice_z

  # ── Stage 1: XY consensus spike filter ────────────────────────────────────
  # A slice is spike_xy only when BOTH cx AND cy are independently flagged
  # (consensus requirement from the reference chunk).
  spike_cx            <- .ellipse_two_pass_spikes(df$center_x, pos,
                                                   params$mad_tol_xy, params$local_k)
  spike_cy            <- .ellipse_two_pass_spikes(df$center_y, pos,
                                                   params$mad_tol_xy, params$local_k)
  df$spike_cx         <- spike_cx
  df$spike_cy         <- spike_cy
  spike_xy            <- spike_cx & spike_cy
  df$spike_xy         <- spike_xy
  df$n_removed_xy     <- sum(spike_xy)

  # ── Stage 2: Semi-axis spike filters ──────────────────────────────────────
  # Exclude XY-flagged rows from axis fits.
  # spike_a and spike_b are detected independently; spike_radii requires BOTH.
  # This matches the reference: spike_radii <- spike_a & spike_b
  exclude_for_r       <- spike_xy
  a_for_fit           <- df$semi_a;  a_for_fit[exclude_for_r] <- NA_real_
  b_for_fit           <- df$semi_b;  b_for_fit[exclude_for_r] <- NA_real_

  spike_a             <- .ellipse_two_pass_spikes(a_for_fit, pos,
                                                   params$mad_tol_r, params$local_k_r) &
                         !exclude_for_r
  spike_b             <- .ellipse_two_pass_spikes(b_for_fit, pos,
                                                   params$mad_tol_r, params$local_k_r) &
                         !exclude_for_r
  spike_radii         <- spike_a & spike_b

  df$spike_a          <- spike_a
  df$spike_b          <- spike_b
  df$spike_radii      <- spike_radii
  df$n_removed_radii  <- sum(spike_radii)

  # ── Composite flag ─────────────────────────────────────────────────────────
  df$spike <- spike_xy | spike_radii

  df
}

#' Apply spike-removal and axis cleaning to a per-tree ellipse list.
#'
#' Uses a flagging-first, filter-second approach:
#'   1. .ellipse_flag_one_tree() annotates each row with spike_xy, spike_a,
#'      spike_b, spike_radii, and the composite spike flag.
#'   2. Rows where spike == TRUE are removed to produce $valid.
#'   3. The full pre-filter data.frame (with flags) is stored as $flagged for
#'      diagnostic inspection.
#'
#' On error for a single tree the pre-filter data.frame is returned with all
#' spike flags set to TRUE, so the tree drops out of $valid downstream.
#'
#' @param raw_ellipses Named list as returned by ransac_ellipse_fit_one() or
#'                     ransac_ellipse_filter_window().
#' @param params       List from ransac_ellipse_clean_params().
#'
#' @return Named list matching raw_ellipses. Each element is a list with:
#'   $flagged        — full data.frame with spike_* flag columns
#'   $valid          — data.frame with spike == FALSE rows only
#'   $ellipses_before, $ellipses_after
#'   $n_removed_xy, $n_removed_radii
#'   $prop_kept
ransac_ellipse_clean <- function(raw_ellipses,
                                 params = ransac_ellipse_clean_params()) {
  lapply(raw_ellipses, function(df) {

    df <- df[!is.na(df$center_x) & !is.na(df$semi_a) & !is.na(df$semi_b), ]
    n_before <- nrow(df)

    # Short-circuit: too few rows to clean — return unflagged, all valid
    if (n_before < 5L) {
      df$spike_cx <- df$spike_cy <- df$spike_xy <- df$spike_a <- df$spike_b <-
        df$spike_radii <- df$spike <- FALSE
      df$n_removed_xy <- df$n_removed_radii <- 0L
      return(list(
        flagged         = df,
        valid           = df,
        ellipses_before = n_before,
        ellipses_after  = n_before,
        n_removed_xy    = 0L,
        n_removed_radii = 0L,
        prop_kept       = NA_real_
      ))
    }

    flagged <- tryCatch(
      .ellipse_flag_one_tree(df, params),
      error = function(e) {
        warning(sprintf("[ransac_ellipse_clean] flagging failed for one tree: %s",
                        conditionMessage(e)))
        # On error: mark everything as spiked so the tree is dropped cleanly
        df$spike_xy <- df$spike_a <- df$spike_b <-
          df$spike_radii <- df$spike <- TRUE
        df$n_removed_xy <- df$n_removed_radii <- 0L
        df
      }
    )

    valid      <- flagged[!flagged$spike, ]
    n_after    <- nrow(valid)

    list(
      flagged         = flagged,
      valid           = valid,
      ellipses_before = n_before,
      ellipses_after  = n_after,
      n_removed_xy    = if ("n_removed_xy"    %in% names(flagged)) flagged$n_removed_xy[1]    else 0L,
      n_removed_radii = if ("n_removed_radii" %in% names(flagged)) flagged$n_removed_radii[1] else 0L,
      prop_kept       = if (n_before > 0L) n_after / n_before else NA_real_
    )
  })
}


# =============================================================================
# SECTION 5 — Point-ID extraction
# =============================================================================

#' Extract per-tree point_id sets from a raw (or windowed) ellipse list.
#'
#' @param raw_ellipses Named list as returned by ransac_ellipse_fit_one().
#' @param type         "inliers"  — only RANSAC inlier point_ids (default)
#'                     "appeared" — all point_ids from any fitted slice window
#'
#' @return Named list of integer vectors.
ransac_ellipse_extract_ids <- function(raw_ellipses,
                                       type = c("inliers", "appeared")) {
  type <- match.arg(type)
  col  <- if (type == "inliers") "inlier_ids" else "slice_point_ids"
  lapply(raw_ellipses, function(raw) {
    if (is.null(raw) || nrow(raw) == 0L) return(integer(0))
    fitted <- raw[!is.na(raw$center_x), ]
    if (nrow(fitted) == 0L) return(integer(0))
    unique(as.integer(unlist(fitted[[col]], use.names = FALSE)))
  })
}


# =============================================================================
# SECTION 6 — High-level orchestrators
# =============================================================================

#' Full RANSAC ellipse pipeline for a single LAS object.
#'
#' Runs: fit → optional window filter → [clean] → extract IDs.
#'
#' @param las          lidR LAS object with $TreeID and $point_id columns.
#' @param tree_ids     Integer vector of TreeIDs. NULL = all non-NA, non-zero.
#' @param fit_params   List from ransac_ellipse_fit_params().
#' @param clean_params List from ransac_ellipse_clean_params().
#' @param window_bot   Lower Z bound for sub-window. NULL = no windowing.
#' @param window_top   Upper Z bound for sub-window. NULL = no windowing.
#' @param clean        Logical. TRUE runs spike removal. FALSE skips cleaning.
#' @param seed         Integer RNG seed.
#' @param label        Character label for progress and console output.
#' @param verbose      Print timing and count summaries.
#'
#' @return Named list:
#'   $raw, $windowed_raw, $clean, $valid, $inlier_ids, $appeared_ids
ransac_ellipse_run <- function(las,
                               tree_ids     = NULL,
                               fit_params   = ransac_ellipse_fit_params(),
                               clean_params = ransac_ellipse_clean_params(),
                               window_bot   = NULL,
                               window_top   = NULL,
                               clean        = TRUE,
                               seed         = 42L,
                               label        = "las",
                               verbose      = TRUE,
                               # Per-tree progress during the fit. Defaults to
                               # `verbose` so existing callers are unaffected;
                               # pass FALSE to silence the bar while keeping the
                               # summary lines.
                               show_progress = verbose) {

  # ── Fit ────────────────────────────────────────────────────────────────────
  if (verbose) {
    bot_lbl <- if (is.null(fit_params$bot)) "adaptive (tree min Z)" else sprintf("%.2f", fit_params$bot)
    top_lbl <- if (is.null(fit_params$top)) "adaptive (tree max Z)" else sprintf("%.2f", fit_params$top)
    cat(sprintf(
      "\n[%s] Fitting RANSAC ellipses (%s\u2013%s m, step=%.3f m) ...\n",
      label, bot_lbl, top_lbl, fit_params$step
    ))
  }
  t0 <- proc.time()["elapsed"]
  # NOTE: seeding is now handled per-tree INSIDE ransac_ellipse_fit_one()
  # (keyed on seed + TreeID), so a single set.seed(seed) call here is no
  # longer used for the fit itself. `seed` is passed straight through so
  # every tree's result is reproducible independent of which other trees
  # are in `las` / what order tree_ids is processed in — this is what makes
  # automatic (all-trees) and manual/CSV (subset) runs agree on shared trees.
  raw <- ransac_ellipse_fit_one(las, tree_ids, fit_params, label,
                                show_progress = show_progress, seed = seed)
  if (verbose) cat(sprintf(
    "  [%s] Fit done: %d trees in %.1f min\n",
    label, length(raw),
    round((proc.time()["elapsed"] - t0) / 60, 1)
  ))

  # ── Optional window filter ──────────────────────────────────────────────────
  apply_window <- !is.null(window_bot) && !is.null(window_top)
  windowed_raw <- if (apply_window) {
    if (verbose) cat(sprintf(
      "  [%s] Filtering to Z window [%.2f, %.2f] m ...\n",
      label, window_bot, window_top
    ))
    w <- ransac_ellipse_filter_window(raw, window_bot, window_top)
    if (verbose) cat(sprintf(
      "  [%s] Ellipses in window: %d\n", label, sum(sapply(w, nrow))
    ))
    w
  } else {
    raw
  }

  # ── Optional cleaning ───────────────────────────────────────────────────────
  # ransac_ellipse_clean() uses a flagging-first approach: it returns $flagged
  # (all rows with spike_* columns) and $valid (rows where spike == FALSE).
  # $valid is used directly downstream; $flagged is available for diagnostics.
  if (clean) {
    if (verbose) cat(sprintf("  [%s] Cleaning ellipses ...\n", label))
    clean_result <- ransac_ellipse_clean(windowed_raw, clean_params)
    valid        <- lapply(clean_result, `[[`, "valid")
    if (verbose) cat(sprintf(
      "  [%s] Valid ellipses after cleaning: %d\n",
      label, sum(sapply(valid, nrow))
    ))
  } else {
    if (verbose) cat(sprintf("  [%s] Cleaning skipped (clean = FALSE).\n", label))
    clean_result <- NULL
    valid        <- lapply(windowed_raw, function(df) {
      df[!is.na(df$center_x) & !is.na(df$semi_a) & !is.na(df$semi_b), , drop = FALSE]
    })
    if (verbose) cat(sprintf(
      "  [%s] Ellipses passed through (no cleaning): %d\n",
      label, sum(sapply(valid, nrow))
    ))
  }

  # ── Extract point IDs (always from the full-range raw list) ────────────────
  inlier_ids   <- ransac_ellipse_extract_ids(raw, "inliers")
  appeared_ids <- ransac_ellipse_extract_ids(raw, "appeared")
  if (verbose) cat(sprintf(
    "  [%s] Trees with RANSAC inlier IDs: %d\n",
    label, sum(sapply(inlier_ids, length) > 0L)
  ))

  list(
    raw          = raw,
    windowed_raw = windowed_raw,
    clean        = clean_result,
    valid        = valid,
    inlier_ids   = inlier_ids,
    appeared_ids = appeared_ids
  )
}


#' Full RANSAC ellipse pipeline for multiple LAS objects.
#'
#' @param las_list     Named list of lidR LAS objects.
#' @param fit_params   Shared fit parameters from ransac_ellipse_fit_params().
#' @param clean_params Shared cleaning parameters from ransac_ellipse_clean_params().
#' @param window_bot   Shared Z sub-window lower bound. NULL = no window.
#' @param window_top   Shared Z sub-window upper bound. NULL = no window.
#' @param seed         Shared RNG seed.
#' @param verbose      Print per-cloud progress.
#'
#' @return Named list matching the names of las_list; each element is the
#'         complete result of ransac_ellipse_run() for that cloud.
ransac_ellipse_run_multi <- function(las_list,
                                     fit_params   = ransac_ellipse_fit_params(),
                                     clean_params = ransac_ellipse_clean_params(),
                                     window_bot   = NULL,
                                     window_top   = NULL,
                                     seed         = 42L,
                                     verbose      = TRUE,
                                     show_progress = verbose) {

  if (is.null(names(las_list)) || any(nchar(names(las_list)) == 0L))
    stop(paste(
      "las_list must be a fully named list.",
      "Example: list(las_1 = my_las_1, las_2 = my_las_2)"
    ))

  lapply(
    setNames(names(las_list), names(las_list)),
    function(nm) {
      ransac_ellipse_run(
        las          = las_list[[nm]],
        tree_ids     = NULL,
        fit_params   = fit_params,
        clean_params = clean_params,
        window_bot   = window_bot,
        window_top   = window_top,
        seed         = seed,
        label        = nm,
        verbose      = verbose,
        show_progress = show_progress
      )
    }
  )
}


# =============================================================================
# SECTION 7 — Diagnostic helpers
# =============================================================================

#' Print a per-tree ellipse cleaning summary table to the console.
#'
#' Reports the two removal stages (XY spike, radius spike) that correspond to
#' the flagging-first approach in ransac_ellipse_clean().
#'
#' @param clean_result Named list as returned by ransac_ellipse_clean(), or the
#'                     $clean element of a ransac_ellipse_run() result.
#' @param label        Character label printed in the section header.
#'
#' @return The summary data.frame, invisibly.
ransac_ellipse_print_clean_summary <- function(clean_result, label = "") {
  hdr <- if (nchar(label) > 0L) sprintf("[%s] ", label) else ""
  cat(sprintf("\n%sPer-tree ellipse cleaning summary:\n", hdr))
  df <- do.call(rbind, lapply(names(clean_result), function(tid) {
    r <- clean_result[[tid]]
    data.frame(
      TreeID          = as.integer(tid),
      ellipses_before = r$ellipses_before,
      n_removed_xy    = r$n_removed_xy,
      n_removed_radii = r$n_removed_radii,
      ellipses_after  = r$ellipses_after,
      prop_kept       = round(r$prop_kept, 3)
    )
  }))
  print(df, row.names = FALSE)
  invisible(df)
}


#' Summarise early-exit usage across all trees in a raw ellipse list.
#'
#' Counts, per tree, how many fitted slices (center_x not NA) converged via
#' the early-exit rule versus running to the full iteration budget.  Only
#' meaningful when early_exit = TRUE was used during fitting; if the
#' early_exit_hit column is absent the function returns NULL invisibly.
#'
#' @param raw_ellipses Named list as returned by ransac_ellipse_fit_one() or
#'                     the $raw element of ransac_ellipse_run().
#' @param label        Character label printed in the section header.
#'
#' @return A data.frame with columns:
#'   tree_id, slices_fitted, slices_early_exit, slices_full_iter,
#'   pct_early_exit — invisibly.  Also prints a summary to the console.
ransac_ellipse_summarise_early_exit <- function(raw_ellipses, label = "") {

  hdr <- if (nchar(label) > 0L) sprintf("[%s] ", label) else ""

  # Check the column exists in at least one tree
  has_col <- any(sapply(raw_ellipses, function(dt)
    !is.null(dt) && "early_exit_hit" %in% names(dt)))

  if (!has_col) {
    cat(sprintf(
      "%sNo early_exit_hit column found — was early_exit = TRUE used?\n", hdr))
    return(invisible(NULL))
  }

  df <- do.call(rbind, lapply(names(raw_ellipses), function(tid) {
    dt      <- raw_ellipses[[tid]]
    if (is.null(dt) || nrow(dt) == 0L)
      return(data.frame(tree_id           = as.integer(tid),
                        slices_fitted     = 0L,
                        slices_early_exit = 0L,
                        slices_full_iter  = 0L,
                        pct_early_exit    = NA_real_))
    fitted  <- dt[!is.na(dt$center_x), ]
    n_fit   <- nrow(fitted)
    n_early <- if (n_fit > 0L) sum(fitted$early_exit_hit, na.rm = TRUE) else 0L
    data.frame(
      tree_id           = as.integer(tid),
      slices_fitted     = n_fit,
      slices_early_exit = n_early,
      slices_full_iter  = n_fit - n_early,
      pct_early_exit    = if (n_fit > 0L) round(100 * n_early / n_fit, 1) else NA_real_
    )
  }))

  # Cloud-level totals
  tot_fit   <- sum(df$slices_fitted)
  tot_early <- sum(df$slices_early_exit)
  tot_full  <- sum(df$slices_full_iter)
  pct_early <- if (tot_fit > 0L) round(100 * tot_early / tot_fit, 1) else NA_real_

  cat(sprintf(
    "\n%sEarly-exit summary (%d trees, %d fitted slices total):\n",
    hdr, nrow(df), tot_fit
  ))
  cat(sprintf(
    "  Early-exit : %d slices (%.1f%%)\n", tot_early, pct_early
  ))
  cat(sprintf(
    "  Full iters : %d slices (%.1f%%)\n", tot_full, 100 - pct_early
  ))
  cat(sprintf("  Per-tree breakdown:\n"))
  print(df, row.names = FALSE)

  invisible(df)
}
