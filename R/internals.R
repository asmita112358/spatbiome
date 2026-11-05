compute_Kcross <- function(X,
                           base_taxa,
                           shift_taxa,
                           r,
                           type = "inhom",
                           correction = "isotropic",
                           sigma_val = NULL) {

  i <- base::as.character(base_taxa)
  j <- base::as.character(shift_taxa)

  if (type == "inhom") {
    Kobj <- spatstat.explore::Kcross.inhom(
      X = X,
      i = i,
      j = j,
      r = r,
      correction = correction,
      sigma = sigma_val
    )
  } else {
    Kobj <- spatstat.explore::Kcross(
      X = X,
      i = i,
      j = j,
      r = r,
      correction = correction
    )
  }

  if (correction == "isotropic") {
    return(Kobj$iso)
  }

  if (correction %in% c("translate", "translation")) {
    return(Kobj$trans)
  }

  if (correction == "border") {
    return(Kobj$border)
  }

  base::stop("Unsupported correction: ", correction)
}


compute_envelope_pval<- function(r, obs, sim_m, direction = "greater") {

  if (is.null(sim_m)) {
    return(NA_real_)
  }

  if (!is.matrix(sim_m)) {
    sim_m <- as.matrix(sim_m)
  }

  keep <- is.finite(r) &
    is.finite(obs) &
    apply(sim_m, 1, function(z) all(is.finite(z)))

  if (!any(keep)) {
    return(NA_real_)
  }

  cset <- GET::create_curve_set(list(
    r = r[keep],
    obs = obs[keep],
    sim_m = sim_m[keep, , drop = FALSE]
  ))

  result <- GET::rank_envelope(
    cset,
    type = "erl",
    alternative = direction
  )

  attr(result, "p")
}
safe_pmclapply <- function(X, FUN, workers, ...) {
  workers <- as.integer(workers)

  if (length(workers) != 1L || is.na(workers) || workers < 1L) {
    stop("`workers` must be a single positive integer.")
  }

  if (workers == 1L) {
    return(lapply(X, FUN, ...))
  }

  if (.Platform$OS.type == "windows") {
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)

    return(tryCatch({
      future::plan(future::multisession, workers = workers)
      future.apply::future_lapply(
        X, FUN, ..., future.seed = TRUE
      )
    }, error = function(e) {
      message("future_lapply failed, falling back to lapply: ", e$message)
      lapply(X, FUN, ...)
    }))
  }

  if (.Platform$OS.type == "unix") {
    return(tryCatch({
      parallel::mclapply(X, FUN, ..., mc.cores = workers)
    }, error = function(e) {
      message("mclapply failed, falling back to lapply: ", e$message)
      lapply(X, FUN, ...)
    }))
  }

  message("Unsupported OS; falling back to lapply.")
  lapply(X, FUN, ...)
}

one_area_permutation_attempt <- function(attempt_id,
                                         points_shift_original,
                                         points_other_original,
                                         original_window,
                                         base_taxa,
                                         shift_taxa,
                                         r,
                                         type,
                                         correction,
                                         sigma_val,
                                         jump_radius) {
  base_chr <- as.character(base_taxa)
  shift_chr <- as.character(shift_taxa)
  out <- list(valid = FALSE)

  shift_vector <- spatstat.random::runifdisc(1, radius = jump_radius)
  shift_x <- shift_vector$x
  shift_y <- shift_vector$y

  points_shifted <- tryCatch(
    {
      spatstat.geom::shift(points_shift_original, vec = c(shift_x, shift_y))
    },
    error = function(e) NULL
  )
  if (is.null(points_shifted)) {
    return(out)
  }

  window_shifted <- tryCatch(
    {
      spatstat.geom::shift(original_window, vec = c(shift_x, shift_y))
    },
    error = function(e) NULL
  )
  if (is.null(window_shifted)) {
    return(out)
  }

  window_reduced <- tryCatch(
    {
      spatstat.geom::intersect.owin(original_window, window_shifted)
    },
    error = function(e) NULL
  )
  if (is.null(window_reduced)) {
    return(out)
  }

  reduced_area <- spatstat.geom::area.owin(window_reduced)
  if (!is.finite(reduced_area) || reduced_area <= 0) {
    return(out)
  }

  pp_combined <- tryCatch(
    {
      spatstat.geom::superimpose(points_other_original, points_shifted, check = FALSE)
    },
    error = function(e) NULL
  )
  if (is.null(pp_combined)) {
    return(out)
  }

  pp_reduced <- tryCatch(
    {
      pp_combined[window_reduced]
    },
    error = function(e) NULL
  )
  if (is.null(pp_reduced) || spatstat.geom::npoints(pp_reduced) == 0) {
    return(out)
  }

  freq_vc <- table(spatstat.geom::marks(pp_reduced))
  n_base_vc <- as.numeric(freq_vc[base_chr])
  n_shift_vc <- as.numeric(freq_vc[shift_chr])
  if (is.na(n_base_vc) || is.na(n_shift_vc) || n_base_vc == 0 || n_shift_vc == 0) {
    return(out)
  }

  K_uncorr <- tryCatch(
    {
      compute_Kcross(
        X = pp_reduced,
        base_taxa = base_taxa,
        shift_taxa = shift_taxa,
        r = r,
        type = type,
        correction = correction,
        sigma_val = sigma_val
      )
    },
    error = function(e) NULL
  )
  if (is.null(K_uncorr) || length(K_uncorr) != length(r) || any(!is.finite(K_uncorr))) {
    return(out)
  }

  list(
    valid = TRUE,
    K_uncorr = K_uncorr,
    area_shift = reduced_area,
    shift_vector = c(shift_x, shift_y)
  )
}

one_area_permutation_attempt <- function(attempt_id,
                                         points_shift_original,
                                         points_other_original,
                                         original_window,
                                         base_taxa,
                                         shift_taxa,
                                         r,
                                         type,
                                         correction,
                                         sigma_val,
                                         jump_radius) {
  base_chr <- base::as.character(base_taxa)
  shift_chr <- base::as.character(shift_taxa)
  out <- base::list(valid = FALSE)

  shift_vector <- spatstat.random::runifdisc(1, radius = jump_radius)
  shift_x <- shift_vector$x
  shift_y <- shift_vector$y

  points_shifted <- base::tryCatch(
    spatstat.geom::shift(
      points_shift_original,
      vec = c(shift_x, shift_y)
    ),
    error = function(e) NULL
  )
  if (base::is.null(points_shifted)) {
    return(out)
  }

  window_shifted <- base::tryCatch(
    spatstat.geom::shift(
      original_window,
      vec = c(shift_x, shift_y)
    ),
    error = function(e) NULL
  )
  if (base::is.null(window_shifted)) {
    return(out)
  }

  window_reduced <- base::tryCatch(
    spatstat.geom::intersect.owin(original_window, window_shifted),
    error = function(e) NULL
  )
  if (base::is.null(window_reduced)) {
    return(out)
  }

  reduced_area <- spatstat.geom::area.owin(window_reduced)
  if (!base::is.finite(reduced_area) || reduced_area <= 0) {
    return(out)
  }

  pp_combined <- base::tryCatch(
    spatstat.geom::superimpose(
      points_other_original,
      points_shifted,
      check = FALSE
    ),
    error = function(e) NULL
  )
  if (base::is.null(pp_combined)) {
    return(out)
  }

  pp_reduced <- base::tryCatch(
    pp_combined[window_reduced],
    error = function(e) NULL
  )
  if (base::is.null(pp_reduced) ||
      spatstat.geom::npoints(pp_reduced) == 0) {
    return(out)
  }

  freq_vc <- base::table(spatstat.geom::marks(pp_reduced))
  n_base_vc <- base::as.numeric(freq_vc[base_chr])
  n_shift_vc <- base::as.numeric(freq_vc[shift_chr])

  if (base::is.na(n_base_vc) || base::is.na(n_shift_vc) ||
      n_base_vc == 0 || n_shift_vc == 0) {
    return(out)
  }

  K_uncorr <- base::tryCatch(
    compute_Kcross(
      X = pp_reduced,
      base_taxa = base_taxa,
      shift_taxa = shift_taxa,
      r = r,
      type = type,
      correction = correction,
      sigma_val = sigma_val
    ),
    error = function(e) NULL
  )

  if (base::is.null(K_uncorr) ||
      base::length(K_uncorr) != base::length(r) ||
      base::any(!base::is.finite(K_uncorr))) {
    return(out)
  }

  base::list(
    valid = TRUE,
    K_uncorr = K_uncorr,
    area_shift = reduced_area,
    shift_vector = c(shift_x, shift_y)
  )
}
