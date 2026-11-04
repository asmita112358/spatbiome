compute_Kcross <- function(X,
                                base_taxa,
                                shift_taxa,
                                r,
                                type = "inhom",
                                correction = "isotropic",
                                sigma_val = NULL) {

  i <- as.character(base_taxa)
  j <- as.character(shift_taxa)

  if (type == "inhom") {
    Kobj <- spatstat::Kcross.inhom(
      X = X,
      i = i,
      j = j,
      r = r,
      correction = correction,
      sigma = sigma_val
    )
  } else {
    Kobj <- spatstat::Kcross(
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

  stop("Unsupported correction: ", correction)
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
