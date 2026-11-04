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

vc_test <- function(K_vc_matrix,
                    obs_K,
                    shift_vectors,
                    bw_shift,
                    r) {

  n_perm <- ncol(K_vc_matrix)

  K_full <- cbind(K_vc_matrix, obs_K)
  K_mean <- rowMeans(K_full, na.rm = TRUE)
  K_centered <- K_full - K_mean

  var_shift <- t(nadaraya_watson(
    X = t(shift_vectors),
    Y = t(K_centered^2),
    X_pred = t(shift_vectors),
    bw = bw_shift
  ))

  var_shift[var_shift <= 1e-10] <- 1e-10

  K_std_shift <- K_centered / sqrt(var_shift)
  K_std_shift[!is.finite(K_std_shift)] <- 0

  CS_shift <- GET::create_curve_set(list(
    r = r,
    obs = K_std_shift[, n_perm + 1],
    sim_m = K_std_shift[, seq_len(n_perm), drop = FALSE]
  ))

  rank_env_obj <- GET::rank_envelope(
    CS_shift,
    type = "erl",
    alternative = "greater"
  )

  list(
    pval_vc_shift = attr(rank_env_obj, "p")
  )
}


nadaraya_watson <- function(X, Y, X_pred = X, bw) {

  n <- nrow(X)
  m <- nrow(X_pred)

  if (is.vector(Y)) {
    Y <- matrix(Y, ncol = 1)
  }

  n_response <- ncol(Y)
  predictions <- matrix(0, nrow = m, ncol = n_response)

  gaussian_kernel <- function(dist_sq, bandwidth) {
    exp(-0.5 * dist_sq / bandwidth^2)
  }

  for (i in seq_len(m)) {
    dist_sq <- rowSums(
      (
        X -
          matrix(
            X_pred[i, ],
            nrow = n,
            ncol = ncol(X),
            byrow = TRUE
          )
      )^2
    )

    weights <- gaussian_kernel(dist_sq, bw)
    sum_weights <- sum(weights)

    if (sum_weights > 1e-10) {
      predictions[i, ] <- as.numeric(crossprod(weights, Y) / sum_weights)
    } else {
      predictions[i, ] <- NA_real_
    }
  }

  if (n_response == 1) {
    predictions <- as.vector(predictions)
  }

  predictions
}

select_bandwidth <- function(X, method = "silverman") {

  n_sim <- ncol(X)
  d <- nrow(X)

  sd_X <- apply(X, 1, sd)
  sd_X[sd_X < 1e-10] <- 1e-10

  sigma_X <- exp(mean(log(sd_X)))

  if (method == "silverman") {
    bw <- sigma_X * (n_sim * (d + 2) / 4)^(-1 / (d + 4))
  } else if (method == "scott") {
    bw <- sigma_X * n_sim^(-1 / (d + 4))
  } else {
    bw <- as.numeric(method)
  }

  bw <- max(bw, 1e-10)

  bw
}
