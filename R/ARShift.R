#' Tests whether two types of points co-cluster
#'
#' Function to implement Area-adjusted Random Shift (ARShift), a non-paramteric test to assess co-clustering of two point patterns. The test is based on the bivariate K-function and uses random shifts of one of the point patterns to generate a null distribution. The test is implemented in a parallelized manner for efficiency.
#' As with bivariate Ripley's K-function, ARShift is directional. ARShift tests whether taxa/cell-type j co-clusters with taxa/cell-type i. This may be different from the co-clustering pattern of taxa/cell-type i with taxa/cell-type j. Therefore, it is recommended to run ARShift based on the scientific hypothesis of interest.
#'
#' @param data The observed point process pattern on which the test will be applied. It must be a multitype point pattern where marks are factor.
#' @param i The type of points(mark) in data which stays fixed. Should be input as a character string.
#' @param j The type of points(mark) in data which is shifted. Should be input as a character string. ARShift tests whether taxa/cell-type j co-clusters with taxa/cell-type i.
#' @param r Numeric vector at which the K-function is evaluated. There is a sensible default.
#' @param rmax Single numeric value indicating the maximum value of r to be considered. There is a sensible default.
#' @param correction A character vector specifying the type of border correction to be used in the K function.
#' Available options are "isotropic", "translate", "translation" and "border". For very large datasets "border" is recommended.
#' @param n_perm Number of random shifts
#' @param k factor by which number of permutations is increased. Defaults to 1.2. For sparse processes k can be increased.
#' @param type The type of K-function to compute, homogneous or inhomogeneous. Can be "hom" or "inhom". Defaults to homogeneous K.
#' @param sigma_val Bandwidth for evaluating inhomogeneous intensity in the inhomogeneous K. has a sensible default.
#' @param lite logical, indicating whether the point process should be randomly thinned for efficiency. Recommended only for very large datasets. Deafults to \code{FALSE}.
#' @param thin_prob numeric scalar between 0 and 1, indicating the thinning probability. Defaults to 0.5. Ignored if \code{lite = FALSE}.
#' @param jump_radius Radius of the uniform disc from which the random shift vector is sampled. There is a sensible default.
#' @param verbose Logical, indicating whether to print progress messages. Defaults to \code{FALSE}.
#' @param workers Number of cores to be used in case of parallel computation. Defaults to 1, serial computation.
#'
#' @return A list containing:
#' \describe{
#'   \item{p}{P-value for the area adjusted random shift test.}
#'   \item{K}{Observed bivariate K-function.}
#'   \item{r}{Radii at which the K-function was evaluated.}
#'   \item{n_attempts}{Total number of shift attempts made. This is equal to \code{ceiling(n_perm * k)}.}
#'   \item{n_valid}{Number of valid random shift attempts. If this is less than \code{n_perm}, increase \code{k} in the function input.}
#'   \item{n_perm}{Targeted number of random shifts.}
#' }
#' @examples
#' win <- spatstat.geom::square(1)
#' a <- spatstat.random::rThomas(8, 0.05, 20, win)
#' b <- spatstat.random::rThomas(7, 0.04, 20, win)
#' spatstat.geom::marks(a) <- factor(
#'   rep(1, spatstat.geom::npoints(a)),
#'   levels = c(1, 2)
#' )
#' spatstat.geom::marks(b) <- factor(
#'   rep(2, spatstat.geom::npoints(b)),
#'   levels = c(1, 2)
#' )
#' X <- spatstat.geom::superimpose(a, b, W = win)
#' spatstat.geom::marks(X) <- factor(
#'   spatstat.geom::marks(X),
#'   levels = c(1, 2)
#' )
#'
#' res <- ARShift(
#'   data = X,
#'   i = "1",
#'   j = "2",
#'   correction = "isotropic",
#'   n_perm = 499,
#'   workers = 5,
#'   type = "hom"
#' )
#' print(res)
#'
#' @author Asmita Roy
#' @export
ARShift <- function(data, i, j, r = NULL, rmax = NULL, correction = "isotropic",
                    n_perm = 999, k = 1.1, type = "hom",
                    sigma_val = NULL, lite = FALSE, thin_prob = 0.5,
                    jump_radius = NULL, workers = 1, verbose = FALSE) {

  log_progress <- function(...) {
    if (base::isTRUE(verbose)) {
      base::message(...)
    }
  }

  log_progress("ARShift: starting.")

  if (lite) {
    log_progress("ARShift: thinning input data.")
    data <- spatstat.random::rthin(data, thin_prob)
  }

  original_window <- spatstat.geom::Window(data)
  win_area <- spatstat.geom::area.owin(original_window)

  if (base::is.null(rmax) && base::is.null(r)) {
    rmax <- 0.3 * spatstat.geom::incircle(original_window)$r
    r <- base::seq(0, rmax, length.out = 50)
  }else if(base::is.null(rmax) && !base::is.null(r)) {
    rmax <- base::max(r)
  }else if(!base::is.null(rmax) && base::is.null(r)){
    r <- base::seq(0, rmax, length.out = 50)
  }else{
    #check whether max(r) conflicts with rmax
    if(base::max(r) > rmax){stop("ARShift: max(r) is greater than rmax. Please check your inputs.")}
  }


  base_chr <- base::as.character(i)
  shift_chr <- base::as.character(j)

  freq_marks <- base::table(spatstat.geom::marks(data))
  n_base <- base::as.numeric(freq_marks[base_chr])
  n_shift <- base::as.numeric(freq_marks[shift_chr])

  if (base::is.na(n_base) || base::is.na(n_shift) ||
      n_base == 0 || n_shift == 0) {
    log_progress("ARShift: one or both requested marks are absent; stopping.")
    return(base::list(
      p = NA_real_,
      K = NA_real_,
      r = r,
      n_attempts = 0,
      n_valid = 0,
      n_perm = n_perm

    ))
  }

  if (base::is.null(jump_radius)) {
    jump_radius <- spatstat.geom::incircle(original_window)$r
  }

  if (type == "inhom" && base::is.null(sigma_val)) {
    log_progress("ARShift: estimating intensity bandwidth.")
    sigma_val <- spatstat.explore::bw.CvL(data)
  }

  log_progress("ARShift: computing observed Kcross.")
  obs_Kcross <- compute_Kcross(
    X = data,
    base_taxa = base_chr,
    shift_taxa = shift_chr,
    r = r,
    type = type,
    correction = correction,
    sigma_val = sigma_val
  )

  points_shift_original <- base::subset(
    data,
    spatstat.geom::marks(data) == shift_chr
  )
  points_other_original <- base::subset(
    data,
    spatstat.geom::marks(data) != shift_chr
  )

  n_attempts <- base::ceiling(n_perm * k)

  log_progress(
    "ARShift: running ", n_attempts, " shift attempts using ",
    workers, if (workers == 1) " worker." else " workers."
  )

  attempt_results <- safe_pmclapply(
    X = base::seq_len(n_attempts),
    FUN = one_area_permutation_attempt,
    workers = workers,
    points_shift_original = points_shift_original,
    points_other_original = points_other_original,
    original_window = original_window,
    base_taxa = i,
    shift_taxa = j,
    r = r,
    type = type,
    correction = correction,
    sigma_val = sigma_val,
    jump_radius = jump_radius
  )

  valid_results <- base::Filter(
    function(z) base::is.list(z) && base::isTRUE(z$valid),
    attempt_results
  )
  n_valid_total <- base::length(valid_results)

  log_progress(
    "ARShift: ", n_valid_total, " of ", n_attempts,
    " shift attempts were valid."
  )

  if (n_valid_total < n_perm) {
    base::warning(
      "Only ", n_valid_total, " valid permutations obtained out of requested ",
      n_perm, ". Consider increasing k. Current k = ", k, ".",
      call. = FALSE
    )
  }

  n_keep <- base::min(n_perm, n_valid_total)

  if (n_keep == 0L) {
    return(base::list(
      p = NA_real_,
      K = obs_Kcross,
      r = r,
      n_attempts = n_attempts,
      n_valid = n_valid_total,
      n_perm = n_perm

      ))
  }

  valid_results <- valid_results[base::seq_len(n_keep)]
  Kcross_uncorr <- base::matrix(
    NA_real_,
    nrow = base::length(r),
    ncol = n_keep
  )
  area_shift <- base::numeric(n_keep)

  for (b in base::seq_len(n_keep)) {
    Kcross_uncorr[, b] <- valid_results[[b]]$K_uncorr
    area_shift[b] <- valid_results[[b]]$area_shift
  }

  log_progress("ARShift: computing T and envelope p-value.")

  meanK <- base::rowMeans(Kcross_uncorr, na.rm = TRUE)
  T_0 <- base::sqrt(win_area) * (obs_Kcross - meanK)
  T_sim <- Kcross_uncorr - meanK
  T_scaled <- base::sweep(T_sim, 2, base::sqrt(area_shift), "*")

  pval_Kcross_area <- compute_envelope_pval(
    r = r,
    obs = T_0,
    sim_m = T_scaled
  )

  log_progress("ARShift: finished.")

  base::list(
    p = pval_Kcross_area,
    K = obs_Kcross,
    r = r,
    n_attempts = n_attempts,
    n_valid = n_valid_total,
    n_perm = n_perm
  )
}
