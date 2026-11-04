#' Function to implement ARShift
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
#' @param lite logical, indicating whether the point process should be randomly thinned for very large datasets. Recommended only for very large datasets. Deafults to \code{FALSE}.
#' @param thin_prob numeric scalar between 0 and 1, indicating the thinning probability. Defaults to 0.5.
#' @param jump_radius Radius of the uniform disc from which the random shift vector is sampled. There is a sensible default.
#' @param workers Number of cores to be used in case of parallel computation. Defaults to 1, serial computation.
#' @return Returns a list with the following:
#' \item{K_fun}{Bivariate K function}
#' \item{r}{numeric vector of radii on which the K-function was evaluated}
#' \item{n_valid}{Number of valid random shift attempts. A random shift attempt is valid if points of both type are present in the intersection window after shift. If n_valid is less than n_perm increasing k is recommended.}
#' \item{k}{ratio of number of attempted permutations and number of required permutations.}
#' @author Asmita Roy
#' @export
ARShift <- function(data, i, j, r, rmax = NULL, correction, n_perm = 999, k = 1.1,
                    type = "hom", sigma_val = NULL, lite = FALSE, thin_prob = 0.5, jump_radius = NULL,
                    workers = 1){
  ##Set parameters
  if (lite) {
    data <- rthin(data, thin_prob)
  }

  original_window <- Window(data)
  win_area <- area.owin(original_window)

  if (is.null(r)) {
    rmax <- 0.3 * incircle(original_window)$r
    r <- seq(0, rmax, length.out = 50)
  }

  n_dist <- length(r)
  base_chr <- as.character(i)
  shift_chr <- as.character(j)

  freq_marks <- table(marks(data))
  n_base <- as.numeric(freq_marks[base_chr])
  n_shift <- as.numeric(freq_marks[shift_chr])


  #Data check
  if (is.na(n_base) || is.na(n_shift) || n_base == 0 || n_shift == 0) {
    return(list(K_fun =  NA_real_, r = r, n_valid = 0, n_attempts = 0, k = k))
  }

  #More parameters

  if (is.null(jump_radius)) {
    jump_radius <- incircle(original_window)$r
  }

  sigma_val <- NULL
  if (type == "inhom") {
    if(is.null(sigma_val)){
      sigma_val <- bw.CvL(data)
    }
   }

  #Compute Kcross for original data
  obs_Kcross <- compute_Kcross(
    X = data,
    base_taxa = base_chr,
    shift_taxa = shift_chr,
    r = r,
    type = type,
    correction = correction,
    sigma_val = sigma_val
  )


}
