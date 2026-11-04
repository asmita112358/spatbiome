#' Function to implement ARShift
#'
#' @param data The observed point process pattern on which the test will be applied. It must be a multitype point pattern where marks are factor.
#' @param i The type of points(mark) in data which stays fixed. Should be input as a character string.
#' @param j The type of points(mark) in data which is shifted. Should be input as a character string. ARShift tests whether taxa/cell-type j co-clusters with taxa/cell-type i.
#' @param r Numeric vector at which the K-function is evaluated. There is a sensible default.
#' @param correction A character vector specifying the type of border correction to be used in the K function.
#' Available options are "isotropic", "translate", "translation" and "border". For very large datasets "border" is recommended.
#' @param n_perm Number of random shifts
#' @param k factor by which number of permutations is increased. Defaults to 1.2. For sparse processes k can be increased.
#' @param type The type of K-function to compute, homogneous or inhomogeneous. Can be "hom" or "inhom". Defaults to homogeneous K.
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
ARShift <- function(data, i, j, r, correction, n_perm = 999, k = 1.1,
                    lite = FALSE, thin_prob = 0.5, jump_radius = NULL,
                    workers = 1){

}
