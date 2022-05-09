#' Halton Sequences
#'
#' @param draw Number of draws for the sequence
#' @param dim Number of dimensions to draw
#' @param prime If \code{NULL}, the function will use prime numbers in ascending order. If specified, should be an atomic vector of length \code{dim} containing prime numbers to be used for each dimension.
#' @param scramble If \code{TRUE}, randomly permute the sequences
#' @param randomise If \code{TRUE}, randomise each draw
#' @param normal If \code{TRUE}, convert the sequence to a standard normal distribution
#' @param discard If provided, simulate an extra number of draws and discard the first few. For example, if 10, \code{draw} + 10 draws will be drawn, and the first 10 discarded.
#'
#' @return
#'
#' @references
#' \insertRef{bhat2003}{rngbox}
#'
#' \insertRef{train2003}{rngbox}
#'
#' \insertRef{okten2012}{rngbox}
#'
#' @export
halton <- function(draw, dim, prime=NULL, scramble=F, randomise=F, normal=F, discard=NULL) {

  if (is.null(prime)) {prime <- rngbox:::prime; pckgPrime <- T} else pckgPrime <- F
  if (length(prime) < dim) stop("Length of prime is smaller than number of dimensions")
  if (!pckgPrime & length(prime) > dim) warning(paste0(length(prime), " prime numbers have been included but only ", dim, " dimensions have been specified. Using the first ", dim, " prime numbers."))

  mat <- matrix(numeric(draw*dim), draw, dim)

  halton_r <- function(draw, seq, interval, t) {
    breaks <- eval(interval)
    seq <- c(seq, rowSums(expand.grid(seq, breaks)))
    if (length(seq) > draw) {
      return(seq[2:(draw+1)])
    } else {
      return(halton_r(draw, seq, interval, t+1))
    }
  }

  for (i in 1:dim) {
    mat[,i] <- halton_r(draw=ifelse(is.null(discard), draw, draw+discard), seq=0, interval=quote(1:(prime[i]-1)/(prime[i]^t)), t=1)
  }

  if (!is.null(discard)) {
    mat[,]
  }

  if (scramble) {
    for (i in 1:dim) {
      mat[,i] <- mat[sample(draw, replace=F),i]
    }
  }

  if (randomise) {
    for (i in 1:dim) {
      mu <- runif(1, 0, 1)
      mat[,i] <- mat[,i] + 1
      mat[,i] <- ifelse(mat[,i] > 1, mat[,i]-1, mat[,i])
    }
  }

  if (normal) {
    for (i in 1:dim) {
      mat[,i] <- stats::qnorm(mat[,i])
    }
  }

  return(mat)

}
