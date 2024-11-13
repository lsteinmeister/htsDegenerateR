#library("SparseM")
#library("Matrix")

#' @import Matrix
#' @import SparseM


# 3 algorithms for forecasting reconciliation through trace minimization
# Only used for BLUF
# Author: Shanika Wickramasuriya
# Paper: Forecasting hierarchical and grouped time series through trace
#        minimization
# All these functions return a reverse reconciled matrix with all ts.

#LU decomposition is fast but sometimes instable. Use QR decomposition if LU decomposition fails
solveLUQR <- function(lhs.l, rhs.l) {
  tryCatch(solve(lhs.l, rhs.l), error=function(cond){

    #browser()
    warning("An error in LU decomposition occurred, the message was the following:\n",
            cond$message, "\n Trying QR decomposition instead...")
    solve(qr(lhs.l), rhs.l)
  })
}

# LU factorization (Matrix pkg)
LU <- function(fcasts, S, weights, allow.changes = FALSE) {
  nts <- nrow(S)
  nbts <- ncol(S)
  nagg <- nts - nbts
  seqagg <- 1L:nagg

  if (!allow.changes) {
    utmat <- cbind2(sparseMatrix(i = seqagg, j = seqagg, x = 1),
                    -1 * S[1L:nagg, ])
  } else {
    # Identifying rows with one 1 element to make the Identity matrix in S
    indx <- rowSums(S)
    idx <- tail(which(indx == 1L), nbts)

    # Permulation vector to rearrange rows of S, rows/col of W and forecasts
    pvec <- c(setdiff(1:nts, idx) , idx)
    S2 <- S[pvec, ]
    weights <- weights[pvec, pvec]
    fcasts <- fcasts[pvec, ]
    utmat <- cbind2(sparseMatrix(i = seqagg, j = seqagg, x = 1),
                    -1 * S2[1L:nagg, ])
  }
  jmat <- sparseMatrix(i = 1L:nbts, j = (nagg + 1L):nts, x = rep(1L, nbts),
                       dims = c(nbts, nts))
  rhs.l <-  methods::as(utmat %*% fcasts, "CsparseMatrix")
  if (is.null(weights)) {
    lhs.l <- utmat %*% t(utmat)
    lhs.l <- (t(lhs.l) + lhs.l)/2
    lin.sol <- solveLUQR(lhs.l, rhs.l)
    p1 <- jmat %*% fcasts - (jmat %*% t(utmat) %*% lin.sol)
  } else {
    lhs.l <- utmat %*% weights %*% t(utmat)
    lhs.l <- (t(lhs.l) + lhs.l)/2
    lin.sol <- solveLUQR(lhs.l, rhs.l)
    p1 <- jmat %*% fcasts - (jmat %*% weights %*% t(utmat) %*% lin.sol)
  }
  if (!allow.changes) {
    comb <- as.matrix(S %*% p1)
  } else {
    comb <- numeric()
    comb[pvec] <- as.matrix(S2 %*% p1)
  }
  return(comb)
}

# Conjugate Gradient (Matrix and RcppEigen pkgs)
CG <- function(fcasts, S, weights, allow.changes = FALSE) {
  nts <- nrow(S)
  nbts <- ncol(S)
  nagg <- nts - nbts
  seqagg <- 1L:nagg
  if (!allow.changes) {
    utmat <- cbind2(Matrix::sparseMatrix(i = seqagg, j = seqagg, x = 1),
                    -1 * S[1L:nagg, ])
  } else {
    # Identifying rows with one 1 element to make the Identity matrix in S
    indx <- rowSums(S)
    idx <- tail(which(indx == 1L), nbts)

    # Permulation vector to rearrange rows of S, rows/col of W and forecasts
    pvec <- c(setdiff(1:nts, idx) , idx)
    S2 <- S[pvec, ]
    weights <- weights[pvec, pvec]
    fcasts <- fcasts[pvec, ]
    utmat <- cbind2(Matrix::sparseMatrix(i = seqagg, j = seqagg, x = 1),
                    -1 * S2[1L:nagg, ])
  }
  jmat <- Matrix::sparseMatrix(i = 1L:nbts, j = (nagg + 1L):nts, x = rep(1L, nbts),
                               dims = c(nbts, nts))
  rhs.l <- as.matrix(utmat %*% fcasts)
  if (is.null(weights)) {
    lhs.l <- utmat %*% t(utmat)
    lin.sol <- as.matrix(cgm_c(lhs.l, rhs.l)) # cgm_c is a C++ function
    p1 <- jmat %*% fcasts - (jmat %*% t(utmat) %*% lin.sol)
  } else {
    lhs.l <- utmat %*% weights %*% t(utmat)
    lin.sol <- as.matrix(cgm_c(lhs.l, rhs.l))
    p1 <- jmat %*% fcasts - (jmat %*% weights %*% t(utmat) %*% lin.sol)
  }
  if (!allow.changes) {
    comb <- as.matrix(S %*% p1)
  } else {
    comb <- numeric()
    comb[pvec] <- as.matrix(S2 %*% p1)
  }
  return(comb)
}
# Cholesky factorization
CHOL <- function(fcasts, S, weights, allow.changes = FALSE) {
  fcasts <- t(stats::na.omit(t(fcasts)))
  nts <- nrow(S)
  nbts <- ncol(S)
  nagg <- nts - nbts
  seqagg <- 1L:nagg
  h = if(is.matrix(fcasts)) nrow(fcasts) else 1
  if (!allow.changes) {
    utmat <- cbind(methods::as(nagg, "matrix.diag.csr"), -1 * S[1L:nagg, ])
  } else {
    # Identifying rows with one 1 element to make the Identity matrix in S
    Sm <- as(S, "dgCMatrix")
    indx <- rowSums(Sm)
    idx <- tail(which(indx == 1L), nbts) # idx are the indeces of the bts

    # Permulation vector to rearrange rows of S, rows/col of W and forecasts
    pvec <- c(setdiff(1:nts, idx) , idx)
    S2 <- S[pvec, ]
    weights <- weights[pvec, pvec]
    fcasts <- fcasts[pvec, ]
    utmat <- cbind(methods::as(nagg, "matrix.diag.csr"), -1 * S2[1L:nagg, ])
  }
  jmat <- methods::new("matrix.csr", ra = rep(1L, nbts), ja = seq((nagg + 1L), nts),
                       ia = 1L:(nbts + 1L), dimension = as.integer(c(nbts, nts)))
  rhs.l <- utmat %*% fcasts
  if (is.null(weights)) {
    lhs.l <- utmat %*% t(utmat)
    lhs.l <- (t(lhs.l) + lhs.l)/2
    lin.sol <- backsolve(chol(lhs.l), rhs.l)
    p1 <- jmat %*% fcasts - (jmat %*% t(utmat) %*% lin.sol)
  } else {
    lhs.l <- utmat %*% weights %*% t(utmat)
    lhs.l <- (t(lhs.l) + lhs.l)/2
    lin.sol <- backsolve(chol(lhs.l), rhs.l) # same as solve as long as default twice = T
    #p1 <- jmat %*% fcasts - kronecker(jmat %*% weights %*% t(utmat) , t(lin.sol))#jmat %*% fcasts - (jmat %*% weights %*% t(utmat) %*% lin.sol)
    p1 <-if(nagg == 1 & h > 1) jmat %*% fcasts - (jmat %*% weights %*% kronecker(t(utmat), t(lin.sol))) else
      jmat %*% fcasts - (jmat %*% weights %*% t(utmat) %*% lin.sol)
  }
  if (!allow.changes) {
    comb <- as.matrix(S %*% p1)
  } else {
    comb <- numeric()
    comb[pvec] <- as.matrix(S2 %*% p1)
  }
  return(comb)
}



is.posdef <- function (x, tol = 1e-08) {
  n <- NROW(x)
  if(n != NCOL(x))
    stop("x is not a square matrix")
  if(sum(c(abs(x - t(x)))) > 1e-8)
    stop("x is not a symmetric matrix")
  eigenvalues <- eigen(x, only.values = TRUE)$values
  eigenvalues[abs(eigenvalues) < tol] <- 0
  all(eigenvalues >= 0)
}

# Shrunk covariance matrix - Schafer and strimmer approach
# adapted to work with NA entries (different TS lengths)
shrink.estim <- function(x, tar, cov.type = "complete.obs", substract.mean = F)
{
  if (is.matrix(x) == TRUE && is.numeric(x) == FALSE)
    stop("The data matrix must be numeric!", call. = FALSE)
  p <- ncol(x)
  n <- nrow(x)
  if(cov.type =="complete.obs" & !substract.mean) x = na.omit(x)
  covm <-if(substract.mean)
    cov(x, use = cov.type)
  else crossprod(x) / n

  corm <- cov2cor(covm)
  xs <- scale(x, center = FALSE, scale = sqrt(diag(covm)))
  v <- if(substract.mean) (1/(n)) * (cov(xs^2, use = cov.type) -
                                       1/n * (cov(xs, use = cov.type))^2) else
                                         (1/(n * (n - 1))) * (crossprod(xs^2) - 1/n * (crossprod(xs))^2)

  diag(v) <- 0
  corapn <- cov2cor(tar)
  d <- (corm - corapn)^2
  lambda <- sum(v)/sum(d)
  lambda <- max(min(lambda, 1), 0)
  shrink.cov <- lambda * tar + (1 - lambda) * covm
  return(list(shrink.cov, c("The shrinkage intensity lambda is:",
                            round(lambda, digits = 4))))
}


lowerD <- function(x, cov.type)
{
  n <- nrow(x)

  return(diag(apply(x, 2, function(z){
    crossprod(na.omit(z))/sum(!is.na(z))} )))
}
#' accuracy.gts
#'
#' @param fcasts forecasts to be evaluated
#' @param actuals actuals to compare the forecasts against
#'
#' @return Averaged error measures across all time series in matrix form.
#' @export
accuracy.gts <- function(fcasts, actuals) {

  x <- actuals
  tspf <- tsp(fcasts)
  tspx <- tsp(x)
  start <- max(tspf[1], tspx[1])
  end <- min(tspf[2], tspx[2])
  start <- min(start, end)
  end <- max(start, end)
  res <- as.matrix(x) - as.matrix(fcasts)


  scale <- colMeans(abs(diff(x, lag = max(1, round(stats::frequency(x))))),
                    na.rm = TRUE)
  q <- sweep(res, 2, scale, "/")

  #mean absolute scales error
  mase <- colMeans(abs(q), na.rm = TRUE)

  pe <- res/x * 100  # percentage error

  # mean error
  me <- colMeans(res, na.rm = TRUE)
  # root mean squared error
  rmse <- sqrt(colMeans(res^2, na.rm = TRUE))
  # mean absolute error
  mae <- colMeans(abs(res), na.rm = TRUE)
  # mean absolute percentage error
  mape <- colMeans(abs(pe), na.rm = TRUE)
  # mean percentage error
  mpe <- colMeans(pe, na.rm = TRUE)

  out <- rbind(me, rmse, mae, mape, mpe)
  rownames(out) <- c("ME", "RMSE", "MAE", "MAPE", "MPE")
  if (exists("mase")) {
    out <- rbind(out, mase)
    rownames(out)[6L] <- "MASE"
  }
  if (exists("fcasts")) {
    colnames(out) <- colnames(fcasts)
  }
  return(out)
}
