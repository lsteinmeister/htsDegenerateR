#' Structural Scaling reconciliation
#'
#' @param fcasts forecasts to be reconciled
#' @param Smat structure matrix representing the hierarchical structure of the hts
#' @param weights use the default for structural scaling and a vector of the residual variances for variance scaling
#'
#' @return reconciled forecasts
#' @export
strucScaling = function(fcasts, Smat, weights = rowSums(Smat)){
  if(is.matrix(weights))
    weights = methods::as(weights, "matrix.csr") else
      weights = methods::as(weights, "matrix.diag.csr")
    Smat = methods::as(Smat, "matrix.csr")
    if(is.matrix(fcasts))
      fcasts <- t(fcasts)

    allf <- t(CHOL(fcasts = fcasts, S = Smat, weights = weights, allow.changes = FALSE))
    colnames(allf) = if(is.matrix(fcasts)) colnames(fcasts) else  names(fcasts)
    allf
}
