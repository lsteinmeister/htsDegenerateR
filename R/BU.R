#' Bottom-up reconciliation
#'
#' @param fcasts forecasts to be reconciled
#' @param S structure matrix representing the hierarchical structure of the hts
#'
#' @return reconciled forecasts
#' @export
BU <- function(fcasts, S) {
  #fcasts <- t(stats::na.omit(t(fcasts)))

  fcasts = if(is.matrix(fcasts)) fcasts else t(fcasts)

  bo_lvl_names = colnames(S)

  up_lvl = which(!rownames(S) %in% bo_lvl_names)

  G= t(S)

  G[,up_lvl] <- 0

  t(S %*% G %*% t(fcasts))
}
