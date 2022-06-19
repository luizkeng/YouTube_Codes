## Rotina adaptada do pacote chemometrics

pca_Var_explicada_by_PC <- function (X, a, center = TRUE, scale = TRUE, plot = TRUE, ...) 
{
  if (a < 1 || a > min(nrow(X) - 1, ncol(X))) 
    stop("Invalid number of components, a")
  X <- as.matrix(scale(X, center = center, scale = scale))
  if (ncol(X) > nrow(X)) {
    e <- eigen(X %*% t(X))
    T <- e$vectors %*% diag(sqrt(e$values))
    P <- t(X) %*% T %*% diag(1/e$values)
  }
  else {
    X_svd <- svd(X)
    P <- X_svd$v
    T <- X %*% P
  }
  varexpl = 1 - apply((X - T[, a] %*% t(P[, a]))^2, 2, 
                      sum)/apply(X^2, 2, sum)
  if (plot) {
    barplot(varexpl, ylab = "Explained variance", ylim = c(0, 
                                                           1), ...)
  }
  list(ExplVar = varexpl)
}

