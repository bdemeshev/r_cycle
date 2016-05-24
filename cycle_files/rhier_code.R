# rhierMnlRwMixture from bayesm code


Data = list(p = 5, lgtdata = choice_data)
Prior = list(ncomp = 1) 
Mcmc = mcmc_pars

# from reference manual:
BayesmConstant.A <- 0.01 # A and a?????? to match ordinary function call
BayesmConstant.nuInc <- 3 
BayesmConstant.RRScaling <- 2.93
BayesmConstant.w <- 0.1
BayesmConstant.a <- 5 # A and a??????
BayesmConstant.keep <- 1
BayesmConstant.nprint <- 100 
  
pandterm <- stop

function (Data, Prior, Mcmc) 
{
  if (missing(Data)) {
    pandterm("Requires Data argument -- list of p,lgtdata, and (possibly) Z")
  }
  if (is.null(Data$p)) {
    pandterm("Requires Data element p (# chce alternatives)")
  }
  p = Data$p
  if (is.null(Data$lgtdata)) {
    pandterm("Requires Data element lgtdata (list of data for each unit)")
  }
  lgtdata = Data$lgtdata
  nlgt = length(lgtdata)
  drawdelta = TRUE
  if (is.null(Data$Z)) {
    cat("Z not specified", fill = TRUE)
    fsh()
    drawdelta = FALSE
  } else {
    if (nrow(Data$Z) != nlgt) {
      pandterm(paste("Nrow(Z) ", nrow(Z), "ne number logits ", 
                     nlgt))
    } else {
      Z = Data$Z
    }
  }
  if (drawdelta) {
    nz = ncol(Z)
    colmeans = apply(Z, 2, mean)
    if (sum(colmeans) > 1e-05) {
      pandterm(paste("Z does not appear to be de-meaned: colmeans= ", 
                     colmeans))
    }
  }
  ypooled = NULL
  Xpooled = NULL
  if (!is.null(lgtdata[[1]]$X)) {
    oldncol = ncol(lgtdata[[1]]$X)
  }
  for (i in 1:nlgt) {
    if (is.null(lgtdata[[i]]$y)) {
      pandterm(paste("Requires element y of lgtdata[[", 
                     i, "]]"))
    }
    if (is.null(lgtdata[[i]]$X)) {
      pandterm(paste("Requires element X of lgtdata[[", 
                     i, "]]"))
    }
    ypooled = c(ypooled, lgtdata[[i]]$y)
    nrowX = nrow(lgtdata[[i]]$X)
    if ((nrowX/p) != length(lgtdata[[i]]$y)) {
      pandterm(paste("nrow(X) ne p*length(yi); exception at unit", 
                     i))
    }
    newncol = ncol(lgtdata[[i]]$X)
    if (newncol != oldncol) {
      pandterm(paste("All X elements must have same # of cols; exception at unit", 
                     i))
    }
    Xpooled = rbind(Xpooled, lgtdata[[i]]$X)
    oldncol = newncol
  }
  nvar = ncol(Xpooled)
  levely = as.numeric(levels(as.factor(ypooled)))
  if (length(levely) != p) {
    pandterm(paste("y takes on ", length(levely), " values -- must be = p"))
  }
  bady = FALSE
  for (i in 1:p) {
    if (levely[i] != i) 
      bady = TRUE
  }
  cat("Table of Y values pooled over all units", fill = TRUE)
  print(table(ypooled))
  if (bady) {
    pandterm("Invalid Y")
  }
  if (missing(Prior)) {
    pandterm("Requires Prior list argument (at least ncomp)")
  }
  if (is.null(Prior$ncomp)) {
    pandterm("Requires Prior element ncomp (num of mixture components)")
  } else {
    ncomp = Prior$ncomp
  }
  if (is.null(Prior$mubar)) {
    mubar = matrix(rep(0, nvar), nrow = 1)
  } else {
    mubar = matrix(Prior$mubar, nrow = 1)
  }
  if (ncol(mubar) != nvar) {
    pandterm(paste("mubar must have ncomp cols, ncol(mubar)= ", 
                   ncol(mubar)))
  }
  if (is.null(Prior$Amu)) {
    Amu = matrix(BayesmConstant.A, ncol = 1)
  } else {
    Amu = matrix(Prior$Amu, ncol = 1)
  }
  if (ncol(Amu) != 1 | nrow(Amu) != 1) {
    pandterm("Am must be a 1 x 1 array")
  }
  if (is.null(Prior$nu)) {
    nu = nvar + BayesmConstant.nuInc
  } else {
    nu = Prior$nu
  }
  if (nu < 1) {
    pandterm("invalid nu value")
  }
  if (is.null(Prior$V)) {
    V = nu * diag(nvar)
  } else {
    V = Prior$V
  }
  if (sum(dim(V) == c(nvar, nvar)) != 2) 
    pandterm("Invalid V in prior")
  if (is.null(Prior$Ad) & drawdelta) {
    Ad = BayesmConstant.A * diag(nvar * nz)
  } else {
    Ad = Prior$Ad
  }
  if (drawdelta) {
    if (ncol(Ad) != nvar * nz | nrow(Ad) != nvar * nz) {
      pandterm("Ad must be nvar*nz x nvar*nz")
    }
  }
  if (is.null(Prior$deltabar) & drawdelta) {
    deltabar = rep(0, nz * nvar)
  } else {
    deltabar = Prior$deltabar
  }
  if (drawdelta) {
    if (length(deltabar) != nz * nvar) {
      pandterm("deltabar must be of length nvar*nz")
    }
  }
  if (is.null(Prior$a)) {
    a = rep(BayesmConstant.a, ncomp)
  } else {
    a = Prior$a
  }
  if (length(a) != ncomp) {
    pandterm("Requires dim(a)= ncomp (no of components)")
  }
  bada = FALSE
  for (i in 1:ncomp) {
    if (a[i] < 0) 
      bada = TRUE
  }
  if (bada) 
    pandterm("invalid values in a vector")
  if (missing(Mcmc)) {
    pandterm("Requires Mcmc list argument")
  } else {
    if (is.null(Mcmc$s)) {
      s = BayesmConstant.RRScaling/sqrt(nvar)
    } else {
      s = Mcmc$s
    }
    if (is.null(Mcmc$w)) {
      w = BayesmConstant.w
    } else {
      w = Mcmc$w
    }
    if (is.null(Mcmc$keep)) {
      keep = BayesmConstant.keep
    } else {
      keep = Mcmc$keep
    }
    if (is.null(Mcmc$R)) {
      pandterm("Requires R argument in Mcmc list")
    } else {
      R = Mcmc$R
    }
    if (is.null(Mcmc$nprint)) {
      nprint = BayesmConstant.nprint
    } else {
      nprint = Mcmc$nprint
    }
    if (nprint < 0) {
      pandterm("nprint must be an integer greater than or equal to 0")
    }
  }
  cat(" ", fill = TRUE)
  cat("Starting MCMC Inference for Hierarchical Logit:", fill = TRUE)
  cat("   Normal Mixture with", ncomp, "components for first stage prior", 
      fill = TRUE)
  cat(paste("  ", p, " alternatives; ", nvar, " variables in X"), 
      fill = TRUE)
  cat(paste("   for ", nlgt, " cross-sectional units"), fill = TRUE)
  cat(" ", fill = TRUE)
  cat("Prior Parms: ", fill = TRUE)
  cat("nu =", nu, fill = TRUE)
  cat("V ", fill = TRUE)
  print(V)
  cat("mubar ", fill = TRUE)
  print(mubar)
  cat("Amu ", fill = TRUE)
  print(Amu)
  cat("a ", fill = TRUE)
  print(a)
  if (drawdelta) {
    cat("deltabar", fill = TRUE)
    print(deltabar)
    cat("Ad", fill = TRUE)
    print(Ad)
  }
  cat(" ", fill = TRUE)
  cat("MCMC Parms: ", fill = TRUE)
  cat(paste("s=", round(s, 3), " w= ", w, " R= ", R, " keep= ", 
            keep, " nprint= ", nprint), fill = TRUE)
  cat("", fill = TRUE)
  oldbetas = matrix(double(nlgt * nvar), ncol = nvar)
  llmnlFract = function(beta, y, X, betapooled, rootH, w, wgt) {
    z = as.vector(rootH %*% (beta - betapooled))
    return((1 - w) * llmnl(beta, y, X) + w * wgt * (-0.5 * 
                                                      (z %*% z)))
  }
  cat("initializing Metropolis candidate densities for ", nlgt, 
      " units ...", fill = TRUE)
  fsh()
  betainit = c(rep(0, nvar))
  out = optim(betainit, llmnl, method = "BFGS", control = list(fnscale = -1, 
                                                               trace = 0, reltol = 1e-06), X = Xpooled, y = ypooled)
  betapooled = out$par
  H = mnlHess(betapooled, ypooled, Xpooled)
  rootH = chol(H)
  for (i in 1:nlgt) {
    wgt = length(lgtdata[[i]]$y)/length(ypooled)
    out = optim(betapooled, llmnlFract, method = "BFGS", 
                control = list(fnscale = -1, trace = 0, reltol = 1e-04), 
                X = lgtdata[[i]]$X, y = lgtdata[[i]]$y, betapooled = betapooled, 
                rootH = rootH, w = w, wgt = wgt)
    if (out$convergence == 0) {
      hess = mnlHess(out$par, lgtdata[[i]]$y, lgtdata[[i]]$X)
      lgtdata[[i]] = c(lgtdata[[i]], list(converge = 1, 
                                          betafmle = out$par, hess = hess))
    } else {
      lgtdata[[i]] = c(lgtdata[[i]], list(converge = 0, 
                                          betafmle = c(rep(0, nvar)), hess = diag(nvar)))
    }
    oldbetas[i, ] = lgtdata[[i]]$betafmle
    if (i%%50 == 0) 
      cat("  completed unit #", i, fill = TRUE)
    fsh()
  }
  ind = NULL
  ninc = floor(nlgt/ncomp)
  for (i in 1:(ncomp - 1)) {
    ind = c(ind, rep(i, ninc))
  }
  if (ncomp != 1) {
    ind = c(ind, rep(ncomp, nlgt - length(ind)))
  } else {
    ind = rep(1, nlgt)
  }
  oldprob = rep(1/ncomp, ncomp)
  if (drawdelta) {
    olddelta = rep(0, nz * nvar)
  } else {
    olddelta = 0
    Z = matrix(0)
    deltabar = 0
    Ad = matrix(0)
  }
  draws = rhierMnlRwMixture_rcpp_loop(lgtdata, Z, deltabar, 
                                      Ad, mubar, Amu, nu, V, s, R, keep, nprint, drawdelta, 
                                      as.matrix(olddelta), a, oldprob, oldbetas, ind)
  if (drawdelta) {
    attributes(draws$Deltadraw)$class = c("bayesm.mat", "mcmc")
    attributes(draws$Deltadraw)$mcpar = c(1, R, keep)
  }
  attributes(draws$betadraw)$class = c("bayesm.hcoef")
  attributes(draws$nmix)$class = "bayesm.nmix"
  return(draws)
}