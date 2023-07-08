# 07/03/2023
# binsreg package, supporting functions
# p: degree of polynomial
# s: number of cts (deriv) constraints
# J: # of bins
# knot: complete knot sequence

genKnot.es <- function(x.min, x.max, J) {
  knot <- seq(x.min, x.max, length.out = J+1)
  return(knot)
}

# Quantile knot list (including xmin and xmax as boundaries)
genKnot.qs <- function(x, J) {
  knot <- quantile(x, seq(0, 1, 1/J), names = F, type=2)
  return(knot)
}

# grid generation
binsreg.grid <- function(knot, ngrid, addmore=F) {
  eval <- cumsum(c(knot[1], rep(diff(knot)/(ngrid+1), each=ngrid+1)))
  eval <- eval[-c(1, length(eval))]
  bin <- isknot <- mid <- NA
  if (addmore) {
     bin <- rep(1:(length(knot)-1), each=ngrid+1)
     bin <- bin[-length(bin)]
     isknot <- rep(c(rep(0, ngrid), 1), length(knot)-1)
     isknot <- isknot[-length(isknot)]
     id.mid <- rep(0, ngrid+1); id.mid[ceiling(ngrid/2)] <- 1
     mid <- rep(id.mid, length(knot)-1)
     mid <- mid[-length(mid)]
  }
  return(list(eval=eval, bin=bin, isknot=isknot, mid=mid))
}

# Generate Design
binsreg.spdes <- function(eval, p, s, knot, deriv) {
  if (s == 0) {
     # mimic STATA irecode
     pos <- findInterval(eval, knot, rightmost.closed = T, left.open = T)
     polyx <- matrix(0, length(eval), p+1)
     h <- diff(knot)
     eval.cen <- (eval-knot[-length(knot)][pos]) / h[pos]
     for (j in (deriv+1):(p+1)) {
       polyx[,j] <- eval.cen^(j-1-deriv) * factorial(j-1)/factorial(j-1-deriv) / h[pos]^deriv
     }
     P <- matrix(sapply(1:(length(knot)-1), function(i) polyx*(pos==i)), nrow = length(eval))
  } else {
    if (length(knot) >= 3) {
       ext.knot <- c(rep(knot[1], p+1), rep(knot[2:(length(knot)-1)], each = p-s+1), rep(knot[length(knot)], p+1))
    } else {
       ext.knot <- c(rep(knot[1], p+1), rep(knot[length(knot)], p+1))
    }
    P <- splineDesign(knots = ext.knot, eval, ord = p+1, derivs = deriv)
  }
  return(P)
}

# check drop, display warning
check.drop <- function(beta, k) {
  if (any(is.na(beta[1:k]))) {
    warning("some X-based variables dropped")
  }
}

# wrapper of vcov and vcovCL
binsreg.vcov <- function(model, type, cluster, is.qreg=FALSE, ...) {
  if (is.qreg) {
    V <- summary.rq(model, se=type, covariance = TRUE, cluster=cluster, ...)$cov
  } else {
    if (type=="const") {
      V <- vcov(model)
    } else {
      V <- vcovCL(model, type=type, cluster=cluster)
    }
  }
  return(V)
}

# internal pred function (model is long regression, NA handled inside)
binsreg.pred <- function(X, model, type="xb", vce, cluster=NULL, deriv=0, wvec=NULL, is.qreg=FALSE, avar=FALSE,...) {
   k <- ncol(X)
   fit <- NA
   if (type == "xb" | type == "all") {
      coef <- model$coeff
      coef[is.na(coef)] <- 0
      if (!is.null(wvec) & deriv==0) {
        fit <- as.vector(X %*% coef[1:k]) + sum(wvec*coef[-(1:k)])
      } else {
        fit <- as.vector(X %*% coef[1:k])
      }
   }
   se <- NA
   if (type == "se" | type == "all") {
      if (is.qreg) {
        vcv <- binsreg.vcov(model, type=vce, cluster=cluster, is.qreg=TRUE, ...)
      } else {
        vcv <- binsreg.vcov(model, type=vce, cluster=cluster)
      }
      if (avar) {
        pos <- !is.na(model$coeff[1:k])
        k.new <- sum(pos)
        vcv <- vcv[1:k.new, 1:k.new]
        se <- sqrt(rowSums((X[, pos, drop=F] %*% vcv) * X[, pos, drop=F]))
      } else {
        if (!is.null(wvec)) {
          if (deriv==0) X <- cbind(X, outer(rep(1,nrow(X)), wvec))
          else          X <- cbind(X, outer(rep(1,nrow(X)), rep(0,length(wvec))))
        }
        pos <- !is.na(model$coeff)
        se <- sqrt(rowSums((X[, pos, drop=F] %*% vcv) * X[, pos, drop=F]))
      }
   }
   return(list(fit=fit, se=se))
}

# square root matrix
lssqrtm <- function(A) {
  decomp <- svd(A)
  decomp$d[decomp$d < 0] <- 0
  rootA  <- decomp$u %*% diag(sqrt(decomp$d)) %*% t(decomp$u)
  return(rootA)
}

# pval, cval simulation (tstat should be 2-col matrix)
binsreg.pval <- function(num, denom, rep, tstat=NULL, side=NULL, alpha, lp=Inf) {
  tvec <- c()
  pval <- NA
  if (!is.null(tstat)) pval <- rep(0, nrow(tstat))
  cval <- NA
  k <- ncol(num)

  for (i in 1:rep) {
    eps <- matrix(rnorm(k, 0, 1), ncol = 1)
    tx  <- (num %*% eps) / denom

    # for critical value
    if (!is.null(side)) {
       if (side == "two") {
          if (is.infinite(lp)) tvec[i] <- max(abs(tx))
          else                 tvec[i] <- mean(abs(tx)^lp)^(1/lp)
       } else if (side == "left") {
          tvec[i] <- max(tx)
       } else if (side == "right") {
          tvec[i] <- min(tx)
       }
    }

    # for p value
    if (!is.null(tstat)) {
       for (j in 1:nrow(tstat)) {
         # 1: left; 2: right; 3: two-sided
         if (tstat[j,2] == 1) {
           pval[j] <- pval[j] + (max(tx) >= tstat[j,1])
         } else if (tstat[j,2] == 2) {
           pval[j] <- pval[j] + (min(tx) <= tstat[j,1])
         } else if (tstat[j,2] == 3) {
           if (is.infinite(lp)) pval[j] <- pval[j] + (max(abs(tx)) >= tstat[j,1])
           else                 pval[j] <- pval[j] + (mean(abs(tx)^lp)^(1/lp) >= tstat[j,1])
         }
       }
    }
  }

  if (!is.null(tstat)) {
    pval <- pval / rep
  }

  if (!is.null(side)) {
    cval <- quantile(tvec, alpha/100, na.rm=T, names = F, type=2)
  }

  return(list(pval=pval, cval=cval))
}

# pval used only by binspwc
binspwc.pval <- function(nummat1, nummat2, denom1, denom2, rep, tstat=NULL, testtype=NULL, lp=Inf) {
  pval <- 0
  k1 <- ncol(nummat1); k2 <- ncol(nummat2)

  for (i in 1:rep) {
    eps1 <- matrix(rnorm(k1, 0, 1), ncol = 1)
    eps2 <- matrix(rnorm(k2, 0, 1), ncol = 1)
    tx  <- (nummat1 %*% eps1 - nummat2 %*% eps2) / sqrt(denom1^2+denom2^2)

    # for p value
    if (testtype == "left") {
      pval <- pval + (max(tx) >= tstat)
    } else if (testtype == "right") {
      pval <- pval + (min(tx) <= tstat)
    } else {
      if (is.infinite(lp)) pval <- pval + (max(abs(tx)) >= tstat)
      else                 pval <- pval + (mean(abs(tx)^lp)^(1/lp) >= tstat)
    }
  }
  pval <- pval / rep
  return(pval)
}

# IMSE V constant
imse.vcons <- function(p, deriv) {
  n <- p + 1
  V <- as.matrix(1/sapply(1:n, function(x) x:(x+n-1)))
  Vderiv <- matrix(0, n, n)
  for (r in (deriv+1):n) {
    for (c in (deriv+1):n) {
      Vderiv[r,c] <- 1/(r+c-1-2*deriv)*
            (factorial(r-1)/factorial(r-1-deriv))*
            (factorial(c-1)/factorial(c-1-deriv))
    }
  }
  vcons <- sum(diag(solve(V, Vderiv)))
  return(vcons)
}

# IMSE B constant
imse.bcons <- function(p, s=0, deriv) {
  ord <- p + 1
  if (s == 0) {
     bcons <- 1 / (2*(ord-deriv) + 1) / factorial(ord-deriv)^2 / choose(2*(ord-deriv), ord-deriv)^2
  } else {
    if (p==0) {
      bernum <- 1/6
    } else if (p==1) {
      bernum <- 1/30
    } else if (p==2) {
      bernum <- 1/42
    } else  if (p==3) {
      bernum <- 1/30
    } else if (p==4) {
      bernum <- 5/66
    } else if (p==5) {
      bernum <- 691/2730
    } else if (p==6) {
      bernum <- 7/6
    } else {
      print("p>6 not allowed.")
      stop()
    }
     bcons <- 1 / factorial(2*(ord-deriv)) * bernum
  }
  return(bcons)
}

bernpoly <- function(x, p) {
  n <- length(x)
  if (p==0) {
    bernx <- rep(1, n)
  } else if (p==1) {
    bernx <- x - 0.5
  } else if (p==2) {
    bernx <- x^2-x+1/6
  } else if (p==3) {
    bernx <- x^3 - 1.5*x^2 + 0.5*x
  } else if (p==4) {
    bernx <- x^4 - 2*x^3 + x^2 - 1/30
  } else if (p==5) {
    bernx <- x^5 - 2.5*x^4 + 5/3*x^3 - 1/6*x
  } else if (p==6) {
    bernx <- x^6 - 3*x^5 + 2.5*x^4 - 0.5*x^2 + 1/42
  } else {
    print("p>6 not allowed.")
    stop()
  }
  return(bernx)
}

# ROT selector
binsregselect.rot <- function(y, x, w, p, s, deriv, es=F, eN, norotnorm=F, qrot=2, den.alpha=0.975, weights=NULL) {
  x <- (x - min(x)) / (max(x) - min(x))
  ord <- p+1
  N <- length(x)

  x.p <- matrix(NA, N, p+qrot+1)
  for (j in 1:(p+qrot+1))  x.p[,j] <- x^(j-1)
  P <- cbind(x.p, w)
  est <- lm(y~P-1, weights=weights)
  beta <- est$coefficients; est <- est$fitted.values

  # variance constant
  s2 <- lm(y^2~P-1, weights=weights)$fitted.values - est^2
  if (norotnorm) {
    fz <- 1
  } else {
    x.summ <- binsreg.summ(x, w=weights, std=T)
    fz <- dnorm(x, x.summ$mu, x.summ$sig)
    # trim density from below
    cutval <- dnorm(qnorm(den.alpha)*x.summ$sig, 0, x.summ$sig)
    fz[fz<cutval] <- cutval
  }
  if (es) {
    s2 <- s2 / fz
  } else {
    s2 <- s2 * (fz^(2*deriv))
  }

  s2 <- binsreg.summ(s2, w=weights, std=F)$mu
  vcons <- imse.vcons(p, deriv)
  imse.v <- vcons * s2

  # bias constant
  bcons <- imse.bcons(p, s=0, deriv)
  mu.m.hat <- 0
  for (j in (p+1):(p+qrot)) {
    mu.m.hat <- mu.m.hat + x^(j-p-1)*beta[j+1]*factorial(j)/factorial(j-p-1)
  }
  mu.m.hat <- mu.m.hat^2
  if (!es) {
    mu.m.hat <- mu.m.hat / (fz^(2*ord-2*deriv))
  }
  imse.b <- bcons * binsreg.summ(mu.m.hat, w=weights, std=F)$mu

  J.rot <- ceiling((imse.b*2*(ord-deriv)/(imse.v*(1+2*deriv)))^(1/(2*ord+1)) * eN^(1/(2*ord+1)))
  return(list(J.rot=J.rot, imse.b=imse.b, imse.v=imse.v))
}

# locate h
binsreg.locate <- function(x, knot, type="all") {
  h <- tl <- NA
  loc.ind <- findInterval(x, knot, rightmost.closed = T)
  if (type=="all"|type=="h") {
    size    <- diff(knot)
    h       <- size[loc.ind]
  }
  if (type=="all"|type=="tl") {
    pos     <- knot[-length(knot)]
    tl      <- pos[loc.ind]
  }
  return(list(h=h, tl=tl))
}

# IMSE V cons
genV <- function(y, x, w, p, s, deriv, knot, vce, cluster=NULL, weights=NULL) {
  B  <- binsreg.spdes(eval=x, p=p, s=s, knot=knot, deriv=0)
  k  <- ncol(B)
  if (deriv>0) {
     basis <- binsreg.spdes(eval=x, p=p, s=s, knot=knot, deriv=deriv)
  } else {
     basis <- B
  }
  P     <- cbind(B, w)
  model  <- lm(y ~ P-1, weights = weights)
  pos <- !is.na(model$coeff[1:k])
  k.new <- sum(pos)
  vcv <- binsreg.vcov(model, type=vce, cluster=cluster)[1:k.new, 1:k.new]
  m.s2   <- binsreg.summ(rowSums((basis[,pos, drop=F] %*% vcv) * basis[,pos, drop=F]), w=weights, std=F)$mu
  return(m.s2)
}

# bias term
bias <- function(x, p, s, deriv, knot) {
  locate <- binsreg.locate(x, knot)
  h  <- locate$h
  tl <- locate$tl
  bern <- bernpoly((x-tl)/h, p+1-deriv) / factorial(p+1-deriv) * (h^(p+1-deriv))
  return(bern)
}

genB <- function(y, x, w, p, s, deriv, knot, weights=NULL) {
  B  <- binsreg.spdes(eval=x, p=p+1, s=s+1, knot=knot, deriv=0)   # use smoothest spline
  k    <- ncol(B)
  P    <- cbind(B, w)
  beta <- lm(y~P-1, weights=weights)$coefficients[1:k]
  pos  <- !is.na(beta)
  basis <- binsreg.spdes(eval=x, p=p+1, s=s+1, knot=knot, deriv=p+1)
  mu.m.fit  <- basis[, pos, drop=F] %*% beta[pos]

  bias.0 <- mu.m.fit * bias(x, p, s, 0, knot)    # proj component, v=0!!!
  B    <- binsreg.spdes(eval=x, p=p, s=s, knot=knot, deriv=0)
  beta <- lm(bias.0~B-1, weights=weights)$coefficients
  pos <- !is.na(beta)
  if (deriv > 0) {
     basis <- binsreg.spdes(eval=x, p=p, s=s, knot=knot, deriv=deriv)
     bias.v <- mu.m.fit * bias(x, p, s, deriv, knot)    # need to recalculate for v>0!!!
  } else {
     basis <- B
     bias.v <- bias.0
  }
  bias.l2 <- bias.v - basis[,pos,drop=F] %*% beta[pos]
  bias.cons <- binsreg.summ(bias.l2^2, w=weights, std=F)$mu

  return(bias.cons)
}

# DPI selector
binsregselect.dpi <- function(y, x, w, p, s, deriv, es=F, vce, cluster=NULL, nbinsrot, weights=NULL) {
  J.rot <- nbinsrot
  x <- (x - min(x))/ (max(x) - min(x))
  ord <- p + 1

  if (es) {
    knot <- genKnot.es(0, 1, J.rot)
  } else {
    knot <- genKnot.qs(x, J.rot)
  }

  # bias constant
  imse.b <- genB(y, x, w, p, s, deriv, knot, weights=weights) * J.rot^(2*(ord-deriv))

  # variance constant
  imse.v <- genV(y, x, w, p, s, deriv, knot, vce, cluster, weights=weights) / (J.rot^(1+2*deriv))

  J.dpi <- ceiling((imse.b*2*(ord-deriv)/((1+2*deriv)*imse.v))^(1/(2*ord+1)))

  return(list(J.dpi=J.dpi, imse.v=imse.v, imse.b=imse.b))
}

# Check local mass points
binsreg.checklocalmass <- function(x, J, es, knot=NULL) {
  if (is.null(knot)) {
    if (es) {
      knot <- genKnot.es(min(x), max(x), J)
    } else {
      knot <- genKnot.qs(x, J)
    }
  }
  pos <- findInterval(x, knot, rightmost.closed = T, left.open = T)
  uniqnum <- sapply(1:J, function(z) length(unique(x[pos==z])))
  return(min(uniqnum))
}

# slightly modified mean fun
binsreg.summ <- function(x, w=NULL, std=F) {
  mu <- sig <- NA
  if (is.null(w)) {
    mu <- mean(x)
    if (std) sig <- sd(x)
  } else {
    mu <- weighted.mean(x, w=w)
    if (std) sig <- sqrt(sum(w*(x-mu)^2)/(sum(w)-1))
  }
  return(list(mu=mu, sig=sig))
}

# extract formula information
binsreg.model.mat <- function(formula=NULL, data=NULL) {
  des <- model.matrix(formula, data=data)
  col.term <- attr(des, "assign")                    # each column comes from which var
  factor.strlist <- names(attr(des, "contrasts"))    # names of factor vars.

  formula.info <- terms.formula(formula)
  term.mat <- attr(formula.info, "factors")          # which term uses which var.
  row.name <- rownames(term.mat)
  factor.index <- row.name %in% factor.strlist       # pos of factor vars in the string list
  usefactor.pos <- which(colSums(term.mat[factor.index,,drop=F])>0)   # index for terms that use factor vars

  factor.colnum <- col.term %in% usefactor.pos      # pos of columns in the design that use factor vars

  if (attr(formula.info, "intercept")==1) {
    des <- des[,-1,drop=F]
    factor.colnum <- factor.colnum[-1]
  }
  return(list(design=des, factor.colnum=factor.colnum))
}

# check if an object is a formula
is.formula <- function(x) {
  ind <- try(is.call(x) && x[[1]] == quote(`~`), silent = T)
  if (inherits(ind, "try-error")) {
    ind <- F
  }
  return(ind)
}
