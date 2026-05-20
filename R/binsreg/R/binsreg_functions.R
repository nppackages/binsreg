# 07/16/2024
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

binsreg.bin.counts <- function(x, knot, nbins) {
  if (nbins == length(knot)) {
    bin <- match(x, knot)
  } else {
    bin <- findInterval(x, knot, rightmost.closed = T, left.open = T)
  }
  return(tabulate(bin, nbins=nbins))
}

binsreg.stata.irecode <- function(x, knot) {
  if (length(knot) <= 2L) return(rep.int(1L, length(x)))

  # Stata's binsreg_irecode passes cutpoints to irecode() through local macro
  # expansion.  On normalized selector knots this matches %18.0g: 16 decimals.
  cuts <- as.numeric(sprintf("%.16f", knot[-c(1L, length(knot))]))
  pos <- findInterval(x, cuts, left.open=TRUE) + 1L
  pmax.int(1L, pmin.int(length(knot)-1L, pos))
}

binsreg.which.min <- function(x) {
  if (all(is.na(x))) return(1L)
  which.min(x)
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
binsreg.spdes <- function(eval, p, s, knot, deriv, pos=NULL) {
  n <- length(eval)
  k <- length(knot)
  if (is.null(pos)) {
    pos <- findInterval(eval, knot, rightmost.closed = T, left.open = T)
  } else {
    pos <- as.integer(pos)
  }

  if (p == 0 && s == 0) {
    P <- matrix(0, n, k-1L)
    if (deriv == 0 && n > 0) P[cbind(seq_len(n), pos)] <- 1
    return(P)
  }

  if (p == 1 && (s == 0 || s == 1) && (deriv == 0 || deriv == 1)) {
    width <- p - s + 1L
    P <- matrix(0, n, (k-2L)*width + p + 1L)
    if (n > 0) {
      h <- knot[pos+1L] - knot[pos]
      w <- (eval - knot[pos]) / h
      if (deriv == 0) {
        local <- cbind(1-w, w)
      } else {
        local <- cbind(-1/h, 1/h)
      }
      col.base <- if (s == 1) pos else 2L*pos - 1L
      P[cbind(seq_len(n), col.base)] <- local[, 1L]
      P[cbind(seq_len(n), col.base+1L)] <- local[, 2L]
    }
    return(P)
  }

  width <- p - s + 1L
  if (k >= 3L) {
    ext.knot <- c(rep(knot[1L], p+1L), rep(knot[2L:(k-1L)], each=width), rep(knot[k], p+1L))
  } else {
    ext.knot <- c(rep(knot[1L], p+1L), rep(knot[k], p+1L))
  }

  ind.lk <- p + 1L + (pos - 1L)*width
  lk <- rk <- matrix(NA_real_, n, p)
  if (p > 0 && n > 0) {
    for (i in seq_len(p)) {
      lk[, p-i+1L] <- ext.knot[ind.lk-i+1L]
      rk[, i] <- ext.knot[ind.lk+i]
    }
  }

  bs <- matrix(1, n, 1L)
  zero <- matrix(0, n, 1L)
  if (p >= 1L) {
    if (p < deriv) {
      bs <- matrix(0, n, p+1L)
    } else if (p > deriv) {
      for (i in seq_len(p-deriv)) {
        tl <- lk[, (p-i+1L):p, drop=FALSE]
        tr <- rk[, seq_len(i), drop=FALSE]
        w <- (eval - tl) / (tr - tl)
        bs <- cbind((1-w)*bs, zero) + cbind(zero, w*bs)
      }
      if (deriv > 0) {
        for (i in (p-deriv+1L):p) {
          tl <- lk[, (p-i+1L):p, drop=FALSE]
          tr <- rk[, seq_len(i), drop=FALSE]
          w <- 1 / (tr - tl)
          bs <- (cbind(zero, w*bs) - cbind(w*bs, zero)) * i
        }
      }
    } else {
      for (i in seq_len(p)) {
        tl <- lk[, (p-i+1L):p, drop=FALSE]
        tr <- rk[, seq_len(i), drop=FALSE]
        w <- 1 / (tr - tl)
        bs <- (cbind(zero, w*bs) - cbind(w*bs, zero)) * i
      }
    }
  }

  P <- matrix(0, n, (k-2L)*width + p + 1L)
  if (n > 0) {
    col.base <- (pos - 1L)*width + 1L
    col.index <- as.vector(outer(col.base, 0:p, `+`))
    row.index <- rep.int(seq_len(n), times=p+1L)
    P[cbind(row.index, col.index)] <- as.vector(bs)
  }
  return(P)
}

# check drop, display warning
check.drop <- function(beta, k) {
  if (any(is.na(beta[1:k]))) {
    warning("some X-based variables dropped")
  }
}

binsreg.rq.formula <- function(y, P, design.name) {
  fit.formula <- reformulate(design.name, response="y", intercept=FALSE)
  fit.env <- new.env(parent=parent.frame())
  fit.env$y <- y
  assign(design.name, P, envir=fit.env)
  environment(fit.formula) <- fit.env
  return(fit.formula)
}

binsreg.rq.complete <- function(model) {
  if (is.null(model$binsreg.design.name)) return(model)
  if (is.null(model$formula)) {
    model$formula <- binsreg.rq.formula(model$y, model$x, model$binsreg.design.name)
  }
  if (is.null(model$terms)) model$terms <- terms(model$formula)
  if (is.null(model$call)) model$call <- call("rq", formula=model$formula, tau=model$tau)
  if (is.null(model$xlevels)) model$xlevels <- list()
  if (is.null(model$contrasts)) model$contrasts <- NULL
  if (is.null(model$rho)) {
    resid <- as.vector(model$residuals)
    rho.weights <- if (is.null(model$weights)) 1 else model$weights
    model$rho <- sum(rho.weights * resid * (model$tau - (resid < 0)))
  }
  if (is.null(model$model)) {
    model.frame <- data.frame(y=model$y, check.names=FALSE)
    model.frame[[model$binsreg.design.name]] <- I(model$x)
    if (!is.null(model$weights)) model.frame[["(weights)"]] <- model$weights
    attr(model.frame, "terms") <- model$terms
    model$model <- model.frame
  }
  return(model)
}

# wrapper of vcov and vcovCL
binsreg.lm.normal.eq.cond.max <- 1e6

binsreg.vcov.fast.lm.supports <- function(type, cluster) {
  if (is.null(type)) return(TRUE)
  if (is.null(cluster)) return(type %in% c("const", "HC1", "HC2", "HC3"))
  return(type == "HC1")
}

binsreg.vcov.fast.lm <- function(model, type, cluster) {
  if (is.null(model$x)) return(NULL)
  if (!(type %in% c("const", "HC1", "HC2", "HC3"))) return(NULL)
  if (!is.null(cluster) && type != "HC1") return(NULL)

  X <- model$x
  rank <- model$rank
  coefficients <- model$coefficients
  nonalias <- !is.na(coefficients)
  if (sum(nonalias) != rank) return(NULL)
  if (type == "const" && rank != ncol(X)) return(NULL)
  df.resid <- model$df.residual
  if (is.null(df.resid) || is.na(df.resid) || df.resid <= 0) return(NULL)

  if (!is.null(model$binsreg.xtx)) {
    if (rank != ncol(X)) return(NULL)
    XtX.inv <- tryCatch(solve(model$binsreg.xtx), error=function(e) NULL)
    if (is.null(XtX.inv)) return(NULL)
    X.vcov <- X
  } else {
    R <- qr.R(model$qr)[seq_len(rank), seq_len(rank), drop=FALSE]
    XtX.inv.pivot <- tryCatch(chol2inv(R), error=function(e) NULL)
    if (is.null(XtX.inv.pivot)) return(NULL)

    XtX.inv.full <- matrix(0, ncol(X), ncol(X))
    pivot <- model$qr$pivot[seq_len(rank)]
    XtX.inv.full[pivot, pivot] <- XtX.inv.pivot
    if (rank == ncol(X)) {
      XtX.inv <- XtX.inv.full
      X.vcov <- X
    } else {
      XtX.inv <- XtX.inv.full[nonalias, nonalias, drop=FALSE]
      X.vcov <- X[, nonalias, drop=FALSE]
    }
  }

  resid <- as.vector(model$residuals)
  weights <- model$weights
  if (is.null(weights)) weights <- rep.int(1, length(resid))
  if (type == "const") {
    return(XtX.inv * sum(weights * resid^2) / df.resid)
  }

  if (is.null(cluster)) {
    scale <- weights^2 * resid^2
    if (type == "HC1") {
      scale <- scale * length(resid) / df.resid
    } else if (type %in% c("HC2", "HC3")) {
      leverage <- weights * rowSums((X.vcov %*% XtX.inv) * X.vcov)
      leverage <- pmin(leverage, 1 - .Machine$double.eps)
      scale <- scale / (1 - leverage)
      if (type == "HC3") scale <- scale / (1 - leverage)
    }
    X.scaled <- X.vcov * sqrt(scale)
    return(XtX.inv %*% crossprod(X.scaled) %*% XtX.inv)
  }

  cluster <- as.data.frame(cluster)
  if (ncol(cluster) != 1L || nrow(cluster) != length(resid) || anyNA(cluster)) return(NULL)
  cluster <- cluster[[1L]]
  G <- if (is.factor(cluster)) length(levels(cluster)) else length(unique(cluster))
  if (G <= 1L) return(NULL)

  estfun <- X.vcov * (weights * resid)
  score <- rowsum(estfun, group=cluster, reorder=FALSE)
  scale <- (G / (G - 1L)) * ((length(resid) - 1L) / df.resid)
  return(XtX.inv %*% crossprod(score) %*% XtX.inv * scale)
}

binsreg.vcov <- function(model, type, cluster, is.qreg=FALSE, ...) {
  if (is.qreg) {
    model <- binsreg.rq.complete(model)
    V <- summary.rq(model, se=type, covariance = TRUE, cluster=cluster, ...)$cov
  } else {
    if (!is.null(model$binsreg.vcov) &&
        identical(model$binsreg.vcov.type, type) &&
        identical(model$binsreg.vcov.has.cluster, !is.null(cluster))) {
      return(model$binsreg.vcov)
    }
    V <- binsreg.vcov.fast.lm(model, type=type, cluster=cluster)
    if (is.null(V)) {
      if (type=="const") {
        V <- vcov(model)
      } else {
        V <- vcovCL(model, type=type, cluster=cluster)
      }
    }
  }
  return(V)
}

binsreg.lm.block.supports <- function(type, cluster) {
  if (!is.null(cluster)) return(!is.null(type) && type == "HC1")
  return(is.null(type) || type %in% c("const", "HC1"))
}

binsreg.block.crossprod <- function(B, w) {
  Bw <- crossprod(B, w)
  rbind(cbind(crossprod(B), Bw), cbind(t(Bw), crossprod(w)))
}

binsreg.fit.lm.design <- function(y, B, w=NULL, weights=NULL, vcov.type=NULL, cluster=NULL) {
  if (is.null(w) || length(w)==0) {
    return(binsreg.fit.lm(y, B, weights=weights, vcov.type=vcov.type, cluster=cluster))
  }
  if (!binsreg.lm.block.supports(vcov.type, cluster)) {
    return(binsreg.fit.lm(y, binsreg.cbind(B, w), weights=weights, vcov.type=vcov.type, cluster=cluster))
  }

  y.vec <- as.vector(y)
  B <- as.matrix(B)
  w <- as.matrix(w)
  n <- length(y.vec)
  k.b <- ncol(B)
  k.w <- ncol(w)
  k <- k.b + k.w
  if (nrow(B) != n || nrow(w) != n || n < k || k <= 0L) {
    return(binsreg.fit.lm(y, binsreg.cbind(B, w), weights=weights, vcov.type=vcov.type, cluster=cluster))
  }

  if (is.null(weights)) {
    fit.B <- B
    fit.w <- w
    fit.y <- y.vec
  } else {
    sqrt.weights <- sqrt(weights)
    fit.B <- B * sqrt.weights
    fit.w <- w * sqrt.weights
    fit.y <- y.vec * sqrt.weights
  }

  XtX <- binsreg.block.crossprod(fit.B, fit.w)
  cond <- tryCatch(kappa(XtX, exact=TRUE), error=function(e) Inf)
  if (!is.finite(cond) || cond > binsreg.lm.normal.eq.cond.max) {
    return(binsreg.fit.lm(y, binsreg.cbind(B, w), weights=weights, vcov.type=vcov.type, cluster=cluster))
  }

  Xty <- c(crossprod(fit.B, fit.y), crossprod(fit.w, fit.y))
  beta <- tryCatch(solve(XtX, Xty), error=function(e) NULL)
  if (is.null(beta)) {
    return(binsreg.fit.lm(y, binsreg.cbind(B, w), weights=weights, vcov.type=vcov.type, cluster=cluster))
  }
  beta <- as.vector(beta)
  coef.names <- c(colnames(B), colnames(w))
  if (length(coef.names)==length(beta) && any(nzchar(coef.names))) names(beta) <- coef.names

  fitted <- as.vector(B %*% beta[seq_len(k.b)] + w %*% beta[k.b + seq_len(k.w)])
  resid <- y.vec - fitted
  model <- list(
    coefficients = beta,
    residuals = resid,
    fitted.values = fitted,
    rank = k,
    df.residual = n - k,
    y = y.vec,
    binsreg.xtx = XtX
  )
  if (!is.null(weights)) model$weights <- weights
  model$terms <- terms(y ~ -1 + P)
  model$call <- match.call()
  class(model) <- "lm"

  if (!is.null(vcov.type) || !is.null(cluster)) {
    XtX.inv <- tryCatch(solve(XtX), error=function(e) NULL)
    if (is.null(XtX.inv)) {
      return(binsreg.fit.lm(y, binsreg.cbind(B, w), weights=weights, vcov.type=vcov.type, cluster=cluster))
    }
    df.resid <- model$df.residual
    if (is.null(cluster)) {
      if (vcov.type == "const") {
        V <- XtX.inv * sum(if (is.null(weights)) resid^2 else weights * resid^2) / df.resid
      } else {
        wt <- if (is.null(weights)) 1 else weights
        scale <- wt^2 * resid^2
        if (vcov.type == "HC1") scale <- scale * n / df.resid
        sqrt.scale <- sqrt(scale)
        meat <- binsreg.block.crossprod(B * sqrt.scale, w * sqrt.scale)
        V <- XtX.inv %*% meat %*% XtX.inv
      }
    } else {
      cluster.data <- as.data.frame(cluster)
      if (ncol(cluster.data) != 1L || nrow(cluster.data) != n || anyNA(cluster.data)) {
        return(binsreg.fit.lm(y, binsreg.cbind(B, w), weights=weights, vcov.type=vcov.type, cluster=cluster))
      }
      cluster.vec <- cluster.data[[1L]]
      G <- if (is.factor(cluster.vec)) length(levels(cluster.vec)) else length(unique(cluster.vec))
      if (G <= 1L) {
        return(binsreg.fit.lm(y, binsreg.cbind(B, w), weights=weights, vcov.type=vcov.type, cluster=cluster))
      }
      wt.resid <- if (is.null(weights)) resid else weights * resid
      score <- cbind(rowsum(B * wt.resid, group=cluster.vec, reorder=FALSE),
                     rowsum(w * wt.resid, group=cluster.vec, reorder=FALSE))
      V <- XtX.inv %*% crossprod(score) %*% XtX.inv *
        (G / (G - 1L)) * ((n - 1L) / df.resid)
    }
    model$binsreg.vcov <- V
    model$binsreg.vcov.type <- vcov.type
    model$binsreg.vcov.has.cluster <- !is.null(cluster)
  }

  return(model)
}

# faster internal lm/glm fitters for already-built numeric design matrices
binsreg.fit.lm <- function(y, P, weights=NULL, vcov.type=NULL, cluster=NULL) {
  if (binsreg.vcov.fast.lm.supports(vcov.type, cluster)) {
    y.vec <- as.vector(y)
    n <- length(y.vec)
    k <- ncol(P)
    if (n >= k && k > 0L) {
      if (is.null(weights)) {
        fit.x <- P
        fit.y <- y.vec
      } else {
        sqrt.weights <- sqrt(weights)
        fit.x <- P * sqrt.weights
        fit.y <- y.vec * sqrt.weights
      }
      XtX <- crossprod(fit.x)
      cond <- tryCatch(kappa(XtX, exact=TRUE), error=function(e) Inf)
      if (is.finite(cond) && cond <= binsreg.lm.normal.eq.cond.max) {
        beta <- tryCatch(solve(XtX, crossprod(fit.x, fit.y)), error=function(e) NULL)
        if (!is.null(beta)) {
          beta <- as.vector(beta)
          names(beta) <- colnames(P)
          fitted <- as.vector(P %*% beta)
          model <- list(
            coefficients = beta,
            residuals = y.vec - fitted,
            fitted.values = fitted,
            rank = k,
            df.residual = n - k,
            x = P,
            y = y.vec,
            binsreg.xtx = XtX
          )
          if (!is.null(weights)) model$weights <- weights
          model$terms <- terms(y ~ -1 + P)
          model$call <- match.call()
          class(model) <- "lm"
          return(model)
        }
      }
    }
  }

  if (is.null(weights)) {
    model <- lm.fit(x=P, y=y)
  } else {
    model <- lm.wfit(x=P, y=y, w=weights)
  }
  model$terms <- terms(y ~ -1 + P)
  model$call <- match.call()
  model$x <- P
  model$y <- y
  if (!is.null(weights)) model$weights <- weights
  class(model) <- "lm"
  return(model)
}

binsreg.fit.glm <- function(y, P, family, weights=NULL, ...) {
  dots <- list(...)
  dot.names <- names(dots)
  allowed <- setdiff(names(formals(glm.fit)), c("x", "y", "weights", "family"))
  if (length(dots) && (is.null(dot.names) || any(dot.names == "") || any(!dot.names %in% allowed))) {
    return(glm(y ~ P - 1, family=family, weights=weights, ...))
  }
  model <- do.call(glm.fit, c(list(x=P, y=y, weights=weights, family=family), dots))
  model$terms <- terms(y ~ -1 + P)
  model$call <- match.call()
  model$x <- P
  model$y <- y
  model$prior.weights <- if (is.null(weights)) rep.int(1, length(y)) else weights
  model$contrasts <- NULL
  model$xlevels <- list()
  model$formula <- y ~ -1 + P
  class(model) <- c("glm", "lm")
  return(model)
}

binsreg.fit.rq <- function(y, P, tau, weights=NULL, qregopt=NULL, design.name="P") {
  if (is.null(qregopt)) qregopt <- list()

  opt.names <- names(qregopt)
  reserved <- c("formula", "data", "subset", "na.action", "model", "contrasts",
                "x", "y", "tau", "weights")
  if (length(tau)!=1L || (length(qregopt) &&
      (is.null(opt.names) || any(opt.names=="") || any(opt.names %in% reserved)))) {
    fit.formula <- binsreg.rq.formula(y, P, design.name)
    return(do.call(rq, c(list(formula=fit.formula, tau=tau, weights=weights), qregopt)))
  }

  if (is.null(weights)) {
    model <- if (length(qregopt)==0L) {
      rq.fit(P, y, tau=tau)
    } else {
      do.call(rq.fit, c(list(x=P, y=y, tau=tau), qregopt))
    }
  } else {
    model <- if (length(qregopt)==0L) {
      rq.wfit(P, y, tau=tau, weights=weights)
    } else {
      do.call(rq.wfit, c(list(x=P, y=y, tau=tau, weights=weights), qregopt))
    }
  }

  coef.names <- colnames(P)
  if (is.null(coef.names)) {
    coef.names <- if (ncol(P)==1L) design.name else paste0(design.name, seq_len(ncol(P)))
  } else {
    coef.names <- paste0(design.name, coef.names)
  }
  if (length(model$coefficients)==length(coef.names)) names(model$coefficients) <- coef.names

  model$x <- P
  model$y <- y
  if (!is.null(weights)) model$weights <- weights
  model$binsreg.design.name <- design.name
  model$tau <- tau
  model$method <- if (is.null(qregopt$method)) "br" else qregopt$method

  class(model) <- "rq"
  return(model)
}

binsreg.cbind <- function(x, w) {
  if (is.null(w) || length(w)==0) x else cbind(x, w)
}

binsregselect.cache.env <- new.env(parent=emptyenv())
binsregselect.cache.env$entries <- list()

binsregselect.eval <- function(args) {
  arg.names <- names(args)
  if (is.null(arg.names)) arg.names <- rep("", length(args))
  unnamed <- !nzchar(arg.names)
  if (any(unnamed)) {
    arg.names[unnamed] <- names(formals(binsregselect))[seq_along(args)][unnamed]
  }
  names(args) <- arg.names

  eval.env <- list2env(args, parent=parent.frame())
  call.args <- lapply(arg.names, as.name)
  names(call.args) <- arg.names
  eval(as.call(c(list(quote(binsregselect)), call.args)), envir=eval.env)
}

binsregselect.cached <- function(...) {
  args <- list(...)
  # randcut consumes RNG, so cache only deterministic internal selections.
  if (!is.null(args$randcut)) {
    return(binsregselect.eval(args))
  }

  key <- tryCatch(serialize(args, NULL, version=2), error=function(e) NULL)
  if (is.null(key)) {
    return(binsregselect.eval(args))
  }

  entries <- binsregselect.cache.env$entries
  for (entry in entries) {
    if (identical(entry$key, key)) {
      for (msg in entry$warnings) warning(msg, call.=FALSE)
      return(entry$value)
    }
  }

  warnings <- character()
  out <- withCallingHandlers(binsregselect.eval(args),
                             warning=function(w) warnings <<- c(warnings, conditionMessage(w)))
  entries <- c(list(list(key=key, value=out, warnings=warnings)), entries)
  if (length(entries) > 8) entries <- entries[seq_len(8)]
  binsregselect.cache.env$entries <- entries
  return(out)
}

# internal pred function (model is long regression, NA handled inside)
binsreg.pred <- function(X, model, type="xb", vce, cluster=NULL, deriv=0, wvec=NULL, is.qreg=FALSE, avar=FALSE, vcv=NULL,...) {
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
      if (is.null(vcv)) {
        if (is.qreg) {
          vcv <- binsreg.vcov(model, type=vce, cluster=cluster, is.qreg=TRUE, ...)
        } else {
          vcv <- binsreg.vcov(model, type=vce, cluster=cluster)
        }
      }
      if (avar) {
        pos <- !is.na(model$coeff[1:k])
        k.new <- sum(pos)
        vcv <- vcv[1:k.new, 1:k.new, drop=F]
        X.pos <- X[, pos, drop=F]
        se <- sqrt(rowSums((X.pos %*% vcv) * X.pos))
      } else {
        if (!is.null(wvec) && deriv != 0) {
          pos <- !is.na(model$coeff[1:k])
          k.new <- sum(pos)
          vcv <- vcv[1:k.new, 1:k.new, drop=F]
        } else {
          if (!is.null(wvec)) {
            X <- cbind(X, matrix(wvec, nrow=nrow(X), ncol=length(wvec), byrow=TRUE))
          }
          pos <- !is.na(model$coeff)
        }
        X.pos <- X[, pos, drop=F]
        se <- sqrt(rowSums((X.pos %*% vcv) * X.pos))
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
  tvec <- if (!is.null(side)) numeric(rep) else NULL
  pval <- NA
  if (!is.null(tstat)) pval <- numeric(nrow(tstat))
  cval <- NA
  k <- ncol(num)

  chunk.size <- min(rep, max(1L, floor(5e6 / max(1L, nrow(num)))))
  start <- 1L
  while (start <= rep) {
    chunk <- min(chunk.size, rep - start + 1L)
    eps <- matrix(rnorm(k * chunk, 0, 1), nrow = k)
    tx <- (num %*% eps) / denom

    max.tx <- matrixStats::colMaxs(tx)
    min.tx <- matrixStats::colMins(tx)
    abs.tx <- if (is.infinite(lp)) {
      matrixStats::colMaxs(abs(tx))
    } else {
      colMeans(abs(tx)^lp)^(1/lp)
    }

    if (!is.null(side)) {
      idx <- start:(start + chunk - 1L)
      if (side == "two") {
        tvec[idx] <- abs.tx
      } else if (side == "left") {
        tvec[idx] <- max.tx
      } else if (side == "right") {
        tvec[idx] <- min.tx
      }
    }

    if (!is.null(tstat)) {
      for (j in seq_len(nrow(tstat))) {
        # 1: left; 2: right; 3: two-sided
        if (tstat[j,2] == 1) {
          pval[j] <- pval[j] + sum(max.tx >= tstat[j,1])
        } else if (tstat[j,2] == 2) {
          pval[j] <- pval[j] + sum(min.tx <= tstat[j,1])
        } else if (tstat[j,2] == 3) {
          pval[j] <- pval[j] + sum(abs.tx >= tstat[j,1])
        }
      }
    }
    start <- start + chunk
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
binspwc.pval <- function(nummat1, nummat2, denom1, denom2, rep, tstat=NULL, testtype=NULL, lp=Inf, alpha=95) {
  pval <- 0
  tvec <- numeric(rep)
  k1 <- ncol(nummat1); k2 <- ncol(nummat2)
  denom <- sqrt(denom1^2+denom2^2)

  chunk.size <- min(rep, max(1L, floor(5e6 / max(1L, nrow(nummat1)))))
  start <- 1L
  while (start <= rep) {
    chunk <- min(chunk.size, rep - start + 1L)
    eps <- matrix(rnorm((k1 + k2) * chunk, 0, 1), nrow = k1 + k2)
    eps1 <- eps[seq_len(k1), , drop=F]
    eps2 <- eps[k1 + seq_len(k2), , drop=F]
    tx <- (nummat1 %*% eps1 - nummat2 %*% eps2) / denom

    max.tx <- matrixStats::colMaxs(tx)
    min.tx <- matrixStats::colMins(tx)
    abs.tx <- if (is.infinite(lp)) {
      matrixStats::colMaxs(abs(tx))
    } else {
      colMeans(abs(tx)^lp)^(1/lp)
    }

    if (testtype == "left") {
      pval <- pval + sum(max.tx >= tstat)
    } else if (testtype == "right") {
      pval <- pval + sum(min.tx <= tstat)
    } else {
      pval <- pval + sum(abs.tx >= tstat)
    }

    idx <- start:(start + chunk - 1L)
    tvec[idx] <- matrixStats::colMaxs(abs(tx))
    start <- start + chunk
  }
  pval <- pval / rep
  cval.cb <- quantile(tvec, alpha/100, na.rm=T, names = F, type=2)
  return(list(pval=pval, cval.cb=cval.cb))
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
binsregselect.rot <- function(y, x, w, p, s, deriv, es=F, eN, norotnorm=F, qrot=2, den.alpha=0.975, weights=NULL,
                              x.norm=NULL) {
  if (is.null(x.norm)) {
    x <- (x - min(x)) / (max(x) - min(x))
  } else {
    x <- x.norm
  }
  ord <- p+1
  N <- length(x)

  x.p <- matrix(NA, N, p+qrot+1)
  for (j in 1:(p+qrot+1))  x.p[,j] <- x^(j-1)
  P <- binsreg.cbind(x.p, w)
  est <- binsreg.fit.lm(y, P, weights=weights)
  beta <- est$coefficients; est <- est$fitted.values

  # variance constant
  s2 <- binsreg.fit.lm(y^2, P, weights=weights)$fitted.values - est^2
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
binsreg.locate <- function(x, knot, type="all", pos=NULL) {
  h <- tl <- NA
  if (is.null(pos)) {
    loc.ind <- findInterval(x, knot, rightmost.closed = T, left.open = T)
  } else {
    loc.ind <- as.integer(pos)
  }
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
genV <- function(y, x, w, p, s, deriv, knot, vce, cluster=NULL, weights=NULL,
                 B0=NULL, basis.deriv=NULL, bin.pos=NULL) {
  if (is.null(B0)) B0 <- binsreg.spdes(eval=x, p=p, s=s, knot=knot, deriv=0, pos=bin.pos)
  B <- B0
  k  <- ncol(B)
  if (!is.null(basis.deriv)) {
     basis <- basis.deriv
  } else if (deriv>0) {
     basis <- binsreg.spdes(eval=x, p=p, s=s, knot=knot, deriv=deriv, pos=bin.pos)
  } else {
     basis <- B
  }
  P     <- binsreg.cbind(B, w)
  model  <- binsreg.fit.lm(y, P, weights=weights, vcov.type=vce, cluster=cluster)
  pos <- !is.na(model$coeff[1:k])
  k.new <- sum(pos)
  vcv <- binsreg.vcov(model, type=vce, cluster=cluster)[1:k.new, 1:k.new]
  m.s2   <- binsreg.summ(rowSums((basis[,pos, drop=F] %*% vcv) * basis[,pos, drop=F]), w=weights, std=F)$mu
  return(m.s2)
}

# bias term
bias <- function(x, p, s, deriv, knot, pos=NULL) {
  locate <- binsreg.locate(x, knot, pos=pos)
  h  <- locate$h
  tl <- locate$tl
  bern <- bernpoly((x-tl)/h, p+1-deriv) / factorial(p+1-deriv) * (h^(p+1-deriv))
  return(bern)
}

genB <- function(y, x, w, p, s, deriv, knot, weights=NULL,
                 B0=NULL, basis.deriv=NULL, bin.pos=NULL,
                 basis.smooth=NULL, basis.smooth.deriv=NULL) {
  if (is.null(basis.smooth)) {
    B <- binsreg.spdes(eval=x, p=p+1, s=s+1, knot=knot, deriv=0, pos=bin.pos)   # use smoothest spline
  } else {
    B <- basis.smooth
  }
  k    <- ncol(B)
  P    <- binsreg.cbind(B, w)
  beta <- binsreg.fit.lm(y, P, weights=weights)$coefficients[1:k]
  pos  <- !is.na(beta)
  if (is.null(basis.smooth.deriv)) {
    basis <- binsreg.spdes(eval=x, p=p+1, s=s+1, knot=knot, deriv=p+1, pos=bin.pos)
  } else {
    basis <- basis.smooth.deriv
  }
  basis.pos <- basis[, pos, drop=F]
  mu.m.fit  <- basis.pos %*% beta[pos]

  bias.0 <- mu.m.fit * bias(x, p, s, 0, knot, pos=bin.pos)    # proj component, v=0!!!
  if (is.null(B0)) B0 <- binsreg.spdes(eval=x, p=p, s=s, knot=knot, deriv=0, pos=bin.pos)
  B <- B0
  if (p == 0 && s == 0 && deriv == 0 && is.null(weights)) {
    xcat <- if (is.null(bin.pos)) findInterval(x, knot, rightmost.closed=TRUE, left.open=TRUE) else bin.pos
    projbias <- rep(NA_real_, length(x))
    for (j in seq_len(length(knot) - 1L)) {
      selected <- xcat == j
      if (any(selected)) projbias[selected] <- mean(bias.0[selected])
    }
    bias.l2 <- bias.0 - projbias
  } else {
    beta <- binsreg.fit.lm(bias.0, B, weights=weights)$coefficients
    pos <- !is.na(beta)
    if (deriv > 0) {
       if (is.null(basis.deriv)) {
         basis <- binsreg.spdes(eval=x, p=p, s=s, knot=knot, deriv=deriv, pos=bin.pos)
       } else {
         basis <- basis.deriv
       }
       bias.v <- mu.m.fit * bias(x, p, s, deriv, knot, pos=bin.pos)    # need to recalculate for v>0!!!
    } else {
       basis <- B
       bias.v <- bias.0
    }
    basis.pos <- basis[, pos, drop=F]
    bias.l2 <- bias.v - basis.pos %*% beta[pos]
  }
  bias.cons <- binsreg.summ(bias.l2^2, w=weights, std=F)$mu
  return(bias.cons)
}

# DPI selector
binsregselect.dpi <- function(y, x, w, p, s, deriv, es=F, vce, cluster=NULL, nbinsrot, weights=NULL,
                              x.norm=NULL, selector.cache=NULL) {
  J.rot <- nbinsrot
  if (is.null(x.norm)) {
    x <- (x - min(x))/ (max(x) - min(x))
  } else {
    x <- x.norm
  }
  ord <- p + 1

  cache.key <- paste(as.integer(es), format(J.rot, digits=17), sep=":")
  if (!is.null(selector.cache)) {
    knot.key <- paste("knot", cache.key, sep=":")
    bin.key <- paste("bin", cache.key, sep=":")
    if (exists(knot.key, envir=selector.cache, inherits=FALSE)) {
      knot <- get(knot.key, envir=selector.cache, inherits=FALSE)
    } else {
      if (es) {
        knot <- genKnot.es(0, 1, J.rot)
      } else {
        knot <- genKnot.qs(x, J.rot)
      }
      assign(knot.key, knot, envir=selector.cache)
    }
    if (exists(bin.key, envir=selector.cache, inherits=FALSE)) {
      bin.pos <- get(bin.key, envir=selector.cache, inherits=FALSE)
    } else {
      bin.pos <- binsreg.stata.irecode(x, knot)
      assign(bin.key, bin.pos, envir=selector.cache)
    }
  } else {
    if (es) {
      knot <- genKnot.es(0, 1, J.rot)
    } else {
      knot <- genKnot.qs(x, J.rot)
    }
    bin.pos <- binsreg.stata.irecode(x, knot)
  }

  get.basis <- function(pp, ss, dd) {
    if (is.null(selector.cache)) {
      return(binsreg.spdes(eval=x, p=pp, s=ss, knot=knot, deriv=dd, pos=bin.pos))
    }
    basis.key <- paste("basis", cache.key, pp, ss, dd, sep=":")
    if (!exists(basis.key, envir=selector.cache, inherits=FALSE)) {
      assign(basis.key, binsreg.spdes(eval=x, p=pp, s=ss, knot=knot, deriv=dd, pos=bin.pos),
             envir=selector.cache)
    }
    get(basis.key, envir=selector.cache, inherits=FALSE)
  }
  B0 <- get.basis(p, s, 0)
  basis.deriv <- if (deriv>0) get.basis(p, s, deriv) else B0
  basis.smooth <- get.basis(p+1, s+1, 0)
  basis.smooth.deriv <- get.basis(p+1, s+1, p+1)

  # bias constant
  imse.b <- genB(y, x, w, p, s, deriv, knot, weights=weights,
                 B0=B0, basis.deriv=basis.deriv,
                 bin.pos=bin.pos, basis.smooth=basis.smooth,
                 basis.smooth.deriv=basis.smooth.deriv) * J.rot^(2*(ord-deriv))

  # variance constant
  genV_val <- genV(y, x, w, p, s, deriv, knot, vce, cluster,
                   weights=weights, B0=B0, basis.deriv=basis.deriv,
                   bin.pos=bin.pos)
  imse.v <- genV_val / (J.rot^(1+2*deriv))
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
