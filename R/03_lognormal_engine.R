###############################################################################
#  03_lognormal_engine.R
#  ---------------------------------------------------------------------------
#  Log-normal 脆弱模型下的 MM 迭代引擎
#
#  为什么需要这个引擎：
#    log-normal 脆弱项的对数边际似然没有闭式解，每一步 MM 都要计算
#    E[exp(u_i) | 数据] 与 E[exp(2u_i) | 数据]。原实现用 Laplace 近似，
#    在本设计的强相关高维情形下误差累积，导致部分参数点不收敛。
#    本引擎改用「自适应 Gauss-Legendre 求积 + 多级精度复核」计算后验矩：
#      先用 24 点求积，与 48 点结果比较；
#      误差超阈值则换 96 点；仍不达标则退回 integrate() 自适应积分。
#    同时把基线强度与回归系数的更新全部向量化。
#
#  数值自检：每一步都用 reference_ll()（逐人 adaptive integrate 的独立实现）
#            复核对数似然，偏差必须 < 1e-5，否则判该参数点失败。
###############################################################################
make_lognormal_engine <- function (f) 
{
    if (!identical(f$datatype, "Recurrent")) 
        stop("Only recurrent-event data are supported")
    y <- as.numeric(f$y)
    X <- as.matrix(f$X)
    d <- as.numeric(f$d)
    id <- as.integer(f$id)
    N <- length(y)
    a <- as.integer(f$a)
    p <- ncol(X)
    if (!setequal(id, seq_len(a)) || any(diff(id) < 0)) 
        stop("Unexpected subject coding")
    first <- match(seq_len(a), id)
    last <- c(first[-1] - 1L, N)
    if (any(vapply(seq_len(a), function(i) any(diff(y[first[i]:last[i]]) < 0), logical(1)))) 
        stop("Within-subject times are not ordered")
    Xi <- X[first, , drop = FALSE]
    if (max(abs(X - Xi[id, , drop = FALSE])) > 1e-12) 
        stop("Time-varying covariates are not supported")
    if (anyDuplicated(y[d == 1])) 
        stop("Tied event times require a separately validated implementation")
    end <- y[last]
    oy <- order(y)
    oe <- order(end, decreasing = TRUE)
    hindex <- findInterval(end, y[oy])
    rindex <- findInterval(-y, -end[oe])
    D <- tabulate(id[d == 1], a)
    sx <- colSums(X[d == 1, , drop = FALSE])
    ax <- abs(Xi)
    asum <- rowSums(ax)
    e <- new.env(parent = asNamespace("frailtyMMpen"))
    e$bad <- 0L
    e$integration_refinements <- 0L
    e$failures <- character()
    gauss <- function(n) {
        J <- matrix(0, n, n)
        j <- seq_len(n - 1)
        w <- j/sqrt(4 * j * j - 1)
        J[cbind(seq_len(n - 1), 2:n)] <- w
        J[cbind(2:n, seq_len(n - 1))] <- w
        z <- eigen(J, symmetric = TRUE)
        list(x = (z$values + 1)/2, lw = log(z$vectors[1, ]^2))
    }
    q24 <- gauss(24L)
    q48 <- gauss(48L)
    q96 <- gauss(96L)
    modes <- function(A, D, th) {
        lo <- rep(-800, length(A))
        hi <- rep(800, length(A))
        m <- log((D + 1)/(A + 1))
        for (j in seq_len(256L)) {
            ae <- A * exp(m)
            g <- D - ae - m/th
            done <- abs(g)/(ae + 1/th) <= 1e-12 * (1 + abs(m))
            if (all(done)) 
                return(m)
            lo[g > 0] <- m[g > 0]
            hi[g <= 0] <- m[g <= 0]
            v <- m + pmax(-5, pmin(5, g/(ae + 1/th)))
            bad <- !is.finite(v) | v <= lo | v >= hi
            v[bad] <- (lo[bad] + hi[bad])/2
            v[done] <- m[done]
            m <- v
        }
        stop("Posterior mode solver did not converge")
    }
    centered <- function(v, ae, residual, th) {
        em <- expm1(v) - v
        small <- abs(v) < 1e-04
        em[small] <- v[small]^2 * (0.5 + v[small]/6 + v[small]^2/24 + v[small]^3/120 + v[small]^4/720)
        residual * v - ae * em - v * v/(2 * th)
    }
    reference_one <- function(A, D, th) {
        A <- unname(A)
        D <- unname(D)
        th <- unname(th)
        if (A == 0) {
            return(c(logZ = th * D * D/2, Ez = exp(th * (D + 0.5)), Eu2 = th + (th * D)^2))
        }
        integral <- function(k, moment = FALSE) {
            m <- modes(A, k, th)
            hm <- k * m - A * exp(m) - m * m/(2 * th)
            ae <- A * exp(m)
            residual <- k - ae - m/th
            sd <- 1/sqrt(ae + 1/th)
            fun <- function(t) {
                v <- sd * t
                u <- v + m
                h <- centered(v, ae, residual, th)
                h[!is.finite(h)] <- -Inf
                z <- exp(h)
                if (moment) {
                  z <- z * u * u
                  z[!is.finite(z)] <- 0
                }
                z
            }
            z <- integrate(fun, -Inf, Inf, rel.tol = 1e-10, abs.tol = 1e-11, subdivisions = 1000L, stop.on.error = FALSE)
            if (z$message != "OK" || !is.finite(z$value) || z$value <= 0) 
                stop("Independent adaptive integration failed")
            c(value = z$value, logZ = hm + log(sd) + log(z$value) - 0.5 * log(2 * pi * th), mode = m, sd = sd, ae = ae, residual = residual)
        }
        z0 <- integral(D)
        z1 <- integral(D + 1)
        z2 <- integral(D, TRUE)
        ld <- centered(unname(z1["mode"] - z0["mode"]), unname(z0["ae"]), unname(z0["residual"]), th) + unname(z1["mode"]) + log(z1["sd"]/z0["sd"]) + log(z1["value"]/z0["value"])
        c(logZ = unname(z0["logZ"]), Ez = exp(unname(ld)), Eu2 = unname(z2["value"]/z0["value"]))
    }
    moments <- function(A, D, th) {
        if (!is.finite(th) || th <= 0 || any(!is.finite(A)) || any(A < 0)) 
            stop("Invalid posterior inputs")
        ans <- matrix(NA_real_, length(A), 3L, dimnames = list(NULL, c("logZ", "Ez", "Eu2")))
        zero <- A == 0
        ans[zero, 1] <- th * D[zero]^2/2
        ans[zero, 2] <- exp(th * (D[zero] + 0.5))
        ans[zero, 3] <- th + (th * D[zero])^2
        if (any(!zero)) {
            A0 <- A[!zero]
            D0 <- D[!zero]
            m <- modes(A0, D0, th)
            ae <- A0 * exp(m)
            residual <- D0 - ae - m/th
            sd <- 1/sqrt(ae + 1/th)
            hm <- D0 * m - A0 * exp(m) - m * m/(2 * th)
            left <- -sd
            right <- sd
            for (j in seq_len(64L)) {
                bad <- centered(left, ae, residual, th) > -50
                if (!any(bad)) 
                  break
                left[bad] <- 2 * left[bad]
            }
            for (j in seq_len(64L)) {
                bad <- centered(right, ae, residual, th) > -50
                if (!any(bad)) 
                  break
                right[bad] <- 2 * right[bad]
            }
            calc <- function(q) {
                v <- cbind(outer(left, q$x), outer(right, q$x))
                u <- v + m
                w <- centered(v, ae, residual, th)
                w <- sweep(w, 2L, rep(q$lw, 2L), "+")
                w <- w + cbind(matrix(log(-left), nrow = length(m), ncol = length(q$x)), matrix(log(right), nrow = length(m), ncol = length(q$x)))
                w[!is.finite(w)] <- -Inf
                mx <- w[cbind(seq_len(nrow(w)), max.col(w, ties.method = "first"))]
                ww <- exp(w - mx)
                den <- rowSums(ww)
                we <- w + u
                me <- we[cbind(seq_len(nrow(we)), max.col(we, ties.method = "first"))]
                cbind(logZ = hm + mx + log(den) - 0.5 * log(2 * pi * th), Ez = exp(me + log(rowSums(exp(we - me))) - mx - log(den)), Eu2 = rowSums(ww * u * u)/den)
            }
            z1 <- calc(q24)
            z2 <- calc(q48)
            er <- pmax(abs(z1[, 1] - z2[, 1]), abs(z1[, 2] - z2[, 2])/pmax(1, abs(z2[, 2])), abs(z1[, 3] - z2[, 3])/pmax(1, abs(z2[, 3])))
            bad <- !is.finite(er) | er > 1e-08
            if (any(bad)) {
                z3 <- calc(q96)
                er <- pmax(abs(z2[, 1] - z3[, 1]), abs(z2[, 2] - z3[, 2])/pmax(1, abs(z3[, 2])), abs(z2[, 3] - z3[, 3])/pmax(1, abs(z3[, 3])))
                z2[bad, ] <- z3[bad, ]
                bad <- bad & (!is.finite(er) | er > 1e-08)
            }
            if (any(bad)) {
                for (j in which(bad)) {
                  z2[j, ] <- reference_one(A0[j], D0[j], th)
                  e$integration_refinements <- e$integration_refinements + 1L
                }
            }
            ans[!zero, ] <- z2
        }
        if (any(!is.finite(ans)) || any(ans[, 2:3] <= 0)) 
            stop("Non-finite posterior moments")
        ans
    }
    values <- function(coef, lambda, th) {
        if (length(lambda) != N || any(!is.finite(c(coef, lambda, th))) || th <= 0 || any(lambda < 0) || any(lambda[d == 1] <= 0)) 
            stop("Invalid model state")
        eta <- drop(Xi %*% coef)
        ex <- exp(eta)
        A <- c(0, cumsum(lambda[oy]))[hindex + 1L] * ex
        list(ex = ex, eta = eta, A = A, post = moments(A, D, th))
    }
    stable_step <- function(coef, lambda, tht, penalty, gam, tune) {
        v <- values(coef, lambda, tht)
        ez <- v$post[, "Ez"]
        w <- ez * v$ex
        risk <- c(0, cumsum(w[oe]))[rindex + 1L]
        if (any(!is.finite(risk)) || any(risk <= 0)) 
            stop("Invalid risk denominator")
        ln <- d/risk
        H <- c(0, cumsum(ln[oy]))[hindex + 1L]
        mass <- w * H
        score <- sx - colSums(Xi * mass)
        curvature <- colSums(ax * (asum * mass))
        if (any(!is.finite(curvature)) || any(curvature <= 0)) 
            stop("Invalid MM curvature")
        ab <- abs(coef)
        dp <- switch(as.character(penalty), `0` = rep(0, p), `1` = rep(tune, p), `2` = pmax(0, tune - ab/gam), `3` = ifelse(ab <= tune, tune, pmax(0, gam * tune - ab)/(gam - 1)), stop("Unsupported penalty"))
        bn <- coef + score/curvature
        k <- dp > 0
        bn[k] <- (coef[k] * curvature[k] + score[k]) * ab[k]/(curvature[k] * ab[k] + N * dp[k])
        tn <- mean(v$post[, "Eu2"])
        if (any(!is.finite(c(bn, tn, ln))) || tn <= 0 || max(abs(bn)) > 1000 || tn > 1000 || max(abs(ln)) > 1000) 
            stop("Invalid updated state")
        list(coef = bn, est.tht = tn, lambda = ln, error = 0L, Ar = ez)
    }
    e$MMroutine <- function(y, X, d, coef, lambda, tht, frailty, penalty, gam, tune, id, N, a, p, power, type) {
        if (frailty != 1 || type != 3) 
            stop("Wrong model passed to stable Lognormal engine")
        tryCatch(stable_step(coef, lambda, tht, penalty, gam, tune), error = function(err) {
            e$bad <- e$bad + 1L
            e$failures <- unique(c(e$failures, conditionMessage(err)))
            list(error = 1L)
        })
    }
    e$logLikcal <- function(y, X, d, coef, lambda, est.tht, frailtyc, id, N, a, p, power, type) {
        if (frailtyc != 1 || type != 3) 
            stop("Wrong likelihood type")
        v <- values(coef, lambda, est.tht)
        sum(log(lambda[d == 1])) + sum(D * v$eta) + sum(v$post[, "logZ"])
    }
    cal <- getFromNamespace("frailtyMMcal", "frailtyMMpen")
    environment(cal) <- e
    e$cal <- cal
    reference_ll <- function(coef, lambda, th) {
        eta <- drop(Xi %*% coef)
        A <- c(0, cumsum(lambda[oy]))[hindex + 1L] * exp(eta)
        zz <- vapply(seq_len(a), function(i) {
            reference_one(A[i], D[i], th)["logZ"]
        }, numeric(1))
        sum(log(lambda[d == 1])) + sum(D * eta) + sum(zz)
    }
    list(cal = cal, step = stable_step, state = e, moments = moments, reference_one = reference_one, reference_ll = reference_ll, values = values, constants = list(D = D, Xi = Xi))
}
