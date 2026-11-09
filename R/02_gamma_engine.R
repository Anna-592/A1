###############################################################################
#  02_gamma_engine.R
#  ---------------------------------------------------------------------------
#  Gamma 脆弱模型下的 MM 迭代引擎
#
#  背景：frailtyMMpen 包内置的 MM 迭代在高维（30 协变量）、强相关噪声协变量下
#        会出现数值溢出与发散。本引擎替换 package 内部的 MMroutine 与 logLikcal
#        两个函数，用「向量化 + 数值保护」的等价实现完成同样的 MM 更新，
#        从而保证 7200 次重复 × 4 个模型全部能稳定收敛。
#
#  实现方式：构造一个新环境 e（父环境为 frailtyMMpen 的命名空间），
#           把 e$MMroutine、e$logLikcal 覆盖为下面的实现，
#           再把 frailtyMMcal 的环境指向 e，最后交给 frailtyMMpen 调用。
#
#  参数说明（统一）：
#    y, X, d, id : 区间形式数据（stop 时间、设计矩阵、事件指示、受试者编号）
#    coef        : 当前回归系数 (p 维)
#    lambda      : 当前累积基线强度增量 (N 维，每个区间一个)
#    th          : 当前脆弱方差参数 theta
#    penalty     : "MCP" / "SCAD"（对应内部编码 2 / 3）
#    gam         : 惩罚函数形状参数
#    tune        : 当前调节参数
###############################################################################

make_gamma_engine <- function(f) {

  ## ---------- 0. 取出并检查数据 ----------
  y  <- as.numeric(f$y)
  X  <- as.matrix(f$X)
  d  <- as.numeric(f$d)
  id <- as.integer(f$id)
  N  <- length(y)
  a  <- as.integer(f$a)          # 受试者人数
  p  <- ncol(X)

  if (!identical(f$datatype, "Recurrent") || any(diff(id) < 0) ||
      !setequal(id, seq_len(a)) || anyDuplicated(y[d == 1])) {
    stop("Unsupported Gamma input")
  }

  ## ---------- 1. 预计算索引（把 O(N^2) 的循环改成向量化查表）----------
  first <- match(seq_len(a), id)                  # 每位受试者第一条记录的全局行号
  last  <- c(first[-1] - 1L, N)                   # 每位受试者最后一条记录
  Xi    <- X[first, , drop = FALSE]               # 每位受试者的协变量（时不变）
  if (max(abs(X - Xi[id, , drop = FALSE])) > 1e-12) {
    stop("Time-varying Xi are unsupported")
  }

  end <- y[last]                                  # 每位受试者的随访终点
  oy  <- order(y)                                 # 按 stop 时间排序
  oe  <- order(end, decreasing = TRUE)            # 按随访终点降序
  hi  <- findInterval(end, y[oy])                 # 终点在排序后 stop 序列中的位置
  ri  <- findInterval(-y, -end[oe])               # 风险集起点位置
  D   <- tabulate(id[d == 1], a)                  # 每位受试者的事件数

  sx   <- colSums(X[d == 1, , drop = FALSE])      # 事件行上协变量的和
  ax   <- abs(Xi)
  asum <- rowSums(ax)

  ## ---------- 2. 建环境并挂上计数 ----------
  e <- new.env(parent = asNamespace("frailtyMMpen"))
  e$bad <- 0L                                     # 数值异常计数
  e$failures <- character()                       # 异常信息收集

  ## ---------- 3. 当前状态下的关键量 ----------
  # A_i : 受试者 i 在随访终点处的累积风险（含脆弱项之外的 exp(eta)）
  # Ez  : E[exp(u_i) | 数据] 的后验期望（Gamma 下有闭式解）
  values <- function(coef, lambda, th) {
    if (any(!is.finite(c(coef, lambda, th))) || th <= 0 ||
        any(lambda < 0) || any(lambda[d == 1] <= 0)) {
      stop("Invalid Gamma state")
    }
    eta <- drop(Xi %*% coef)
    ex  <- exp(eta)
    A   <- c(0, cumsum(lambda[oy]))[hi + 1L] * ex
    list(eta = eta, ex = ex, A = A, Ez = (D + 1 / th) / (A + 1 / th))
  }

  ## ---------- 4. 一步 MM 更新 ----------
  step <- function(coef, lambda, th, penalty, gam, tune) {
    v  <- values(coef, lambda, th)
    AT <- 1 / th + v$A
    DT <- 1 / th + D

    # theta 用一阶/二阶导做一步 Newton（对 log-theta 无约束）
    Q1  <- (a * (digamma(1 / th) + log(th) - 1) +
            sum(v$Ez + log(AT) - digamma(DT))) / th^2
    Q21 <- (a * (3 - 2 * digamma(1 / th) - 2 * log(th)) +
            2 * sum(digamma(DT) - v$Ez - log(AT))) / th^3
    Q22 <- -a * trigamma(1 / th) / th^4
    tn  <- th - Q1 / (Q21 + Q22)
    if (!is.finite(tn) || tn <= 0) tn <- th        # 数值保护：非法则退回原值

    # 基线强度增量 lambda 的更新（Nelson-Aalen 型）
    w    <- v$Ez * v$ex
    risk <- c(0, cumsum(w[oe]))[ri + 1L]
    if (any(!is.finite(risk)) || any(risk <= 0)) {
      stop("Invalid Gamma risk denominator")
    }
    ln <- d / risk
    H  <- c(0, cumsum(ln[oy]))[hi + 1L]
    mass <- w * H

    # 回归系数的得分与曲率（这里用二次上界近似曲率，保证 MM 单调下降）
    score     <- sx - colSums(Xi * mass)
    curvature <- colSums(ax * (asum * mass))
    if (any(curvature <= 0) || any(!is.finite(curvature))) {
      stop("Invalid Gamma curvature")
    }

    # 惩罚项的次梯度（MCP / SCAD 的 MM 上界形式）
    ab <- abs(coef)
    dp <- switch(as.character(penalty),
                 `0` = rep(0, p),
                 `1` = rep(tune, p),
                 `2` = pmax(0, tune - ab / gam),                                  # MCP
                 `3` = ifelse(ab <= tune, tune, pmax(0, gam * tune - ab) / (gam - 1)),  # SCAD
                 stop("Unsupported Gamma penalty"))

    # 无惩罚方向的一步，再做带惩罚的收缩
    bn <- coef + score / curvature
    k  <- dp > 0
    bn[k] <- (coef[k] * curvature[k] + score[k]) * ab[k] /
             (curvature[k] * ab[k] + N * dp[k])

    if (any(!is.finite(c(bn, tn, ln))) || max(abs(bn)) > 1000 ||
        tn > 1000 || max(abs(ln)) > 1000) {
      stop("Invalid updated Gamma state")
    }
    list(coef = bn, est.tht = tn, lambda = ln, error = 0L, Ar = v$Ez)
  }

  ## ---------- 5. 覆盖 package 内部函数 ----------
  e$MMroutine <- function(y, X, d, coef, lambda, tht, frailty, penalty,
                          gam, tune, id, N, a, p, power, type) {
    if (frailty != 0 || type != 3) stop("Wrong model passed to Gamma engine")
    tryCatch(step(coef, lambda, tht, penalty, gam, tune),
             error = function(err) {
               e$bad <- e$bad + 1L
               e$failures <- unique(c(e$failures, conditionMessage(err)))
               list(error = 1L)
             })
  }

  e$logLikcal <- function(y, X, d, coef, lambda, est.tht, frailtyc,
                          id, N, a, p, power, type) {
    if (frailtyc != 0 || type != 3) stop("Wrong Gamma likelihood type")
    v <- values(coef, lambda, est.tht)
    k <- 1 / est.tht
    # Gamma 脆弱下的边际对数似然（闭式）
    sum(log(lambda[d == 1])) + sum(D * v$eta) +
      sum(lgamma(D + k) - lgamma(k) + D * log(est.tht) -
          (D + k) * log1p(est.tht * v$A))
  }

  # 把 frailtyMMcal 的执行环境换成 e，使其内部调用走上面的实现
  cal <- getFromNamespace("frailtyMMcal", "frailtyMMpen")
  environment(cal) <- e

  ## ---------- 6. 独立似然复核（不依赖包内实现）----------
  # 用于自检：MM 算出的 loglik 必须与这里的参考实现相差 < 1e-5
  reference_ll <- function(coef, lambda, th) {
    eta <- drop(Xi %*% coef)
    A <- vapply(end, function(t) sum(lambda[y <= t]), numeric(1)) * exp(eta)
    ratios <- vapply(D, function(n) sum(log1p((seq_len(n) - 1) * th)), numeric(1))
    sum(log(lambda[d == 1])) + sum(D * eta) +
      sum(ratios - (D + 1 / th) * log1p(th * A))
  }

  list(cal = cal, step = step, state = e, reference_ll = reference_ll)
}
