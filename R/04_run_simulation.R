###############################################################################
#  04_run_simulation.R
#  ---------------------------------------------------------------------------
#  72 个情景 × 重复次数的正式模拟主程序
#
#  情景设计（6 因素组合，共 72 个）：
#      n            ∈ {200, 400}                样本量
#      theta        ∈ {0.5, 1, 2}               真实脆弱方差
#      censor_prop  ∈ {0.15, 0.30}              删失比例
#      rho          ∈ {0, 0.25, 0.75}           噪声协变量 AR(1) 相关强度
#      true_frailty ∈ {gamma, lognormal}        真实脆弱分布
#    2 × 3 × 2 × 3 × 2 = 72 个情景，编号 s01 ~ s72
#
#  每个数据集上拟合 4 个模型（2 种假定的脆弱分布 × 2 种惩罚）：
#      gamma_MCP, gamma_SCAD, lognormal_MCP, lognormal_SCAD
#
#  每个模型沿 27 个调节参数 tune = exp(seq(-5.5, 1, 0.25)) 走一条解路径，
#  按 BIC 选出最优模型，再计算 TPR / FDR / 入选变量个数。
#
#  随机种子：BASE_SEED + sid * 1e5 + rep  （BASE_SEED = 20261004）
#
#  用法（默认只跑 2 个情景 × 2 次重复，先确认环境正常）：
#      source("R/04_run_simulation.R")
#      run_simulation(scenarios = 1:2, reps = 1:2, maxit = 20000)
#
#  完整复现（72 情景 × 100 次 = 7200 次，极耗时，建议多核 + 断点续跑）：
#      run_simulation(scenarios = 1:72, reps = 1:100, maxit = 20000)
###############################################################################

suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(frailtyMMpen))

## ============ 1. 实验网格 ============
build_grid <- function() {
  g <- expand.grid(
    n            = c(200, 400),
    theta        = c(0.5, 1, 2),
    censor_prop  = c(0.15, 0.30),
    rho          = c(0, 0.25, 0.75),
    true_frailty = c("gamma", "lognormal"),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  g$scenario_id <- seq_len(nrow(g))       # 72 行，s01 ~ s72
  g
}

# 四个待拟合模型
build_fit_grid <- function() {
  expand.grid(assumed_frailty = c("gamma", "lognormal"),
              penalty = c("MCP", "SCAD"),
              stringsAsFactors = FALSE)
}

## ============ 2. 生成一次数据 ============
generate_dataset <- function(sid, rep_id, grid) {
  s <- grid[grid$scenario_id == sid, ]
  make_recurrent_data_formal(
    n = s$n, theta = s$theta, censor_prop = s$censor_prop,
    rho = s$rho, true_frailty = s$true_frailty,
    seed = BASE_SEED + sid * 1e5 + rep_id
  )
}

## ============ 3. 捕获 frailtyMMpen 的初始化起点 ============
#  为什么要这一步：
#    MM 迭代需要一个合法的初始值 safe.ini。frailtyMMpen 在真正的迭代开始前
#    会先做一次内部初始化，本函数通过替换环境的方式把那次调用的
#    「拟合输入 f」和「初始值 ini」拦下来，后面 4 个模型共用同一个起点，
#    保证 4 个模型的比较是 same-data / same-start 的公平比较。
capture_initialization <- function(dat) {
  w <- new.env(parent = asNamespace("frailtyMMpen"))
  state <- new.env(parent = emptyenv())
  state$bad <- 0L

  raw_step <- getFromNamespace("MMroutine", "frailtyMMpen")
  w$MMroutine <- function(...) {
    z <- raw_step(...)
    if (isTRUE(as.logical(z$error)) ||
        any(!is.finite(c(z$coef, z$lambda, z$est.tht)))) {
      state$bad <- state$bad + 1L
    }
    z
  }

  cal <- getFromNamespace("frailtyMMcal", "frailtyMMpen")
  environment(cal) <- w
  w$frailtyMMcal <- function(y, X, d, N, a, id, ...) {
    args <- list(...)
    if (!is.null(args$penalty)) {
      # 一旦进入带 penalty 的正式拟合，立刻把输入与起点截获并中断
      state$value <- list(
        f = list(y = y, X = X, d = d, N = N, a = a,
                 id = id + 1, coefname = colnames(X),
                 datatype = "Recurrent"),
        ini = args$safe.ini
      )
      stop(structure(list(message = "captured", call = NULL),
                     class = c("initialization_captured", "error", "condition")))
    }
    do.call(cal, c(list(y = y, X = X, d = d, N = N, a = a, id = id), args))
  }

  pub <- getFromNamespace("frailtyMMpen", "frailtyMMpen")
  environment(pub) <- w
  tryCatch(pub(make_formula(), dat, frailty = "lognormal", penalty = "MCP"),
           initialization_captured = function(e) NULL)

  if (state$bad != 0L || is.null(state$value) ||
      any(!is.finite(unlist(state$value$ini)))) {
    stop("初始化阶段出现数值异常。")
  }
  state$value
}

## ============ 4. 沿一条调节参数路径拟合一个模型 ============
run_one_model <- function(f, ini, family, pen, tunes, maxit, tol) {
  p <- ncol(f$X)

  # 选对应的 MM 引擎
  engine <- if (family == "gamma") make_gamma_engine(f) else make_lognormal_engine(f)

  rows <- list(); cols <- list(); last <- ini
  for (j in seq_along(tunes)) {
    before <- engine$state$bad

    # 以当前解为热启动，在下一个 tune 上继续迭代
    cur <- engine$cal(
      f$y, as.matrix(f$X), f$d, f$N, f$a, f$id - 1,
      coef.ini = last$coef, est.tht.ini = last$est.tht,
      lambda.ini = last$lambda, safe.ini = ini,
      frailty = if (family == "gamma") "Gamma" else "LogN",
      penalty = pen, tune = tunes[j],
      maxit = maxit, threshold = tol, type = 3L, SQS1 = 0L
    )
    bad <- engine$state$bad - before

    # 数值体检：全部有限、theta 合法、无内部异常、未触及迭代上限、已收敛
    ok <- all(is.finite(c(cur$coef, cur$lambda, cur$est.tht, cur$likelihood))) &&
      cur$est.tht > 0 && all(cur$lambda >= 0) && bad == 0L &&
      cur$iter < maxit &&
      (all(cur$coef == 0) || isTRUE(cur$convergence <= tol * (p + 1)))

    # BIC：-2 loglik + max(1, log log(p+1)) * (非零系数个数 + 1) * log(受试者数)
    rows[[j]] <- data.frame(
      index = j, tune = tunes[j],
      selected = sum(abs(cur$coef) > 1e-8),
      theta = cur$est.tht, loglik = cur$likelihood,
      BIC = -2 * cur$likelihood +
        max(1, log(log(p + 1))) * (sum(abs(cur$coef) > 1e-6) + 1) * log(f$a),
      iterations = cur$iter, convergence = cur$convergence,
      internal_errors = bad, passes_screen = ok
    )
    cols[[j]] <- cur$coef
    last <- list(coef = cur$coef, est.tht = cur$est.tht, lambda = cur$lambda)

    if (sum(abs(cur$coef)) < 1e-6) break   # 全部系数为 0，路径自然结束
  }

  list(path = do.call(rbind, rows), B = do.call(cbind, cols), last = last, engine = engine)
}

## ============ 5. 一次重复：四个模型 ============
simulate_once <- function(sid, rep_id, grid, fit_grid, tunes, maxit, tol) {
  dat <- generate_dataset(sid, rep_id, grid)
  cap <- capture_initialization(dat)         # 4 个模型共用同一数据与同一初始值
  ps  <- attr(dat, "patient_summary")

  meta <- data.frame(
    scenario_id = sid, rep_id = rep_id,
    seed = BASE_SEED + sid * 1e5 + rep_id,
    n = grid$n[grid$scenario_id == sid],
    theta = grid$theta[grid$scenario_id == sid],
    censor_prop = grid$censor_prop[grid$scenario_id == sid],
    rho = grid$rho[grid$scenario_id == sid],
    true_frailty = grid$true_frailty[grid$scenario_id == sid],
    n_rows = nrow(dat), n_events = sum(dat$event),
    mean_censor_time = mean(ps$censor_time),
    realized_followup_loss = 1 - mean(ps$censor_time) / T_MAX,
    stringsAsFactors = FALSE
  )

  out <- list()
  for (i in seq_len(nrow(fit_grid))) {
    fg <- fit_grid[i, ]
    z  <- run_one_model(cap$f, cap$ini, fg$assumed_frailty, fg$penalty,
                        tunes, maxit, tol)
    path <- z$path
    k    <- which.min(path$BIC)              # BIC 最小处为最终模型
    co   <- z$B[, k]; names(co) <- rownames(z$B)
    selected <- names(co)[abs(co) > 1e-8]

    # 独立似然复核：与 reference_ll 的偏差必须 < 1e-5
    delta <- path$loglik[k] - z$engine$reference_ll(co, z$last$lambda, z$last$est.tht)

    out[[i]] <- cbind(meta, data.frame(
      assumed_frailty = fg$assumed_frailty, penalty = fg$penalty,
      TPR = mean(true_names %in% selected),
      FDR = sum(noise_names %in% selected) / max(1, length(selected)),
      n_selected = length(selected),
      selected_vars = paste(selected, collapse = ","),
      BIC = path$BIC[k], tune = path$tune[k],
      LL_difference = delta,
      all_points_pass = all(path$passes_screen),
      stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, out)
}

## ============ 6. 主循环 ============
run_simulation <- function(scenarios = 1:2, reps = 1:2,
                           maxit = 20000L, tol = 1e-5,
                           out_dir = "results", verbose = TRUE) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  grid     <- build_grid()
  fit_grid <- build_fit_grid()
  tunes    <- exp(seq(-5.5, 1, 0.25))      # 27 个调节参数
  do.call(RNGkind, list("Mersenne-Twister", "Inversion", "Rejection"))

  res <- list()
  for (sid in scenarios) {
    for (r in reps) {
      if (verbose) {
        message(sprintf("情景 s%02d/72（n=%d, theta=%g, censor=%g, rho=%g, true=%s）；重复 %d；种子=%d",
                        sid, grid$n[sid], grid$theta[sid], grid$censor_prop[sid],
                        grid$rho[sid], grid$true_frailty[sid], r,
                        BASE_SEED + sid * 1e5 + r))
      }
      one <- simulate_once(sid, r, grid, fit_grid, tunes, maxit, tol)
      res[[length(res) + 1L]] <- one
      # 每完成一次重复就落盘，中断后不丢已算结果
      write.csv(do.call(rbind, res), file.path(out_dir, "simulation_results.csv"),
                row.names = FALSE)
    }
  }
  invisible(do.call(rbind, res))
}
