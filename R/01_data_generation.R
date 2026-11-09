###############################################################################
#  01_data_generation.R
#  ---------------------------------------------------------------------------
#  脆弱分布错设下复发事件生存数据的 Plasmode 模拟 —— 数据生成部分
#
#  设计思路（Plasmode）：
#    协变量轮廓取自 survival::bladder2（85 位膀胱癌患者的真实基线协变量），
#    有放回抽样得到 n 位受试者；噪声协变量用 AR(1) 结构生成；
#    事件过程用非齐次泊松过程模拟，强度 λ_i(t) = LAMBDA0 * u_i * exp(X beta)，
#    其中 u_i 为脆弱项（真实分布可为 gamma 或 lognormal），
#    观测窗口 [0, censor_time_i]，censor_time_i ~ U((1-2c) T_MAX, T_MAX)。
#
#  本文件只定义函数与全局常数，不执行任何模拟。
###############################################################################

## ============ 一、全局常数（与正式实验完全一致，不可随意修改）============
# 基础随机种子：第 (sid, rep) 次重复的种子 = BASE_SEED + sid * 1e5 + rep
BASE_SEED <- 20261004

# 三个真实效应参数（对应 bladder2 的 rx / number / size 标准化值）
BETA_TRUE <- c(rx = 1.2, number_z = -1, size_z = 0.8)

# 基线事件强度
LAMBDA0 <- 0.6

# 最长随访时间（年）
T_MAX <- 3.5

# 噪声（无关）协变量个数：真实模型只含 3 个信号变量，其余 27 个系数为 0
N_NOISE <- 27

# 全部 30 个协变量名：前 3 个是信号变量，后 27 个是噪声变量
x_names    <- c("rx", "number_z", "size_z", paste0("z", seq_len(N_NOISE)))
true_names <- c("rx", "number_z", "size_z")
noise_names <- setdiff(x_names, true_names)

## ============ 二、协变量轮廓：从 bladder2 有放回抽样 ============
# 输入 n：需要的受试者人数
# 输出：n 行 3 列（rx, number_z, size_z）外加一列 id
make_profiles <- function(n) {
  d <- survival::bladder2                       # 真实数据源
  ids <- unique(d$id)
  # 每位患者只取一行基线协变量
  p <- d[match(ids, d$id), c("id", "rx", "number", "size")]
  if (nrow(p) != 85) {
    stop("bladder2 中预期有 85 位病人，实际为 ", nrow(p), "。")
  }
  p$rx <- as.integer(p$rx == 2)                 # rx 处理组指示变量转为 0/1
  p$number_z <- as.numeric(scale(p$number))     # 肿瘤数量标准化
  p$size_z <- as.numeric(scale(p$size))         # 肿瘤大小标准化
  idx <- sample.int(nrow(p), size = n, replace = TRUE)  # 有放回抽样
  ans <- p[idx, c("rx", "number_z", "size_z")]
  rownames(ans) <- NULL
  ans$id <- seq_len(n)                          # 重新编号 1..n
  ans
}

## ============ 三、脆弱项抽样 ============
# family = "gamma"     : Gamma(shape = 1/theta, scale = theta)，均值 1、方差 theta
# family = "lognormal" : 对数正态，均值 1、方差 theta
#                        （sig2 = log(1+theta)，meanlog = -sig2/2）
draw_frailty <- function(n, family, theta) {
  if (family == "gamma") {
    return(rgamma(n, shape = 1 / theta, scale = theta))
  }
  if (family == "lognormal") {
    sig2 <- log1p(theta)
    return(rlnorm(n, meanlog = -sig2 / 2, sdlog = sqrt(sig2)))
  }
  stop("未知的真实 frailty 分布：", family)
}

## ============ 四、噪声协变量的 AR(1) 结构 ============
# 生成 n 行 p 列的标准正态矩阵，列与列之间带一阶自相关：
#   x[, j] = rho * x[, j-1] + sqrt(1 - rho^2) * e[, j]
# rho 越大，相邻噪声协变量相关性越强，用于考察高维相关设计下的变量选择
draw_ar1 <- function(n, p, rho) {
  e <- matrix(rnorm(n * p), nrow = n, ncol = p)
  x <- e
  if (p > 1) {
    for (j in 2:p) {
      x[, j] <- rho * x[, j - 1] + sqrt(1 - rho^2) * e[, j]
    }
  }
  x
}

## ============ 五、生成一次复发事件数据集 ============
# 输入：n 受试者数，theta 脆弱方差，censor_prop 删失比例，
#       rho 噪声协变量自相关，true_frailty 真实脆弱分布，seed 随机种子
# 输出：长表（每位受试者若干行，区间形式 start-stop-event），
#       属性 patient_summary 记录每位受试者的删失时间与事件数
make_recurrent_data_formal <- function(n, theta, censor_prop, rho, true_frailty, seed) {
  set.seed(seed)                                            # 固定本次种子

  prof <- make_profiles(n)                                   # 真实协变量轮廓
  z <- draw_ar1(n, N_NOISE, rho)                             # 噪声协变量
  colnames(z) <- noise_names
  x <- cbind(as.matrix(prof[, true_names]), z)               # 拼成 30 列设计矩阵

  beta <- c(BETA_TRUE, rep(0, N_NOISE))                      # 真实系数（后27个为0）

  frailty <- draw_frailty(n, true_frailty, theta)            # 真实脆弱项 u_i

  # 删失时间：删失比例 censor_prop 决定观测窗口的下界
  c_min <- (1 - 2 * censor_prop) * T_MAX
  censor_time <- runif(n, min = c_min, max = T_MAX)          # 每人的随访终点

  # 个体事件强度：λ_i = LAMBDA0 * u_i * exp(X_i beta)
  rate <- LAMBDA0 * frailty * exp(drop(x %*% beta))
  if (any(!is.finite(rate)) || any(rate <= 0)) {
    stop("生成了无效的事件强度。")
  }

  rows <- vector("list", n)
  n_events <- integer(n)
  for (i in seq_len(n)) {
    # 观测窗口内的事件个数 ~ Poisson(rate_i * censor_time_i)
    k <- rpois(1, lambda = rate[i] * censor_time[i])
    n_events[i] <- k
    # 事件发生时刻在 [0, censor_time_i] 上均匀且有序
    event_times <- if (k > 0) sort(runif(k, min = 0, max = censor_time[i])) else numeric(0)
    # 拆成区间形式 (start, stop]，最后一个区间为 (t_k, censor_time] 且 event = 0
    di <- data.frame(id = i,
                     start = c(0, event_times),
                     stop  = c(event_times, censor_time[i]),
                     event = c(rep(1, k), 0))
    for (nm in x_names) di[[nm]] <- x[i, nm]                 # 协变量按人复制
    rows[[i]] <- di
  }

  ans <- do.call(rbind, rows)
  rownames(ans) <- NULL
  attr(ans, "patient_summary") <- data.frame(
    id = seq_len(n), censor_time = censor_time, n_events = n_events)
  ans
}

## ============ 六、拟合公式 ============
# 复发事件专用的计数过程形式：Surv(start, stop, event) + cluster(id)
make_formula <- function() {
  as.formula(paste("Surv(start, stop, event) ~",
                   paste(c(x_names, "cluster(id)"), collapse = " + ")))
}
