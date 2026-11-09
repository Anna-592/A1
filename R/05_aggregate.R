###############################################################################
#  05_aggregate.R
#  ---------------------------------------------------------------------------
#  汇总 simulation_results.csv，输出与论文表格一致的统计量
#
#  指标定义：
#    TPR           真信号变量被选中的平均比例（理想值 1）
#    FDP           FDR 的实际发生值 = 假阳性个数 / 总入选个数
#                  （sum(noise_names %in% selected) / max(1, n_selected) 的均值）
#    size          平均入选变量个数（理想值 3，恰好是 3 个真实信号）
#    exactRec%     4 个模型中「恰好选中 3 个变量且全部正确」的比例
#
#  另输出配对比较：同一数据集上「错设脆弱分布」减「正确脆弱分布」的差异，
#  这是本研究核心结论的直接证据。
###############################################################################

res <- read.csv("results/simulation_results.csv", stringsAsFactors = FALSE)

## ============ 1. 主表：按 真实分布 × 假定分布 × 惩罚 汇总 ============
agg <- aggregate(
  cbind(TPR, FDR, n_selected, exact = as.integer(n_selected == 3)) ~
    true_frailty + assumed_frailty + penalty,
  data = res,
  FUN = mean
)
names(agg)[names(agg) == "n_selected"] <- "size"
names(agg)[names(agg) == "exact"] <- "exactRec"
agg$exactRec <- agg$exactRec * 100

cat("\n=== 主表：TPR / FDP / 入选个数 / 精确恢复率 ===\n")
print(agg, row.names = FALSE, digits = 4)

## ============ 2. 配对比较：错设 减 正确 ============
# 同一 (scenario_id, rep_id, penalty) 上，两种假定的差分
wide <- reshape(
  res[, c("scenario_id", "rep_id", "true_frailty", "assumed_frailty",
          "penalty", "TPR", "FDR", "n_selected")],
  idvar = c("scenario_id", "rep_id", "penalty", "true_frailty"),
  timevar = "assumed_frailty", direction = "wide"
)
# 只保留「真实分布」与「两种假定」都齐全的记录
ok <- complete.cases(wide)
wide <- wide[ok, ]

# 若真实为 gamma，则正确假定是 gamma、错设为 lognormal；反之亦然
wide$dTPR   <- ifelse(wide$true_frailty == "gamma",
                      wide$TPR.lognormal - wide$TPR.gamma,
                      wide$TPR.gamma - wide$TPR.lognormal)
wide$dFDP   <- ifelse(wide$true_frailty == "gamma",
                      wide$FDR.lognormal - wide$FDR.gamma,
                      wide$FDR.gamma - wide$FDR.lognormal)
wide$dSize  <- ifelse(wide$true_frailty == "gamma",
                      wide$n_selected.lognormal - wide$n_selected.gamma,
                      wide$n_selected.gamma - wide$n_selected.lognormal)

paired <- aggregate(cbind(dTPR, dFDP, dSize) ~ true_frailty + penalty,
                    data = wide, FUN = mean)
names(paired)[2] <- "penalty"
paired$n_pairs <- as.integer(table(paste(wide$true_frailty, wide$penalty))[
  paste(paired$true_frailty, paired$penalty)])

cat("\n=== 配对比较：错设 减 正确（同数据集） ===\n")
print(paired, row.names = FALSE, digits = 4)

## ============ 3. 按情景细看（可选） ============
by_scenario <- aggregate(
  cbind(TPR, FDR, n_selected) ~ scenario_id + true_frailty + assumed_frailty +
    penalty + n + theta + censor_prop + rho,
  data = res, FUN = mean
)
write.csv(by_scenario, "results/aggregate_by_scenario.csv", row.names = FALSE)
write.csv(agg, "results/aggregate_main.csv", row.names = FALSE)
write.csv(paired, "results/aggregate_paired.csv", row.names = FALSE)

cat("\n已输出：results/aggregate_main.csv, aggregate_paired.csv, aggregate_by_scenario.csv\n")
