#' 去除季节循环
#'
#' 以月日（MMDD）为季节标记，计算每个标记对应的气候态均值，再从观测序列中扣除，
#' 得到去季节序列（anomaly）。常用于消除气温、降水等要素的年循环。
#'
#' @param y 观测序列。
#' @param t_mmdd 与 `y` 等长的月日标记（MMDD，如 `101` 表示 1 月 1 日）。
#'
#' @return 列表：
#' - `mu`：各月日的气候态均值，以 `t_mmdd` 为名；
#' - `anomaly`：去季节序列，即 `y - mu`。
#'
#' @examples
#' y <- c(1, 2, 3, 2, 3, 4)
#' t_mmdd <- c(101, 101, 101, 102, 102, 102)
#' RemoveSeasonal(y, t_mmdd)
#' @encoding UTF-8
#' @export
RemoveSeasonal <- function(y, t_mmdd) {
  n <- length(y)
  if (length(t_mmdd) != n) {
    stop("`y` and `t_mmdd` must have the same length", call. = FALSE)
  }

  inds <- sort(unique(t_mmdd))
  mu <- vapply(inds, function(i) mean(y[t_mmdd == i], na.rm = TRUE), numeric(1))
  names(mu) <- inds
  anomaly <- y - unname(mu[match(t_mmdd, inds)])
  list(mu = mu, anomaly = anomaly)
}

#' 多断点分段回归（同斜率、异截距）
#'
#' 在 `Ips` 指定的位置把序列截成若干连续段，各段共享同一斜率、各自拥有截距，
#' 即“同斜率、不同截距”的分段回归。
#'
#' @param Y 响应值（观测序列）。
#' @param X 自变量（通常为时间序列）。内部中心化为 `X - mean(X)`，
#'   因此截距表示 `X` 均值处的拟合值。
#' @param Ips 各段末尾的位置，末位必须等于序列长度。例如 `c(100, 200)`
#'   表示在 100 处有一个断点、200 为序列末尾。
#' @param ... 传给 [lm_solve()]。
#'
#' @details 设计矩阵为 `[1, X - mean(X), 段指示变量]`：第 i 个指示变量在第 i
#'   段取 1、其余为 0。因此 `coef` 的前两项是第一段的截距与公共斜率，
#'   其后各项是第 i 段相对第一段的截距差。
#'
#' @return 列表，含 `coef`、`fitted`、`resi`、`SSE`，见 [lm_solve()]。
#'
#' @examples
#' set.seed(1)
#' x <- seq_len(120)
#' y <- 3 + 0.1 * x + ifelse(x > 60, 2, 0) + rnorm(120, sd = 0.1)
#' LS_MultiBreaks(y, x, c(60, 120))$coef
#' @encoding UTF-8
#' @export
LS_MultiBreaks <- function(Y, X, Ips, ...) {
  n <- length(Y)
  if (length(X) != n) {
    stop("`Y` and `X` must have the same length", call. = FALSE)
  }
  if (any(Ips < 1) || any(diff(Ips) <= 0) || Ips[length(Ips)] != n) {
    stop("`Ips` must be strictly increasing, start at 1 or later, and end at `length(Y)`", call. = FALSE)
  }

  n_breaks <- length(Ips) - 1
  design <- cbind(1, X - mean(X))
  for (i in seq_len(n_breaks)) {
    segment <- numeric(n)
    segment[(Ips[i] + 1):Ips[i + 1]] <- 1
    design <- cbind(design, segment)
  }
  lm_solve(design, Y, ...)
}

#' 单断点分段回归（同斜率、异截距）
#'
#' 在 `Ic` 处把序列分为两段，两段共享同一斜率、各有截距。
#' 是 [LS_MultiBreaks()] 在只有一个断点时的特例。
#'
#' @param Y 响应值（观测序列）。
#' @param X 自变量（通常为时间序列）。内部中心化为 `X - mean(X)`。
#' @param Ic 断点位置：`[1, Ic]` 为第一段，`[Ic + 1, n]` 为第二段；
#'   取 `NA` 时退化为普通线性回归。
#' @param ... 传给 [lm_solve()]。
#'
#' @return 列表，含 `coef`、`fitted`、`resi`、`SSE`，见 [lm_solve()]。
#'
#' @examples
#' x <- seq_len(120)
#' y <- 3 + 0.1 * x + ifelse(x > 60, 2, 0)
#' LS_OneBreak(y, x, 60)$coef
#' @encoding UTF-8
#' @export
LS_OneBreak <- function(Y, X, Ic, ...) {
  n <- length(Y)
  if (length(X) != n) {
    stop("`Y` and `X` must have the same length", call. = FALSE)
  }
  if (!is.na(Ic) && (Ic < 1 || Ic >= n)) {
    stop("`Ic` must be within [1, length(Y) - 1]", call. = FALSE)
  }

  design <- cbind(1, X - mean(X))
  if (!is.na(Ic)) {
    design <- cbind(design, c(rep(0, Ic), rep(1, n - Ic)))
  }
  lm_solve(design, Y, ...)
}

#' 最小二乘求解
#'
#' 对设计矩阵 `design` 与响应 `response` 做最小二乘拟合，是
#' [LS_OneBreak()] 与 [LS_MultiBreaks()] 共用的求解器。
#'
#' @param design 设计矩阵，行数等于 `length(response)`。
#' @param response 响应向量，不能含缺失值。
#' @param only.SSE 是否只返回残差平方和。
#' @param ... 未使用，保留以兼容调用方。
#'
#' @return 列表：
#' - `coef`：系数向量；
#' - `fitted`：拟合值；
#' - `resi`：残差；
#' - `SSE`：残差平方和。
#'
#' 当 `only.SSE = TRUE` 时只返回 `SSE`。
#'
#' @encoding UTF-8
#' @keywords internal
lm_solve <- function(design, response, only.SSE = FALSE, ...) {
  if (anyNA(design) || anyNA(response)) {
    stop("`design` and `response` must not contain missing values", call. = FALSE)
  }

  fit <- .lm.fit(design, as.matrix(response))
  SSE <- sum(fit$residuals^2)
  if (only.SSE) {
    return(list(SSE = SSE))
  }
  list(
    coef = as.vector(fit$coefficients),
    fitted = as.vector(design %*% fit$coefficients),
    resi = as.vector(fit$residuals),
    SSE = SSE
  )
}
