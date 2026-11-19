# 地转风（geostrophic wind）----------------------------------------------------
# 自由大气中，水平气压梯度力与地转偏向力平衡时的风。
#
#   水平运动方程（f 平面、无摩擦、定常）：
#     -f v = -(1 / rho) dp/dx
#      f u = -(1 / rho) dp/dy
#
#   解得：
#     u = -(1 / (rho f)) dp/dy
#     v =  (1 / (rho f)) dp/dx
#
#   矢量形式：V_g = 1 / (rho f) * k x grad(p)

#' 科里奥利参数
#'
#' 地转偏向力的大小与风速成正比，比例系数即科里奥利参数 f = 2Ω sinφ。
#' 北半球 f > 0，南半球 f < 0；赤道处 f = 0，地转平衡不成立。
#'
#' @param lat 纬度（度），北半球为正、南半球为负。
#' @param omega 地球自转角速度（rad/s），默认 `7.292115e-5`。
#'
#' @return 科里奥利参数（1/s），与 `lat` 等长。
#'
#' @examples
#' CoriolisParameter(c(0, 30, 45, 60))
#' @encoding UTF-8
#' @export
CoriolisParameter <- function(lat, omega = 7.292115e-5) {
  if (!is.numeric(lat)) {
    if (!(is.logical(lat) && all(is.na(lat)))) {
      stop("`lat` must be numeric", call. = FALSE)
    }
    lat <- as.numeric(lat)
  }
  if (any(abs(lat) > 90, na.rm = TRUE)) {
    stop("`lat` must be within [-90, 90] degrees", call. = FALSE)
  }
  2 * omega * sin(lat * pi / 180)
}

#' 地转风分量
#'
#' 在忽略摩擦与加速度的自由大气中，水平气压梯度力与地转偏向力平衡：
#' \deqn{-fv = -\frac{1}{\rho}\frac{\partial p}{\partial x}, \quad
#'       fu = -\frac{1}{\rho}\frac{\partial p}{\partial y}}{%
#'       -f*v = -(1/rho) dp/dx,  f*u = -(1/rho) dp/dy}
#' 由此得到地转风的东西分量 `u` 与南北分量 `v`。
#'
#' @param dpdx 纬向（东西向）气压梯度，向东为正，单位 Pa/m。
#' @param dpdy 经向（南北向）气压梯度，向北为正，单位 Pa/m。
#' @param f 科里奥利参数（1/s），见 [CoriolisParameter()]。
#' @param rho 空气密度（kg/m^3），默认 `1.2`。
#'
#' @return 数据框，含 `u`（向东为正）与 `v`（向北为正）两个分量，单位 m/s。
#'
#' @examples
#' f <- CoriolisParameter(45)
#' # 1 hPa / 100 km 的纬向气压梯度
#' GeostrophicWind(dpdx = 1e-3, dpdy = 0, f = f)
#' @encoding UTF-8
#' @export
GeostrophicWind <- function(dpdx, dpdy, f, rho = 1.2) {
  if (any(f == 0, na.rm = TRUE)) {
    stop("geostrophic wind is undefined where `f = 0` (the equator)", call. = FALSE)
  }
  if (any(rho <= 0, na.rm = TRUE)) {
    stop("`rho` must be positive", call. = FALSE)
  }
  data.frame(
    u = -dpdy / (rho * f),
    v =  dpdx / (rho * f)
  )
}

#' 风速
#'
#' @param u 风的东西分量（m/s），向东为正。
#' @param v 风的南北分量（m/s），向北为正。
#'
#' @return 风速（m/s），即 \eqn{\sqrt{u^2 + v^2}}。
#'
#' @examples
#' WindSpeed(c(3, 0), c(4, 5))
#' @encoding UTF-8
#' @export
WindSpeed <- function(u, v) {
  sqrt(u^2 + v^2)
}

#' 风向
#'
#' 气象上的风向指风的来向：以正北为 0°、顺时针为正，取值范围 `[0, 360)`。
#' 静风（`u = v = 0`）时风向无定义，返回 `NA`。
#'
#' @param u 风的东西分量（m/s），向东为正。
#' @param v 风的南北分量（m/s），向北为正。
#'
#' @return 风向（度），与 `u`、`v` 等长。
#'
#' @examples
#' # 西风、南风、东风
#' WindDirection(c(1, 0, -1), c(0, 1, 0))
#' @encoding UTF-8
#' @export
WindDirection <- function(u, v) {
  direction <- (270 - atan2(v, u) * 180 / pi) %% 360
  direction[which(WindSpeed(u, v) == 0)] <- NA_real_
  direction
}
