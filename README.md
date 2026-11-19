
# hydromet

<!-- badges: start -->
[![R-CMD-check](https://github.com/kavel-2026/hydromet/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/kavel-2026/hydromet/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/gh/kavel-2026/hydromet/branch/master/graph/badge.svg)](https://app.codecov.io/gh/kavel-2026/hydromet)
<!-- badges: end -->

`hydromet` 是水文气象学课程作业的 R 包，提供三类函数：

- 去季节循环：[`RemoveSeasonal()`](R/Component.R)；
- 同斜率、异截距的分段回归：[`LS_OneBreak()`](R/Component.R)、[`LS_MultiBreaks()`](R/Component.R)；
- 地转风：[`CoriolisParameter()`](R/geostrophic.R)、[`GeostrophicWind()`](R/geostrophic.R)、
  [`WindSpeed()`](R/geostrophic.R)、[`WindDirection()`](R/geostrophic.R)。

理论的推导与用法说明见 [`docs/note_01.typ`](docs/note_01.typ)。

## 安装

``` r
remotes::install_github("kavel-2026/hydromet")
```

## 示例

``` r
library(hydromet)

# 去掉年循环
RemoveSeasonal(c(1, 3, 10, 12, 100), c(101, 101, 102, 102, 103))

# 分段回归：在 60 处有一个断点
x <- seq_len(120)
y <- 3 + 0.1 * x + ifelse(x > 60, 2, 0)
LS_OneBreak(y, x, 60)$coef

# 地转风：45°N，1 hPa / 100 km 的纬向气压梯度
f <- CoriolisParameter(45)
GeostrophicWind(dpdx = 1e-3, dpdy = 0, f = f)
```
