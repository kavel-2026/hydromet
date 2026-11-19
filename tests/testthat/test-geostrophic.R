test_that("CoriolisParameter() follows 2 * omega * sin(lat)", {
  omega <- 7.292115e-5

  expect_equal(CoriolisParameter(0), 0)
  expect_equal(CoriolisParameter(90), 2 * omega)
  expect_equal(CoriolisParameter(-90), -2 * omega)
  expect_equal(CoriolisParameter(30), 2 * omega * sin(pi / 6))
  expect_equal(CoriolisParameter(c(45, -45)), c(1, -1) * 2 * omega * sin(pi / 4))
})

test_that("CoriolisParameter() handles missing and invalid input", {
  expect_true(is.na(CoriolisParameter(NA)))
  expect_error(CoriolisParameter(91), "lat")
  expect_error(CoriolisParameter("45"), "numeric")
})

test_that("GeostrophicWind() solves the geostrophic balance", {
  f <- 1e-4
  rho <- 1.2
  dpdx <- 2e-3
  dpdy <- -3e-3
  wind <- GeostrophicWind(dpdx, dpdy, f, rho)

  expect_equal(wind$u, -dpdy / (rho * f))
  expect_equal(wind$v, dpdx / (rho * f))

  # 平衡方程：f v = dpdx / rho，f u = -dpdy / rho
  expect_equal(f * wind$v, dpdx / rho)
  expect_equal(f * wind$u, -dpdy / rho)

  # 风与水平气压梯度垂直
  expect_equal(wind$u * dpdx + wind$v * dpdy, 0)
})

test_that("GeostrophicWind() reverses with the hemisphere", {
  north <- GeostrophicWind(1e-3, 0, 1e-4)
  south <- GeostrophicWind(1e-3, 0, -1e-4)

  expect_equal(south$v, -north$v)
  expect_equal(north$u, 0)
})

test_that("GeostrophicWind() is vectorised", {
  f <- CoriolisParameter(c(30, 45, 60))
  wind <- GeostrophicWind(rep(1e-3, 3), 0, f)

  expect_equal(nrow(wind), 3)
  expect_equal(wind$v, 1e-3 / (1.2 * f))
})

test_that("GeostrophicWind() rejects f = 0 and invalid density", {
  expect_error(GeostrophicWind(1e-3, 0, 0), "equator")
  expect_error(GeostrophicWind(1e-3, 0, 1e-4, rho = 0), "rho")
})

test_that("WindSpeed() computes the magnitude of the wind vector", {
  expect_equal(WindSpeed(3, 4), 5)
  expect_equal(WindSpeed(c(0, 1), c(0, 1)), c(0, sqrt(2)))
})

test_that("WindDirection() is meteorological, i.e. the direction winds blow from", {
  expect_equal(WindDirection(1, 0), 270) # 西风
  expect_equal(WindDirection(0, 1), 180) # 南风
  expect_equal(WindDirection(-1, 0), 90) # 东风
  expect_equal(WindDirection(0, -1), 0) # 北风
  expect_equal(WindDirection(1, 1), 225) # 西南风
})

test_that("WindDirection() stays in [0, 360) and is NA for calm winds", {
  direction <- WindDirection(c(-1, 1, 0), c(-1, 1, 1))

  expect_true(all(direction >= 0 & direction < 360))
  expect_true(is.na(WindDirection(0, 0)))
})
