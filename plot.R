source("gill.R")

drawvec <- function(x, y, u, v, length = 0.05, uvref = NA, ...) {
  nx <- length(x)
  ny <- length(y)
  usr <- par("usr")
  d <- min((usr[2] - usr[1]) / nx, (usr[4] - usr[3]) / ny)
  if (is.na(uvref)) uvref <- max(abs(c(u, v)))
  d <- d / uvref
  u <- 0.5 * u * d
  v <- 0.5 * v * d
  X <- rep(x, times = length(y))
  Y <- rep(y, each  = length(x))
  tol <- 1e-8
  suppressWarnings(
    arrows(X - u, Y - v, X + u, Y + v, length = length, ...)
  )
}

#res <- symmetric(lon, lat)
res <- antisymmetric(lon, lat)
print(sapply(res, function(x) summary(as.numeric(x))))
wlk <- walker(lon, lev)
print(sapply(wlk, function(x) summary(as.numeric(x))))

oldpar <- par(mfrow = c(3, 1), mar = c(4, 4, 2, 1), oma = c(3, 0, 0,0), mgp = c(2, 1, 0))

plot.new()
plot.window(c(-10, 15), c(-4, 4))

.filled.contour(lon, lat, res$w,
  levels = c(min(res$w), -0.1, 0.1, max(res$w)),
  col = c("lightblue", NA, "lightpink"))
clev_p <- pretty(res$p, 10)
if (any(clev_p %in% 0)) clev_p <- clev_p[-which(clev_p == 0)]
contour(lon, lat, res$p,
  levels = clev_p, add = TRUE)
ivec <- seq(1, nx, 5)
jvec <- seq(1, ny, 5)
drawvec(lon[ivec], lat[jvec], res$u[ivec, jvec], res$v[ivec, jvec],
  uvref = 1, angle = 30, length = 0.02)

axis(1)
axis(2)
title("vector: (u, v), shade: w, contour: p", "", "longitude", "latitude")
box()

plot.new()
plot.window(c(-10, 15), c(0, pi))

clev_kai <- pretty(wlk$kai, 10)
if (any(clev_kai %in% 0)) clev_kai <- clev_kai[-which(clev_kai == 0)]
lty_kai <- rep(1, length(clev_kai))
lty_kai[clev_kai < 0] <- 2
contour(lon, lev, wlk$kai, levels = clev_kai, lty = lty_kai,
        add = TRUE)

axis(1)
axis(2, c(0, pi / 2, pi), c("0", "D/2", "D"))
title("meridionally averaged streamfunction", "", "longitude", "height")
box()

plot.new()
plot.window(c(-10, 15), c(-10, 0))

lines(lon, wlk$p)

axis(1)
axis(2)
title("meridionally averaged pressure", "", "longitude", "pressure")
box()

par(oldpar)
