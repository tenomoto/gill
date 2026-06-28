source("gill.R")

oldpar <- par(mfrow = c(3, 1), mar = c(4, 4, 2, 1), oma = c(3, 0, 0,0), mgp = c(2, 1, 0))

res <- hadley_symmetric(lat, lev) 
#res <- hadley_antisymmetric(lat, lev) 
print(sapply(res, function(x) summary(as.numeric(x))))

plot.new()
plot.window(c(-4, 4), c(0, pi))

clev_u <- pretty(res$u, 10)
if (any(clev_u %in% 0)) clev_u <- clev_u[-which(clev_u == 0)]
lty_u <- rep(1, length(clev_u))
lty_u[clev_u < 0] <- 2

contour(lat, lev, res$u, levels = clev_u, lty = lty_u, add = TRUE)

axis(1)
axis(2, c(0, pi / 2, pi), c("0", "D/2", "D"))
title("Hadley circulation: zonally averaged zonal wind", "", "latitude", "height")
box()

plot.new()
plot.window(c(-4, 4), c(0, pi))

clev_kai <- pretty(res$kai, 10)
if (any(clev_kai %in% 0)) clev_kai <- clev_kai[-which(clev_kai == 0)]
lty_kai <- rep(1, length(clev_kai))
lty_kai[clev_kai < 0] <- 2

contour(lat, lev, res$kai, levels = clev_kai, lty = lty_kai, add = TRUE)

axis(1)
axis(2, c(0, pi / 2, pi), c("0", "D/2", "D"))
title("Hadley circulation: zonally averaged streamfunction", "", "latitude", "height")
box()

plot.new()
plot.window(c(-4, 4), range(res$p))

lines(lat, res$p)

axis(1)
axis(2)
title("Hadley circulation: zonally averaged pressure", "", "latitude", "pressure")
box()

par(oldpar)
