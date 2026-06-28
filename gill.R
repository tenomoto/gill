L <- 2.0
k <- pi / (2.0 * L)
eps <- 0.1
nx <- 251
ny <-  81
nz <-  30
xmin <- -10.0
ymin <-  -4.0
ymax <-   4.0
xmax <-  15.0
zmin <-   0.0
zmax <-    pi
I <- 4.0 * L / pi

lon <- seq(xmin, xmax, length.out = nx)
lat <- seq(ymin, ymax, length.out = ny)
lev <- seq(zmin, zmax, length.out = nz)

q0 <- function(x) {
  ekr <- 1.0 / (eps^2 + k^2)
  kx <- k * x
  ifelse(x <= L, 
    ifelse(x <= -L, 0.0,
      -ekr * (eps * cos(kx) + k * (sin(kx) + exp(-eps * (x + L))))
    ), -ekr * k * (1.0 + exp(-2.0 * eps * L)) * exp(eps * (L - x))
  )
}

qn <- function(n, x) {
  nn <- 2 * n[1] + 1
  ekr <- 1.0 / ((nn * eps)^2 + k^2)
  kx <- k * x
  ifelse(x <= L,
    ifelse (x <= -L,
      -ekr * k * (1.0 + exp(-2 * nn * eps * L)) * exp(nn * eps * (x + L)),
       ekr * (-nn * eps * cos(kx) + k * (sin(kx) - exp(nn * eps * (x - L))))
    ), 0.0
  )
}

ffunc <- function(x) {
  ifelse(abs(x) < L, cos(k * x), 0.0)
}

symmetric <- function(x, y) {
  q0x <- q0(x)
  q2x <- qn(1, x)
  Fx <- ffunc(x)

  y2 <- y^2
  ey2 <- exp(-0.25 * y2)

  p <- 0.5 * (outer(q0x, ey2) + outer(q2x, (1 + y2) * ey2))
  u <- 0.5 * (outer(q0x, ey2) + outer(q2x, (y2 - 3) * ey2))
  v <- outer(Fx + 4 * eps * q2x, y * ey2)
  w <- outer(0.5 * eps * q0x + Fx, ey2) + 0.5 * eps * outer(q2x, (1 + y2) * ey2)

  list(p = p, u = u, v = v, w = w)
}

antisymmetric <- function(x, y) {
  q3x <- qn(2, x)
  Fx <- ffunc(x)

  y2 <- y^2
  y3 <- y2 * y
  ey2 <- exp(-0.25 * y2)

  p <- 0.5 * outer(q3x, y3 * ey2)
  u <- 0.5 * outer(q3x, (y3 - 6 * y) * ey2)
  v <- outer(6 * eps * q3x,  (y2 - 1) * ey2) + outer(Fx, y2 * ey2)
  w <- 0.5 * eps * outer(q3x, y3 * ey2) + outer(Fx,  y * ey2)

  list(p = p, u = u, v = v, w = w)
}

# kai: clockwise + (meteorology, oceanography as opposed to fluid dynamics)
hadley_symmetric <- function(y, z) {
  y2 <- y^2
  eyy <- I * exp(-0.25 * y2)

  p <- -(4 + y2) / (6 * eps) * eyy

  u <- outer(-y2 / (6 * eps) * eyy, cos(z))
  kai <- outer(y / 3 * eyy, sin(z))

  list(p = p, u = u, kai = kai)
}

hadley_antisymmetric <- function(y, z) {
  y2 <- y^2
  y3 <- y2 * y
  eyy <- I * exp(-0.25 * y2)

  p <- -y3 / (10 * eps) * eyy
  u <- outer((6 * y - y3) / (10 * eps) * eyy, cos(z))
  kai <- outer((y2 - 6) / 5 * eyy, sin(z))

  list(p = p, u = u, kai = kai)
}

walker <- function(x, z) {
  q0x <- q0(x)
  q2x <- qn(1, x)

  p <- (q0x + 3 * q2x) * sqrt(pi)
  kai <- outer((q2x - q0x) * sqrt(pi), sin(z))

  list(p = p, kai = kai)
}

