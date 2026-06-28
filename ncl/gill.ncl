pi = acos(-1.0)
pir = 1.0/pi
L = 2.0
k = pi / (2.0 * L)
eps = 0.1
nx = 251
ny =  81
nz =  30
xmin = -10.0
ymin =  -4.0
ymax =   4.0
xmax =  15.0
zmin =   0.0
zmax =    pi
I = 4.0 * L * pir

function lon()
begin
  return fspan(xmin, xmax, nx)
end

function lat()
begin
  return fspan(ymin, ymax, ny)
end

function lev()
begin
  return fspan(zmin, zmax, nz)
end

function q0(x)
local ekr, kx
begin
  ekr = 1.0 / (eps * eps + k * k)
  kx = k * x
  if (x .le. L) then
    if (x .le. -L) then
      return 0.0
    else
      return -ekr * (eps * cos(kx) + k * (sin(kx) + exp(-eps * (x + L))))
    end if
  else
    return -ekr * k * (1.0 + exp(-2.0 * eps * L)) * exp(eps * (L - x))
  end if
end

function qn(n,x)
local ekr, kx
begin
  nn = 2 * n + 1
  ekr = 1.0 / (nn * nn * eps * eps + k * k)
  kx = k * x
  if (x .le. L) then
    if (x .le. -L) then
      return -ekr * k * (1.0 + exp(-2 * nn * eps * L)) * exp(nn * eps * (x + L))
    else
      return ekr * (-nn * eps * cos(kx) + k * (sin(kx) - exp(nn * eps * (x - L))))
    end if
  else
    return 0.0
  end if
end

function F(x)
begin
  if (fabs(x) .lt. L) then
    return cos(k * x)
  else
    return 0.0
  end if
end

procedure symmetric(p, u, v, w, x, y)
local i, j, q0x, q2x
begin
  do j = 0, ny - 1
    do i = 0, nx - 1
      q0x = q0(x(i))
      q2x = qn(1, x(i))
      Fx = F(x(i))
      ey2 = exp(-0.25 * y(j)^2)
      p(j, i) = 0.5 * (q0x + q2x * (1.0 + y(j) * y(j))) * ey2
      u(j, i) = 0.5 * (q0x + q2x * (y(j) * y(j) - 3.0)) * ey2
      v(j, i) = (Fx + 4.0 * eps * q2x) * y(j) * ey2
      w(j, i) = (0.5 * eps * q0x + Fx + 0.5 * eps * q2x * (1.0 + y(j) * y(j))) * ey2
    end do
  end do
end

procedure antisymmetric(p, u, v, w, x, y)
local i, j, q3x, y3, y2, ey2
begin
  do j = 0, ny - 1
    y2 = y(j) * y(j)
    y3 = y2 * y(j)
    do i = 0, nx - 1
      q3x = qn(2,x(i))
      Fx = F(x(i))
      ey2 = exp(-0.25 * y(j)^2)
      p(j, i) = 0.5 * q3x * y3 * ey2
      u(j, i) = 0.5 * q3x *(y3 - 6 * y(j)) * ey2
      v(j, i) = (6.0 * eps * q3x * (y2 - 1.0) + Fx * y2) * ey2
      w(j, i) = (0.5 * eps * q3x * y3 + Fx * y(j)) * ey2
    end do
  end do
end

;procedure hadley_symmetric(p, u, v, w, y, z)
; kai: clockwise + (meteorology, oceanography as opposed to fluid dynamics)
procedure hadley_symmetric(p, u, kai, y, z)
local j, k, yy, eyy
begin
  do j = 0, ny-1
    yy = y(j) * y(j)
    eyy = I * exp(-0.25 * yy)
    p(j) = -(4.0 + yy) / (6.0 * eps) * eyy
    do k = 0, nz-1
      u(k, j) = -yy / (6.0 * eps) * eyy * cos(z(k))
;      v(k, j) = -y(j) / 3.0 * eyy * cos(z(k))
;      w(k, j) = -(2.0 - yy) / 6.0 * eyy * sin(z(k))
      kai(k, j) = y(j) / 3.0 * eyy * sin(z(k))
    end do
  end do
end

procedure hadley_antisymmetric(p, u, kai, y, z)
local j, k, yy, eyy
begin
  do j = 0, ny-1
    yy = y(j) * y(j)
    yyy = yy * y(j)
    eyy = I * exp(-0.25 * yy)
    p(j) = -yyy / (10.0 * eps) * eyy
    do k = 0, nz-1
      u(k, j) = (6.0 * y(j) - yyy) / (10.0 * eps) * eyy * cos(z(k))
      kai(k, j) = (yy - 6.0) / 5.0 * eyy * sin(z(k))
    end do
  end do
end

procedure walker(p, kai, x, z)
local i, k, q0x, q2x, sqrtpi
begin
  sqrtpi = sqrt(pi)
  do i = 0, nx-1
    q0x = q0(x(i))
    q2x = qn(1, x(i))
    p(i) = (q0x + 3.0 * q2x) * sqrtpi
    do k = 0, nz-1
      kai(k, i) = (q2x - q0x) * sqrtpi * sin(z(k))
    end do
  end do
end
