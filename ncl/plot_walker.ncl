load "gill.ncl"
begin
  x = lon()
  z = lev()
  p = new((/nx/), "float")
  p!0 = "x"
  p&x = x
  kai = new((/nz, nx/), "float")
  kai!0 = "z"
  kai&z = z
  kai!1 = "x"
  kai&x = x

  wks = gsn_open_wks("x11", "walker")

  resk = True
  resk@vpWidthF = 0.8
  resk@vpHeightF = 0.4
  resk@vpXF = 0.1
  resk@vpYF = 0.8
  resp = resk
  resp@trXMinF = x(0)
  resp@trXMaxF = x(nx-1)

  walker(p, kai, x, z)

  resk@tmYLMode = "Explicit"
  resk@tmYLValues = (/z(0), 0.5 * pi, z(nz-1)/)
  resk@tmYLLabels = (/"0.0", "D/2", "D"/)
  resk@cnLineThicknessF = 3.0
  resk@gsnContourNegLineDashPattern = 2
  resk@gsnContourZeroLineThicknessF = 0.0

  plot = gsn_csm_contour(wks, kai, resk)
  plot = gsn_csm_xy(wks, x, p, resp)

end
