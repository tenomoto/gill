load "gill.ncl"
;lsymmetric = True
lsymmetric = False
lboth = False
;lboth = True
begin
  y = lat()
  z = lev()
  p = new((/ny/), "float")
  p!0 = "y"
  p&y = y
  u = new((/nz, ny/), "float")
  u!0 = "z"
  u&z = z
  u!1 = "y"
  u&y = y
  kai = u

  wks = gsn_open_wks("x11", "hadley")

  resu = True
  resu@vpWidthF = 0.8
  resu@vpHeightF = 0.4
  resu@vpXF = 0.1
  resu@vpYF = 0.7
  resu@gsnDraw = False
  resu@gsnFrame = False
  resp = resu

;  hadley_symmetric(p, u, v, w, y, z)
  if (lsymmetric.or.lboth) then
    hadley_symmetric(p, u, kai, y, z)
  else
    hadley_antisymmetric(p, u, kai, y, z)
  end if
  if (lboth) then
    pa = p
    ua = u
    kaia = u
    hadley_antisymmetric(pa, ua, kaia, y, z)
    p = p + pa
    u = u + ua
    kai = kai + kaia
  end if

  resu@tmYLMode = "Explicit"
  resu@tmYLValues = (/z(0), 0.5 * pi, z(nz-1)/)
  resu@tmYLLabels = (/"0.0", "D/2", "D"/)
  resu@cnLineThicknessF = 3.0
  resu@gsnContourNegLineDashPattern = 2
  resu@gsnContourZeroLineThicknessF = 0.0
  resk = resu

  plot = new(3, "graphic")
  plot(0) = gsn_csm_contour(wks, u, resu)
  plot(1) = gsn_csm_contour(wks, kai, resk)
  plot(2) = gsn_csm_xy(wks, y, p, resp)

  pres = True
  pres@gsnPanelFigureStrings = (/"(a)", "(b)", "(c)"/)
  pres@gsnPanelFigureStringsPerimOn = False

  gsn_panel(wks, plot, (/3, 1/), pres)
end
