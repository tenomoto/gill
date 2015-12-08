lsymmetric = True
;lsymmetric = False
lboth = False
;lboth = True
load "gill.ncl"

begin
  x = lon()
  y = lat()
  p = new((/ny, nx/), "float")
  p!0 = "y"
  p&y = y
  p!1 = "x"
  p&x = x
  u = p
  v = p
  w = p
  if (lsymmetric.or.lboth) then
    symmetric(p, u, v, w, x, y)
  else
    antisymmetric(p, u, v, w, x, y)
  end if
  if (lboth) then
    pa = p
    ua = p
    va = p
    wa = p
    antisymmetric(pa, ua, va, wa, x, y)
    p = p + pa
    u = u + ua
    v = v + va
    w = w + wa
  end if

  wks = gsn_open_wks("x11", "gill")
;  gsn_define_colormap(wks, "BlAqGrWh2YeOrReVi22")
  gsn_define_colormap(wks, "gs16")

  resp = True
  resp@vpWidthF = 0.8
  resp@vpHeightF = 0.4
  resp@vpXF = 0.1
  resp@vpYF = 0.7
  resp@gsnDraw = False
  resp@gsnFrame = False
  resv = resp
  resp@cnLevelSelectionMode = "ExplicitLevels"
  resp@cnInfoLabelOn = False
  resw = resp

;  resw@cnLevels = (/      -0.3, -0.1, 0.1, 0.3, 0.6, 0.9/)
;  resw@cnFillColors = (/ 6,    4,   -1,  8,  9,  10,  11/)
  resw@cnLevels = (/      -0.3, -0.1, 0.1, 0.3, 0.6, 0.9/)
  resw@cnFillColors = (/ 4,    4,   -1,  9,  9,  9,  9/)
  resw@cnFillOn = True
;  resw@cnLinesOn = False
  resw@lbLabelFontHeightF = 0.02
  plotw = gsn_csm_contour(wks, w, resw)

  resp@cnLevels = (/-1.5, -1.2, -0.9, -0.6, -0.3, 0.3, 0.6, 0.9, 1.2, 1.5/)
  resp@cnFillOn = False
  resp@cnLineThicknessF = 3
  resp@cnLineLabelBackgroundColor = -1
  plotp = gsn_csm_contour(wks, p, resp)

  resv@vcRefAnnoOrthogonalPosF = -1.45
  resv@vcRefMagnitudeF = 1.0
  resv@vcRefLengthF = 0.02
; Causes Segmentation Fault
;  resv@vcGlyphStyle = "CurlyVector"
  resv@vcGlyphStyle = "FillArrow"
  resv@vcMinDistanceF = 0.02
  resv@vcRefAnnoString2On = False
;  resv@vcLineArrowThicknessF = 2.0

;  res@cnLevelSelectionMode = "AutomaticLevels"
;  res@gsnContourNegLineDashPattern = 2
;  plot = gsn_csm_contour(wks, u, res)
;  plot = gsn_csm_contour(wks, v, res)
  plot = gsn_csm_vector(wks, u, v, resv)
  overlay(plotw, plotp)
  overlay(plotw, plot)
  
  resp = True

  draw(plotw)
  frame(wks)
end 
