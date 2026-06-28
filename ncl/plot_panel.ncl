lsymmetric = True
;lsymmetric = False
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
  if (lsymmetric) then
    symmetric(p, u, v, w, x, y)
  else
    antisymmetric(p, u, v, w, x, y)
  end if

  wks = gsn_open_wks("x11", "gill_wp")

  resp = True
  resp@vpWidthF = 0.8
  resp@vpHeightF = 0.4
  resp@vpXF = 0.1
  resp@vpYF = 0.7
  resp@gsnDraw = False
  resp@gsnFrame = False
  resv = resp
  resp@cnLineThicknessF = 3.0
  resp@cnLevelSelectionMode = "ExplicitLevels"
  resp@cnInfoLabelOn = False
  resw = resp

;  resw@gsnContourZeroLineThicknessF = 2
;  resw@cnLevels = (/-0.6, 0.3, -0.1, 0.0, 0.1, 0.3, 0.6/)
  resw@gsnContourNegLineDashPattern = 2
  resw@cnLevels = (/-0.6, 0.3, -0.1, 0.1, 0.3, 0.6/)
  plotw = gsn_csm_contour(wks, w, resw)

  resp@gsnContourNegLineDashPattern = 0
  resp@cnLevels = (/-1.5, -1.2, -0.9, -0.6, -0.3, 0.3, 0.6, 0.9, 1.2, 1.5/)
  resp@cnFillOn = False
  plotp = gsn_csm_contour(wks, p, resp)

  resv@vcRefAnnoOrthogonalPosF = -1.15
  resv@vcRefMagnitudeF = 1.0
  resv@vcRefLengthF = 0.02
  resv@vcMinDistanceF = 0.02
  resv@vcRefAnnoString2On = False
;  resv@vcMinMagnitudeF = 0.3

  plot = new(2, "graphic")

  plot(0) = gsn_csm_vector(wks, u, v, resv)
  overlay(plot(0), plotw)

  plot(1) = gsn_csm_vector(wks, u, v, resv)
  overlay(plot(1), plotp)
  
  pres = True

  gsn_panel(wks, plot, (/2, 1/), pres)
end 
