pro settings
  !except=0
end

function n, var
  return, float(n_elements(var)-1)
end



function n1, var,offset
  if isa(offset) eq 1 then return, float(n_elements(var))+offset
  if isa(offset) eq 0 then return, float(n_elements(var))
end



function xLog, aIn,bIn,count
    a=float(min([aIn,bIn]))
    b=float(max([aIn,bIn]))
    vals=FINDGEN(count)*(ALOG10(b)-ALOG10(a))/(count-1)+ALOG10(a)
    vals=10.^vals
    return, vals
end  



function med, var
  varSort=sort(var)
  varSorted=var[varSort]
return,median(varSorted,/double,/even)
end



function q1, var
  varSort=sort(var)
  varSorted=var[varSort]
  return,varSorted[n_elements(var)*.25]
end



function q3, var
  varSort=sort(var)
  varSorted=var[varSort]
  return,varSorted[n_elements(var)*.75]
end



function f, var
  return, var[where(finite(var) eq 1)]
end



function nan, var
  return, replicate(!values.d_nan,var)
end



function srt, var
  varSort=sort(var)  
  return, var[varSort]
end



function hist, var, binsize=binsize, minval=minval  
  if isa(minval) eq 1 then var=var[where(var ge minval)]
  
  sortvar=sort(var)
  var=var[sortvar]

  if isa(binsize) eq 0 then begin
    h1=histogram(var,locations=inds)
  endif else begin
    h1=histogram(var,binsize=binsize,locations=inds)
  endelse
  
  means=[]
  for i=0,n(inds)-1 do begin
    means=[means,mean([[inds[i]],[inds[i+1]]])]
  endfor
  
  b1=barplot(dindgen(n1(h1)),h1/n1(var),dimensions=[1800,1000],histogram=1,margin=[130,100,50,20],/device,font_size=22)
  b1.xrange=[0,n1(h1)]
  b1.xmajor=n1(h1)+1
  b1.ymajor=8
  b1.xtickname=[string(inds,format='(f0.1)'),string(max(inds)+(inds[2]-inds[1]),format='(f0.1)'),' ',' ']
  b1.xtext_orientation=90
  b1.xticklen=0
  return,b1
end  




function lp, xVar,yVar,op

  if isa(op) eq 0 then begin
    if max(xVar) ne 0 then begin
      p1=plot(xvar,yvar,dimensions=[1600,1200],margin=[150,100,25,25],/device)    
    endif else begin
      p1=plot(yvar,dimensions=[1600,1200],margin=[150,100,25,25],/device)
    endelse
  endif else begin
    if max(xVar) ne 0 then begin
      p1=plot(xvar,yvar,dimensions=[1600,1200],margin=[150,100,25,25],/device,/overplot)
    endif else begin
      p1=plot(yvar,dimensions=[1600,1200],margin=[150,100,25,25],/device,/overplot)
    endelse
  endelse  
  
  p1.font_size=22
  p1.thick=2
  
  return,p1
end





;-----Kevin Ivory from the Max Plank Institute-----
FUNCTION CIRCLE, xcenter, ycenter, radius
  points = (2 * !PI / 99.0) * FINDGEN(100)
  x = xcenter + radius * COS(points )
  y = ycenter + radius * SIN(points )
  RETURN, TRANSPOSE([[x],[y]])
END




;-----Baard Krane (bard.krane@fys.uio.no) at the University of Oslo------
FUNCTION Inside, x, y, px, py

  ;  x - The x coordinate of the point.
  ;  y - The y coordinate of the point.
  ; px - The x coordinates of the polygon.
  ; py - The y coordinates of the polygon.
  ;
  ; The return value of the function is 1 if the point is inside the
  ; polygon and 0 if it is outside the polygon.

  sx = Size(px)
  sy = Size(py)
  IF (sx[0] EQ 1) THEN NX=sx[1] ELSE RETURN, -1    ; Error if px not a vector
  IF (sy[0] EQ 1) THEN NY=sy[1] ELSE RETURN, -1    ; Error if py not a vector
  IF (NX EQ NY) THEN N = NX ELSE RETURN, -1        ; Incompatible dimensions

  tmp_px = [px, px[0]]                             ; Close Polygon in x
  tmp_py = [py, py[0]]                             ; Close Polygon in y

  i = indgen(N)                                    ; Counter (0:NX-1)
  ip = indgen(N)+1                                 ; Counter (1:nx)

  X1 = tmp_px(i)  - x
  Y1 = tmp_py(i)  - y
  X2 = tmp_px(ip) - x
  Y2 = tmp_py(ip) - y

  dp = X1*X2 + Y1*Y2                               ; Dot-product
  cp = X1*Y2 - Y1*X2                               ; Cross-product
  theta = Atan(cp,dp)

  IF (Abs(Total(theta)) GT !PI) THEN RETURN, 1 ELSE RETURN, 0
END