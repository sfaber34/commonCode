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
  var=var[where(finite(var) eq 1)]
  varSort=sort(var)
  varSorted=var[varSort]
  varSorted=varSorted[0:ceil(.5*n(var))]
  return,mean(varSorted)
end



function q3, var
  var=var[where(finite(var) eq 1)]
  varSort=sort(var)
  varSorted=var[varSort]
  varSorted=varSorted[ceil(.5*n1(var)):ceil(n(var))]
  return,mean(varSorted)
end



function q95, var
  var=var[where(finite(var) eq 1)]
  varSort=sort(var)
  varSorted=var[varSort]
  varSorted=varSorted[ceil(n1(var)*.95)]
  return,mean(varSorted)
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



function hist, var, size=size, min=min
  split=0
  if isa(size) eq 1 and isa(min) eq 0 then split=1
  if isa(size) eq 1 and isa(min) eq 1 then split=2

  case split of
    0: h1=HISTOGRAM(var,LOCATIONS=locs)
    1: h1=HISTOGRAM(var,binsize=size,LOCATIONS=locs)
    2: h1=HISTOGRAM(var,binsize=size,min=min,LOCATIONS=locs)
  endcase


  p1=BARPLOT(locs, h1/n1(var),histogram=1,dimension=[1000,800])
  return,p1
end




function cdphist, var, small=small
  
  ;binEdgesB=[1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,10.5,14.,17.,19.,21.,23.,25.,27.,29.,31.,33.,35.,37.,39.,41.,43.,45.,47.,49.,50.]
  binEdges=[indgen(7,start=2),indgen(21, start=10,increment=2)]
  
  widths=shift(binEdges,-1)-binedges
  widths=widths[0:n1(widths,-2)]
  
  freq=var 
  var=freq/total(freq)
  
  x=indgen(27)
  if isa(small) eq 0 then begin
    b1=barplot(x,var,dimensions=[1000,800],histogram=1,margin=[130,70,50,60],/device,font_size=22)
    b1.xtickfont_size=19
  endif else begin
    b1=barplot(x,var,dimensions=[600,500],histogram=1,margin=[70,50,30,40],/device,font_size=14) 
  endelse
  
  b1.ymajor=10
  b1.xrange=[0,27]
  b1.xtickname=[string(binEdges)]
  b1.xtext_orientation=90
  b1.xminor=0
  b1.xtickfont_style=1
  
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




;-----Harris Geospatial-----
function WHEREB, Array_expression, Count, XIND=xind, YIND=yind,ZIND=zind
  ; works for 1, 2 or 3 dimensional arrays
  ;
  ; Returns the 1D indices (same as WHERE)
  ;
  ; ARGUMENTS
  ;  - same as WHERE (see WHERE)
  ;
  ; KEYWORDS
  ; - Optionally returns X, Y, Z locations through:
  ;
  ; XIND: Output keyword, array of locations along the first dimension
  ; YIND: Output keyword, array of locations along the second dimension (if present)
  ; ZIND: Output keyword, array of locations along the third dimension (if present)
  ;
  ; If no matches where found, then XIND returns -1
  ;
  index_array=where(Array_expression, Count)
  dims=size(Array_expression,/dim)
  xind=index_array mod dims[0]
  case n_elements(dims) of
    2: yind=index_array / dims[0]
    3: begin
      yind=index_array / dims[0] mod dims[1]
      zind=index_array / dims[0] / dims[1]
    end
    else:
  endcase
  return, index_array
end




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




pro match, a, b, suba, subb, COUNT = count, SORT = sort, epsilon=epsilon
  ;+
  ; NAME:
  ;       MATCH
  ; PURPOSE:
  ;       Routine to match values in two vectors.
  ;
  ; CALLING SEQUENCE:
  ;       match, a, b, suba, subb, [ COUNT =, /SORT, EPSILON =  ]
  ;
  ; INPUTS:
  ;       a,b - two vectors to match elements, numeric or string data types
  ;
  ; OUTPUTS:
  ;       suba - subscripts of elements in vector a with a match
  ;               in vector b
  ;       subb - subscripts of the positions of the elements in
  ;               vector b with matchs in vector a.
  ;
  ;       suba and subb are ordered such that a[suba] equals b[subb]
  ;
  ; OPTIONAL INPUT KEYWORD:
  ;       /SORT - By default, MATCH uses two different algorithm: (1) the
  ;               /REVERSE_INDICES keyword to HISTOGRAM is used for integer data,
  ;               while (2) a sorting algorithm is used for non-integer data.  The
  ;               histogram algorithm is usually faster, except when the input
  ;               vectors are sparse and contain very large numbers, possibly
  ;               causing memory problems.   Use the /SORT keyword to always use
  ;               the sort algorithm.
  ;       epsilon - if values are within epsilon, they are considered equal. Used only
  ;               only for non-integer matching.  Note that input vectors should
  ;               be unique to within epsilon to provide one-to-one mapping..
  ;               Default=0.
  ;
  ; OPTIONAL KEYWORD OUTPUT:
  ;       COUNT - set to the number of matches, integer scalar
  ;
  ; SIDE EFFECTS:
  ;       The obsolete system variable !ERR is set to the number of matches;
  ;       however, the use !ERR is deprecated in favor of the COUNT keyword
  ;
  ; RESTRICTIONS:
  ;       The vectors a and b should not have duplicate values within them.
  ;       You can use rem_dup function to remove duplicate values
  ;       in a vector
  ;
  ; EXAMPLE:
  ;       If a = [3,5,7,9,11]   & b = [5,6,7,8,9,10]
  ;       then
  ;               IDL> match, a, b, suba, subb, COUNT = count
  ;
  ;       will give suba = [1,2,3], subb = [0,2,4],  COUNT = 3
  ;       and       a[suba] = b[subb] = [5,7,9]
  ;
  ;
  ; METHOD:
  ;       For non-integer data types, the two input vectors are combined and
  ;       sorted and the consecutive equal elements are identified.   For integer
  ;       data types, the /REVERSE_INDICES keyword to HISTOGRAM of each array
  ;       is used to identify where the two arrays have elements in common.
  ; HISTORY:
  ;       D. Lindler  Mar. 1986.
  ;       Fixed "indgen" call for very large arrays   W. Landsman  Sep 1991
  ;       Added COUNT keyword    W. Landsman   Sep. 1992
  ;       Fixed case where single element array supplied   W. Landsman Aug 95
  ;       Use a HISTOGRAM algorithm for integer vector inputs for improved
  ;             performance                W. Landsman         March 2000
  ;       Work again for strings           W. Landsman         April 2000
  ;       Use size(/type)                  W. Landsman         December 2002
  ;       Work for scalar integer input    W. Landsman         June 2003
  ;       Assume since V5.4, use COMPLEMENT to WHERE() W. Landsman Apr 2006
  ;       Added epsilon keyword            Kim Tolbert         March 14, 2008
  ;-
  ;-------------------------------------------------------------------------
  On_error,2
  compile_opt idl2

  if N_elements(epsilon) EQ 0 then epsilon = 0

  if N_params() LT 3 then begin
    print,'Syntax - match, a, b, suba, subb, [ COUNT =, EPSILON=, /SORT]'
    print,'    a,b -- input vectors for which to match elements'
    print,'    suba,subb -- output subscript vectors of matched elements'
    return
  endif

  da = size(a,/type) & db =size(b,/type)
  if keyword_set(sort) then hist = 0b else $
    hist = (( da LE 3 ) || (da GE 12)) &&  ((db LE 3) || (db GE 12 ))

  if ~hist then begin           ;Non-integer calculation

    na = N_elements(a)              ;number of elements in a
    nb = N_elements(b)             ;number of elements in b

    ; Check for a single element array

    if (na EQ 1) || (nb EQ 1) then begin
      if (nb GT 1) then begin
        subb = where(b EQ a[0], nw)
        if (nw GT 0) then suba = replicate(0,nw) else suba = [-1]
      endif else begin
        suba = where(a EQ b[0], nw)
        if (nw GT 0) then subb = replicate(0,nw) else subb = [-1]
      endelse
      count = nw
      return
    endif

    c = [ a, b ]                   ;combined list of a and b
    ind = [ lindgen(na), lindgen(nb) ]       ;combined list of indices
    vec = [ bytarr(na), replicate(1b,nb) ]  ;flag of which vector in  combined
    ;list   0 - a   1 - b

    ; sort combined list

    sub = sort(c)
    c = c[sub]
    ind = ind[sub]
    vec = vec[sub]

    ; find duplicates in sorted combined list

    n = na + nb                            ;total elements in c
    if epsilon eq 0. then $
      firstdup = where( (c EQ shift(c,-1)) and (vec NE shift(vec,-1)), Count ) $
    else $
      firstdup = where( (abs(c - shift(c,-1)) lt epsilon) and (vec NE shift(vec,-1)), Count )

    if Count EQ 0 then begin               ;any found?
      suba = lonarr(1)-1
      subb = lonarr(1)-1
      return
    end

    dup = lonarr( Count*2 )                     ;both duplicate values
    even = lindgen( N_elements(firstdup))*2     ;Changed to LINDGEN 6-Sep-1991
    dup[even] = firstdup
    dup[even+1] = firstdup+1
    ind = ind[dup]                         ;indices of duplicates
    vec = vec[dup]                         ;vector id of duplicates
    subb = ind[ where( vec, complement = vzero) ]             ;b subscripts
    suba = ind[ vzero]

  endif else begin             ;Integer calculation using histogram.

    minab = min(a, MAX=maxa) > min(b, MAX=maxb) ;Only need intersection of ranges
    maxab = maxa < maxb

    ;If either set is empty, or their ranges don't intersect:
    ;  result = NULL (which is denoted by integer = -1)
    !ERR = -1
    suba = -1
    subb = -1
    COUNT = 0L
    if (maxab lt minab) || (maxab lt 0) then return

    ha = histogram([a], MIN=minab, MAX=maxab, reverse_indices=reva)
    hb = histogram([b], MIN=minab, MAX=maxab, reverse_indices=revb)

    r = where((ha ne 0) and (hb ne 0), count)
    if count gt 0 then begin
      suba = reva[reva[r]]
      subb = revb[revb[r]]
    endif
  endelse

  return

end






;+
; NAME:
;       MATCH2
; PURPOSE:
;       Routine to cross-match values in two vectors (including non-matches)
; EXPLANATION:
;       MATCH2 reports matching elements of two arrays.

;       This procedure *appears* similar to MATCH of the IDL astronomy
;       library.  However, this routine is quite different in that it
;       reports an index value for each element of the input arrays.
;       In other words, while MATCH reports the *existence* of
;       matching elements in each array, MATCH2 reports explicitly
;       *which* elements match.
;
;       Furthermore, while MATCH reports only unique matching
;       elements, MATCH2 will always report a cross-match for every
;       element in each array, even if it is a repeat.
;
;       In cases where no match was found, an index of -1 is
;       reported.
;
; CALLING SEQUENCE:
;       match2, a, b, suba, subb
;
; INPUTS:
;       a,b - two vectors to match elements, numeric or string data
;             types.  (See below for RESTRICTIONS on A and B)
;
;
; OUTPUTS:
;       suba - vector with same number of elements as A, such that
;              A EQ B[SUBA], except non-matches which are indicated
;              by SUBA EQ -1
;       subb - vector with same number of elements as B, such that
;              B EQ A[SUBB], except non-matches which are indicated
;              by SUBB EQ -1
;
;
; RESTRICTIONS:
;
;       The vectors A and B are allowed to have duplicates in them,
;       but for matching purposes, only the first one found will
;       be reported.
;
;       If A and B are string arrays, then non-printable ASCII values
;       1B and 2B will confuse the algorithm.  Don't use these
;       non-printable characters in strings.
;
; EXAMPLE:
;      A = [0,7,14,23,24,30]
;      B = [7,8,14,25,14]
;      IDL> match2, a, b, suba, subb
;     --> suba = [ -1 ,  0,  4,  -1, -1, -1 ]
;     (indicates that A[1] matches B[1] and A[3] matches B[2])
;     --> subb = [  1 , -1,  2,  -1,  2 ]
;     (indicates that B[1] matches A[1] and B[2] matches A[3])
;
;  Compare to the results of the original MATCH procedure,
;
;      IDL> match, a, b, suba, subb
;     --> suba = [  1,  3]
;  (indicates that A[1] and A[3] match elements in B, but not which ones)
;     --> subb = [  1,  2]
;  (indicates that B[1] and B[2] match elements in A, but not which ones)
;
; MODIFICATION HISTORY
;   Derived from the IDL Astronomy Library MATCH, 14 Feb 2007
;   Updated documentation, 17 Jul 2007
;   More updated documentation (example), 03 Sep 2007
;   Bug fix for string arrays with numerical contents; the subset
;   string is now 1B and 2B; this is now documented, 2014-10-20 CM
;
;
;-
;-------------------------------------------------------------------------
pro match2, a, b, suba, subb

  On_error,2
  compile_opt idl2

  if N_params() LT 3 then begin
    print,'Syntax - match2, a, b, suba, subb'
    print,'    a,b -- input vectors for which to match elements'
    print,'    suba,subb -- match index lists'
    return
  endif

  da = size(a,/type) & db =size(b,/type)

  na = N_elements(a)              ;number of elements in a
  nb = N_elements(b)             ;number of elements in b
  suba = lonarr(na)-1 & subb = lonarr(nb)-1

  ; Check for a single element array

  if (na EQ 1) or (nb EQ 1) then begin
    if (nb GT 1) then begin
      wh = where(b EQ a[0], nw)
      if nw GT 0 then begin
        subb[wh] = 0L
        suba[0]  = wh[0]
      endif
    endif else begin
      wh = where(a EQ b[0], nw)
      if nw GT 0 then begin
        suba[wh] = 0L
        subb[0]  = wh[0]
      endif
    endelse
    return
  endif

  c = [ a, b ]                   ;combined list of a and b
  ind = [ lindgen(na), lindgen(nb) ]       ;combined list of indices
  vec = [ intarr(na), replicate(1,nb) ]  ;flag of which vector in  combined
  ;list   0 - a   1 - b

  ; sort combined list

  if da EQ 7 OR db EQ 7 then begin
    vecstr = [string(1b), string(2b)]
    ;; String sort (w/ double key)
    sub = sort(c+vecstr[vec])
  endif else begin
    ;; Number sort (w/ double key)
    eps = (machar(/double)).eps
    sub = sort(double(c)*(1d + vec*eps))
  endelse

  c = c[sub]
  ind = ind[sub]
  vec = vec[sub]

  n = na + nb                    ;total elements in c
  wh = where( c[1:*] NE c, ct)
  if ct EQ 0 then begin
    whfirst = [0]
    whlast  = [n-1]
  endif else begin
    whfirst = [0, wh+1]
    whlast  = [wh, n-1]
  endelse

  vec0 = vec[whfirst]
  vec1 = vec[whlast]
  ;; 0 = present in A but not B
  ;; 1 = can't occur (since the array was sorted on 'VEC')
  ;; 2 = present in both
  ;; 3 = present in B but not A
  matchtype = vec0 + vec1*2

  nm = n_elements(matchtype)
  mm = ind*0L & wa = mm & wb = mm
  for i = 0, nm-1 do begin
    mm[whfirst[i]:whlast[i]] = matchtype[i]
    wa[whfirst[i]:whlast[i]] = ind[whfirst[i]]
    wb[whfirst[i]:whlast[i]] = ind[whlast[i]]
  endfor

  suba = lonarr(na)-1 & subb = lonarr(nb)-1

  wh = where(mm EQ 2 AND vec EQ 0, ct)
  if ct GT 0 then suba[ind[wh]] = wb[wh]
  wh = where(mm EQ 2 AND vec EQ 1, ct)
  if ct GT 0 then subb[ind[wh]] = wa[wh]

  return
end