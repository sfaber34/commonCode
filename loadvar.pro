function loadvar,varname,ncid=ncid,filename=filename,cal=cal, double_prec=double_prec, $
  scale = scale, offset = offset, coef = coef, adcoef=adcoef, span=span, attr = attr, help=help

  ;+
  ;
  ;
  ;  $Id: loadvar.pro,v 1.1.1.1 2005/05/03 17:48:13 french Exp $
  ;  MODIFIED 2015/03/26:  added ability to read and return variable attributes 'Span' and 'AdditionalCoefficients'
  ;
  ;
  ; A generic function to load a variable from a netCDF file,
  ;
  ; The function will also look for standard variable attributes (scale_factor &
  ;   add_offset) and non-standard attributes (cal_coef) and apply them if
  ; they exist.
  ;
  ; Added in 2015, if user provides a string for the keyword 'attr', then function
  ;   will search if that variable attrib exists, and if it does return the value of
  ;   that attribute rather than the value of the variable. If the attrib does not exist
  ;   it returns -1
  ;
  ; The function is written to be backwards compatible with old version of
  ; of ARL aircraft data which contain calibration coefficients written
  ; with 'cal_pwr' and 'coefn' where n is 0 or 1
  ;
  ; The function returns the values for varname
  ;
  ; One can either pass the NAME of a CLOSED netCDF file through the 'filename'
  ; keyword, or the netCDF FILE ID of an OPEN netCDF file through the 'ncid'
  ; keyword. If using 'filename', the file will be closed by loadvar.
  ; If using 'ncid', the file will remain open.
  ;
  ;
  if KEYWORD_SET(help) then begin
    print,''
    print,'Result = loadvar(varname, [ncid=ncid | filename=filename])'
    print,''
    print,'ARGUMENTS'
    print,'varname  : name of netCDF variable'
    print,'ncid   : netCDF file ID of OPEN netCDF file'
    print,'filename : name (& path) of CLOSED netCDF file'
    print,' NOTE: call requires either ncid or filename KEYWORD'
    print,''
    print,'KEYWORDS'
    print,'attr    : set this keyword to the name of a variable attribute'
    print,'   if set and the attribute exists function returns value of attribute'
    print,'   if set and the attribute doesnot exist function returns -1'
    print,' '
    print,'cal    : set this keyword to apply calibration coefficients'
    print,'   if they exist (ignored if they do not exist)'
    print,'double_prec: if set, result is returned as double precision'
    print,'   otherwise returned as float'
    print,'scale    : name of variable that will contain, upon exit, value'
    print,'   of scale factor'
    print,'offset   : name of variable that will contain, upon exit, value'
    print,'   of add_offset'
    print,'coef   : name of variable that will contain, upon exit, value(s)'
    print,'   of calibration coefficient(s)'
    RETURN, -1
  endif
  ;-

  IF NOT isa(ncid) THEN $
    close_file = 1

  IF KEYWORD_SET(filename) THEN $
    ncid = NCDF_OPEN(filename)


  IF NOT isa(ncid) THEN BEGIN
    PRINT, "loadvar: no filename or valid file ID given"
    RETURN,-1
  ENDIF

  IF (KEYWORD_SET(double_prec)) THEN dbl = 1 ELSE dbl = 0

  varid = NCDF_VARID(ncid, varname)


   vardata = NCDF_VARINQ(ncid, varid)


  scale = 1. & offset = 0.
  FOR i = 1, vardata.natts DO BEGIN
    attname = NCDF_ATTNAME(ncid, varid, i-1)
    ; print,attname
    IF (attname eq 'add_offset') THEN BEGIN
      NCDF_ATTGET, ncid, varid, 'add_offset', offset
      IF (dbl) THEN offset = double(offset)
    ENDIF
    IF (attname eq 'scale_factor') THEN BEGIN
      NCDF_ATTGET, ncid, varid, 'scale_factor', scale
      IF (dbl) THEN scale = double(scale)
    ENDIF
    IF ( (KEYWORD_SET(cal)) AND (attname EQ 'cal_coef') ) THEN BEGIN
      NCDF_ATTGET, ncid, varid, 'cal_coef', cal_coef
      IF (dbl) THEN cal_coef = double(cal_coef)
    ENDIF

    ;the following is for special attributes contained in some netcdf variables in UWKA files

    IF (  (attname EQ 'AdditionalCoefficients') ) THEN $
      NCDF_ATTGET, ncid, varid, 'AdditionalCoefficients', adcoef
    IF (  (attname EQ 'Span') ) THEN $
      NCDF_ATTGET, ncid, varid, 'Span', span

    ; the following is to be back compatible with old ncr files
    IF ( (KEYWORD_SET(cal)) AND (attname EQ 'cal_pwr') ) THEN $
      NCDF_ATTGET, ncid, varid, 'cal_pwr', cal_pwr
    IF ( (KEYWORD_SET(cal)) AND (attname EQ 'coef_0') ) THEN BEGIN
      NCDF_ATTGET, ncid, varid, 'coef_0', coef_0
      IF (dbl) THEN coef_0 = double(coef_0)
    ENDIF
    IF ( (KEYWORD_SET(cal)) AND (attname EQ 'coef_1') ) THEN BEGIN
      NCDF_ATTGET, ncid, varid, 'coef_1', coef_1
      IF (dbl) THEN coef_1 = double(coef_1)
    ENDIF

  ENDFOR


  ;check to see if keyword 'attr' is present -- if it is, get the attribute and return it
  ;if the keyword is not present, then continue getting the variable data and return that

  IF isa(attr) THEN BEGIN
    NCDF_ATTGET, ncid, varid, attr, temp
    RETURN, temp
  ENDIF ELSE BEGIN
    NCDF_VARGET, ncid, varid, temp_u
  ENDELSE

  temp = temp_u*scale + offset

  ; the following is to be back compatible with old ncr files
  IF isa(cal_pwr) THEN BEGIN
    cal_coef = FLTARR(2)
    cal_coef[0] = coef_0
    cal_coef[1] = coef_1
  ENDIF

  engin = 0
  IF ( KEYWORD_SET(cal) AND (N_ELEMENTS(cal_coef) gt 0) ) THEN BEGIN
    coef = cal_coef
    FOR i = 0, N_ELEMENTS(cal_coef)-1 DO $
      engin = cal_coef[i]*(temp^i)+engin
    temp = engin
  ENDIF

  IF isa(close_file) THEN $
    NCDF_CLOSE, ncid

  RETURN, temp

end