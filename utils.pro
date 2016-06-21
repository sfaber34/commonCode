function n, var
  return, float(n_elements(var)-1)
end

function n1, var
  return, float(n_elements(var))
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



function hist, var, binsize
  
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
  
  b1=barplot(dindgen(n1(h1)),h1/n1(var),dimensions=[1800,1000],histogram=1,margin=[50,100,50,50],/device)
  b1.xrange=[0,n1(h1)]
  b1.xmajor=n1(h1)+1
  b1.xtickname=[string(inds),string(max(var)),' ',' ']
  b1.xtext_orientation=90
  b1.xticklen=0
  return,b1
end  




function lp, xVar,yVar
  if xVar eq 0 then begin
    p1=plot(yvar,dimensions=[1600,1200],thick=2,margin=50,/device)
  endif else begin
    p1=plot(xvar,yvar,dimensions=[1600,1200],thick=2,margin=50,/device)
  endelse
  
  return,p1
end