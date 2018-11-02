FUNCTION deseason, month, year, ratio

;this function will perfrom the deseasonal process for an entire data set
;input: month, year, ratio
;output: momth, year, deseasonal data

residual = fltarr(n_elements(ratio))
holder = fltarr(n_elements(ratio))
b = fltarr(13)
for j = 1, 12 do begin
	a = where(month eq j, count)
	if count gt 0 then begin
		b[j] = mean(ratio[a], /NAN)
		;the holder array contains the monthly averages that
		;is subtracted out to get the residual
		holder[a] = b[j]
	endif
endfor
residual = ratio - holder
return, residual

end
