FUNCTION neg999toNAN, array

;convert -999 values to NAN values
neg999 = where(array eq -999, count)
if (count gt 0) then begin
	array[neg999] = !VALUES.F_NAN
endif

return, array
end



PRO time_series_ihr_BR_CP

compile_opt idl2

;this program plots the annual mean and annual IHR of the Barrow and Cape Grim sites

;read Barrow's annual means (output of time_series_barrow)
infile = '/home/excluded-from-backup/ethane/IDL/temp_file/annual_mean_barrow_NoHeaders.dat'
openr, lun, infile, /get_lun

n_line = file_lines(infile)

c1 = 0.0
c2 = 0.0
c3 = 0.0
c4 = 0.0
file_arr = fltarr(4, n_line)

for i = 0, n_line - 1 do begin
	readf, lun, c1, c2, c3, c4
	file_arr[0, i] = c1
	file_arr[1, i] = c2
	file_arr[2, i] = c3
	file_arr[3, i] = c4
endfor

free_lun, lun

;change the -999 values to NAN values
x = where(file_arr[1, *] eq -999)
file_arr[1, x] = !VALUES.F_NAN
file_arr[2, x] = !VALUES.F_NAN
file_arr[3, x] = !VALUES.F_NAN


noaa_idx = where(file_arr[0, *] eq 1)
uci_idx = where(file_arr[0, *] eq 2)
ogi_idx = where(file_arr[0, *] eq 3)


br_noaa = { year: (file_arr[1, noaa_idx]), $
			avg: (file_arr[2, noaa_idx]), $
			err: (file_arr[3, noaa_idx]) }

br_uci = { year: (file_arr[1, uci_idx]), $
			avg: (file_arr[2, uci_idx]), $
			err: (file_arr[3, uci_idx]) }

br_ogi = { year: (file_arr[1, ogi_idx]), $
			avg: (file_arr[2, ogi_idx]), $
			err: (file_arr[3, ogi_idx]) }


;read Cape Grim's annual means (output of time_series_capegrim)
infile = '/home/excluded-from-backup/ethane/IDL/temp_file/annual_mean_capegrim_NoHeaders.dat'
openr, lun, infile, /get_lun

n_line = file_lines(infile)

c1 = 0.0
c2 = 0.0
c3 = 0.0
c4 = 0.0
file_arr = fltarr(4, n_line)

for i = 0, n_line - 1 do begin
	readf, lun, c1, c2, c3, c4
	file_arr[0, i] = c1
	file_arr[1, i] = c2
	file_arr[2, i] = c3
	file_arr[3, i] = c4
endfor

free_lun, lun

;change the -999 values to NAN values
x = where(file_arr[1, *] eq -999)
file_arr[1, x] = !VALUES.F_NAN
file_arr[2, x] = !VALUES.F_NAN
file_arr[3, x] = !VALUES.F_NAN

noaa_idx = where(file_arr[0, *] eq 1)
uci_idx = where(file_arr[0, *] eq 2)
ogi_idx = where(file_arr[0, *] eq 3)

cg_noaa = { year: (file_arr[1, noaa_idx]), $
			avg: (file_arr[2, noaa_idx]), $
			err: (file_arr[3, noaa_idx]) }

cg_uci = { year: (file_arr[1, uci_idx]), $
			avg: (file_arr[2, uci_idx]), $
			err: (file_arr[3, uci_idx]) }

cg_ogi = { year: (file_arr[1, ogi_idx]), $
			avg: (file_arr[2, ogi_idx]), $
			err: (file_arr[3, ogi_idx]) }
			
print, 'Make sure the years of Barrow UCI match the years of Cape Grim band UCI, or else the IHR is nor reliable'
if (n_elements(br_uci.year) ne n_elements(cg_uci.year)) then $
	print, 'Error!!! The size of the Barrow and Cape Grim UCI does not match, IHR is not reliable'
;the OGI and UCI data happen to be the same size and with the years match up, so IHR
;can be calculated easily. NOAA data don't have the same size so algorithm to get IHR
;for NOAA will be different from OGI and UCI
ihr_uci = br_uci.avg/cg_uci.avg
ihr_ogi = br_ogi.avg/cg_ogi.avg
;calculate the standard error of the IHR
ihr_ogi_err = sqrt( (br_ogi.err/cg_ogi.avg)^2 + (br_ogi.avg*cg_ogi.err/cg_ogi.avg^2)^2 )
ihr_uci_err = sqrt( (br_uci.err/cg_uci.avg)^2 + (br_uci.avg*cg_uci.err/cg_uci.avg^2)^2 )

ihr_noaa = fltarr(n_elements(br_noaa.year))
ihr_noaa_err = fltarr(n_elements(br_noaa.year))
;calculate the interhemispheric ratio for NOAA
for i = 0, n_elements(br_noaa.year)-1 do begin
	x = where(cg_noaa.year[*] eq br_noaa.year[i], count)
	if (count gt 0) then begin
		ihr_noaa[i] = br_noaa.avg[i] / cg_noaa.avg[x]
		ihr_noaa_err[i] = sqrt( (br_noaa.err[i]/cg_noaa.avg[x])^2 + (br_noaa.avg[i]*cg_noaa.err[x]/cg_noaa.avg[x]^2)^2 )
	endif else begin
		ihr_noaa[i] = !VALUES.F_NAN
		ihr_noaa_err[i] = !VALUES.F_NAN
	endelse
endfor

;write the ihr and ihr error array out to a file to plot along with the Barrow Cape Grim 
;sampling sensitivity study
infile = '/home/excluded-from-backup/ethane/IDL/temp_file/time_series_ihr_BR_CP.dat'
openw, lun, infile, /get_lun

for i = 0, n_elements(br_noaa.year)-1 do begin
	printf, lun, '1', br_noaa.year[i], ihr_noaa[i], ihr_noaa_err[i]
endfor
 
for i = 0, n_elements(br_uci.year)-1 do begin
	printf, lun, '2', br_uci.year[i], ihr_uci[i], ihr_uci_err[i]
endfor

for i = 0, n_elements(br_ogi.year)-1 do begin
	printf, lun, '3', br_ogi.year[i], ihr_ogi[i], ihr_ogi_err[i]
endfor

free_lun, lun
;plotting procedure
;set up plot

open_device, /ps, /color, file='temp.eps', margin=0.05, xsize = 10.0, ysize = 7.5
!x.thick=1
!y.thick=1
!p.font=1
!p.thick =0.5

multiplot, /default    ; resets multiplot settings
multiplot, [1,3], ygap=0.002, xgap=0;  sets up multiplot 

;plot barrow
cgPlot, br_noaa.year, br_noaa.avg, xrange = [1982, 2016], xticklen = 1, xgridstyle = 1, $
	xticks = 17, /nodata, yrange = [1000,2200], ytitle = 'Mixing ratio(pptv)', $
	title = 'Time series of Barrow and Cape Grim along with IHR. UCI pulls data from lat 38 to 46 south'
cgPlot, br_noaa.year, br_noaa.avg, /overplot, err_yhigh = br_noaa.err, err_ylow = br_noaa.err, $
	psym = 5, color = 'steelblue'
cgPlot, br_uci.year, br_uci.avg, /overplot, err_yhigh = br_uci.err, err_ylow = br_uci.err, $
	psym = 2, color = 'forest green'
cgPlot, br_ogi.year, br_ogi.avg, /overplot, err_yhigh = br_ogi.err, err_ylow = br_ogi.err, $
	psym = 4, color = 'red'
	
multiplot, /doyaxis, /doxaxis

;plot Cape Grim
cgPlot, cg_noaa.year, cg_noaa.avg, /nodata, xrange = [1982, 2016], xticklen = 1, xgridstyle = 1, $
	xticks = 17, XTickformat='(A1)', yrange = [150,350], ytitle = 'Mixing ratio(pptv)'
cgPlot, cg_noaa.year, cg_noaa.avg, /overplot, err_yhigh = cg_noaa.err, err_ylow = cg_noaa.err, $
	psym = 5, color = 'steelblue'
cgPlot, cg_uci.year, cg_uci.avg, /overplot, err_yhigh = cg_uci.err, err_ylow = cg_uci.err, $
	psym = 2, color = 'forest green'
cgPlot, cg_ogi.year, cg_ogi.avg, /overplot, err_yhigh = cg_ogi.err, err_ylow = cg_ogi.err, $
	psym = 4, color = 'red'

multiplot, /doyaxis, /doxaxis

	
cgPlot, br_uci.year, ihr_uci, /nodata, xtitle = 'Years', ytitle = 'IHR', $
	xrange = [1982, 2016], yrange = [4, 9], $
	xticklen = 1, xgridstyle = 1, xticks = 17
cgPlot, br_uci.year, ihr_uci, /overplot, psym = 2, color = 'forest green', err_yhigh = ihr_uci_err, $
	err_ylow = ihr_uci_err
cgPlot, br_ogi.year, ihr_ogi, /overplot, psym = 4, color = 'red', err_yhigh = ihr_ogi_err, $
	err_ylow = ihr_ogi_err
cgPlot, br_noaa.year, ihr_noaa, /overplot, psym = 5, color = 'steelblue', err_ylow = ihr_noaa_err, $
	err_yhigh = ihr_noaa_err

AL_Legend, ['OGI', 'UCI', 'NOAA'], psym = [4, 2, 5], linestyle = [0, 0, 0], box = 1, $
	position = [1988, 8.7], color = ['red', 'forest green', 'steelblue'], $
	background_color = 'rose'
close_device

spawn, 'gv temp.eps'

multiplot, /default
end
