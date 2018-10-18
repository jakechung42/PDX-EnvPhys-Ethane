FUNCTION check_season, avail_month
;this function performs the season checking algorithm to make sure that 
;a year has all the season, if any season is missing, the function returns 0

;per the new edition that calculate annual means using the monthly averages, this
;function is not used in this program anymore. 
hasWinter = 0
hasSummer = 0
hasSpring = 0
hasFall = 0
hasAllSeasons = 0
;check Winter
lg_arr = avail_month ge 1 and avail_month le 3
for p = 0, n_elements(lg_arr)-1 do begin
	if lg_arr[p] eq 1 then hasWinter = 1
endfor
;check Spring
lg_arr = avail_month ge 4 and avail_month le 6
for p = 0, n_elements(lg_arr)-1 do begin
	if lg_arr[p] eq 1 then hasSpring = 1
endfor
;check Summer
lg_arr = avail_month ge 7 and avail_month le 9
for p = 0, n_elements(lg_arr)-1 do begin
	if lg_arr[p] eq 1 then hasSummer = 1
endfor	
;check Fall
lg_arr = avail_month ge 10 and avail_month le 12
for p = 0, n_elements(lg_arr)-1 do begin
	if lg_arr[p] eq 1 then hasFall = 1
endfor	
	
if (hasWinter eq 0 or hasSpring eq 0 or hasSummer eq 0 or hasFall eq 0) then hasAllSeasons = 0 $
	else hasAllSeasons = 1
	
return, hasAllSeasons

end



FUNCTION NAN2neg999, array
;this function changes NAN values to -999 so the output file can be read by another program.
notNAN = finite(array)
x = where(notNAN eq 0, count)
if (count gt 0) then begin 
	array[x] = -999
endif

return, array

end



FUNCTION annual_mean, ratio_vec, month_vec
;this function calculates the annual average of a site by calculating the 
;average of each month and then calculate the average of the year using
;the results from the all months. If a month is missing, the program will
;produce a NAN value 
;the input is the ratio vector and the corresponding month vector 
;the output is the annual mean and the annual standard error
;NotEnoughData is a variable to verify if the year has all 12 months
NotEnoughData = 0
annual_avg = !VALUES.F_NAN
annual_err = !VALUES.F_NAN

;get months that are available
hasMonth = month_vec[sort(month_vec)]
hasMonth = hasMonth[uniq(hasMonth)]

;calculate the annual mean and annual standard error
if n_elements(hasMonth) ne 12 then begin
	NotEnoughData = 1
endif else begin
	month_avg = fltarr(12)
	month_err = fltarr(12)
	for i = 1, 12 do begin
		x = where(month_vec eq i, count0)
		month_avg[i - 1] = mean(ratio_vec[x])
		;calculate standard error of each month. If a month only has one value
		;the stddev output will be NAN. 
		month_err[i - 1] = stddev(ratio_vec[x], /DOUBLE)/sqrt(count0)
	endfor 
	
	;find which month that has only one value
	notNAN = finite(month_err)
	NAN_idx = where(notNAN eq 0, count1)
	if count1 gt 0 then begin
		;the standard error for the month that only has one value will be the 
		;average of the standard error of the entire year.
		month_err[NAN_idx] = mean(month_err, /NAN)
	endif
	annual_avg = mean(month_avg)
	annual_err = 1/12.0 * sqrt(total(month_err^2))
endelse 
out = [annual_avg, annual_err]
return, out

end



FUNCTION annual_mean_uci, ratio_vec, month_vec
;this function is ONLY used for the UCI data set because the UCI data set 
;only has data in March, June, September and December, so it slightly different 
;from the function annual_mean
;this function calculates the annual average of a site by calculating the 
;average of each month and then calculate the average of the year using
;the results from the all months. If a month is missing, the program will
;produce a NAN value 
;the input is the ratio vector and the corresponding month vector 
;the output is the annual mean and the annual standard error
;EnoughData is a variable to verify if the year has all 12 months
EnoughData = 1
annual_avg = !VALUES.F_NAN
annual_err = !VALUES.F_NAN

;get months that are available
hasMonth = month_vec[sort(month_vec)]
hasMonth = hasMonth[uniq(hasMonth)]
month = [3, 6, 9, 12]

;check to see if all the months have data
if n_elements(hasMonth) eq 4 then begin 
	for j = 0, 3 do begin
		if (hasMonth[j] ne month[j]) then begin
			EnoughData = 0
			break
		endif
	endfor	
endif else EnoughData = 0

;calculate the annual mean and annual standard error
if EnoughData eq 1 then begin
	month_avg = fltarr(4)
	month_err = fltarr(4)
	for i = 0, 3 do begin
		x = where(month_vec eq month[i], count0)
		month_avg[i] = mean(ratio_vec[x])
		;calculate standard error of each month. If a month only has one value
		;the stddev output will be NAN. 
		month_err[i] = stddev(ratio_vec[x], /DOUBLE)/sqrt(count0)
	endfor 
	
	;find which month that has only one value
	notNAN = finite(month_err)
	NAN_idx = where(notNAN eq 0, count1)
	if count1 gt 0 then begin
		;the standard error for the month that only has one value will be the 
		;average of the standard error of the entire year.
		month_err[NAN_idx] = mean(month_err, /NAN)
	endif
	annual_avg = mean(month_avg)
	annual_err = 1/4.0 * sqrt(total(month_err^2))
endif 
out = [annual_avg, annual_err]
return, out

end





PRO time_series_barrow

;this program plots the annual mean of the Barrow data from all 3 networks. 
;there's a program already for this task but that programs doesn't use the
;new functions made during Summber 2018 and it doesn't have the 1984 - 1996 data.
;So to be consistant with the Cape Grim time series program, this program
;will be built similar to it.

compile_opt idl2

;read in the 1984 1996 UCI data. uci_new is the data structure of the 1984 - 1996 UCI data
uci_new = read_uci8496()

;read in NOAA data
noaa = read_noaa()

;read in UCI data
uci = read_uci()

;read in OGI data
ogi = read_ogi()

;prepare the data
;get Barrow(br) data from the OGI and NOAA data
;br data for NOAA is BRW
;for OGI, the coordinate of br is longitude -156.5, latitude 71.16

br_noaa_idx = where(strmatch(noaa.site[*], 'BRW', /FOLD_CASE) EQ 1)
;reconstrut the strucure for Barrow NOAA
;the time tag in the br_noaa structure is decimal year
br_noaa = { month: noaa.month[br_noaa_idx], $
            year: noaa.year[br_noaa_idx], $
            ratio: noaa.ratio[br_noaa_idx], $
            time: noaa.year[br_noaa_idx] + noaa.month[br_noaa_idx]/13, $
            deseason: fltarr(n_elements(br_noaa_idx))}

print, 'Barrow NOAA: '
help, br_noaa
;OGI:
br_ogi_idx = where(ogi.lat eq 71.16 and ogi.lon eq -156.5)
;reconstrut the strucure for Brrow OGI
;the time tag in the br_ogi structure is decimal year
br_ogi = { month: ogi.month[br_ogi_idx], $
            year: ogi.year[br_ogi_idx], $
            ratio: ogi.ratio[br_ogi_idx], $
            time: ogi.year[br_ogi_idx] + ogi.month[br_ogi_idx]/13, $
            deseason: fltarr(n_elements(br_ogi_idx))}

print, 'Barrow OGI: '
help, br_ogi
;get UCI data
br_uci_idx = where(uci.lat eq 72 and uci.lon eq -157.5)
br_uci_new_idx = where(uci_new.lat eq 72 and uci_new.lon eq -157.5)
;reconstrut the strucure for Brrow UCI
;the time tag in the br_uci structure is decimal year
temp_month = [uci.month[br_uci_idx], uci_new.month[br_uci_new_idx]]
temp_year = [uci.year[br_uci_idx], uci_new.year[br_uci_new_idx]]
temp_ratio = [uci.ratio[br_uci_idx], uci_new.ratio[br_uci_new_idx]]
temp_time = [uci.year[br_uci_idx] + uci.month[br_uci_idx]/13, uci_new.year[br_uci_new_idx] + uci_new.month[br_uci_new_idx]/13]

br_uci = { month: temp_month, $
            year: temp_year, $
            ratio: temp_ratio, $
            time: temp_time, $
            deseason: fltarr(n_elements(br_uci_idx) + n_elements(br_uci_new_idx))}

print, 'Barrow UCI: '
help, br_uci
;deseasonalize the data
br_noaa.deseason = deseason(br_noaa.month, br_noaa.year, br_noaa.ratio)
br_ogi.deseason = deseason(br_ogi.month, br_ogi.year, br_ogi.ratio)
br_uci.deseason = deseason(br_uci.month, br_uci.year, br_uci.ratio)


;filter the data for NOAA
;data that are out of the 3 sigma bound are set as NAN value
stddev_noaa = std_dev_fgauss(br_noaa.deseason)
outlier_idx = where(br_noaa.deseason gt stddev_noaa*3 or br_noaa.deseason lt -stddev_noaa*3, count)
if count gt 0 then begin
  br_noaa.ratio[outlier_idx] = !VALUES.F_NAN
  br_noaa.month[outlier_idx] = !VALUES.F_NAN
  br_noaa.year[outlier_idx] = !VALUES.F_NAN
  br_noaa.time[outlier_idx] = !VALUES.F_NAN
  br_noaa.deseason[outlier_idx] = !VALUES.F_NAN
endif

print, 'Removed ', count, ' data points out of the Barrow NOAA data set'
;filter the data for OGI
stddev_ogi = std_dev_fgauss(br_ogi.deseason)
outlier_idx = where(br_ogi.deseason gt stddev_ogi*3 or br_ogi.deseason lt -stddev_ogi*3, count)
if count gt 0 then begin
  br_ogi.ratio[outlier_idx] = !VALUES.F_NAN
  br_ogi.month[outlier_idx] = !VALUES.F_NAN
  br_ogi.year[outlier_idx] = !VALUES.F_NAN
  br_ogi.time[outlier_idx] = !VALUES.F_NAN
  br_ogi.deseason[outlier_idx] = !VALUES.F_NAN
endif

print, 'Removed ', count, ' data points out of the Barrow OGI data set'
;filter the data for UCI
stddev_uci = std_dev_fgauss(br_uci.deseason)
outlier_idx = where(br_uci.deseason gt stddev_uci*3 or br_uci.deseason lt -stddev_uci*3, count)
if count gt 0 then begin
  br_uci.ratio[outlier_idx] = !VALUES.F_NAN
  br_uci.month[outlier_idx] = !VALUES.F_NAN
  br_uci.year[outlier_idx] = !VALUES.F_NAN
  br_uci.time[outlier_idx] = !VALUES.F_NAN
  br_uci.deseason[outlier_idx] = !VALUES.F_NAN
endif

print, 'Removed ', count, ' data points out of the Barrow UCI data set'



;make a quick plot of the overlapping period of NOAA and UCI data
;set up plot
; cgDisplay, 1344, 756
; cgPlot, br_uci.time, br_uci.ratio, /nodata, xtitle = 'Time', ytitle = 'Mixing ratio (pptv)', $
	; title = 'Overlapping period between UCI and NOAA at Barrow', xrange = [2004, 2010], $
	; yrange = [0,3500]
; cgPlot, br_uci.time, br_uci.ratio, /overplot, psym = 2, color = 'forest green', symsize = 1.5
; cgPlot, br_noaa.time, br_noaa.ratio, /overplot, psym = 4, color = 'steelblue', symsize = 1.5

;make another quick plot of the overlapping period of UCI and OGI
; cgDisplay, 1344, 756
; cgPlot, br_ogi.time, br_ogi.ratio, /nodata, xtitle = 'Time', ytitle = 'Mixing ratio (pptv)', $
	; title = 'Overlapping period between OGI and UCI at Barrow', xrange = [1983, 1987], $
	; yrange = [0,3500]
; cgPlot, br_uci.time, br_uci.ratio, /overplot, psym = 2, color = 'forest green', symsize = 1.5
; cgPlot, br_ogi.time, br_ogi.ratio, /overplot, psym = 4, color = 'steelblue', symsize = 1.5



;-----br_noaa calculations-----
;find how many years and which years the NOAA network has data for
temp_year_noaa = br_noaa.year[sort(br_noaa.year)]
;since there are NAN values in the array, the sort and uniq function
;does not deal with NAN value very satisfactory, need some preliminary 
;algorithm to clean up the NAN values. If all the NAN indexes are removed
;at this step, all subsequent steps should not have problems with NAN 
notNAN = finite(temp_year_noaa)
x = where(notNAN eq 0, count)
if (count gt 0) then begin
	;remove all the rows with NAN value
	temp_year_noaa = RemoveRows(rotate(temp_year_noaa, 1), x)
endif

year_noaa = temp_year_noaa[uniq(temp_year_noaa)]
avg_noaa = fltarr(n_elements(year_noaa));contains the annual mean for NOAA
stderr_noaa = fltarr(n_elements(year_noaa));contains the annual error for NOAA

;write the annual data out to a file
infile = '/home/excluded-from-backup/ethane/IDL/temp_file/annual_mean_barrow.dat'
openw, lun, infile, /get_lun


printf, lun, 'Year|Annual Average Mixing Ratio|Std Error|Number of Samples'
printf, lun, '-----------'
printf, lun, '>>>>>>>>>>NOAA report'
for i = 0, n_elements(year_noaa) - 1 do begin
	x = where(br_noaa.year[*] eq year_noaa[i], count)
	
	;print out the available months of a year
	flt_month = br_noaa.month[x]
	flt_month = flt_month[sort(flt_month)]
	
	;run the function to calculate the annual mean and annual error
	temp = annual_mean(br_noaa.ratio[x], br_noaa.month[x])
	
	;pull the mean and the standard error
	avg_noaa[i] = temp[0]
	stderr_noaa[i] = temp[1]
	printf, lun, year_noaa[i], avg_noaa[i], stderr_noaa[i], count
	printf, lun, 'Months with data: ', flt_month[uniq(flt_month)]
	printf, lun, '----------------------------'
	;change the years with no data to NAN 
	if finite(avg_noaa[i]) eq 0 then begin
		year_noaa[i] = !VALUES.F_NAN
	endif

endfor
printf, lun, '>>>>>END NOAA REPORT<<<<<'


;-----br_ogi-----
;find how many years and which years the OGI network has data for
temp_year_ogi = br_ogi.year[sort(br_ogi.year)]

;since there are NAN values in the array, the sort and uniq function
;does not deal with NAN value very satisfactory, need some preliminary 
;algorithm to clean up the NAN values. If all the NAN indexes are removed
;at this step, all subsequent steps should not have problems with NAN 
notNAN = finite(temp_year_ogi)
x = where(notNAN eq 0, count)
if (count gt 0) then begin
	;remove all the rows with NAN value
	temp_year_ogi = RemoveRows(rotate(temp_year_ogi, 1), x)
endif

year_ogi = temp_year_ogi[uniq(temp_year_ogi)]
avg_ogi = fltarr(n_elements(year_ogi));contains the annual mean for OGI
stderr_ogi = fltarr(n_elements(year_ogi));contains the annual error for OGI


printf, lun, '-----------'
printf, lun, '>>>>>>>>>OGI report'
for i = 0, n_elements(year_ogi) - 1 do begin
	x = where(br_ogi.year[*] eq year_ogi[i], count)
	
	;print out the available months of a year
	flt_month = br_ogi.month[x]
	flt_month = flt_month[sort(flt_month)]
	
	;run the function to get average and the std error
	temp = annual_mean(br_ogi.ratio[x], br_ogi.month[x])
		
	;pull out the average and the std error from the output.
	avg_ogi[i] = temp[0]
	stderr_ogi[i] = temp[1]
	printf, lun, year_ogi[i], avg_ogi[i], stderr_ogi[i], count
	printf, lun, 'Months with data: ', flt_month[uniq(flt_month)]
	printf, lun, '----------------------------'

	;change the years with no data to NAN 
	if finite(avg_ogi[i]) eq 0 then begin
		year_ogi[i] = !VALUES.F_NAN
	endif
	
endfor
printf, lun, '>>>>>END OGI REPORT<<<<<'



;-----uci-----
;find how many years and which years the uci1 network has data for
temp_year_uci = br_uci.year[sort(br_uci.year)]

;since there are NAN values in the array, the sort and uniq function
;does not deal with NAN value very satisfactory, need some preliminary 
;algorithm to clean up the NAN values. If all the NAN indexes are removed
;at this step, all subsequent steps should not have problems with NAN 
notNAN = finite(temp_year_uci)
x = where(notNAN eq 0, count)
if (count gt 0) then begin
	;remove all the rows with NAN value
	temp_year_uci = RemoveRows(rotate(temp_year_uci, 1), x)
endif

year_uci = temp_year_uci[uniq(temp_year_uci)]
avg_uci = fltarr(n_elements(year_uci));contains the annual mean for OGI
stderr_uci = fltarr(n_elements(year_uci));contains the annual error for OGI


printf, lun, '-----------'
printf, lun, '>>>>>>>>>UCI report'
for i = 0, n_elements(year_uci) - 1 do begin
	x = where(br_uci.year[*] eq year_uci[i], count)
	
	;print out the available months of a year
	flt_month = br_uci.month[x]
	flt_month = flt_month[sort(flt_month)]
	
	;run the annual mean function to get the average and std error
	temp = annual_mean_uci(br_uci.ratio[x], br_uci.month[x])

	
	avg_uci[i] = temp[0]
	stderr_uci[i] = temp[1]
	printf, lun, year_uci[i], avg_uci[i], stderr_uci[i], count
	printf, lun, 'Months with data: ', flt_month[uniq(flt_month)]
	printf, lun, '----------------------------'

	;change the years with no data to NAN 
	if finite(avg_uci[i]) eq 0 then begin
		year_uci[i] = !VALUES.F_NAN
	endif
endfor
printf, lun, '>>>>>END UCI REPORT<<<<<'


free_lun, lun


;plotting procedure
open_device, /ps, /color, file='temp.eps', margin=0.05, xsize = 10.0, ysize = 7.5
!x.thick=1
!y.thick=1
!p.font=1
!p.thick =0.5


;set up plot
cgPlot, year_noaa, avg_noaa, /nodata, xrange = [1980, 2016], xticklen = 1, $
	xgridstyle = 1, yrange = [1000,2200], xticks = 18, $
	xtitle = 'Time', ytitle = 'Mixing Ratio (pptv)', $
	title = 'Time series of Barrow'
cgPlot, year_noaa, avg_noaa, /overplot, psym = 4, color = 'forest green'
cgPlot, year_noaa, avg_noaa, /overplot, color = 'forest green', $
	err_ylow = stderr_noaa, err_yhigh = stderr_noaa
cgPlot, year_ogi, avg_ogi, /overplot, psym = 5, color = 'steelblue'
cgPlot, year_ogi, avg_ogi, /overplot, color = 'steelblue', $
	err_ylow = stderr_ogi, err_yhigh = stderr_ogi
cgPlot, year_uci, avg_uci, /overplot, psym = 2, color = 'tomato'
cgPlot, year_uci, avg_uci, /overplot, color = 'tomato', $
	err_ylow = stderr_uci, err_yhigh = stderr_uci


cgLegend, SymColors = ['forest green', 'tomato', 'steelblue'], $
	PSyms = [4, 2, 5], Symsize = 1.5, Location = [0.7, 0.85], $
	titles = ['NOAA', 'UCI', 'OGI'], $
	/Box, /Background, BG_Color = 'rose', /center_sym, length = 0
close_device

spawn, 'gv temp.eps'




;write to an output file to make IHR in a different program

;change NAN values to -999 in order to read by other programs
year_noaa = NAN2neg999(year_noaa)
avg_noaa = NAN2neg999(avg_noaa)
stderr_noaa = NAN2neg999(stderr_noaa)

year_ogi = NAN2neg999(year_ogi)
avg_ogi = NAN2neg999(avg_ogi)
stderr_ogi = NAN2neg999(stderr_ogi)

year_uci = NAN2neg999(year_uci)
avg_uci = NAN2neg999(avg_uci)
stderr_uci = NAN2neg999(stderr_uci)

;write to file
infile = '/home/excluded-from-backup/ethane/IDL/temp_file/annual_mean_barrow_NoHeaders.dat'
openw, lun, infile, /get_lun

for i = 0, n_elements(year_noaa)-1 do begin
	printf, lun, '1', year_noaa[i], avg_noaa[i], stderr_noaa[i]
endfor

for i = 0, n_elements(year_uci)-1 do begin
	printf, lun, '2', year_uci[i], avg_uci[i], stderr_uci[i]
endfor

for i = 0, n_elements(year_ogi)-1 do begin
	printf, lun, '3', year_ogi[i], avg_ogi[i], stderr_ogi[i]
endfor

free_lun, lun

end