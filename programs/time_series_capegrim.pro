FUNCTION check_season, avail_month
;this function performs the season checking algorithm to make sure that 
;a year has all the season, if any season is missing, the function returns 0
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





PRO time_series_capegrim

;this program plots the annual mean of the Cape Grim data for NOAA and OGI data sets
;For UCI, it uses the latitude bands -46 to -38 and treat all the sites in that band
;as just one single site.

compile_opt idl2


;read in the 1984 1996 UCI data. uci_new is the data structure of the 1984 - 1996 UCI data
uci_new = read_uci8496()
print, 'UCI 1984 - 1996 data: '
help, uci_new, /str

;read in NOAA data
noaa = read_noaa()
print, 'NOAA data: '
help, noaa, /str

;read in UCI data
uci = read_uci()
print, 'UCI data: '
help, uci, /str

;read in OGI data
ogi = read_ogi()
print, 'OGI data: '
help, ogi, /str

;prepare the data
;get Cape Grim (cg) data from the OGI and NOAA data
;cg data for NOAA is CGO
;for OGI, the coordinate of cg is longitude 145, latitude -42
;NOAA:
cg_noaa_idx = where(strmatch(noaa.site[*], 'CGO', /FOLD_CASE) EQ 1)
;reconstrut the strucure for Cape Grim NOAA
;the time tag in the cg_noaa structure is decimal year
cg_noaa = { month: noaa.month[cg_noaa_idx], $
            year: noaa.year[cg_noaa_idx], $
            ratio: noaa.ratio[cg_noaa_idx], $
            time: noaa.year[cg_noaa_idx] + noaa.month[cg_noaa_idx]/13, $
            deseason: fltarr(n_elements(cg_noaa_idx))}


;OGI:
cg_ogi_idx = where(ogi.lat eq -42 and ogi.lon eq 145)
;reconstrut the strucure for Cape Grim OGI
;the time tag in the cg_ogi structure is decimal year
cg_ogi = { month: ogi.month[cg_ogi_idx], $
            year: ogi.year[cg_ogi_idx], $
            ratio: ogi.ratio[cg_ogi_idx], $
            time: ogi.year[cg_ogi_idx] + ogi.month[cg_ogi_idx]/13, $
            deseason: fltarr(n_elements(cg_ogi_idx))}

;process the UCI data, isolate the sites on latitude band -40 to -44
;this algorithm will treat all the data between latitude -40 and -44 as 1 site
;ignoring longitude
;change the up_bound and low_bound values to change the latitude band boundaries.
up_bound1 = -38
low_bound1 = -46
cg_band1_uci_idx = where(uci.lat ge low_bound1 and uci.lat le up_bound1)
cg_band1_uci = fltarr(4, n_elements(cg_band1_uci_idx))
cg_band1_uci[0,*] = uci.ratio[cg_band1_uci_idx]
cg_band1_uci[1,*] = uci.month[cg_band1_uci_idx]
cg_band1_uci[2,*] = uci.year[cg_band1_uci_idx]
cg_band1_uci[3,*] = uci.month[cg_band1_uci_idx]/13 + uci.year[cg_band1_uci_idx]

;process UCI new data (1984-1996)
cg_band1_new_uci_idx = where(uci_new.lat ge low_bound1 and uci_new.lat le up_bound1)
cg_band1_new_uci = fltarr(4, n_elements(cg_band1_new_uci_idx))
cg_band1_new_uci[0,*] = uci_new.ratio[cg_band1_new_uci_idx]
cg_band1_new_uci[1,*] = uci_new.month[cg_band1_new_uci_idx]
cg_band1_new_uci[2,*] = uci_new.year[cg_band1_new_uci_idx]
cg_band1_new_uci[3,*] = uci_new.month[cg_band1_new_uci_idx]/13 + uci_new.year[cg_band1_new_uci_idx]

;combine the uci 8496 data and the other UCI data into one data structure  called cg_uci
temp1 = [rotate(cg_band1_uci[1,*], 1), rotate(cg_band1_new_uci[1,*], 1)]
temp2 = [rotate(cg_band1_uci[2,*], 1), rotate(cg_band1_new_uci[2,*], 1)]
temp0 = [rotate(cg_band1_uci[0,*], 1), rotate(cg_band1_new_uci[0,*], 1)]
temp3 = [rotate(cg_band1_uci[3,*], 1), rotate(cg_band1_new_uci[3,*], 1)]

cg_uci = {month: temp1, $
          year: temp2, $
          ratio: temp0, $
          time: temp3, $
          deseason: fltarr(n_elements(cg_band1_new_uci_idx)+n_elements(cg_band1_uci_idx)) }

;deseasonalize the data
cg_noaa.deseason = deseason(cg_noaa.month, cg_noaa.year, cg_noaa.ratio)
cg_ogi.deseason = deseason(cg_ogi.month, cg_ogi.year, cg_ogi.ratio)
cg_uci.deseason = deseason(cg_uci.month, cg_uci.year, cg_uci.ratio)

;plot the pre-filtered data and the post-filtered data to compare
;temp0 - temp5 store the pre-filtered deseasonal data
temp0 = cg_uci.deseason
temp1 = cg_uci.time

temp2 = cg_noaa.deseason
temp3 = cg_noaa.time

temp4 = cg_ogi.deseason
temp5 = cg_ogi.time

;filter the data for NOAA
;data that are out of the 3 sigma bound are set as NAN value
stddev_noaa = std_dev_fgauss(cg_noaa.deseason)
outlier_idx = where(cg_noaa.deseason gt stddev_noaa*3 or cg_noaa.deseason lt -stddev_noaa*3, count)
if count gt 0 then begin
  cg_noaa.ratio[outlier_idx] = !VALUES.F_NAN
  cg_noaa.month[outlier_idx] = !VALUES.F_NAN
  cg_noaa.year[outlier_idx] = !VALUES.F_NAN
  cg_noaa.time[outlier_idx] = !VALUES.F_NAN
  cg_noaa.deseason[outlier_idx] = !VALUES.F_NAN
  print, 'Removed ', count, ' data points out of the Cape Grim NOAA data set'
endif

;filter the data for OGI
stddev_ogi = std_dev_fgauss(cg_ogi.deseason)
outlier_idx = where(cg_ogi.deseason gt stddev_ogi*3 or cg_ogi.deseason lt -stddev_ogi*3, count)
if count gt 0 then begin
  cg_ogi.ratio[outlier_idx] = !VALUES.F_NAN
  cg_ogi.month[outlier_idx] = !VALUES.F_NAN
  cg_ogi.year[outlier_idx] = !VALUES.F_NAN
  cg_ogi.time[outlier_idx] = !VALUES.F_NAN
  cg_ogi.deseason[outlier_idx] = !VALUES.F_NAN
  print, 'Removed ', count, ' data points out of the Cape Grim OGI data set'
endif

;filter the data for UCI
stddev_uci = std_dev_fgauss(cg_uci.deseason)
outlier_idx = where(cg_uci.deseason gt stddev_uci*3 or cg_uci.deseason lt -stddev_uci*3, count)
if count gt 0 then begin
  cg_uci.ratio[outlier_idx] = !VALUES.F_NAN
  cg_uci.month[outlier_idx] = !VALUES.F_NAN
  cg_uci.year[outlier_idx] = !VALUES.F_NAN
  cg_uci.time[outlier_idx] = !VALUES.F_NAN
  cg_uci.deseason[outlier_idx] = !VALUES.F_NAN
  print, 'Removed ', count, ' data points out of the', low_bound1, ' to ', up_bound1, ' lat band 1 UCI data set'
endif

;write out the data to a textfile to manually validate the means
; infile = '/home/excluded-from-backup/ethane/IDL/temp_file/capegrim_data.dat'
; openw, lun, infile, /get_lun

; printf, lun, 'NOAA'
; for i = 0, n_elements(cg_noaa.time)-1 do begin
	; printf, lun, cg_noaa.month[i], cg_noaa.year[i], cg_noaa.ratio[i]
; endfor

; printf, lun, 'UCI'
; for i = 0, n_elements(cg_uci.time)-1 do begin
	; printf, lun, cg_uci.month[i], cg_uci.year[i], cg_uci.ratio[i]
; endfor

; free_lun, lun


;uncomment this if need to see the deseasoned data.
;plot the pre-filtered data and the post-filtered data to compare
; cgDisplay, 960, 540, /free, /window
; cgPlot, cg_noaa.time, cg_noaa.ratio, /nodata, ytitle = 'Mixing ratio (pptv)', $
  ; xrange = [1980, 2017], yrange = [-500,500], xtitle = 'Time', $
  ; title = 'Cape Grim Deseasonal Data Pre-filtered vs. Post-filtered (3 sigma)', $
  ; xticklen = 1, xgridstyle = 1, ygridstyle = 1, xticks = 35
; cgPlot, temp1, temp0, /overplot, color = 'red', psym = 9 ;UCI
; cgPlot, temp3, temp2, /overplot, color = 'blue', psym = 9 ;NOAA
; cgPlot, temp5, temp4, /overplot, color = 'blue', psym = 9 ;OGI
; cgPlot, cg_uci.time, cg_uci.deseason, color = 'forest green', psym = 2, /overplot
; cgPlot, cg_noaa.time, cg_noaa.deseason, color = 'forest green', psym = 2, /overplot
; cgPlot, cg_ogi.time, cg_ogi.deseason, color = 'forest green', psym = 2, /overplot


;+
;even though the data has been removed, the dimensions of the structure doesn't
;change. This could be because of how structure works where you cannot change the
;dimensions of a structure without rebuilding it.
;-

;this section calculates the annual means for all data sets separately
;i.e.: during the overlapping period, the calculation will not merge
;the different data networks
;this section will also pull an extra latitude band of the UCI data
;and plot it with the data to compare its sensitivity to latitudes.

;pulling the new lat band 
low_bound2 = -46
up_bound2 = -38

cg_band2_uci_idx = where(uci.lat ge low_bound2 and uci.lat le up_bound2)
cg_band2_uci = fltarr(4, n_elements(cg_band2_uci_idx))
cg_band2_uci[0,*] = uci.ratio[cg_band2_uci_idx]
cg_band2_uci[1,*] = uci.month[cg_band2_uci_idx]
cg_band2_uci[2,*] = uci.year[cg_band2_uci_idx]
cg_band2_uci[3,*] = uci.month[cg_band2_uci_idx]/13 + uci.year[cg_band2_uci_idx]

;process UCI new data (1984-1996)
cg_band2_new_uci_idx = where(uci_new.lat ge low_bound2 and uci_new.lat le up_bound2)
cg_band2_new_uci = fltarr(4, n_elements(cg_band2_new_uci_idx))
cg_band2_new_uci[0,*] = uci_new.ratio[cg_band2_new_uci_idx]
cg_band2_new_uci[1,*] = uci_new.month[cg_band2_new_uci_idx]
cg_band2_new_uci[2,*] = uci_new.year[cg_band2_new_uci_idx]
cg_band2_new_uci[3,*] = uci_new.month[cg_band2_new_uci_idx]/13 + uci_new.year[cg_band2_new_uci_idx]

;combine the uci 8496 data and the other UCI data into one data structure  called cg_uci2
temp1 = [rotate(cg_band2_uci[1,*], 1), rotate(cg_band2_new_uci[1,*], 1)]
temp2 = [rotate(cg_band2_uci[2,*], 1), rotate(cg_band2_new_uci[2,*], 1)]
temp0 = [rotate(cg_band2_uci[0,*], 1), rotate(cg_band2_new_uci[0,*], 1)]
temp3 = [rotate(cg_band2_uci[3,*], 1), rotate(cg_band2_new_uci[3,*], 1)]

cg_uci2 = {month: temp1, $
          year: temp2, $
          ratio: temp0, $
          time: temp3, $
          deseason: fltarr(n_elements(cg_band2_new_uci_idx)+n_elements(cg_band2_uci_idx))}

;deseasonalize the data
cg_uci2.deseason = deseason(cg_uci2.month, cg_uci2.year, cg_uci2.ratio)

;filter the data for the second UCI band
stddev_uci = std_dev_fgauss(cg_uci2.deseason)
outlier_idx = where(cg_uci2.deseason gt stddev_uci*3 or cg_uci2.deseason lt -stddev_uci*3, count)
if count gt 0 then begin
  cg_uci2.ratio[outlier_idx] = !VALUES.F_NAN
  cg_uci2.month[outlier_idx] = !VALUES.F_NAN
  cg_uci2.year[outlier_idx] = !VALUES.F_NAN
  cg_uci2.time[outlier_idx] = !VALUES.F_NAN
  cg_uci2.deseason[outlier_idx] = !VALUES.F_NAN
  print, 'Removed ', count, ' data points out of the', low_bound2, ' to ', up_bound2, ' lat band 2 UCI data set'
endif


;make a quick plot of the overlapping period of NOAA and UCI data
;set up plot
cgDisplay, 1344, 756
cgPlot, cg_uci.time, cg_uci.ratio, /nodata, xtitle = 'Time', ytitle = 'Mixing ratio (pptv)', $
	title = 'Overlapping period between UCI and NOAA at Cape Grim', xrange = [2004, 2010], $
	yrange = [0,500]
cgPlot, cg_uci.time, cg_uci.ratio, /overplot, psym = 2, color = 'forest green', symsize = 1.5
cgPlot, cg_noaa.time, cg_noaa.ratio, /overplot, psym = 4, color = 'steelblue', symsize = 1.5
cgPlot, cg_uci2.time, cg_uci2.ratio, /overplot, psym = 5, color = 'tomato', symsize = 1.5

;there are 4 arrays that need to calculate the annual mean for:
;cg_noaa, cg_uci, cg_uci2, cg_ogi

;-----cg_noaa calculations-----
;find how many years and which years the NOAA network has data for
temp_year_noaa = cg_noaa.year[sort(cg_noaa.year)]
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
infile = '/home/excluded-from-backup/ethane/IDL/temp_file/annual_mean_cape_grim_SECTION_B'
openw, lun, infile, /get_lun


printf, lun, 'Year|Annual Average Mixing Ratio|Std Error|Number of Samples'
printf, lun, '-----------'
printf, lun, '>>>>>>>>>>NOAA report'
for i = 0, n_elements(year_noaa) - 1 do begin
	x = where(cg_noaa.year[*] eq year_noaa[i], count)
	
	;print out the available months of a year
	flt_month = cg_noaa.month[x]
	flt_month = flt_month[sort(flt_month)]
	
	;run the function to calculate the annual mean and annual error
	temp = annual_mean(cg_noaa.ratio[x], cg_noaa.month[x])
	
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

;-----cg_ogi-----
;find how many years and which years the OGI network has data for
temp_year_ogi = cg_ogi.year[sort(cg_ogi.year)]

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
	x = where(cg_ogi.year[*] eq year_ogi[i], count)
	
	;print out the available months of a year
	flt_month = cg_ogi.month[x]
	flt_month = flt_month[sort(flt_month)]
	
	;run the function to get average and the std error
	temp = annual_mean(cg_ogi.ratio[x], cg_ogi.month[x])
	
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


;-----uci band 1-----
;find how many years and which years the uci1 network has data for
temp_year_uci1 = cg_uci.year[sort(cg_uci.year)]

;since there are NAN values in the array, the sort and uniq function
;does not deal with NAN value very satisfactory, need some preliminary 
;algorithm to clean up the NAN values. If all the NAN indexes are removed
;at this step, all subsequent steps should not have problems with NAN 
notNAN = finite(temp_year_uci1)
x = where(notNAN eq 0, count)
if (count gt 0) then begin
	;remove all the rows with NAN value
	temp_year_uci1 = RemoveRows(rotate(temp_year_uci1, 1), x)
endif

year_uci1 = temp_year_uci1[uniq(temp_year_uci1)]
avg_uci1 = fltarr(n_elements(year_uci1));contains the annual mean for OGI
stderr_uci1 = fltarr(n_elements(year_uci1));contains the annual error for OGI


printf, lun, '-----------'
printf, lun, '>>>>>>>>>UCI band 1 report'
for i = 0, n_elements(year_uci1) - 1 do begin
	x = where(cg_uci.year[*] eq year_uci1[i], count)
	
	;print out the available months of a year
	flt_month = cg_uci.month[x]
	flt_month = flt_month[sort(flt_month)]

	;run the annual mean function to get the average and std error
	temp = annual_mean_uci(cg_uci.ratio[x], cg_uci.month[x])
	
	;pull the annual mean and the std error out from the temp variable
	avg_uci1[i] = temp[0]
	stderr_uci1[i] = temp[1]
	printf, lun, year_uci1[i], avg_uci1[i], stderr_uci1[i], count
	printf, lun, 'Months with data: ', flt_month[uniq(flt_month)]
	printf, lun, '----------------------------'
	
	;change the years with no data to NAN 
	if finite(avg_uci1[i]) eq 0 then begin
		year_uci1[i] = !VALUES.F_NAN
	endif
endfor
printf, lun, '>>>>>END UCI BAND 1 REPORT<<<<<'

;-----uci band 2-----
;find how many years and which years the uci1 network has data for
temp_year_uci2 = cg_uci2.year[sort(cg_uci2.year)]

;since there are NAN values in the array, the sort and uniq function
;does not deal with NAN value very satisfactory, need some preliminary 
;algorithm to clean up the NAN values. If all the NAN indexes are removed
;at this step, all subsequent steps should not have problems with NAN 
notNAN = finite(temp_year_uci2)
x = where(notNAN eq 0, count)
if (count gt 0) then begin
	;remove all the rows with NAN value
	temp_year_uci2 = RemoveRows(rotate(temp_year_uci2, 1), x)
endif

year_uci2 = temp_year_uci2[uniq(temp_year_uci2)]
avg_uci2 = fltarr(n_elements(year_uci2));contains the annual mean for OGI
stderr_uci2 = fltarr(n_elements(year_uci2));contains the annual error for OGI


printf, lun, '-----------'
printf, lun, '>>>>>>>>>UCI band 2 report'
for i = 0, n_elements(year_uci2) - 1 do begin
	x = where(cg_uci2.year[*] eq year_uci2[i], count)
	
	;print out the available months of a year
	flt_month = cg_uci2.month[x]
	flt_month = flt_month[sort(flt_month)]
	
	;run the annual mean function to get the average and std error
	temp = annual_mean_uci(cg_uci2.ratio[x], cg_uci2.month[x])
	
	avg_uci2[i] = temp[0]
	stderr_uci2[i] = temp[1]
	printf, lun, year_uci2[i], avg_uci2[i], stderr_uci2[i], count
	printf, lun, 'Months with data: ', flt_month[uniq(flt_month)]
	printf, lun, '----------------------------'

	;change the years with no data to NAN 
	if finite(avg_uci2[i]) eq 0 then begin
		year_uci2[i] = !VALUES.F_NAN
	endif	
endfor
printf, lun, '>>>>>END UCI BAND 2 REPORT<<<<<'

free_lun, lun

; plotting procedure
open_device, /ps, /color, file='temp.eps', margin=0.05, xsize = 10.0, ysize = 7.5
!x.thick=1
!y.thick=1
!p.font=1
!p.thick =0.5


;set up plot
cgPlot, year_noaa, avg_noaa, /nodata, xrange = [1980, 2016], xticklen = 1, $
	xgridstyle = 1, yrange = [100,400], xticks = 18, $
	xtitle = 'Time', ytitle = 'Mixing Ratio (pptv)', $
	title = 'Time series of Cape Grim 3 sigma filter'
cgPlot, year_noaa, avg_noaa, /overplot, psym = 4, color = 'forest green'
cgPlot, year_noaa, avg_noaa, /overplot, color = 'forest green', $
	err_ylow = stderr_noaa, err_yhigh = stderr_noaa
cgPlot, year_ogi, avg_ogi, /overplot, psym = 5, color = 'steelblue'
cgPlot, year_ogi, avg_ogi, /overplot, color = 'steelblue', $
	err_ylow = stderr_ogi, err_yhigh = stderr_ogi
cgPlot, year_uci1, avg_uci1, /overplot, psym = 2, color = 'tomato'
cgPlot, year_uci1, avg_uci1, /overplot, color = 'tomato', $
	err_ylow = stderr_uci1, err_yhigh = stderr_uci1
cgPlot, year_uci2, avg_uci2, /overplot, psym = 2, color = 'violet'
cgPlot, year_uci2, avg_uci2, /overplot, color = 'violet', $
	err_ylow = stderr_uci2, err_yhigh = stderr_uci2

cgLegend, SymColors = ['forest green', 'tomato', 'violet', 'steelblue'], $
	PSyms = [4, 2, 2, 5], Symsize = 1.5, Location = [0.7, 0.85], $
	titles = ['NOAA', 'UCI 44 to 40', 'UCI 46 to 38', 'OGI'], $
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

year_uci1 = NAN2neg999(year_uci1)
avg_uci1 = NAN2neg999(avg_uci1)
stderr_uci1 = NAN2neg999(stderr_uci1)

;write to file for time_series_ihr_BR_CP.pro to read
infile = '/home/excluded-from-backup/ethane/IDL/temp_file/annual_mean_capegrim_NoHeaders.dat'
openw, lun, infile, /get_lun

for i = 0, n_elements(year_noaa)-1 do begin
	printf, lun, '1', year_noaa[i], avg_noaa[i], stderr_noaa[i]
endfor

for i = 0, n_elements(year_uci1)-1 do begin
	printf, lun, '2', year_uci1[i], avg_uci1[i], stderr_uci1[i]
endfor

for i = 0, n_elements(year_ogi)-1 do begin
	printf, lun, '3', year_ogi[i], avg_ogi[i], stderr_ogi[i]
endfor

free_lun, lun



end
