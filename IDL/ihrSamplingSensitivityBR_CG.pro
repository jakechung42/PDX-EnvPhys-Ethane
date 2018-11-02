FUNCTION removeMonth, year, month, ratio, time

;this program will only use the months 3, 6, 9, 12 of NOAA and OGI. The following
;section will remove the irrelevant months. This will require deconstructing 
;the structure and rebuilding it.

rmIndex = where(month eq 1 or $
			month eq 2 or $
			month eq 4 or $
			month eq 5 or $
			month eq 7 or $
			month eq 8 or $
			month eq 10 or $
			month eq 11, count)

year = RemoveRows(Rotate(year, 1), rmIndex)
month = RemoveRows(Rotate(month, 1), rmIndex)
ratio = RemoveRows(Rotate(ratio, 1), rmIndex)
time = RemoveRows(Rotate(time, 1), rmIndex)

year = Rotate(year, 3)
month = Rotate(month, 3)
ratio = Rotate(ratio, 3)
time = Rotate(time, 3)

out = {ratio: ratio, $
		year: year, $
		month: month, $
		time: time, $
		deseason: fltarr(n_elements(ratio))}

return, out

end

;======================================================================
FUNCTION rmOutliers, ratio, year, month, time, deseason

;this function removes the outliers using the function std_dev_fgauss
;outliers are replace by NaN values

;change variable bound to change how many standard deviation
bound = 3

std_dev = std_dev_fgauss(deseason)

;get index
outlier_idx = where(deseason gt std_dev*bound or deseason lt (std_dev*(-1)*bound), count)

;remove the outliers
if count gt 0 then begin
	ratio[outlier_idx] = !VALUES.F_NAN
	year[outlier_idx] = !VALUES.F_NAN
	month[outlier_idx] = !VALUES.F_NAN
	time[outlier_idx] = !VALUES.F_NAN
	deseason[outlier_idx] = !VALUES.F_NAN
endif
;reconstruct the structure 
out = {ratio: ratio, $
		year: year, $
		month: month, $
		time: time, $
		deseason: deseason}
		
return, out

end

;======================================================================
FUNCTION annual_mean, ratio_vec, month_vec

;this function calculates the annual mean of the data
;despite the fact that the name of this function is annual mean, it only
;calculate the mean of one year when the function is called. So to actually
;calculate monthly means, this function would have to be looped through all the years

;if there's not enough data for a year, the function will return NaN
EnoughData = 1
annual_avg = !VALUES.F_NAN
annual_err = !VALUES.F_NAN

;get months that are available
hasMonth = month_vec[sort(month_vec)]
hasMonth = hasMonth[uniq(hasMonth)]
month = [3, 6, 9, 12]

;check to see if all data are available in all 4 months
if n_elements(hasMonth) eq 4 then begin 
	for j = 0, 3 do begin
		if (hasMonth[j] ne month[j]) then begin
			EnoughData = 0
			break
		endif
	endfor	
endif else EnoughData = 0

;calculate the annual mean and annual standard error
if EnoughData eq 1 then begin ;only proceed if have data in 3 6 9 and 12
	month_avg = fltarr(4) ;create array to store monthly averages of a year
	month_err = fltarr(4) ;create array to store monthly error of a year
	for i = 0, 3 do begin
		x = where(month_vec eq month[i], count0) ;loop through each month
		month_avg[i] = mean(ratio_vec[x]) ;take the average of all the month in one year
		;calculate standard error of each month. If a month only has one value
		;the stddev output will be NAN. The NAN value will be taken care of later in the script
		month_err[i] = stddev(ratio_vec[x], /DOUBLE)/sqrt(count0)
	endfor 
	
	;find which month that has only one value by looking at the NAN values in the monthly error
	notNAN = finite(month_err)
	NAN_idx = where(notNAN eq 0, count1)
	if count1 gt 0 then begin
		;the standard error for the month that only has one value will be the 
		;average of the standard error of the entire year.
		print, 'Number of months with one data point: ', count1 ;print out how many months in a year
		;with only 1 data point
		month_err[NAN_idx] = mean(month_err, /NAN)
	endif
	annual_avg = mean(month_avg) ;the average of a year is the average of the monthly averages
	annual_err = 1/4.0 * sqrt(total(month_err^2)) ;propagation of error to determine the year error
endif 
out = [annual_avg, annual_err] ;combine the average and error for output
return, out

end

;======================================================================
FUNCTION rmNaN, inputArr

;this function remove NaN values from a vector 
;first locate the NaN values in the vector
notNaN_idx = finite(inputArr) ;0 is NaN, 1 is finite
out = []
for i = 0, n_elements(inputArr)-1 do begin
	if notNaN_idx[i] eq 1 then begin
		temp = inputArr[i]
		out = [out, temp]
	endif
endfor

return, out

end

;======================================================================
FUNCTION reqNetwork, inputArr, n
;just a quick function for moving the NOAA, UCI and OGI data out
;from the input array
;input n is either 1 for NOAA, 2 for UCI or 3 for OGI

;fetch the index for the network corresponding to n
x = where(inputArr[0, *] eq n)
;fetch the data using the fetched index
out = inputArr[1 : 3, x]
;return the output 
return, out

end

;======================================================================
FUNCTION neg999toNAN, array

;convert -999 values to NAN values
neg999 = where(array eq -999, count)
if (count gt 0) then begin
	array[neg999] = !VALUES.F_NAN
endif

return, array
end
	
;======================================================================
PRO ihrSamplingSensitivityBR_CG

;this program looks at the sensitivity of the IHR of Barrow and Cape Grim when 
;the data is calculated using the same months as UCI for NOAA and OGI (3, 6, 9, 12)
;for consistency. The 4-month results will be plotted to compare with the 12-month
;result to look at how sensitive is IHR when the different months are used. 

;ihrSamplingSensitivity is similar to this program but it uses global data

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

;OGI:
br_ogi_idx = where(ogi.lat eq 71.16 and ogi.lon eq -156.5)
;reconstrut the strucure for Brrow OGI
;the time tag in the br_ogi structure is decimal year
br_ogi = { month: ogi.month[br_ogi_idx], $
            year: ogi.year[br_ogi_idx], $
            ratio: ogi.ratio[br_ogi_idx], $
            time: ogi.year[br_ogi_idx] + ogi.month[br_ogi_idx]/13, $
            deseason: fltarr(n_elements(br_ogi_idx))}

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

;get Cape Grim data from NOAA, OGI and UCI networks
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
          deseason: fltarr(n_elements(cg_band1_new_uci_idx)+n_elements(cg_band1_uci_idx))}

;remove irrelavent months from NOAA and OGI data.
br_noaa = removeMonth(br_noaa.year, br_noaa.month, br_noaa.ratio, br_noaa.time)
br_ogi = removeMonth(br_ogi.year, br_ogi.month, br_ogi.ratio, br_ogi.time)

cg_noaa = removeMonth(cg_noaa.year, cg_noaa.month, cg_noaa.ratio, cg_noaa.time)
cg_ogi = removeMonth(cg_ogi.year, cg_ogi.month, cg_ogi.ratio, cg_ogi.time)

print, 'Barrow NOAA:'
help, br_noaa, /str
print, 'Barrow UCI:'
help, br_uci, /str
print, 'Barrow OGI: '
help, br_ogi, /str
print, 'Cape Grim NOAA: '
help, cg_noaa, /str
print, 'Cape Grim UCI: '
help, cg_uci, /str
print, 'Cape Grim OGI: '
help, cg_ogi, /str

;deseaon the data and store it in the array
br_noaa.deseason = deseason(br_noaa.month, br_noaa.year, br_noaa.ratio)
br_ogi.deseason = deseason(br_ogi.month, br_ogi.year, br_ogi.ratio)
br_uci.deseason = deseason(br_uci.month, br_uci.year, br_uci.ratio)
cg_noaa.deseason = deseason(cg_noaa.month, cg_noaa.year, cg_noaa.ratio)
cg_ogi.deseason = deseason(cg_ogi.month, cg_ogi.year, cg_ogi.ratio)
cg_uci.deseason = deseason(cg_uci.month, cg_uci.year, cg_uci.ratio)

;remove outliers 
br_noaa = rmOutliers(br_noaa.ratio, br_noaa.year, br_noaa.month, br_noaa.time, br_noaa.deseason)
br_uci = rmOutliers(br_uci.ratio, br_uci.year, br_uci.month, br_uci.time, br_uci.deseason)
br_ogi = rmOutliers(br_ogi.ratio, br_ogi.year, br_ogi.month, br_ogi.time, br_ogi.deseason)
cg_noaa = rmOutliers(cg_noaa.ratio, cg_noaa.year, cg_noaa.month, cg_noaa.time, cg_noaa.deseason)
cg_uci = rmOutliers(cg_uci.ratio, cg_uci.year, cg_uci.month, cg_uci.time, cg_uci.deseason)
cg_ogi = rmOutliers(cg_ogi.ratio, cg_ogi.year, cg_ogi.month, cg_ogi.time, cg_ogi.deseason)

;print out the data array for infomation sake
print, '>>>>Barrow NOAA:'
help, br_noaa, /str
print, '>>>>Barrow UCI:'
help, br_uci, /str
print, '>>>>Barrow OGI: '
help, br_ogi, /str
print, '>>>>Cape Grim NOAA: '
help, cg_noaa, /str
print, '>>>>Cape Grim UCI: '
help, cg_uci, /str
print, '>>>>Cape Grim OGI: '
help, cg_ogi, /str

;calculate the annual mean by looping through each year of the data array for each
;data set
;create a time array that contains the years of a data set
temp = br_noaa.year[sort(br_noaa.year)]
noaaBR_year = temp[uniq(temp)]
temp = br_ogi.year[sort(br_ogi.year)]
ogiBR_year = temp[uniq(temp)]
temp = br_uci.year[sort(br_uci.year)]
uciBR_year = temp[uniq(temp)]
temp = cg_ogi.year[sort(cg_ogi.year)]
ogiCG_year = temp[uniq(temp)]
temp = cg_noaa.year[sort(cg_noaa.year)]
noaaCG_year = temp[uniq(temp)]
temp = cg_uci.year[sort(cg_uci.year)]
uciCG_year = temp[uniq(temp)]

;remove NaN values from the year arrays, this NaN stuff is such a headache at sometime
noaaBR_year = rmNaN(noaaBR_year)
ogiBR_year = rmNaN(ogiBR_year)
uciBR_year = rmNaN(uciBR_year)
ogiCG_year = rmNaN(ogiCG_year)
noaaCG_year = rmNaN(noaaCG_year)
uciCG_year = rmNaN(uciCG_year)

;create arrays to store the annual mean and error for different data sets
noaaBR_mean = fltarr(2, n_elements(noaaBR_year))
ogiBR_mean = fltarr(2, n_elements(ogiBR_year))
uciBR_mean = fltarr(2, n_elements(uciBR_year))
noaaCG_mean = fltarr(2, n_elements(noaaCG_year))
ogiCG_mean = fltarr(2, n_elements(ogiCG_year))
uciCG_mean = fltarr(2, n_elements(uciCG_year))

;starting looping through each year 
;Barrow NOAA
print, 'Barrow NOAA: mixing ratio and standard error'
for i = 0, n_elements(noaaBR_year)-1 do begin
	print, noaaBR_year[i]
	x = where(br_noaa.year eq noaaBR_year[i], count)
	noaaBR_mean[*, i]= annual_mean(br_noaa.ratio[x], br_noaa.month[x])
	print, ' ', noaaBR_mean[*, i]
endfor 
;Barrow UCI
print, 'Barrow UCI: mixing ratio and standard error'
for i = 0, n_elements(uciBR_year)-1 do begin
	print, uciBR_year[i]
	x = where(br_uci.year eq uciBR_year[i], count)
	uciBR_mean[*, i]= annual_mean(br_uci.ratio[x], br_uci.month[x])
	print, ' ', uciBR_mean[*, i]
endfor 
;Barrow OGI
print, 'Barrow OGI: mixing ratio and standard error'
for i = 0, n_elements(ogiBR_year)-1 do begin
	print, ogiBR_year[i]
	x = where(br_ogi.year eq ogiBR_year[i], count)
	ogiBR_mean[*, i]= annual_mean(br_ogi.ratio[x], br_ogi.month[x])
	print, ' ', ogiBR_mean[*, i]
endfor 
;Cape Grim NOAA
print, 'Cape Grim NOAA: mixing ratio and standard error'
for i = 0, n_elements(noaaCG_year)-1 do begin
	print, noaaCG_year[i]
	x = where(cg_noaa.year eq noaaCG_year[i], count)
	noaaCG_mean[*, i]= annual_mean(cg_noaa.ratio[x], cg_noaa.month[x])
	print, ' ', noaaCG_mean[*, i]
endfor 
;Cape Grim UCI
print, 'Cape Grim UCI: mixing ratio and standard error'
for i = 0, n_elements(uciCG_year)-1 do begin
	print, uciCG_year[i]
	x = where(cg_uci.year eq uciCG_year[i], count)
	uciCG_mean[*, i]= annual_mean(cg_uci.ratio[x], cg_uci.month[x])
	print, ' ', uciCG_mean[*, i]
endfor 
;Cape Grim OGI
print, 'Cape Grim OGI: mixing ratio and standard error'
for i = 0, n_elements(ogiCG_year)-1 do begin
	print, ogiCG_year[i]
	x = where(cg_ogi.year eq ogiCG_year[i], count)
	ogiCG_mean[*, i]= annual_mean(cg_ogi.ratio[x], cg_ogi.month[x])
	print, ' ', ogiCG_mean[*, i]
endfor 

;in order to calculate IHR, the Barrow and Cape Grim data array of a data network
;must match, the following algorithm calculates the IHR and make sure that the 
;arrays match
;resize the NOAA Barrow and CG
noaaYear = [noaaCG_year, noaaBR_year]
noaaYear = noaaYear[sort(noaaYear)]
noaaYear = noaaYear[uniq(noaaYear)]
noaaIHR = fltarr(n_elements(noaaYear))
noaaIHR_err = fltarr(n_elements(noaaYear))
for i = 0, n_elements(noaaYear)-1 do begin
	a = where(noaaBR_year eq noaaYear[i], count0)
	b = where(noaaCG_year eq noaaYear[i], count1)
	if count0 && count1 then begin 
		noaaIHR[i] = noaaBR_mean[0, a]/noaaCG_mean[0, b]
		noaaIHR_err[i] = sqrt((noaaBR_mean[1, a]/noaaCG_mean[0, b])^2 + (noaaBR_mean[0, a]*noaaCG_mean[1, b]/noaaCG_mean[0, b]^2)^2)
	endif else begin
		noaaIHR[i] = !VALUES.F_NAN
		noaaIHR_err[i] = !VALUES.F_NAN
	endelse
endfor 

;resize the UCI Barrow and CG
uciYear = [uciCG_year, uciBR_year]
uciYear = uciYear[sort(uciYear)]
uciYear = uciYear[uniq(uciYear)]
uciIHR = fltarr(n_elements(uciYear))
uciIHR_err = fltarr(n_elements(uciYear))
for i = 0, n_elements(uciYear)-1 do begin
	a = where(uciBR_year eq uciYear[i], count0)
	b = where(uciCG_year eq uciYear[i], count1)
	if count0 && count1 then begin 
		uciIHR[i] = uciBR_mean[0, a]/uciCG_mean[0, b]
		uciIHR_err[i] = sqrt((uciBR_mean[1, a]/uciCG_mean[0, b])^2 + (uciBR_mean[0, a]*uciCG_mean[1, b]/uciCG_mean[0, b]^2)^2)
	endif else begin
		uciIHR[i] = !VALUES.F_NAN
		uciIHR_err[i] = !VALUES.F_NAN
	endelse
endfor 


;reseize for OGI
;resize the UCI Barrow and CG
ogiYear = [ogiCG_year, ogiBR_year]
ogiYear = ogiYear[sort(ogiYear)]
ogiYear = ogiYear[uniq(ogiYear)]
ogiIHR = fltarr(n_elements(ogiYear))
ogiIHR_err = fltarr(n_elements(ogiYear))
for i = 0, n_elements(ogiYear)-1 do begin
	a = where(ogiBR_year eq ogiYear[i], count0)
	b = where(ogiCG_year eq ogiYear[i], count1)
	if count0 && count1 then begin 
		ogiIHR[i] = ogiBR_mean[0, a]/ogiCG_mean[0, b]
		ogiIHR_err[i] = sqrt((ogiBR_mean[1, a]/ogiCG_mean[0, b])^2 + (ogiBR_mean[0, a]*ogiCG_mean[1, b]/ogiCG_mean[0, b]^2)^2)
	endif else begin
		ogiIHR[i] = !VALUES.F_NAN
		ogiIHR_err[i] = !VALUES.F_NAN
	endelse
endfor 

;read the exported data the Barrow and Cape Grim IHR plot to plot with the sensitivity plot
;specigy the directory
infile1 = '/home/excluded-from-backup/ethane/IDL/temp_file/time_series_ihr_BR_CP.dat'
;create an array to store the data
inputIHR = fltarr(4, file_lines(infile1))
;open to read
openr, lun1, infile1, /get_lun
;read the entire file, and put it in inputIHR variable
readf, lun1, inputIHR
;release the logical parser 
free_lun, lun1

;read the annual average of Barrow and Cape Grim
;specify directorey
infile = '/home/excluded-from-backup/ethane/IDL/temp_file/annual_mean_barrow_NoHeaders.dat'
;create an array to store the data
inputBR = fltarr(4, file_lines(infile))
;open to read
openr, lun, infile, /get_lun
;read the fine and store the data in inputBR
readf, lun, inputBR
;release
free_lun, lun

;read the annual average of Cape Grim
infile = '/home/excluded-from-backup/ethane/IDL/temp_file/annual_mean_capegrim_NoHeaders.dat'
;create an array to store the data
inputCG = fltarr(4, file_lines(infile))
;open to read
openr, lun, infile, /get_lun
;read the fine and store the data in inputCG
readf, lun, inputCG
;release
free_lun, lun

;convert negative 999 in the inputBR and inputCG data to NaN
inputBR = neg999toNAN(inputBR)
inputCG = neg999toNAN(inputCG)

;seperate the NOAA, UCI and OGI data out from each input array
fullNoaaIHR = reqNetwork(inputIHR, 1)
fullUciIHR = reqNetwork(inputIHR, 2)
fullOgiIHR = reqNetwork(inputIHR, 3)

fullNoaaBR = reqNetwork(inputBR, 1)
fullUciBR = reqNetwork(inputBR, 2)
fullOgiBR = reqNetwork(inputBR, 3)

fullNoaaCG = reqNetwork(inputCG, 1)
fullUciCG = reqNetwork(inputCG, 2)
fullOgiCG = reqNetwork(inputCG, 3)

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
cgPlot, noaaBR_year, xrange = [1982, 2016], xticklen = 1, xgridstyle = 1, $
	xticks = 17, /nodata, yrange = [1000,2400], ytitle = 'Mixing ratio(pptv)', $
	title = 'Time series of Barrow and Cape Grim along with IHR (Mar, Jun, Sep, Dec). UCI pulls data from lat 38 to 46 south'
cgPlot, noaaBR_year, noaaBR_mean[0, *], /overplot, err_yhigh = noaaBR_mean[1, *], err_ylow = noaaBR_mean[1, *], $
	psym = 5, color = 'steelblue'
cgPlot, uciBR_year, uciBR_mean[0, *], /overplot, err_yhigh = uciBR_mean[1, *], err_ylow = uciBR_mean[1, *], $
	psym = 2, color = 'forest green'
cgPlot, ogiBR_year, ogiBR_mean[0, *], /overplot, err_yhigh = ogiBR_mean[1, *], err_ylow = ogiBR_mean[1, *], $
	psym = 4, color = 'red'
	
;plot full-sampled barrow
cgPlot, fullNoaaBR[0, *], fullNoaaBR[1, *], /overplot, err_yhigh = fullNoaaBR[2, *], err_ylow = fullNoaaBR[2, *], $
	psym = 5, color = 'blue violet'
cgPlot, fullUciBR[0, *], fullUciBR[1, *], /overplot, err_yhigh = fullUciBR[2, *], err_ylow = fullUciBR[2, *], $
	psym = 2, color = 'forest green'
cgPlot, fullOgiBR[0, *], fullOgiBR[1, *], /overplot, err_yhigh = fullOgiBR[2, *], err_ylow = fullOgiBR[2, *], $
	psym = 4, color = 'violet'
multiplot, /doyaxis, /doxaxis

;plot Cape Grim
cgPlot, noaaCG_year, /nodata, xrange = [1982, 2016], xticklen = 1, xgridstyle = 1, $
	xticks = 17, XTickformat='(A1)', yrange = [150,350], ytitle = 'Mixing ratio(pptv)'
cgPlot, noaaCG_year, noaaCG_mean[0, *], /overplot, err_yhigh = noaaCG_mean[1, *], err_ylow = noaaCG_mean[1, *], $
	psym = 5, color = 'steelblue'
cgPlot, uciCG_year, uciCG_mean[0, *], /overplot, err_yhigh = uciCG_mean[1, *], err_ylow = uciCG_mean[1, *], $
	psym = 2, color = 'forest green'
cgPlot, ogiCG_year, ogiCG_mean[0, *], /overplot, err_yhigh = ogiCG_mean[1, *], err_ylow = ogiCG_mean[1, *], $
	psym = 4, color = 'red'
	
;plot full-sampled Cape Grim
cgPlot, fullNoaaCG[0, *], fullNoaaCG[1, *], /overplot, err_yhigh = fullNoaaCG[2, *], err_ylow = fullNoaaCG[2, *], $
	psym = 5, color = 'blue violet'
cgPlot, fullUciCG[0, *], fullUciCG[1, *], /overplot, err_yhigh = fullUciCG[2, *], err_ylow = fullUciCG[2, *], $
	psym = 2, color = 'forest green'
cgPlot, fullOgiCG[0, *], fullOgiCG[1, *], /overplot, err_yhigh = fullOgiCG[2, *], err_ylow = fullOgiCG[2, *], $
	psym = 4, color = 'violet'
multiplot, /doyaxis, /doxaxis

;plot IHR
cgPlot, noaaYear, /nodata, xtitle = 'Years', ytitle = 'IHR', $
	xrange = [1982, 2016], yrange = [4, 11], $
	xticklen = 1, xgridstyle = 1, xticks = 17
cgPlot, noaaYear, noaaIHR, /overplot, psym = 5, color = 'steelblue', err_yhigh = noaaIHR_err, $
	err_ylow = noaaIHR_err
cgPlot, uciYear, uciIHR, /overplot, psym = 2, color = 'forest green', err_yhigh = uciIHR_err, $
	err_ylow = uciIHR_err
cgPlot, ogiYear, ogiIHR, /overplot, psym = 4, color = 'red', err_ylow = ogiIHR_err, $
	err_yhigh = ogiIHR_err

;plot full-sampled IHR
cgPlot, fullNoaaIHR[0, *], fullNoaaIHR[1, *], /overplot, psym = 5, color = 'blue violet', err_yhigh = fullNoaaIHR[2, *], $
	err_ylow = fullNoaaIHR[2, *]
cgPlot, fullUciIHR[0, *], fullUciIHR[1, *], /overplot, psym = 2, color = 'forest green', err_yhigh = fullUciIHR[2, *], $
	err_ylow = fullUciIHR[2, *]
cgPlot, fullOgiIHR[0, *], fullOgiIHR[1, *], /overplot, psym = 4, color = 'violet', err_ylow = fullOgiIHR[2, *], $
	err_yhigh = fullOgiIHR[2, *]


AL_Legend, ['4-month OGI', '4-month UCI', '4-month NOAA', 'full OGI', 'full NOAA'], $
	psym = [4, 2, 5, 4, 5], linestyle = [0, 0, 0, 0, 0], box = 1, $
	position = [1991, 10.5], color = ['red', 'forest green', 'steelblue', 'violet', 'blue violet'], $
	background_color = 'rose'
close_device

spawn, 'gv temp.eps'

multiplot, /default


end



