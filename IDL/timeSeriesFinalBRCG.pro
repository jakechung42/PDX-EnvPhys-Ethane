;this program plot the annual mean time series of Barrow and Cape Grim from all
;3 data networks and compare it with the simulated data. 
;This program is built based on ihrSamplingSensitivityBR_CG.pro and simDataAnalysisBR_CG.pro

;>>>>>>>>>>[declare subfunctions]<<<<<<<<<<
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

;this function calculates the annual mean of the data by first calculating the
;mnea of the average of each month to ensure that the data is not biased by the months with more 
;data
;despite the fact that the name of this function is annual mean, it only
;calculate the mean of one year when the function is called. So to actually
;calculate annual means, this function would have to be looped through all the years

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
FUNCTION annualSimMean, inputArr, year

;the output of GEOS-Chem is a 3 dimensional array with time, lon, lat
;corresponding to x, y, z. The time dimension for our simulations 
;is monthly average. We need annual average calculated from March, 
;June, September and December.
;This program takes the 3D-array and calculate the annual average

;make sure that the data has 12 months for every year
modulo = n_elements(inputArr(*, 0, 0)) mod n_elements(year)
if modulo ne 0 then begin
	print, 'Do not have enough months, stop the program!!!'
	stop	;if there's not enough months, stop the program
endif
;start calculating the annual average
;first, need to pull out only the months Mar, Jun, Sep, Dec from the array
;shortArr will contain the months of Mar, Jun, Sep, and Dec only
shortArr = fltarr(n_elements(inputArr[*, 0, 0])/3, n_elements(inputArr[0, *, 0]), n_elements(inputArr[0, 0, *]))
for n = 0, n_elements(shortArr[*, 0, 0])-1 do begin
	shortArr[n, *, *] = inputArr[3 * n + 2, *, *]
endfor
;make an array to store the annual average 
out = fltarr(n_elements(year), n_elements(shortArr[0, *, 0]), n_elements(shortArr[0, 0, *]))
for i = 0, n_elements(year)-1 do begin
	;take every 4 indeces out to calculate the annual mean
	tempArr = shortArr[4 * i : 4 * i + 3, *, *] 
	out[i, *, *] = mean(tempArr, 1)/2 * 1000 ;adjust to pptv
endfor
return, out ;return value

end

;======================================================================
FUNCTION getTimeSeriesSite, masterArr, lat, lon

;this function takes the input as the three dimensional simulated annual mean array and uses
;the input lon and lat to return the time series of the input site

;employ the ctm_index function to obtain the index of the input lon lat in the array
CTM_INDEX, CTM_TYPE('GEOS1', RESOLUTION= 2), i, j, CENTER= [lat, lon], /non_interactive
out = masterArr[*, i-1, j-1] ;get the annual mean time series of the input coordinate.

return, out

end

;======================================================================
FUNCTION getSimRawData, filename

;this function take the input as the path to the simulated data and read the 
;data then output the 3D data array

;read in the file
ctm_get_data, datainfo, 'IJ-AVG-$', tracer=1, filename=filename, tau0=tau0
;Get MODELINFO and GRIDINFO structures, xmid, ymid hold lon/lat centers
getmodelandgridinfo, datainfo[0], modelinfo, gridinfo
;nt is the number of data points in the simulation.
nt = n_elements(datainfo)
tempSimArr = fltarr(nt, 144, 91)
;sim_yymmdd stores the entire time element of the input simulation
sim_ymd = tau2yymmdd(tau0, /GEOS1)
;this section call out the data blocks of the simulated data and store it in simArr
;simArr is a 3-D array with the following attributes:
;simArr[ *, 0, 0 ]: the time dimension
;simArr[ 0, *, 0 ]: longitudes
;simArr[ 0, 0, * ]: latitudes
for i = 0, nt - 1 do begin
	data= CTM_EXTRACT(*(datainfo[i].data), modelinfo= modelinfo, $
	gridinfo= gridinfo, lat= [-90, 90], lon= [-180, 180], alrange = [1, 2], $
	average= 4)
	for k = 0, n_elements(data[* ,0 ]) - 1 do begin
		for j = 0, n_elements(data[0 ,* ]) - 1 do begin
			tempSimArr[i, k, j] = data[k, j]
		endfor
	endfor
endfor
;clean up memory
CTM_CLEANUP 

return, tempSimArr

end

;======================================================================
FUNCTION rmAbsSim, input

;this short function removes the absolute value to show only the trend of the simulated data
;the input variable should be a mixing ratio vector
avgVar = mean(input)
out = input[*] - avgVar

return, out 

end
;>>>>>>>>>>[main program]<<<<<<<<<<
PRO timeSeriesFinalBRCG

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

;remove months from NOAA and OGI data. 
;to be consistent with the UCI data, NOAA and OGI networks will only use data from Mar, Jun, Sep, Dec
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

;since we don't care much about the absolute value of ethane ratio but more about the trend, the following section
;takes the average of the data and subtract it from each data point to compare the trend of the observed
;and simulated data.
allBR = [rotate(noaaBR_mean[0, *], 3), rotate(uciBR_mean[0, *], 3), rotate(ogiBR_mean[0, *], 3)]
allCG = [rotate(noaaCG_mean[0, *], 3), rotate(uciCG_mean[0, *], 3), rotate(ogiCG_mean[0, *], 3)]
all_IHR = [noaaIHR, uciIHR, ogiIHR]
;get the average for all data 
allBR_mean = mean(allBR, /nan)
allCG_mean = mean(allCG, /nan)
all_IHR_mean = mean(all_IHR, /nan)
;remove absolute value and just show the trend
noaaBR_mean[0, *] = noaaBR_mean[0, *] - allBR_mean
uciBR_mean[0, *] = uciBR_mean[0, *] - allBR_mean
ogiBR_mean[0, *] = ogiBR_mean[0, *] - allBR_mean
noaaCG_mean[0, *] = noaaCG_mean[0, *] - allCG_mean
uciCG_mean[0, *] = uciCG_mean[0, *] - allCG_mean
ogiCG_mean[0, *] = ogiCG_mean[0, *] - allCG_mean
;same process for IHR values
noaaIHR = noaaIHR - all_IHR_mean
uciIHR = uciIHR - all_IHR_mean
ogiIHR = ogiIHR - all_IHR_mean
print, noaaBR_mean
;>>>>>>>>>>>[Process Simulated data]<<<<<<<<<<
;sim_title = 'PSU emissions scaled to Xiao et al over 1996-2003' Sce E
filename1 = "/home/excluded-from-backup/data/C2H6/trac_avg.PSUSF_1981_2015.bpch"

;Sce D
filename2 = "/home/excluded-from-backup/data/C2H6/trac_avg.AydinSF_1981_2015.bpch"
;sim_title = 'Constant default base emissions'  Sce A
filename3 = "/home/excluded-from-backup/data/C2H6/trac_avg.spinup_GFED4_MAVG_1981_2015.198101010000"
;sim_title = 'Aydin et al. no normalization to Xiao et al.'  Sce C
filename4 = "/home/excluded-from-backup/data/C2H6/trac_avg.AydinAbsSF_1981_2015.bpch"
;sim_title = 'Simpson et al. no normalization to Xiao et al.' Sce B
filename5 = "/home/excluded-from-backup/data/C2H6/trac_avg.SimpsonSF_1981_2015.bpch"
;sim_title = 'Unscaled PSU emission, MER3BB, MER18FF' Sce F
filename6 = "/home/excluded-from-backup/data/C2H6/trac_avg.PSU_MER3BB_MER18FF_1981_2015.bpch"


;read in the file
ctm_get_data, datainfo, 'IJ-AVG-$', tracer=1, filename=filename1, tau0=tau0
;sim_yymmdd stores the entire time element of the input simulation
sim_ymd = tau2yymmdd(tau0, /GEOS1)
;get the years of the simulated data, should be the same for all simulation
simYear = sim_ymd.year[sort(sim_ymd.year)]
simYear = simYear[uniq(simYear)]

;get raw data for first sim and store in sim1 to sim 6
sim1 = getSimRawData(filename1)
sim2 = getSimRawData(filename2)
sim3 = getSimRawData(filename3)
sim4 = getSimRawData(filename4)
sim5 = getSimRawData(filename5)
sim6 = getSimRawData(filename6)
;sim1 is a 3D array, [year, longitudes, latitudes]
;dimenstion for sim1 is [34, 144, 91] 34 is 34 years, 144 and 91 are the index for lon and lat 
;since the simulation is a 2x2.5 degree grid.
sim1 = annualSimMean(sim1, simYear)
sim2 = annualSimMean(sim2, simYear)
sim3 = annualSimMean(sim3, simYear)
sim4 = annualSimMean(sim4, simYear)
sim5 = annualSimMean(sim5, simYear)
sim6 = annualSimMean(sim6, simYear)

;obtain the annual mean time series for Barrow
barrow1 = getTimeSeriesSite(sim1, 71.31, -156.6)
barrow2 = getTimeSeriesSite(sim2, 71.31, -156.6)
barrow3 = getTimeSeriesSite(sim3, 71.31, -156.6)
barrow4 = getTimeSeriesSite(sim4, 71.31, -156.6)
barrow5 = getTimeSeriesSite(sim5, 71.31, -156.6)
barrow6 = getTimeSeriesSite(sim6, 71.31, -156.6)
;obtain the annual mean time series for Cape Grim 
capegrim1 = getTimeSeriesSite(sim1, -40.66, 144.66)
capegrim2 = getTimeSeriesSite(sim2, -40.66, 144.66)
capegrim3 = getTimeSeriesSite(sim3, -40.66, 144.66)
capegrim4 = getTimeSeriesSite(sim4, -40.66, 144.66)
capegrim5 = getTimeSeriesSite(sim5, -40.66, 144.66)
capegrim6 = getTimeSeriesSite(sim6, -40.66, 144.66)
;calculate simulated ihr
sim1ihr = barrow1/capegrim1
sim2ihr = barrow2/capegrim2
sim3ihr = barrow3/capegrim3
sim4ihr = barrow4/capegrim4
sim5ihr = barrow5/capegrim5
sim6ihr = barrow6/capegrim6
;remove absolute values and just show the trend
barrow1 = rmAbsSim(barrow1)
barrow2 = rmAbsSim(barrow2)
barrow3 = rmAbsSim(barrow3)
barrow4 = rmAbsSim(barrow4)
barrow5 = rmAbsSim(barrow5)
barrow6 = rmAbsSim(barrow6)
;same for cape grim
capegrim1 = rmAbsSim(capegrim1)
capegrim2 = rmAbsSim(capegrim2)
capegrim3 = rmAbsSim(capegrim3)
capegrim4 = rmAbsSim(capegrim4)
capegrim5 = rmAbsSim(capegrim5)
capegrim6 = rmAbsSim(capegrim6)
;same for ihr
sim1ihr = rmAbsSim(sim1ihr)
sim2ihr = rmAbsSim(sim2ihr)
sim3ihr = rmAbsSim(sim3ihr)
sim4ihr = rmAbsSim(sim4ihr)
sim5ihr = rmAbsSim(sim5ihr)
sim6ihr = rmAbsSim(sim6ihr)
;>>>>>>>>>>>[Plotting Procedure]<<<<<<<<<<
;set up plot
open_device, /ps, /color, file='temp.eps', margin=0.05, xsize = 10.0, ysize = 7.5
!x.thick=1
!y.thick=1
!p.font=1
!p.thick =0.5
multiplot, /default    ; resets multiplot settings
multiplot, [1,3], ygap=0.002, xgap=0;  sets up multiplot 

;plot barrow
cgPlot, noaaBR_year, xrange = [1982, 2015], xticklen = 1, xgridstyle = 1, $
	xticks = 17, /nodata, yrange = [-300,700], ytitle = 'Mixing ratio(pptv)', $
	title = 'Time series of Barrow and Cape Grim along with IHR (Mar, Jun, Sep, Dec). UCI pulls data from lat 38 to 46 south'
cgPlot, noaaBR_year, noaaBR_mean[0, *], /overplot, err_yhigh = noaaBR_mean[1, *], err_ylow = noaaBR_mean[1, *], $
	psym = 5, color = 'black'
cgPlot, uciBR_year, uciBR_mean[0, *], /overplot, err_yhigh = uciBR_mean[1, *], err_ylow = uciBR_mean[1, *], $
	psym = 2, color = 'black'
cgPlot, ogiBR_year, ogiBR_mean[0, *], /overplot, err_yhigh = ogiBR_mean[1, *], err_ylow = ogiBR_mean[1, *], $
	psym = 4, color = 'black'
;plot Simulated Barrow
cgPlot, simYear, barrow1, /overplot, linestyle = 2, color = 'red', thick = 1
cgPlot, simYear, barrow2, /overplot, linestyle = 0, color = 'steelblue', thick = 1
cgPlot, simYear, barrow3, /overplot, linestyle = 0, color = 'violet', thick = 1
cgPlot, simYear, barrow4, /overplot, linestyle = 2, color = 'steelblue', thick = 1
cgPlot, simYear, barrow5, /overplot, linestyle = 0, color = 'forest green', thick = 1
cgPlot, simYear, barrow6, /overplot, linestyle = 0, color = 'red', thick = 1

multiplot, /doyaxis, /doxaxis

;plot Cape Grim 
cgPlot, noaaCG_year, /nodata, xrange = [1982, 2015], xticklen = 1, xgridstyle = 1, $
	xticks = 17, XTickformat='(A1)', yrange = [-100,100], ytitle = 'Mixing ratio(pptv)'
cgPlot, noaaCG_year, noaaCG_mean[0, *], /overplot, err_yhigh = noaaCG_mean[1, *], err_ylow = noaaCG_mean[1, *], $
	psym = 5, color = 'black'
cgPlot, uciCG_year, uciCG_mean[0, *], /overplot, err_yhigh = uciCG_mean[1, *], err_ylow = uciCG_mean[1, *], $
	psym = 2, color = 'black'
cgPlot, ogiCG_year, ogiCG_mean[0, *], /overplot, err_yhigh = ogiCG_mean[1, *], err_ylow = ogiCG_mean[1, *], $
	psym = 4, color = 'black'
;plot Simulated Cape Grim
cgPlot, simYear, capegrim1, /overplot, linestyle = 2, color = 'red', thick = 1
cgPlot, simYear, capegrim2, /overplot, linestyle = 0, color = 'steelblue', thick = 1
cgPlot, simYear, capegrim3, /overplot, linestyle = 0, color = 'violet', thick = 1
cgPlot, simYear, capegrim4, /overplot, linestyle = 2, color = 'steelblue', thick = 1
cgPlot, simYear, capegrim5, /overplot, linestyle = 0, color = 'forest green', thick = 1
cgPlot, simYear, capegrim6, /overplot, linestyle = 0, color = 'red', thick = 1

multiplot, /doyaxis, /doxaxis

;plot IHR 
cgPlot, noaaYear, /nodata, xtitle = 'Years', ytitle = 'IHR', $
	xrange = [1982, 2015], yrange = [-3, 3], $
	xticklen = 1, xgridstyle = 1, xticks = 17
cgPlot, noaaYear, noaaIHR, /overplot, psym = 5, color = 'black', err_yhigh = noaaIHR_err, $
	err_ylow = noaaIHR_err
cgPlot, uciYear, uciIHR, /overplot, psym = 2, color = 'black', err_yhigh = uciIHR_err, $
	err_ylow = uciIHR_err
cgPlot, ogiYear, ogiIHR, /overplot, psym = 4, color = 'black', err_ylow = ogiIHR_err, $
	err_yhigh = ogiIHR_err

;plot simulated ihr
cgPlot, simYear, sim1ihr, /overplot, linestyle = 2, color = 'red', thick = 1
cgPlot, simYear, sim2ihr, /overplot, linestyle = 0, color = 'steelblue', thick = 1
cgPlot, simYear, sim3ihr, /overplot, linestyle = 0, color = 'violet', thick = 1
cgPlot, simYear, sim4ihr, /overplot, linestyle = 2, color = 'steelblue', thick = 1
cgPlot, simYear, sim5ihr, /overplot, linestyle = 0, color = 'forest green', thick = 1
cgPlot, simYear, sim6ihr, /overplot, linestyle = 0, color = 'red', thick = 1

AL_Legend, ['Sce A', 'Sce B', 'Sce C', 'Sce D', 'Sce E', 'Sce F'], $
	psym = [0, 0, 0, 0, 0, 0], linestyle = [0, 0, 2, 0, 2, 0], box = 0, $
	position = [1992, 3], color = ['violet', 'forest green', 'steelblue', 'steelblue', 'red', 'red'], $
	charsize = 1
close_device

spawn, 'gv temp.eps'

multiplot, /default


end
