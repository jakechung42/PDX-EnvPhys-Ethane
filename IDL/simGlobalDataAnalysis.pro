;this program analyzes the simulated data on the global scale to compare with with the
;global observed data. To be consistent with the analysis, the sim data will use calculations
;that are similar to the observed data calculations except for the error calculations,
;propagation of error, 3-sigma filter

FUNCTION MakeLatBin, bin_bound, lat, month, year, ratio

;this function takes in the above inputs, then run through an algorithm to bin the latitudes
;as specicfied by the bin_bound input. Latitude bands that is not specified in the bin_bound 
;array will be discarded.
;band id will be set as the mid latitude of that lat band.
;the longitudes data is irrelevant in this calculation. So we can ignore the longitudes.
;we are going to make a new array and convert the latitudes to their appropriate bands

bin_lat = []
bin_month = []
bin_year = []
bin_ratio = []
;algorithm to move the latitude into bins
;n_elements(bin_bound)-2 because we don't need it to run to the end. It just how the algorithm works
for i = 0, n_elements(bin_bound)-2 do begin
	x = where(lat ge bin_bound[i] and lat lt bin_bound[i + 1], count)
	if (count gt 0) then begin
		temp0 = (bin_bound[i] + bin_bound[i + 1])/2
		temp1 = fltarr(n_elements(x))
		temp1[*] = temp0
		bin_lat = [bin_lat, temp1]
		bin_month = [bin_month, month[x]]
		bin_year = [bin_year,year[x]]
		bin_ratio = [bin_ratio, ratio[x]]
	endif
endfor

bin_out = fltarr(4, n_elements(bin_ratio))
bin_out[0, *] = bin_lat
bin_out[1, *] = bin_month
bin_out[2, *] = bin_year
bin_out[3, *] = bin_ratio

return, bin_out
end


;======================================================================
FUNCTION DeseasonBin, bin_mid, data_array

;this function uses the mid latitude of each band to loop through and deseasonlize the entire lat band
;the input data array should have 4 columns that is the output of the MakeLatBin function

DeseasonData = []

for i = 0, n_elements(bin_mid)-1 do begin
	x = where(data_array[0, *] eq bin_mid[i])
	temp = deseason(data_array[1, x], data_array[2, x], data_array[3, x])
	temp = Rotate(temp, 3)
	DeseasonData = [DeseasonData, temp]
endfor

return, DeseasonData

end
	
	
;======================================================================
FUNCTION rmOutliers, bin_mid, DeseasonData, data_array

;this function removes data that is more than 3 sigmas
;data_array input should be a matrix with 4 columns 

;change the varable bound to change how many sigmas to remove the data
bound = 3

for i = 0, n_elements(bin_mid) - 1 do begin
	print, ' bin mid lat: ', bin_mid[i]
	x = where(data_array[0, *] eq bin_mid[i], count1)
	print, 'n bin: ', count1
	std_dev = std_dev_fgauss(DeseasonData[x])
	outlier = where(DeseasonData[x] gt std_dev*bound or DeseasonData[x] lt (std_dev*(-1)*bound), count2)
	print, 'n outlier: ', count2
	DeseasonData[outlier] = !VALUES.F_NAN
	data_array[0, outlier] = !VALUES.F_NAN
	data_array[1, outlier] = !VALUES.F_NAN	
	data_array[2, outlier] = !VALUES.F_NAN	
	data_array[3, outlier] = !VALUES.F_NAN	
	;just change to float to get the percent data removed
	count1 = float(count1)
	count2 = float(count2)
	print, 'percent removed: ', count2/count1*100, '%'
endfor

;rotate to combine the deseason data to the main data array
DeseasonData = Rotate(DeseasonData, 1)
data_array = [data_array, DeseasonData]

;remove NAN values complately out of the arrays because some functions don't work 
;well with NAN
notNAN = finite(data_array[0, *])
notNAN_idx = where(notNAN eq 0)

data_array = RemoveRows(data_array, notNAN_idx)

return, data_array

end

;======================================================================
FUNCTION takeAnnualMean, bin_mid, input_data


;this function will calculate the annual mean of each latitude band.
;the annual mean is calculated using the months 3, 6, 9, 12 averages of a year
;the uncertainty of that year is calculated using propagation of error of the 
;months that contribute to that year.
;this function is complicated so be careful before changing anything
;the input data should be a matrix with 5 columns, which is the output of the 
;function rmOutliers.
;the output is a matrix with 4 columns that contains the latitude band, the year, 
;the average of that year and the standard error of that year
;calculating the annual mean of each latitude band by calculating the
;averages of each month of a year first.
annual_lat = []
annual_year = []
annual_ratio = []
annual_error = []
monthALL4 = [3, 6, 9, 12]
for i = 0, n_elements(bin_mid)-1 do begin
	print, 'Latitude: ', bin_mid[i]
	x = where(input_data[0, *] eq bin_mid[i], count0)
	;temp_arr0 contains one single latitude band data
	temp_arr0 = input_data[*, x]
	temp1 = temp_arr0[2, sort(temp_arr0[2, *])]
	avaiYear = temp1[uniq(temp1)]
	;loop through data for calculation for each year
	for j = 0, n_elements(avaiYear)-1 do begin
		;temp_arr1 contains a single year of a single latitude band
		y = where(temp_arr0[2, *] eq avaiYear[j], count1)
		temp_arr1 = temp_arr0[*, y]
		;pulling available months
		temp1 = temp_arr1[1, sort(temp_arr1[1, *])]
		avaiMonth = temp1[uniq(temp1)]
		;check to see if all 4 months are available(3, 6, 9, 12)
		for r = 0, n_elements(monthALL4)-1 do begin
			haveAllMonth = 1
			p = where(avaiMonth eq monthALL4[r], count4)
			if (count4 ne 1) then begin
				haveAllMonth = 0
				break
			endif
		endfor 
		;if all 4 months are available
		month_avg = fltarr(4)
		month_stderr = fltarr(4)
		if haveAllMonth then begin	
			;loop through to each month
			for k = 0, n_elements(monthALL4)-1 do begin
				z = where(temp_arr1[1, *] eq monthALL4[k], count2)
				;check to see if there's enough data for each month
				;shouldn't be a problem but gonna convert count2 to float to take the standrd error just to be sure.
				count2 = double(count2)
				if (count2 gt 1) then begin
					;third column is ratio
					month_avg[k] = mean(temp_arr1[3, z])
					;fourth column is deseasoned data
					month_stderr[k] = stddev(temp_arr1[4, z])/sqrt(count2)
				endif else begin
					;if only have 1 value for a month
					print, 'Year: ', avaiYear[j], ', month: ', monthALL4[k], ' only has one value'
					month_avg[k] = temp_arr1[3, z]
					month_stderr[k] = !VALUES.F_NAN
				endelse
			endfor
			;now for the months with only 1 data point, replace the stderr
			;of that month with the average standard error of all the months
			;of that year
			notNAN = finite(month_stderr)
			NAN_idx = where(notNAN eq 0, count3)
			if (count3 gt 0) then begin
				month_stderr[NAN_idx] = mean(month_stderr, /NAN)
			endif
		endif else begin
			;if not all 4 months are available, scrap the entire year
			print, 'Year: ', avaiYear[j], ' do not have enough months, remove year.'
			print, 'Available months of that year: ', avaiMonth
			month_avg[*] = !VALUES.F_NAN
			month_stderr[*] = !VALUES.F_NAN
		endelse
		annual_lat = [annual_lat, bin_mid[i]]
		annual_year = [annual_year, avaiYear[j]]
		annual_ratio = [annual_ratio, mean(month_avg)]
		;propagation of error from each month

		year_error = 1/4.0 * sqrt(total(month_stderr^2))
		annual_error = [annual_error, year_error]
	endfor
endfor		

out_arr = fltarr(4, n_elements(annual_lat))

;for i = 0, n_elements(annual_lat)-1 do begin
;	print, annual_lat[i], annual_year[i], annual_ratio[i], annual_error[i]
;endfor

out_arr[0, *] = annual_lat
out_arr[1, *] = annual_year
out_arr[2, *] = annual_ratio
out_arr[3, *] = annual_error

return, out_arr

end

;======================================================================
FUNCTION takeSimAnnualMean, bin_mid, input_data

;this function is for the simulated data. It is similar to the function
;takeAnnualMean but without the error calculation
;this function will calculate the annual mean of each latitude band.
;the annual mean is calculated using the months 3, 6, 9, 12 averages of a year
;this function is complicated so be careful before changing anything
;the output is a matrix with 3 columns that contains the latitude band, the year, 
;the simulated average of that year 
;calculating the annual mean of each latitude band by calculating the
;averages of each month of a year first.
annual_lat = []
annual_year = []
annual_ratio = []
monthALL4 = [3, 6, 9, 12]
for i = 0, n_elements(bin_mid)-1 do begin
	print, 'Latitude: ', bin_mid[i]
	x = where(input_data[0, *] eq bin_mid[i], count0)
	;temp_arr0 contains one single latitude band data
	temp_arr0 = input_data[*, x]
	temp1 = temp_arr0[2, sort(temp_arr0[2, *])]
	avaiYear = temp1[uniq(temp1)]
	;loop through data for calculation for each year
	for j = 0, n_elements(avaiYear)-1 do begin
		;temp_arr1 contains a single year of a single latitude band
		y = where(temp_arr0[2, *] eq avaiYear[j], count1)
		temp_arr1 = temp_arr0[*, y]
		;pulling available months
		temp1 = temp_arr1[1, sort(temp_arr1[1, *])]
		avaiMonth = temp1[uniq(temp1)]
		;check to see if all 4 months are available(3, 6, 9, 12)
		for r = 0, n_elements(monthALL4)-1 do begin
			haveAllMonth = 1
			p = where(avaiMonth eq monthALL4[r], count4)
			if (count4 ne 1) then begin
				haveAllMonth = 0
				break
			endif
		endfor 
		;if all 4 months are available
		month_avg = fltarr(4)
		if haveAllMonth then begin	
			;loop through to each month
			for k = 0, n_elements(monthALL4)-1 do begin
				z = where(temp_arr1[1, *] eq monthALL4[k], count2)
				month_avg[k] = mean(temp_arr1[3, z])
			endfor
		endif else begin
			;if not all 4 months are available, scrap the entire year
			print, 'Year: ', avaiYear[j], ' do not have enough months, remove year.'
			print, 'Available months of that year: ', avaiMonth
			month_avg[*] = !VALUES.F_NAN
		endelse
		annual_lat = [annual_lat, bin_mid[i]]
		annual_year = [annual_year, avaiYear[j]]
		annual_ratio = [annual_ratio, mean(month_avg)]
	endfor
endfor		

out_arr = fltarr(4, n_elements(annual_lat))

out_arr[0, *] = annual_lat
out_arr[1, *] = annual_year
out_arr[2, *] = annual_ratio

return, out_arr

end
;======================================================================
FUNCTION getNorth, input_array

;this function get the northern hemisphere latitude band of an annual data array 
;that is an output of the function takeAnnualMean
x = where(input_array[0, *] gt 0)

out_nor = input_array[*, x]

return, out_nor

end

;======================================================================
FUNCTION getSouth, input_array

;this function get the southern hemisphere latitude band of an annual data array 
;that is an output of the function takeAnnualMean
x = where(input_array[0, *] lt 0)

out_sou = input_array[*, x]

return, out_sou

end

;======================================================================
FUNCTION fitYear, allYear, annual_data

;this function performs a fit year algorithm. Because the northern and the 
;southern hemisphere don't gaurantee to have the same years so this algorithm
;put the data into a year set to prevent the wrong year data being calculated.

;this function uses the input allYear to get all the year that a single 
;data set has (noaa.year, ogi.year, uci.year), the annual data is the data of
;one single latitude that will be fitted with all the years that a data set has

temp = allYear[sort(allYear)]
year_series = temp[uniq(temp)]

;sort by year
idx_sort = sort(annual_data[1, *])
annual_data[*, *] = annual_data[*, idx_sort]


out_array = fltarr(4, n_elements(year_series))
;performing the fit year algorithm
for i = 0, n_elements(year_series)-1 do begin
	x = where(annual_data[1, *] eq year_series[i], count)

	if (count eq 0) then begin
		out_array[0, i] = annual_data[0, 0]
		out_array[1, i] = year_series[i]
		out_array[2, i] = !VALUES.F_NAN
		out_array[3, i] = !VALUES.F_NAN
	endif else begin
		out_array[0, i] = annual_data[0, 0]
		out_array[1, i] = year_series[i]
		out_array[2, i] = annual_data[2, x]
		out_array[3, i] = annual_data[3, x]
	endelse
endfor

return, out_array

end

;======================================================================
FUNCTION sind, angle

;this function takes sine value in degrees

angle_rad = angle*!PI/180

out = sin(angle_rad)

return, out

end

;======================================================================
FUNCTION weightedMean, bin_mid, bin_bound, hem_annual

;this function will calculate the temporal hemispheric mean with the latitudinal weights 
;the input data should be the annual data of just one single hemisphere
temp = hem_annual[1, sort(hem_annual[1, *])]
year_series = temp[uniq(temp)]
out = fltarr(3, n_elements(year_series))

weight = fltarr(n_elements(bin_mid))

for j = 0, n_elements(bin_mid)-1 do begin
	;check hemisphere to make the correct weighting fomular.
	if (bin_mid[j] lt 0) then begin
		weight[j] = abs((sind(bin_bound[j+1]) - sind(bin_bound[j]))/sind(bin_bound[0]))
	endif else if (bin_mid[j] gt 0) then begin
		weight[j] = abs((sind(bin_bound[j+1]) - sind(bin_bound[j]))/sind(bin_bound[n_elements(bin_bound)-1]))
	endif
endfor

;replace the latitude values with the weights
for i = 0, n_elements(hem_annual[0, *])-1 do begin
	x = where(hem_annual[0, i] eq bin_mid[*])
	hem_annual[0, i] = weight[x]
endfor
		
;calculate the temporal hemispheric mean
for i = 0, n_elements(year_series)-1 do begin
	x = where(hem_annual[1, *] eq year_series[i])
	temp_arr = hem_annual[*, x]
	sum = 0
	error = 0
	for j = 0, n_elements(temp_arr[0, *])-1 do begin
		sum = sum + temp_arr[2, j]*temp_arr[0, j]
		error = error + (temp_arr[3, j]^2)*(temp_arr[0, j]^2)
	endfor
	out[0, i] = year_series[i]
	out[1, i] = sum
	out[2, i] = sqrt(error)

end

return, out

end
;======================================================================
FUNCTION removeMonth, year, month, ratio, lat, lon

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
			
print, 'Number of Irrelevant months: ', count

year = RemoveRows(year, rmIndex)
month = RemoveRows(month, rmIndex)
ratio = RemoveRows(ratio, rmIndex)
lat = RemoveRows(lat, rmIndex)
lon = RemoveRows(lon, rmIndex)

out = {ratio: ratio, $
		lat: lat, $
		lon: lon, $
		year: year, $
		month: month}

return, out

end

;======================================================================
FUNCTION normData, input

;this short function normalizes the data to compare the trend
avgVar = mean(input, /NAN)
out = input[*] - avgVar

return, out

end

;======================================================================
FUNCTION calSimGlobal, path, obsArr, bin_bound, bin_mid

;This function calculates the simulated global averages and IHR from the input
;parameters for only one data network.
;description of each variable
;path: the path directory for the bpch file 
;obsArr: the observed data structure
;bin_bound: the vector that contains the latitude bin boundaries 
;bin_mid: the middle latitude of the latitude bands

;build the sim data
simArr = buildSimData(obsArr, path)
;make the latitude bins 
sim_bin = MakeLatBin(bin_bound, simArr.lat, simArr.month, simArr.year, simArr.ratio)
;calculated the annual sim average
simAnnualMean = takeSimAnnualMean(bin_mid, sim_bin)
;allign the data temporally using the fitYeat function
out = []
for q = 0, n_elements(bin_mid)-1 do begin
	x = where(simAnnualMean[0, *] eq bin_mid[q])
	temp_array = simAnnualMean[*, x]
	temp_out = fitYear(simArr.year, temp_array)
	temp_out = Rotate(temp_out, 3)
	out = [out, temp_out]
endfor 
simAnnualMean = Rotate(out, 1)

;separate the northern and southern hemispheres
souSim = getSouth(simAnnualMean)
norSim = getNorth(simAnnualMean)
;calculate hemispheric mean 
souAnnualMean = weightedMean(bin_mid, bin_bound, souSim)
norAnnualMean = weightedMean(bin_mid, bin_bound, norSim)
;calculate IHR 
annual_IHR = norAnnualMean[1, *]/souAnnualMean[1, *]
;combine everything into a matrix to output from the function
out = [norAnnualMean[0, *], norAnnualMean[1, *], souAnnualMean[1, *], annual_IHR]

return, out

end

;======================================================================
;start main program
PRO simGlobalDataAnalysis, $
	normalized = normalized

compile_opt idl2

;some important arrays:
;nor_bin contains the boundaries for the latitude bands of the Northern Hemisphere
;sou_bin contains the boundaries for the latitude bands of the Southern Hemisphere
;use negative values when addressing the southern latitudes.
;The arrays that specifies the latitude bands need to stay with latitude 0

print, 'Initiate Program========================================================='
;read in the 1984 1996 UCI data. uci_new is the data structure of the 1984 - 1996 UCI data
uci_new = read_uci8496()
;read in NOAA data
noaa = read_noaa()
;read in UCI data
uci = read_uci()
;read in OGI data
ogi = read_ogi()
;combine the UCI data and the UCI_NEW data sets
temp = {ratio: [uci.ratio, uci_new.ratio], $
		lat: [uci.lat, uci_new.lat], $
		lon: [uci.lon, uci_new.lon], $
		year: [uci.year, uci_new.year], $
		month: [uci.month, uci_new.month]}

uci = temp
;remove irrelevant months
;need to rotate the OGI array so the removeMonth function could work
tempYear = rotate(ogi.year, 1)
tempMonth = rotate(ogi.month, 1)
tempRatio = rotate(ogi.ratio, 1)
tempLat = rotate(ogi.lat, 1)
tempLon = rotate(ogi.lon, 1)
noaa = removeMonth(noaa.year, noaa.month, noaa.ratio, noaa.lat, noaa.lon)
ogi = removeMonth(tempYear, tempMonth, tempRatio, tempLat, tempLon)

;change bin_bound to change the latitude boundaries. Order from low latitude to high latitude
;southern latitudes are negative
bin_bound = [-50, -30, 0, 30, 50, 75]
bin_bound = float(bin_bound)
;number of bands for each hemisphere
n_bin = n_elements(bin_bound)-1 

print, 'Latitude boundaries: ', bin_bound

bin_noaa = MakeLatBin(bin_bound, noaa.lat, noaa.month, noaa.year, noaa.ratio)
bin_uci = MakeLatBin(bin_bound, uci.lat, uci.month, uci.year, uci.ratio)
bin_ogi = MakeLatBin(bin_bound, ogi.lat, ogi.month, ogi.year, ogi.ratio)
;make an array of the mid lattitude of each bin
bin_mid = fltarr(n_elements(bin_bound)-1)
for i = 0, n_elements(bin_bound)-2 do begin
	bin_mid[i] = (bin_bound[i] + bin_bound[i + 1])/2
endfor
;deseason each lat bin
uci_deseason = DeseasonBin(bin_mid, bin_uci)
ogi_deseason = DeseasonBin(bin_mid, bin_ogi)
noaa_deseason = DeseasonBin(bin_mid, bin_noaa)
;at this point, the deseason arrays are a vector but the bin arrays are 
;matrices. Just be aware!!!

;remove data with more than 3 sigmas
print, '-----NOAA: '
bin_noaa = rmOutliers(bin_mid, noaa_deseason, bin_noaa)
print, '-----UCI: '
bin_uci = rmOutliers(bin_mid, uci_deseason, bin_uci)
print, '-----OGI: '
bin_ogi = rmOutliers(bin_mid, ogi_deseason, bin_ogi)
;calculate annual mean
annual_noaa = takeAnnualMean(bin_mid, bin_noaa)
annual_ogi = takeAnnualMean(bin_mid, bin_ogi)
annual_uci = takeAnnualMean(bin_mid, bin_uci)
;the three most important arrays are annual_uci, annual_noaa, annual_ogi
;these three arrays contain the annual mean of each latitude band of each data set

;the next step is to make the weighting algorithm to calculate the hemispheric 
;annual mean and to put in the correct propagation of error for the annual
;mean uncertainty

;the problem with this is that we can never be sure if the starting year and 
;the ending year of each latitude band are the same also that we cannot be sure,
;if the sampling is continuous every year. 
;To solve those problems, we match the data to a known year series, whichever year
;the data is missing, we can just assign a NAN value, and since the year series
;is known, we also know the starting and ending year.
;the following section fits the years for each latitude and reconstruct the annual 
;array
;NOAA array
out = []
for q = 0, n_elements(bin_mid)-1 do begin
	x = where(annual_noaa[0, *] eq bin_mid[q])
	temp_array = annual_noaa[*, x]
	temp_out = fitYear(noaa.year, temp_array)
	temp_out = Rotate(temp_out, 3)
	out = [out, temp_out]
endfor 

annual_noaa = Rotate(out, 1)

;UCI array
out = []
for q = 0, n_elements(bin_mid)-1 do begin
	x = where(annual_uci[0, *] eq bin_mid[q])
	temp_array = annual_uci[*, x]
	temp_out = fitYear(uci.year, temp_array)
	temp_out = Rotate(temp_out, 3)
	out = [out, temp_out]
endfor 

annual_uci = Rotate(out, 1)

;OGI array
out = []
for q = 0, n_elements(bin_mid)-1 do begin
	x = where(annual_ogi[0, *] eq bin_mid[q])
	temp_array = annual_ogi[*, x]
	temp_out = fitYear(ogi.year, temp_array)
	temp_out = Rotate(temp_out, 3)
	out = [out, temp_out]
endfor 

annual_ogi = Rotate(out, 1)

;break the annual array into 2 smaller arrays contains the northern hemisphere data 
;and southern hemisphere data
sou_uci = getSouth(annual_uci)
nor_uci = getNorth(annual_uci)
sou_noaa = getSouth(annual_noaa)
nor_noaa = getNorth(annual_noaa)
sou_ogi = getSouth(annual_ogi)
nor_ogi = getNorth(annual_ogi) 
;calculate the northern hemispheric mean with weights and everything
annual_sou_noaa = weightedMean(bin_mid, bin_bound, sou_noaa)
annual_nor_noaa = weightedMean(bin_mid, bin_bound, nor_noaa)

annual_sou_uci = weightedMean(bin_mid, bin_bound, sou_uci)
annual_nor_uci = weightedMean(bin_mid, bin_bound, nor_uci)

annual_sou_ogi = weightedMean(bin_mid, bin_bound, sou_ogi)
annual_nor_ogi = weightedMean(bin_mid, bin_bound, nor_ogi)

annual_ihr_noaa = annual_nor_noaa[1, *]/annual_sou_noaa[1, *]
annual_ihr_uci = annual_nor_uci[1, *]/annual_sou_uci[1, *]
annual_ihr_ogi = annual_nor_ogi[1, *]/annual_sou_ogi[1, *]
;calculating the uncertainty for the IHR
annual_ihr_noaa_err = sqrt((annual_nor_noaa[2, *]/annual_sou_noaa[1, *])^2 + $
	(annual_nor_noaa[1, *]*annual_sou_noaa[2, *]/annual_sou_noaa[1, *]^2)^2)
annual_ihr_uci_err = sqrt((annual_nor_uci[2, *]/annual_sou_uci[1, *])^2 + $
	(annual_nor_uci[1, *]*annual_sou_uci[2, *]/annual_sou_uci[1, *]^2)^2)
annual_ihr_ogi_err = sqrt((annual_nor_ogi[2, *]/annual_sou_ogi[1, *])^2 + $
	(annual_nor_ogi[1, *]*annual_sou_ogi[2, *]/annual_sou_ogi[1, *]^2)^2)
	
;getting the simulated data 

;Sce A: default emissions based on Xiao et al. (2008)
SceA_path = "/home/excluded-from-backup/data/C2H6/trac_avg.spinup_GFED4_MAVG_1981_2015.198101010000"

;Sce B: emissions based on budget from Simpson et al. (2012)
SceB_path = "/home/excluded-from-backup/data/C2H6/trac_avg.SimpsonSF_1981_2015.bpch"

;Sce C: Firn C2H6 based on Aydin et al. (2011)
SceC_path = "/home/excluded-from-backup/data/C2H6/trac_avg.AydinSF_1981_2015.bpch"

;Sce D: same as C but normailized to A
SceD_path = "/home/excluded-from-backup/data/C2H6/trac_avg.AydinAbsSF_1981_2015.bpch"

;Sce E: Rice et al. (2016) emission scenario normailized to A
SceE_path = "/home/excluded-from-backup/data/C2H6/trac_avg.PSUSF_1981_2015.bpch"

;Sce F: same as E but normalized to B 
SceF_path = "/home/excluded-from-backup/data/C2H6/trac_avg.PSU_MER3BB_MER18FF_1981_2015.bpch"

;calling function calSimGlobal to calculate IHR and sim average from simData for all networks and scenarios
;output of the sim arrays are year, nor avg, sou avg, IHR

;Sce A:
SceA_noaa = calSimGlobal(SceA_path, noaa, bin_bound, bin_mid)
SceA_uci = calSimGlobal(SceA_path, uci, bin_bound, bin_mid)
SceA_ogi = calSimGlobal(SceA_path, ogi, bin_bound, bin_mid)

;Sce B:
SceB_noaa = calSimGlobal(SceB_path, noaa, bin_bound, bin_mid)
SceB_uci = calSimGlobal(SceB_path, uci, bin_bound, bin_mid)
SceB_ogi = calSimGlobal(SceB_path, ogi, bin_bound, bin_mid)

;Sce C:
SceC_noaa = calSimGlobal(SceC_path, noaa, bin_bound, bin_mid)
SceC_uci = calSimGlobal(SceC_path, uci, bin_bound, bin_mid)
SceC_ogi = calSimGlobal(SceC_path, ogi, bin_bound, bin_mid)

;Sce D:
SceD_noaa = calSimGlobal(SceD_path, noaa, bin_bound, bin_mid)
SceD_uci = calSimGlobal(SceD_path, uci, bin_bound, bin_mid)
SceD_ogi = calSimGlobal(SceD_path, ogi, bin_bound, bin_mid)

;Sce E:
SceE_noaa = calSimGlobal(SceE_path, noaa, bin_bound, bin_mid)
SceE_uci = calSimGlobal(SceE_path, uci, bin_bound, bin_mid)
SceE_ogi = calSimGlobal(SceE_path, ogi, bin_bound, bin_mid)

;Sce F:
SceF_noaa = calSimGlobal(SceF_path, noaa, bin_bound, bin_mid)
SceF_uci = calSimGlobal(SceF_path, uci, bin_bound, bin_mid)
SceF_ogi = calSimGlobal(SceF_path, ogi, bin_bound, bin_mid)

;export the simulated data of Scenario A for sensitivity study in simDataSensitivity.pro
infile = "/home/excluded-from-backup/ethane/IDL/temp_file/SceF_SenStudy_allNetworks.dat"

;prepare data for exporting 
;give each data network an ID (1, 2, 3) so that when read can easily seperate the different network from 1 array
one = fltarr(n_elements(SceF_ogi[0, *]))
one[*] = 1
one = rotate(one, 1) 
outSceF_ogi = [one, SceF_ogi]

two = fltarr(n_elements(SceF_uci[0, *]))
two[*] = 2
two = rotate(two, 1) 
outSceF_uci = [two, SceF_uci]

three = fltarr(n_elements(SceF_noaa[0, *]))
three[*] = 3
three = rotate(three, 1) 
outSceF_noaa = [three, SceF_noaa]
openw, lun1, infile, /get_lun
printf, lun1, outSceF_ogi
printf, lun1, outSceF_uci
printf, lun1, outSceF_noaa

free_lun, lun1

;if normalized is set
if KEYWORD_SET(normalized) then begin
;plotting procedure
;set up plot
open_device, /ps, /color, file='temp.eps', margin=0.05, xsize = 10.0, ysize = 7.5
!x.thick=1
!y.thick=1
!p.font=1
!p.thick =0.5

multiplot, /default    ; resets multiplot settings
multiplot, [1,3], ygap=0.002, xgap=0;  sets up multiplot 

;plot northern hemisphere
;yrange for normalized -200, 300
cgPlot, annual_nor_noaa[0, *], annual_nor_noaa[1, *], xrange = [1982, 2016], xticklen = 1, xgridstyle = 1, $
	xticks = 17, /nodata, yrange = [-200,300], ytitle = 'Mixing ratio(pptv)', $
	title = 'Time series of global Ethane NOAA, OGI sampled in Mar, Jun, Sep, Dec compared with full data'
cgPlot, annual_nor_noaa[0, *], normData(annual_nor_noaa[1, *]), /overplot, psym = 5, color = 'black', $
	err_yhigh = annual_nor_noaa[2, *], err_ylow = annual_nor_noaa[2, *]
cgPlot, annual_nor_uci[0, *], normData(annual_nor_uci[1, *]), /overplot, psym = 2, color = 'black', $
	err_yhigh = annual_nor_uci[2, *], err_ylow = annual_nor_uci[2, *]
cgPlot, annual_nor_ogi[0, *], normData(annual_nor_ogi[1, *]), /overplot, psym = 4, color = 'black', $
	err_yhigh = annual_nor_ogi[2, *], err_ylow = annual_nor_ogi[2, *]
	
;plot the sim data
;Sce A:
cgPlot, SceA_noaa[0, *], normData(SceA_noaa[1, *]), /overplot, psym = 5, color = 'blue violet', $
	symsize = 0.75
cgPlot, SceA_noaa[0, *], normData(SceA_noaa[1, *]), /overplot, linestyle = 0, color = 'blue violet' ;add line
cgPlot, SceA_uci[0, *], normData(SceA_uci[1, *]), /overplot, psym = 2, color = 'blue violet', $
	symsize = 0.75
cgPlot, SceA_uci[0, *], normData(SceA_uci[1, *]), /overplot, linestyle = 0, color = 'blue violet';add line
cgPlot, SceA_ogi[0, *], normData(SceA_ogi[1, *]), /overplot, psym = 4, color = 'blue violet', $
	symsize = 0.75
cgPlot, SceA_ogi[0, *], normData(SceA_ogi[1, *]), /overplot, linestyle = 0, color = 'blue violet';add line

;Sce B:
cgPlot, SceB_noaa[0, *], normData(SceB_noaa[1, *]), /overplot, psym = 5, color = 'blue', $
	symsize = 0.75
cgPlot, SceB_noaa[0, *], normData(SceB_noaa[1, *]), /overplot, linestyle = 0, color = 'blue' ;add line
cgPlot, SceB_uci[0, *], normData(SceB_uci[1, *]), /overplot, psym = 2, color = 'blue', $
	symsize = 0.75
cgPlot, SceB_uci[0, *], normData(SceB_uci[1, *]), /overplot, linestyle = 0, color = 'blue';add line
cgPlot, SceB_ogi[0, *], normData(SceB_ogi[1, *]), /overplot, psym = 4, color = 'blue', $
	symsize = 0.75
cgPlot, SceB_ogi[0, *], normData(SceB_ogi[1, *]), /overplot, linestyle = 0, color = 'blue';add line

;Sce C:
cgPlot, SceC_noaa[0, *], normData(SceC_noaa[1, *]), /overplot, psym = 5, color = 'forest green', $
	symsize = 0.75
cgPlot, SceC_noaa[0, *], normData(SceC_noaa[1, *]), /overplot, linestyle = 0, color = 'forest green' ;add line
cgPlot, SceC_uci[0, *], normData(SceC_uci[1, *]), /overplot, psym = 2, color = 'forest green', $
	symsize = 0.75
cgPlot, SceC_uci[0, *], normData(SceC_uci[1, *]), /overplot, linestyle = 0, color = 'forest green';add line
cgPlot, SceC_ogi[0, *], normData(SceC_ogi[1, *]), /overplot, psym = 4, color = 'forest green', $
	symsize = 0.75
cgPlot, SceC_ogi[0, *], normData(SceC_ogi[1, *]), /overplot, linestyle = 0, color = 'forest green';add line

;Sce D:
cgPlot, SceD_noaa[0, *], normData(SceD_noaa[1, *]), /overplot, psym = 5, color = 'lime green', $
	symsize = 0.75
cgPlot, SceD_noaa[0, *], normData(SceD_noaa[1, *]), /overplot, linestyle = 0, color = 'lime green' ;add line
cgPlot, SceD_uci[0, *], normData(SceD_uci[1, *]), /overplot, psym = 2, color = 'lime green', $
	symsize = 0.75
cgPlot, SceD_uci[0, *], normData(SceD_uci[1, *]), /overplot, linestyle = 0, color = 'lime green';add line
cgPlot, SceD_ogi[0, *], normData(SceD_ogi[1, *]), /overplot, psym = 4, color = 'lime green', $
	symsize = 0.75
cgPlot, SceD_ogi[0, *], normData(SceD_ogi[1, *]), /overplot, linestyle = 0, color = 'lime green';add line

;Sce E:
cgPlot, SceE_noaa[0, *], normData(SceE_noaa[1, *]), /overplot, psym = 5, color = 'red', $
	symsize = 0.75
cgPlot, SceE_noaa[0, *], normData(SceE_noaa[1, *]), /overplot, linestyle = 0, color = 'red' ;add line
cgPlot, SceE_uci[0, *], normData(SceE_uci[1, *]), /overplot, psym = 2, color = 'red', $
	symsize = 0.75
cgPlot, SceE_uci[0, *], normData(SceE_uci[1, *]), /overplot, linestyle = 0, color = 'red';add line
cgPlot, SceE_ogi[0, *], normData(SceE_ogi[1, *]), /overplot, psym = 4, color = 'red', $
	symsize = 0.75
cgPlot, SceE_ogi[0, *], normData(SceE_ogi[1, *]), /overplot, linestyle = 0, color = 'red';add line

;Sce F:
cgPlot, SceF_noaa[0, *], normData(SceF_noaa[1, *]), /overplot, psym = 5, color = 'Magenta', $
	symsize = 0.75
cgPlot, SceF_noaa[0, *], normData(SceF_noaa[1, *]), /overplot, linestyle = 0, color = 'Magenta' ;add line
cgPlot, SceF_uci[0, *], normData(SceF_uci[1, *]), /overplot, psym = 2, color = 'Magenta', $
	symsize = 0.75
cgPlot, SceF_uci[0, *], normData(SceF_uci[1, *]), /overplot, linestyle = 0, color = 'Magenta';add line
cgPlot, SceF_ogi[0, *], normData(SceF_ogi[1, *]), /overplot, psym = 4, color = 'Magenta', $
	symsize = 0.75
cgPlot, SceF_ogi[0, *], normData(SceF_ogi[1, *]), /overplot, linestyle = 0, color = 'Magenta';add line
	
multiplot, /doyaxis, /doxaxis

;plot southern hemisphere
;yrange for normalized -200, 300
cgPlot, annual_sou_noaa[0, *], annual_sou_noaa[1, *], /nodata, xrange = [1982, 2016], $
	xticklen = 1, xgridstyle = 1, xticks = 17, XTickformat='(A1)', yrange = [-200,300], $
	ytitle = 'Mixing ratio(pptv)'
cgPlot, annual_sou_noaa[0, *], normData(annual_sou_noaa[1, *]), /overplot, psym = 5, color = 'black', $
	err_yhigh = annual_sou_noaa[2, *], err_ylow = annual_sou_noaa[2, *]
cgPlot, annual_sou_uci[0, *], normData(annual_sou_uci[1, *]), /overplot, psym = 2, color = 'black', $
	err_yhigh = annual_sou_uci[2, *], err_ylow = annual_sou_uci[2, *]
cgPlot, annual_sou_ogi[0, *], normData(annual_sou_ogi[1, *]), /overplot, psym = 4, color = 'black', $
	err_yhigh = annual_sou_ogi[2, *], err_ylow = annual_sou_ogi[2, *]
	
;plot the sim data
;Sce A:
cgPlot, SceA_noaa[0, *], normData(SceA_noaa[2, *]), /overplot, psym = 5, color = 'blue violet', $
	symsize = 0.75
cgPlot, SceA_noaa[0, *], normData(SceA_noaa[2, *]), /overplot, linestyle = 0, color = 'blue violet' ;add line
cgPlot, SceA_uci[0, *], normData(SceA_uci[2, *]), /overplot, psym = 2, color = 'blue violet', $
	symsize = 0.75
cgPlot, SceA_uci[0, *], normData(SceA_uci[2, *]), /overplot, linestyle = 0, color = 'blue violet';add line
cgPlot, SceA_ogi[0, *], normData(SceA_ogi[2, *]), /overplot, psym = 4, color = 'blue violet', $
	symsize = 0.75
cgPlot, SceA_ogi[0, *], normData(SceA_ogi[2, *]), /overplot, linestyle = 0, color = 'blue violet';add line

;Sce B:
cgPlot, SceB_noaa[0, *], normData(SceB_noaa[2, *]), /overplot, psym = 5, color = 'blue', $
	symsize = 0.75
cgPlot, SceB_noaa[0, *], normData(SceB_noaa[2, *]), /overplot, linestyle = 0, color = 'blue' ;add line
cgPlot, SceB_uci[0, *], normData(SceB_uci[2, *]), /overplot, psym = 2, color = 'blue', $
	symsize = 0.75
cgPlot, SceB_uci[0, *], normData(SceB_uci[2, *]), /overplot, linestyle = 0, color = 'blue';add line
cgPlot, SceB_ogi[0, *], normData(SceB_ogi[2, *]), /overplot, psym = 4, color = 'blue', $
	symsize = 0.75
cgPlot, SceB_ogi[0, *], normData(SceB_ogi[2, *]), /overplot, linestyle = 0, color = 'blue';add line

;Sce C:
cgPlot, SceC_noaa[0, *], normData(SceC_noaa[2, *]), /overplot, psym = 5, color = 'forest green', $
	symsize = 0.75
cgPlot, SceC_noaa[0, *], normData(SceC_noaa[2, *]), /overplot, linestyle = 0, color = 'forest green' ;add line
cgPlot, SceC_uci[0, *], normData(SceC_uci[2, *]), /overplot, psym = 2, color = 'forest green', $
	symsize = 0.75
cgPlot, SceC_uci[0, *], normData(SceC_uci[2, *]), /overplot, linestyle = 0, color = 'forest green';add line
cgPlot, SceC_ogi[0, *], normData(SceC_ogi[2, *]), /overplot, psym = 4, color = 'forest green', $
	symsize = 0.75
cgPlot, SceC_ogi[0, *], normData(SceC_ogi[2, *]), /overplot, linestyle = 0, color = 'forest green';add line

;Sce D:
cgPlot, SceD_noaa[0, *], normData(SceD_noaa[2, *]), /overplot, psym = 5, color = 'lime green', $
	symsize = 0.75
cgPlot, SceD_noaa[0, *], normData(SceD_noaa[2, *]), /overplot, linestyle = 0, color = 'lime green' ;add line
cgPlot, SceD_uci[0, *], normData(SceD_uci[2, *]), /overplot, psym = 2, color = 'lime green', $
	symsize = 0.75
cgPlot, SceD_uci[0, *], normData(SceD_uci[2, *]), /overplot, linestyle = 0, color = 'lime green';add line
cgPlot, SceD_ogi[0, *], normData(SceD_ogi[2, *]), /overplot, psym = 4, color = 'lime green', $
	symsize = 0.75
cgPlot, SceD_ogi[0, *], normData(SceD_ogi[2, *]), /overplot, linestyle = 0, color = 'lime green';add line

;Sce E:
cgPlot, SceE_noaa[0, *], normData(SceE_noaa[2, *]), /overplot, psym = 5, color = 'red', $
	symsize = 0.75
cgPlot, SceE_noaa[0, *], normData(SceE_noaa[2, *]), /overplot, linestyle = 0, color = 'red' ;add line
cgPlot, SceE_uci[0, *], normData(SceE_uci[2, *]), /overplot, psym = 2, color = 'red', $
	symsize = 0.75
cgPlot, SceE_uci[0, *], normData(SceE_uci[2, *]), /overplot, linestyle = 0, color = 'red';add line
cgPlot, SceE_ogi[0, *], normData(SceE_ogi[2, *]), /overplot, psym = 4, color = 'red', $
	symsize = 0.75
cgPlot, SceE_ogi[0, *], normData(SceE_ogi[2, *]), /overplot, linestyle = 0, color = 'red';add line

;Sce F:
cgPlot, SceF_noaa[0, *], normData(SceF_noaa[2, *]), /overplot, psym = 5, color = 'Magenta', $
	symsize = 0.75
cgPlot, SceF_noaa[0, *], normData(SceF_noaa[2, *]), /overplot, linestyle = 0, color = 'Magenta' ;add line
cgPlot, SceF_uci[0, *], normData(SceF_uci[2, *]), /overplot, psym = 2, color = 'Magenta', $
	symsize = 0.75
cgPlot, SceF_uci[0, *], normData(SceF_uci[2, *]), /overplot, linestyle = 0, color = 'Magenta';add line
cgPlot, SceF_ogi[0, *], normData(SceF_ogi[2, *]), /overplot, psym = 4, color = 'Magenta', $
	symsize = 0.75
cgPlot, SceF_ogi[0, *], normData(SceF_ogi[2, *]), /overplot, linestyle = 0, color = 'Magenta';add line

multiplot, /doyaxis, /doxaxis

;plot IHR
;yrange for non-normalized -1.5, 1
cgPlot, annual_nor_noaa[1, *], annual_ihr_noaa, /nodata, xtitle = 'Years', ytitle = 'IHR', $
	xrange = [1982, 2016], yrange = [-1.5, 1], $
	xticklen = 1, xgridstyle = 1, xticks = 17
cgPlot, annual_sou_noaa[0, *], normData(annual_ihr_noaa), /overplot, psym = 5, color = 'black', $
	err_yhigh = annual_ihr_noaa_err, err_ylow = annual_ihr_noaa_err
cgPlot, annual_sou_uci[0, *], normData(annual_ihr_uci), /overplot, psym = 2, color = 'black', $
	err_yhigh = annual_ihr_uci_err, err_ylow =annual_ihr_uci_err
cgPlot, annual_sou_ogi[0, *], normData(annual_ihr_ogi), /overplot, psym = 4, color = 'black', $
	err_yhigh = annual_ihr_ogi_err, err_ylow = annual_ihr_ogi_err

;plot the sim data
;Sce A:
cgPlot, SceA_noaa[0, *], normData(SceA_noaa[3, *]), /overplot, psym = 5, color = 'blue violet', $
	symsize = 0.75
cgPlot, SceA_noaa[0, *], normData(SceA_noaa[3, *]), /overplot, linestyle = 0, color = 'blue violet' ;add line
cgPlot, SceA_uci[0, *], normData(SceA_uci[3, *]), /overplot, psym = 2, color = 'blue violet', $
	symsize = 0.75
cgPlot, SceA_uci[0, *], normData(SceA_uci[3, *]), /overplot, linestyle = 0, color = 'blue violet';add line
cgPlot, SceA_ogi[0, *], normData(SceA_ogi[3, *]), /overplot, psym = 4, color = 'blue violet', $
	symsize = 0.75
cgPlot, SceA_ogi[0, *], normData(SceA_ogi[3, *]), /overplot, linestyle = 0, color = 'blue violet';add line

;Sce B:
cgPlot, SceB_noaa[0, *], normData(SceB_noaa[3, *]), /overplot, psym = 5, color = 'blue', $
	symsize = 0.75
cgPlot, SceB_noaa[0, *], normData(SceB_noaa[3, *]), /overplot, linestyle = 0, color = 'blue' ;add line
cgPlot, SceB_uci[0, *], normData(SceB_uci[3, *]), /overplot, psym = 2, color = 'blue', $
	symsize = 0.75
cgPlot, SceB_uci[0, *], normData(SceB_uci[3, *]), /overplot, linestyle = 0, color = 'blue';add line
cgPlot, SceB_ogi[0, *], normData(SceB_ogi[3, *]), /overplot, psym = 4, color = 'blue', $
	symsize = 0.75
cgPlot, SceB_ogi[0, *], normData(SceB_ogi[3, *]), /overplot, linestyle = 0, color = 'blue';add line

;Sce C:
cgPlot, SceC_noaa[0, *], normData(SceC_noaa[3, *]), /overplot, psym = 5, color = 'forest green', $
	symsize = 0.75
cgPlot, SceC_noaa[0, *], normData(SceC_noaa[3, *]), /overplot, linestyle = 0, color = 'forest green' ;add line
cgPlot, SceC_uci[0, *], normData(SceC_uci[3, *]), /overplot, psym = 2, color = 'forest green', $
	symsize = 0.75
cgPlot, SceC_uci[0, *], normData(SceC_uci[3, *]), /overplot, linestyle = 0, color = 'forest green';add line
cgPlot, SceC_ogi[0, *], normData(SceC_ogi[3, *]), /overplot, psym = 4, color = 'forest green', $
	symsize = 0.75
cgPlot, SceC_ogi[0, *], normData(SceC_ogi[3, *]), /overplot, linestyle = 0, color = 'forest green';add line

;Sce D:
cgPlot, SceD_noaa[0, *], normData(SceD_noaa[3, *]), /overplot, psym = 5, color = 'lime green', $
	symsize = 0.75
cgPlot, SceD_noaa[0, *], normData(SceD_noaa[3, *]), /overplot, linestyle = 0, color = 'lime green' ;add line
cgPlot, SceD_uci[0, *], normData(SceD_uci[3, *]), /overplot, psym = 2, color = 'lime green', $
	symsize = 0.75
cgPlot, SceD_uci[0, *], normData(SceD_uci[3, *]), /overplot, linestyle = 0, color = 'lime green';add line
cgPlot, SceD_ogi[0, *], normData(SceD_ogi[3, *]), /overplot, psym = 4, color = 'lime green', $
	symsize = 0.75
cgPlot, SceD_ogi[0, *], normData(SceD_ogi[3, *]), /overplot, linestyle = 0, color = 'lime green';add line

;Sce E:
cgPlot, SceE_noaa[0, *], normData(SceE_noaa[3, *]), /overplot, psym = 5, color = 'red', $
	symsize = 0.75
cgPlot, SceE_noaa[0, *], normData(SceE_noaa[3, *]), /overplot, linestyle = 0, color = 'red' ;add line
cgPlot, SceE_uci[0, *], normData(SceE_uci[3, *]), /overplot, psym = 2, color = 'red', $
	symsize = 0.75
cgPlot, SceE_uci[0, *], normData(SceE_uci[3, *]), /overplot, linestyle = 0, color = 'red';add line
cgPlot, SceE_ogi[0, *], normData(SceE_ogi[3, *]), /overplot, psym = 4, color = 'red', $
	symsize = 0.75
cgPlot, SceE_ogi[0, *], normData(SceE_ogi[3, *]), /overplot, linestyle = 0, color = 'red';add line

;Sce F:
cgPlot, SceF_noaa[0, *], normData(SceF_noaa[3, *]), /overplot, psym = 5, color = 'Magenta', $
	symsize = 0.75
cgPlot, SceF_noaa[0, *], normData(SceF_noaa[3, *]), /overplot, linestyle = 0, color = 'Magenta' ;add line
cgPlot, SceF_uci[0, *], normData(SceF_uci[3, *]), /overplot, psym = 2, color = 'Magenta', $
	symsize = 0.75
cgPlot, SceF_uci[0, *], normData(SceF_uci[3, *]), /overplot, linestyle = 0, color = 'Magenta';add line
cgPlot, SceF_ogi[0, *], normData(SceF_ogi[3, *]), /overplot, psym = 4, color = 'Magenta', $
	symsize = 0.75
cgPlot, SceF_ogi[0, *], normData(SceF_ogi[3, *]), /overplot, linestyle = 0, color = 'Magenta';add line


cgLegend, title = ['Sce A', 'Sce B', 'Sce C', 'Sce D', 'Sce E', 'Sce F'], $
	psym = [3, 3, 3, 3, 3, 3], linestyle = [0, 0, 0, 0, 0, 0], $
	color = ['blue violet', 'blue', 'forest green', 'lime green', 'red', 'Magenta'], $
	location = [1988, 1.3], charsize = 0.82, /data, vspace = 1

; cgLegend, title = ['NOAA', 'UCI', 'OGI'], $
	; psym = [5, 2, 4], $
	; color = ['black', 'black', 'black'], $
	; location = [1988, 0.5], charsize = 0.8, $
	; length = 0.0, /data, vspace = 1
close_device

spawn, 'gv temp.eps'

multiplot, /default

;otherwise, no normalization
endif else begin

;plotting procedure
;set up plot
open_device, /ps, /color, file='temp.eps', margin=0.05, xsize = 10.0, ysize = 7.5
!x.thick=1
!y.thick=1
!p.font=1
!p.thick =0.5

multiplot, /default    ; resets multiplot settings
multiplot, [1,3], ygap=0.002, xgap=0;  sets up multiplot 

;plot northern hemisphere
;yrange for normalized -200, 300
cgPlot, annual_nor_noaa[0, *], annual_nor_noaa[1, *], xrange = [1982, 2016], xticklen = 1, xgridstyle = 1, $
	xticks = 17, /nodata, yrange = [500,1500], ytitle = 'Mixing ratio(pptv)', $
	title = 'Time series of global Ethane NOAA, OGI sampled in Mar, Jun, Sep, Dec compared with full data'
cgPlot, annual_nor_noaa[0, *], (annual_nor_noaa[1, *]), /overplot, psym = 5, color = 'black', $
	err_yhigh = annual_nor_noaa[2, *], err_ylow = annual_nor_noaa[2, *]
cgPlot, annual_nor_uci[0, *], (annual_nor_uci[1, *]), /overplot, psym = 2, color = 'black', $
	err_yhigh = annual_nor_uci[2, *], err_ylow = annual_nor_uci[2, *]
cgPlot, annual_nor_ogi[0, *], (annual_nor_ogi[1, *]), /overplot, psym = 4, color = 'black', $
	err_yhigh = annual_nor_ogi[2, *], err_ylow = annual_nor_ogi[2, *]
	
;plot the sim data
;Sce A:
cgPlot, SceA_noaa[0, *], (SceA_noaa[1, *]), /overplot, psym = 5, color = 'blue violet', $
	symsize = 0.75
cgPlot, SceA_noaa[0, *], (SceA_noaa[1, *]), /overplot, linestyle = 0, color = 'blue violet' ;add line
cgPlot, SceA_uci[0, *], (SceA_uci[1, *]), /overplot, psym = 2, color = 'blue violet', $
	symsize = 0.75
cgPlot, SceA_uci[0, *], (SceA_uci[1, *]), /overplot, linestyle = 0, color = 'blue violet';add line
cgPlot, SceA_ogi[0, *], (SceA_ogi[1, *]), /overplot, psym = 4, color = 'blue violet', $
	symsize = 0.75
cgPlot, SceA_ogi[0, *], (SceA_ogi[1, *]), /overplot, linestyle = 0, color = 'blue violet';add line

;Sce B:
cgPlot, SceB_noaa[0, *], (SceB_noaa[1, *]), /overplot, psym = 5, color = 'blue', $
	symsize = 0.75
cgPlot, SceB_noaa[0, *], (SceB_noaa[1, *]), /overplot, linestyle = 0, color = 'blue' ;add line
cgPlot, SceB_uci[0, *], (SceB_uci[1, *]), /overplot, psym = 2, color = 'blue', $
	symsize = 0.75
cgPlot, SceB_uci[0, *], (SceB_uci[1, *]), /overplot, linestyle = 0, color = 'blue';add line
cgPlot, SceB_ogi[0, *], (SceB_ogi[1, *]), /overplot, psym = 4, color = 'blue', $
	symsize = 0.75
cgPlot, SceB_ogi[0, *], (SceB_ogi[1, *]), /overplot, linestyle = 0, color = 'blue';add line

;Sce C:
cgPlot, SceC_noaa[0, *], (SceC_noaa[1, *]), /overplot, psym = 5, color = 'forest green', $
	symsize = 0.75
cgPlot, SceC_noaa[0, *], (SceC_noaa[1, *]), /overplot, linestyle = 0, color = 'forest green' ;add line
cgPlot, SceC_uci[0, *], (SceC_uci[1, *]), /overplot, psym = 2, color = 'forest green', $
	symsize = 0.75
cgPlot, SceC_uci[0, *], (SceC_uci[1, *]), /overplot, linestyle = 0, color = 'forest green';add line
cgPlot, SceC_ogi[0, *], (SceC_ogi[1, *]), /overplot, psym = 4, color = 'forest green', $
	symsize = 0.75
cgPlot, SceC_ogi[0, *], (SceC_ogi[1, *]), /overplot, linestyle = 0, color = 'forest green';add line

;Sce D:
cgPlot, SceD_noaa[0, *], (SceD_noaa[1, *]), /overplot, psym = 5, color = 'lime green', $
	symsize = 0.75
cgPlot, SceD_noaa[0, *], (SceD_noaa[1, *]), /overplot, linestyle = 0, color = 'lime green' ;add line
cgPlot, SceD_uci[0, *], (SceD_uci[1, *]), /overplot, psym = 2, color = 'lime green', $
	symsize = 0.75
cgPlot, SceD_uci[0, *], (SceD_uci[1, *]), /overplot, linestyle = 0, color = 'lime green';add line
cgPlot, SceD_ogi[0, *], (SceD_ogi[1, *]), /overplot, psym = 4, color = 'lime green', $
	symsize = 0.75
cgPlot, SceD_ogi[0, *], (SceD_ogi[1, *]), /overplot, linestyle = 0, color = 'lime green';add line

;Sce E:
cgPlot, SceE_noaa[0, *], (SceE_noaa[1, *]), /overplot, psym = 5, color = 'red', $
	symsize = 0.75
cgPlot, SceE_noaa[0, *], (SceE_noaa[1, *]), /overplot, linestyle = 0, color = 'red' ;add line
cgPlot, SceE_uci[0, *], (SceE_uci[1, *]), /overplot, psym = 2, color = 'red', $
	symsize = 0.75
cgPlot, SceE_uci[0, *], (SceE_uci[1, *]), /overplot, linestyle = 0, color = 'red';add line
cgPlot, SceE_ogi[0, *], (SceE_ogi[1, *]), /overplot, psym = 4, color = 'red', $
	symsize = 0.75
cgPlot, SceE_ogi[0, *], (SceE_ogi[1, *]), /overplot, linestyle = 0, color = 'red';add line

;Sce F:
cgPlot, SceF_noaa[0, *], (SceF_noaa[1, *]), /overplot, psym = 5, color = 'Magenta', $
	symsize = 0.75
cgPlot, SceF_noaa[0, *], (SceF_noaa[1, *]), /overplot, linestyle = 0, color = 'Magenta' ;add line
cgPlot, SceF_uci[0, *], (SceF_uci[1, *]), /overplot, psym = 2, color = 'Magenta', $
	symsize = 0.75
cgPlot, SceF_uci[0, *], (SceF_uci[1, *]), /overplot, linestyle = 0, color = 'Magenta';add line
cgPlot, SceF_ogi[0, *], (SceF_ogi[1, *]), /overplot, psym = 4, color = 'Magenta', $
	symsize = 0.75
cgPlot, SceF_ogi[0, *], (SceF_ogi[1, *]), /overplot, linestyle = 0, color = 'Magenta';add line
	
multiplot, /doyaxis, /doxaxis

;plot southern hemisphere
;yrange for normalized -200, 300
cgPlot, annual_sou_noaa[0, *], annual_sou_noaa[1, *], /nodata, xrange = [1982, 2016], $
	xticklen = 1, xgridstyle = 1, xticks = 17, XTickformat='(A1)', yrange = [200,800], $
	ytitle = 'Mixing ratio(pptv)'
cgPlot, annual_sou_noaa[0, *], (annual_sou_noaa[1, *]), /overplot, psym = 5, color = 'black', $
	err_yhigh = annual_sou_noaa[2, *], err_ylow = annual_sou_noaa[2, *]
cgPlot, annual_sou_uci[0, *], (annual_sou_uci[1, *]), /overplot, psym = 2, color = 'black', $
	err_yhigh = annual_sou_uci[2, *], err_ylow = annual_sou_uci[2, *]
cgPlot, annual_sou_ogi[0, *], (annual_sou_ogi[1, *]), /overplot, psym = 4, color = 'black', $
	err_yhigh = annual_sou_ogi[2, *], err_ylow = annual_sou_ogi[2, *]
	
;plot the sim data
;Sce A:
cgPlot, SceA_noaa[0, *], (SceA_noaa[2, *]), /overplot, psym = 5, color = 'blue violet', $
	symsize = 0.75
cgPlot, SceA_noaa[0, *], (SceA_noaa[2, *]), /overplot, linestyle = 0, color = 'blue violet' ;add line
cgPlot, SceA_uci[0, *], (SceA_uci[2, *]), /overplot, psym = 2, color = 'blue violet', $
	symsize = 0.75
cgPlot, SceA_uci[0, *], (SceA_uci[2, *]), /overplot, linestyle = 0, color = 'blue violet';add line
cgPlot, SceA_ogi[0, *], (SceA_ogi[2, *]), /overplot, psym = 4, color = 'blue violet', $
	symsize = 0.75
cgPlot, SceA_ogi[0, *], (SceA_ogi[2, *]), /overplot, linestyle = 0, color = 'blue violet';add line

;Sce B:
cgPlot, SceB_noaa[0, *], (SceB_noaa[2, *]), /overplot, psym = 5, color = 'blue', $
	symsize = 0.75
cgPlot, SceB_noaa[0, *], (SceB_noaa[2, *]), /overplot, linestyle = 0, color = 'blue' ;add line
cgPlot, SceB_uci[0, *], (SceB_uci[2, *]), /overplot, psym = 2, color = 'blue', $
	symsize = 0.75
cgPlot, SceB_uci[0, *], (SceB_uci[2, *]), /overplot, linestyle = 0, color = 'blue';add line
cgPlot, SceB_ogi[0, *], (SceB_ogi[2, *]), /overplot, psym = 4, color = 'blue', $
	symsize = 0.75
cgPlot, SceB_ogi[0, *], (SceB_ogi[2, *]), /overplot, linestyle = 0, color = 'blue';add line

;Sce C:
cgPlot, SceC_noaa[0, *], (SceC_noaa[2, *]), /overplot, psym = 5, color = 'forest green', $
	symsize = 0.75
cgPlot, SceC_noaa[0, *], (SceC_noaa[2, *]), /overplot, linestyle = 0, color = 'forest green' ;add line
cgPlot, SceC_uci[0, *], (SceC_uci[2, *]), /overplot, psym = 2, color = 'forest green', $
	symsize = 0.75
cgPlot, SceC_uci[0, *], (SceC_uci[2, *]), /overplot, linestyle = 0, color = 'forest green';add line
cgPlot, SceC_ogi[0, *], (SceC_ogi[2, *]), /overplot, psym = 4, color = 'forest green', $
	symsize = 0.75
cgPlot, SceC_ogi[0, *], (SceC_ogi[2, *]), /overplot, linestyle = 0, color = 'forest green';add line

;Sce D:
cgPlot, SceD_noaa[0, *], (SceD_noaa[2, *]), /overplot, psym = 5, color = 'lime green', $
	symsize = 0.75
cgPlot, SceD_noaa[0, *], (SceD_noaa[2, *]), /overplot, linestyle = 0, color = 'lime green' ;add line
cgPlot, SceD_uci[0, *], (SceD_uci[2, *]), /overplot, psym = 2, color = 'lime green', $
	symsize = 0.75
cgPlot, SceD_uci[0, *], (SceD_uci[2, *]), /overplot, linestyle = 0, color = 'lime green';add line
cgPlot, SceD_ogi[0, *], (SceD_ogi[2, *]), /overplot, psym = 4, color = 'lime green', $
	symsize = 0.75
cgPlot, SceD_ogi[0, *], (SceD_ogi[2, *]), /overplot, linestyle = 0, color = 'lime green';add line

;Sce E:
cgPlot, SceE_noaa[0, *], (SceE_noaa[2, *]), /overplot, psym = 5, color = 'red', $
	symsize = 0.75
cgPlot, SceE_noaa[0, *], (SceE_noaa[2, *]), /overplot, linestyle = 0, color = 'red' ;add line
cgPlot, SceE_uci[0, *], (SceE_uci[2, *]), /overplot, psym = 2, color = 'red', $
	symsize = 0.75
cgPlot, SceE_uci[0, *], (SceE_uci[2, *]), /overplot, linestyle = 0, color = 'red';add line
cgPlot, SceE_ogi[0, *], (SceE_ogi[2, *]), /overplot, psym = 4, color = 'red', $
	symsize = 0.75
cgPlot, SceE_ogi[0, *], (SceE_ogi[2, *]), /overplot, linestyle = 0, color = 'red';add line

;Sce F:
cgPlot, SceF_noaa[0, *], (SceF_noaa[2, *]), /overplot, psym = 5, color = 'Magenta', $
	symsize = 0.75
cgPlot, SceF_noaa[0, *], (SceF_noaa[2, *]), /overplot, linestyle = 0, color = 'Magenta' ;add line
cgPlot, SceF_uci[0, *], (SceF_uci[2, *]), /overplot, psym = 2, color = 'Magenta', $
	symsize = 0.75
cgPlot, SceF_uci[0, *], (SceF_uci[2, *]), /overplot, linestyle = 0, color = 'Magenta';add line
cgPlot, SceF_ogi[0, *], (SceF_ogi[2, *]), /overplot, psym = 4, color = 'Magenta', $
	symsize = 0.75
cgPlot, SceF_ogi[0, *], (SceF_ogi[2, *]), /overplot, linestyle = 0, color = 'Magenta';add line

multiplot, /doyaxis, /doxaxis

;plot IHR
;yrange for non-normalized -1.5, 1
cgPlot, annual_nor_noaa[1, *], annual_ihr_noaa, /nodata, xtitle = 'Years', ytitle = 'IHR', $
	xrange = [1982, 2016], yrange = [1, 6], $
	xticklen = 1, xgridstyle = 1, xticks = 17
cgPlot, annual_sou_noaa[0, *], (annual_ihr_noaa), /overplot, psym = 5, color = 'black', $
	err_yhigh = annual_ihr_noaa_err, err_ylow = annual_ihr_noaa_err
cgPlot, annual_sou_uci[0, *], (annual_ihr_uci), /overplot, psym = 2, color = 'black', $
	err_yhigh = annual_ihr_uci_err, err_ylow =annual_ihr_uci_err
cgPlot, annual_sou_ogi[0, *], (annual_ihr_ogi), /overplot, psym = 4, color = 'black', $
	err_yhigh = annual_ihr_ogi_err, err_ylow = annual_ihr_ogi_err

;plot the sim data
;Sce A:
cgPlot, SceA_noaa[0, *], (SceA_noaa[3, *]), /overplot, psym = 5, color = 'blue violet', $
	symsize = 0.75
cgPlot, SceA_noaa[0, *], (SceA_noaa[3, *]), /overplot, linestyle = 0, color = 'blue violet' ;add line
cgPlot, SceA_uci[0, *], (SceA_uci[3, *]), /overplot, psym = 2, color = 'blue violet', $
	symsize = 0.75
cgPlot, SceA_uci[0, *], (SceA_uci[3, *]), /overplot, linestyle = 0, color = 'blue violet';add line
cgPlot, SceA_ogi[0, *], (SceA_ogi[3, *]), /overplot, psym = 4, color = 'blue violet', $
	symsize = 0.75
cgPlot, SceA_ogi[0, *], (SceA_ogi[3, *]), /overplot, linestyle = 0, color = 'blue violet';add line

;Sce B:
cgPlot, SceB_noaa[0, *], (SceB_noaa[3, *]), /overplot, psym = 5, color = 'blue', $
	symsize = 0.75
cgPlot, SceB_noaa[0, *], (SceB_noaa[3, *]), /overplot, linestyle = 0, color = 'blue' ;add line
cgPlot, SceB_uci[0, *], (SceB_uci[3, *]), /overplot, psym = 2, color = 'blue', $
	symsize = 0.75
cgPlot, SceB_uci[0, *], (SceB_uci[3, *]), /overplot, linestyle = 0, color = 'blue';add line
cgPlot, SceB_ogi[0, *], (SceB_ogi[3, *]), /overplot, psym = 4, color = 'blue', $
	symsize = 0.75
cgPlot, SceB_ogi[0, *], (SceB_ogi[3, *]), /overplot, linestyle = 0, color = 'blue';add line

;Sce C:
cgPlot, SceC_noaa[0, *], (SceC_noaa[3, *]), /overplot, psym = 5, color = 'forest green', $
	symsize = 0.75
cgPlot, SceC_noaa[0, *], (SceC_noaa[3, *]), /overplot, linestyle = 0, color = 'forest green' ;add line
cgPlot, SceC_uci[0, *], (SceC_uci[3, *]), /overplot, psym = 2, color = 'forest green', $
	symsize = 0.75
cgPlot, SceC_uci[0, *], (SceC_uci[3, *]), /overplot, linestyle = 0, color = 'forest green';add line
cgPlot, SceC_ogi[0, *], (SceC_ogi[3, *]), /overplot, psym = 4, color = 'forest green', $
	symsize = 0.75
cgPlot, SceC_ogi[0, *], (SceC_ogi[3, *]), /overplot, linestyle = 0, color = 'forest green';add line

;Sce D:
cgPlot, SceD_noaa[0, *], (SceD_noaa[3, *]), /overplot, psym = 5, color = 'lime green', $
	symsize = 0.75
cgPlot, SceD_noaa[0, *], (SceD_noaa[3, *]), /overplot, linestyle = 0, color = 'lime green' ;add line
cgPlot, SceD_uci[0, *], (SceD_uci[3, *]), /overplot, psym = 2, color = 'lime green', $
	symsize = 0.75
cgPlot, SceD_uci[0, *], (SceD_uci[3, *]), /overplot, linestyle = 0, color = 'lime green';add line
cgPlot, SceD_ogi[0, *], (SceD_ogi[3, *]), /overplot, psym = 4, color = 'lime green', $
	symsize = 0.75
cgPlot, SceD_ogi[0, *], (SceD_ogi[3, *]), /overplot, linestyle = 0, color = 'lime green';add line

;Sce E:
cgPlot, SceE_noaa[0, *], (SceE_noaa[3, *]), /overplot, psym = 5, color = 'red', $
	symsize = 0.75
cgPlot, SceE_noaa[0, *], (SceE_noaa[3, *]), /overplot, linestyle = 0, color = 'red' ;add line
cgPlot, SceE_uci[0, *], (SceE_uci[3, *]), /overplot, psym = 2, color = 'red', $
	symsize = 0.75
cgPlot, SceE_uci[0, *], (SceE_uci[3, *]), /overplot, linestyle = 0, color = 'red';add line
cgPlot, SceE_ogi[0, *], (SceE_ogi[3, *]), /overplot, psym = 4, color = 'red', $
	symsize = 0.75
cgPlot, SceE_ogi[0, *], (SceE_ogi[3, *]), /overplot, linestyle = 0, color = 'red';add line

;Sce F:
cgPlot, SceF_noaa[0, *], (SceF_noaa[3, *]), /overplot, psym = 5, color = 'Magenta', $
	symsize = 0.75
cgPlot, SceF_noaa[0, *], (SceF_noaa[3, *]), /overplot, linestyle = 0, color = 'Magenta' ;add line
cgPlot, SceF_uci[0, *], (SceF_uci[3, *]), /overplot, psym = 2, color = 'Magenta', $
	symsize = 0.75
cgPlot, SceF_uci[0, *], (SceF_uci[3, *]), /overplot, linestyle = 0, color = 'Magenta';add line
cgPlot, SceF_ogi[0, *], (SceF_ogi[3, *]), /overplot, psym = 4, color = 'Magenta', $
	symsize = 0.75
cgPlot, SceF_ogi[0, *], (SceF_ogi[3, *]), /overplot, linestyle = 0, color = 'Magenta';add line


cgLegend, title = ['Sce A', 'Sce B', 'Sce C', 'Sce D', 'Sce E', 'Sce F'], $
	psym = [3, 3, 3, 3, 3, 3], linestyle = [0, 0, 0, 0, 0, 0], $
	color = ['blue violet', 'blue', 'forest green', 'lime green', 'red', 'Magenta'], $
	location = [1988, 5.9], charsize = 0.82, /data, vspace = 1

; cgLegend, title = ['NOAA', 'UCI', 'OGI'], $
	; psym = [5, 2, 4], $
	; color = ['black', 'black', 'black'], $
	; location = [1988, 0.5], charsize = 0.8, $
	; length = 0.0, /data, vspace = 1
close_device

spawn, 'gv temp.eps'

multiplot, /default

endelse


end