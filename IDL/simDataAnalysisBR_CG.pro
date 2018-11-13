;this program analyzes the simulated data from GEOS-Chem using the same 
;algorithm as for ihrSamplingSensitivityBR_CP.pro



FUNCTION annualMean, inputArr, year

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
help, shortArr
out = fltarr(n_elements(year), n_elements(shortArr[0, *, 0]), n_elements(shortArr[0, 0, *]))
for i = 0, n_elements(year)-1 do begin
	;take every 4 indeces out to calculate the annual mean
	tempArr = shortArr[4 * i : 4 * i + 3, *, *] 
	out[i, *, *] = mean(tempArr, 1)/2 * 1000 ;adjust to pptv
endfor
return, out ;return value

end

FUNCTION getTimeSeriesSite, masterArr, lat, lon

;this function takes the input as the three dimensional simulated annual mean array and uses
;the input lon and lat to return the time series of the input site

;employ the ctm_index function to obtain the index of the input lon lat in the array
CTM_INDEX, CTM_TYPE('GEOS1', RESOLUTION= 2), i, j, CENTER= [lat, lon], /non_interactive
out = masterArr[*, i-1, j-1] ;get the annual mean time series of the input coordinate.

return, out

end


PRO simDataAnalysisBR_CG 

compile_opt idl2

;sim_title = 'PSU emissions scaled to Xiao et al over 1996-2003'
filename1 = "/home/excluded-from-backup/data/C2H6/trac_avg.PSUSF_1981_2015.bpch"

;read in the file
ctm_get_data, datainfo, 'IJ-AVG-$', tracer=1, filename=filename1, tau0=tau0

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
;simYear contains the years that the simulation has
simYear = sim_ymd.year[sort(sim_ymd.year)]
simYear = simYear[uniq(simYear)]
;sim1 is a 3D array, [year, longitudes, latitudes]
sim1 = annualMean(tempSimArr, simYear)
;dimenstion for sim1 is [34, 144, 91] 34 is 34 years, 144 and 91 are the index for lon and lat 
;since the simulation is a 2x2.5 degree grid.
;obtain index of the Barrow site, variable i1, j1 contain the index of Barrow
barrow = getTimeSeriesSite(sim1, 71.31, -156.6)
;obtain index of the Cape Grim site
capegrim = getTimeSeriesSite(sim1, -40.66, 144.66)





end
