FUNCTION annualMean, inputArr, year

;the output of GEOS-Chem is a 3 dimensional array with time, lon, lat
;corresponding to x, y, z. The time dimension for our simulations 
;is monthly average. We need annual average. 
;This program takes the 3D-array and calculate the annual average

;make sure that the data has 12 months for every year
modulo = n_elements(inputArr(*, 0, 0)) mod n_elements(year)
if modulo ne 0 then begin
	print, 'Do not have enough months, stop the program!!!'
	stop		;if there's not enough months, stop the program
endif

;start calculating the annual average
;make an array to store the annual average 
out = fltarr(n_elements(year), n_elements(inputArr[0, *, 0]), n_elements(inputArr[0, 0, *]))

for i = 0, n_elements(year)-1 do begin
	;take every 12 indeces out to calculate the annual mean
	tempArr = inputArr[i * 2 : 12 * i + 11, *, *] 
	out[i, *, *] = mean(tempArr, 1)
endfor

return, out	

end





PRO simDataAnalysisBR_CG 

compile_opt idl2

;this program analyzes the simulated data from GEOS-Chem using the same 
;algorithm as for ihrSamplingSensitivityBR_CP.pro

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

print, simYear 

sim1 = annualMean(tempSimArr, simYear)

help, sim1
end
