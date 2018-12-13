;this program studies the sensitivity of the simulated data when sampled with different methods
;refer to issue #10 of the Github repo

PRO simDataSensitivity

compile_opt idl2


;specify the latitude bins
bin_bound = [-50, -30, 0, 30, 50, 75]

;read in sim data first. 
;use default emission scenario
ems_dir = "/home/excluded-from-backup/data/C2H6/trac_avg.spinup_GFED4_MAVG_1981_2015.198101010000"
print, ems_dir

;extract the 3D array that contains the global mixing ratio
;read in the file
ctm_get_data, datainfo, 'IJ-AVG-$', tracer=1, filename=ems_dir, tau0=tau0
;Get MODELINFO and GRIDINFO structures, xmid, ymid hold lon/lat centers
getmodelandgridinfo, datainfo[0], modelinfo, gridinfo
;nt is the number of data points in the simulation.
nt = n_elements(datainfo)
simArr = fltarr(nt, 144, 91)
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
			simArr[i, k, j] = data[k, j]/2 * 1000
		endfor
	endfor
endfor
;clean up memory
CTM_CLEANUP 

;simYear contains the years that the simulation has
simYear = sim_ymd.year[sort(sim_ymd.year)]
simYear = simYear[uniq(simYear)]

;calculate the time series of the simulated data using complete global data
;first, need to calculate the annual average of the sim data
;there are 34 years in the simulated data 
;create another 3D array that will contain the annual mean of the simulated data
annualSimGlobal = fltarr(34, 144, 91)
for i = 0, nt - 1 do begin
	if (i eq 407) then break
	annualSimGlobal[i/12, *, *] = mean(simArr[i : i + 11, *, *], 1)
	i = i + 10 ;jump 10 indices, not 11 because the for loop will add 1
endfor

	
end
