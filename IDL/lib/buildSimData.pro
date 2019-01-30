;this function belongs to the global simulated data and observed data comparison cluster. 
;To reduce the complexity of the cluster this program is kept seperated from the analysis 
;itself.
;This function builds the simulated data set that matches exactly to the observed data 
;temporally and spatially, so a consistent comparison with the observed data can be
;conducted.
;The input for this function is the observed data structure that needed to be modeled
;with sim data, and the path directory to the .bpch output from GEOS-Chem

FUNCTION rmMonth, struct

;this function removes the months that are not 3, 6, 9, 12 from the array
idx = where(struct.month eq 3 or $
			struct.month eq 6 or $
			struct.month eq 9 or $
			struct.month eq 12)
;reconstruct the structure
out = {ratio: struct.ratio[idx], $
		lat: struct.lat[idx], $
		lon: struct.lon[idx], $
		year: struct.year[idx], $
		month: struct.month[idx]}
;return a new array with only month 3, 6, 9, 12
return, out

end
;======================================================================
;main function starts
FUNCTION buildSimData, obsData, simPath

;Some networks have 2015 data, remove it from the data
;array so this function would not throw an error
lt2015 = where(obsData.year lt 2015)
obsData = {ratio: obsData.ratio[lt2015], $
		lat: obsData.lat[lt2015], $
		lon: obsData.lon[lt2015], $
		year: obsData.year[lt2015], $
		month: obsData.month[lt2015]}
;remove the excess months and only keep Mar, June, Sep, Dec
obsData = rmMonth(obsData)
;read in the file from the input directory
ctm_get_data, datainfo, 'IJ-AVG-$', tracer=1, filename=simPath, tau0=tau0
;Get MODELINFO and GRIDINFO structures, xmid, ymid hold lon/lat centers
getmodelandgridinfo, datainfo[0], modelinfo, gridinfo
;nt is the number of months in the simulation.
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
			tempSimArr[i, k, j] = data[k, j]/2 * 1000
		endfor
	endfor
endfor
;clean up memory
CTM_CLEANUP 
;simYear contains the years that the simulation has
simYear = sim_ymd.year[sort(sim_ymd.year)]
simYear = simYear[uniq(simYear)]

;begin retrieving sim data
simRatio = fltarr(n_elements(obsData.ratio)) 
pulledTime = []
pulledLat = []
pulledLon = []
for i = 0, n_elements(obsData.ratio)-1 do begin 
;loop through the entire data set to pull the sim data based on the information
;of the observed data
	lat = obsData.lat[i]
	lon = obsData.lon[i]
	year = obsData.year[i]
	month = obsData.month[i]
	;use the lat and lon to get the indices for spatial data
	CTM_INDEX, CTM_TYPE('GEOS1', RESOLUTION= 2), IDX_lon, IDX_lat, CENTER= [lat, lon], /non_interactive
	;calculate the index for the time dimension
	IDX_time = (year - 1981)*12 + month - 1
	;to prevent duplicate months, build an array that stores the indecies and time 
	;that data have been retrieved, if try to pull data out from an index-set, set
	;pulled-value to be NaN
	calledBefore = where(IDX_time eq pulledTime and IDX_lat eq pulledLat and IDX_lon eq pulledLon, count)
	if count eq 0 then begin ;if not called before
		;add to the list of called values
		pulledTime = [pulledTime, IDX_time]
		pulledLat = [pulledLat, IDX_lat]
		pulledLon = [pulledLon, IDX_lon]
		;pull sim mixing ratio
		simRatio[i] = tempSimArr[IDX_time, IDX_lon-1, IDX_lat-1]
	endif else begin
		;if called before set NaN
		simRatio[i] = !VALUES.F_NAN
	endelse
endfor 
;reconstruct the sim data structure and remove the NaN values
notNaN = finite(simRatio)
notNaN_idx = where(notNaN eq 1)
;remove NaN values from array
simRatio = simRatio[notNaN_idx]
;reconstruct the struct
simDataStruct = {ratio: simRatio, $
				lat: obsData.lat[notNaN_idx], $
				lon: obsData.lon[notNaN_idx], $
				year: obsData.year[notNaN_idx], $
				month: obsData.month[notNaN_idx]}

;return the simulated data structure that modeled from the input observed data structure
return, simDataStruct	
	
end
