;this program belongs to the global simulated data and observed data comparison cluster. 
;To reduce the complexity of the cluster this program is kept seperated from the analysis 
;itself.
;This program builds the simulated data set that matches exactly to the observed data 
;temporally and spatially, so a consistent comparison with the observed data can be
;conducted.

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

PRO buildSimData

compile_opt idl2 

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
;remove months that are not Mar, Jun, Sep, or Dec
uci = rmMonth(uci)
noaa = rmMonth(noaa)
ogi = rmMonth(ogi)


;read in the full sim data 
;sim_title = 'PSU emissions scaled to Xiao et al over 1996-2003'
filename1 = "/home/excluded-from-backup/data/C2H6/trac_avg.PSUSF_1981_2015.bpch"
;read in the file
ctm_get_data, datainfo, 'IJ-AVG-$', tracer=1, filename=filename1, tau0=tau0
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
			tempSimArr[i, k, j] = data[k, j]
		endfor
	endfor
endfor
;clean up memory
CTM_CLEANUP 
;simYear contains the years that the simulation has
simYear = sim_ymd.year[sort(sim_ymd.year)]
simYear = simYear[uniq(simYear)]

simRatio = fltarr(n_elements(uci.ratio))
for i = 0, n_elements(uci.ratio)-1 do begin 
;loop through the entire data set to pull the sim data based on the information
;of the observed data
	lat = uci.lat[i]
	lon = uci.lon[i]
	year = uci.year[i]
	month = uci.month[i]
	;use the lat and lon to get the indices for spacial data
	CTM_INDEX, CTM_TYPE('GEOS1', RESOLUTION= 2), IDX_lon, IDX_lat, CENTER= [lat, lon], /non_interactive
	print, lat, lon, IDX_lat, IDX_lon
	;calculate the index for the time dimension
	IDX_time = (year - 1981)*12 + month - 1
	simRatio[i] = tempSimArr[IDX_time, IDX_lon-1, IDX_lat-1]
endfor 


	
	
end
