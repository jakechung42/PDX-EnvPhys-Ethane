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
globalavg= fltarr(nt)
tempSimArr = fltarr( nt, 144, 91 )

;sim_yymmdd stores the entire time element of the input simulation
sim_yymmdd= tau2yymmdd(tau0, /GEOS1)

;this section call out the data blocks of the simulated data and store it in simArr
;simArr is a 3-D array with the following attributes:
;simArr[ *, 0, 0 ]: the time dimension
;simArr[ 0, *, 0 ]: longitudes
;simArr[ 0, 0, * ]: latitudes

for i = 0, nt - 1 do begin
	data= CTM_EXTRACT( *(datainfo[ i ].data), modelinfo= modelinfo, $
	gridinfo= gridinfo, lat= [-90, 90], lon= [-180, 180], alrange = [1, 2], $
	average= 4 )
	for k = 0, n_elements( data[ * , 0 ] ) - 1 do begin
		for j = 0, n_elements( data[ 0 , * ] ) - 1 do begin
			tempSimArr[ i , k , j ] = data[ k , j ]
		endfor
	endfor
endfor

help, tempSimArr


end
