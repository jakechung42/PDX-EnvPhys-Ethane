PRO LON_DIS

;08/17/2017	-	Jake Chung
;	This program plot the longitudinal distribution at the Samoa, South Pole, 
;Summit Greenland and Barrow AK latitudes using the default emission data.
;08/18/2017 - 	Jake Chung
;	Adding Cape Grim site

compile_opt idl2

;default base emission scenario<<<<<<<<<<<
filename = "/home/excluded-from-backup/data/C2H6/trac_avg.spinup_GFED4_MAVG_1981_2015.198101010000"

;read in the file
ctm_get_data, datainfo, 'IJ-AVG-$', tracer=1, filename=filename, tau0=tau0

;Get MODELINFO and GRIDINFO structures, xmid, ymid hold lon/lat centers
getmodelandgridinfo, datainfo[0], modelinfo, gridinfo

;nt is the number of data points in the simulation.
nt = n_elements(datainfo)
globalavg= fltarr(nt)

;sim_yymmdd stores the entire time element of the input simulation
sim_yymmdd= tau2yymmdd(tau0, /GEOS1)

;this section call out the data blocks of the simulated data and store it in simArr
;simArr is a 3-D array with the following attributes:
;simArr[ *, 0, 0 ]: the time dimension
;simArr[ 0, *, 0 ]: longitudes
;simArr[ 0, 0, * ]: latitudes
simArr = fltarr( nt, 144, 91 )

for i = 0, nt - 1 do begin
	data= CTM_EXTRACT( *(datainfo[ i ].data), modelinfo= modelinfo, $
	gridinfo= gridinfo, lat= [-90, 90], lon= [-180, 180], alrange = [1, 2], $
	average= 4 )
	for k = 0, n_elements( data[ * , 0 ] ) - 1 do begin
		for j = 0, n_elements( data[ 0 , * ] ) - 1 do begin
			simArr[ i , k , j ] = data[ k , j ]
		endfor
	endfor
endfor

;extract the data at the Samoa latitude
avg_samoa = fltarr(144)
avg_samoa = ( mean( simArr[ * , * , 38], 1) / 2 )* 1000

;extract the data at the South Pole Latitude
avg_southpole = fltarr(144)
avg_southpole = ( mean( simArr[ * , * , 0], 1) / 2 )* 1000

;extract the data at the Barrow latitude
avg_barrow = fltarr(144)
avg_barrow = ( mean( simArr[ * , * , 81], 1) / 2 )* 1000

;extract the data at the WAIS-D latitude
avg_wais = fltarr(144)
avg_wais = ( mean( simArr[ * , * , 5], 1) / 2 )* 1000

;extract the data at the Cape Grim latitude
avg_cpgrim = fltarr(144)
avg_cpgrim = ( mean( simArr[ * , * , 25], 1) / 2 )* 1000


longitudes = fltarr(144)
longitudes = findgen(144)*2.5 - 178.75

open_device, /ps, /color, file='temp.eps', /swaplandscape, margin=0.05, xsize = 10.0, ysize = 7.5
!x.thick=5
!y.thick=5
!p.font=0
!p.thick =5

multiplot, /default    ; resets multiplot settings
multiplot, [1,2], ygap=0.002, xgap=0;  sets up multiplot 

;Northern Hemisphere
cgplot, longitudes, avg_samoa, /nodata, xrange = [ -180 , 180 ], yrange = [ 1200 , 1800 ], charsize = 1.3,	$
	ytitle = 'NH mixing ratio (pptv)'
cgplot, longitudes, avg_barrow, thick = 3, /overplot, color = 'steelblue'

AL_Legend, ['Barrow/Summit', 'Samoa', 'South Pole', 'WAIS-D', 'Cape Grim'], $
	color = ['steelblue', 'red', 'forest green', 'brown', 'black'], $
	linestyle = [0,0,0,0,0], position = [50, 1500], charsize = 1, box = 0
	
multiplot, /doyaxis, /doxaxis
;Southern Hemisphere
cgplot, longitudes, avg_samoa, /nodata, xrange = [ -180 , 180 ], yrange = [ 200 , 700 ], charsize = 1.3,	$
	ytitle = 'SH mixing ratio (pptv)', xtitle = 'Longitudes'
	
cgplot, longitudes, avg_samoa, thick = 3, /overplot, color = 'red'
cgplot, longitudes, avg_southpole, thick = 3, /overplot, color = 'forest green'
cgplot, longitudes, avg_wais, thick = 3, /overplot, color = 'brown'
cgplot, longitudes, avg_cpgrim, thick = 3, /overplot, color = 'black'

close_device
spawn, 'gv temp.eps'
 
end
