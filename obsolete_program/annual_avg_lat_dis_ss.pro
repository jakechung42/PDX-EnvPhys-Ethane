PRO ANNUAL_AVG_LAT_DIS_SS

;This program is very similar to annual_avg_lat_dis_noaa but this program deal with the Simpson et al. data

compile_opt idl2		;Set compiler option
;=====Start reading in Simpson et al. data files=====
;get the files
infile1= "/home/jakechung/IDL/myIDL/data_set.dat"
infile2= "/home/jakechung/IDL/myIDL/location_set.dat"

;get number of lines in the data file
n1_ss = file_lines(infile1)
n2_ss = file_lines(infile2)

print, 'The input has ', n1_ss, '   data lines with', n2_ss, '  locations'

;create an empty array to store the data that was read in lon, lat, avg, date and locations
unix_date_ss = fltarr(n1_ss)
lat1_ss = fltarr(n1_ss)
lon1_ss = fltarr(n1_ss)
avg_ss = fltarr(n1_ss)
location_ss = strarr(n2_ss)

;FLOAT input variable
location0= ' '
date0= 0.0
lat0= 0.0
lon0= 0.0
avg0= 0.0

;open input files
openr, iunit1, infile1, /get_lun
openr, iunit2, infile2, /get_lun

;read float data of the data file
for i = 0, n1_ss-1 do begin
	;read one line from the file
	readf, iunit1, date0, lat0, lon0, avg0
	;store the values
	unix_date_ss[ i ] = date0
	lat1_ss[ i ] = lat0
	lon1_ss[ i ] = lon0
	avg_ss[ i ] = avg0
endfor

;read float data of the location file
for i = 0, n2_ss-1 do begin

	readf, iunit2, location0
	location_ss[ i ] = location0
	
endfor

;close input file
free_lun, iunit1, iunit2

;Because time in the Simpson et al. data set denoted time as Unix time, these lines convert Unix time to Geos_chem time and stored in tau_time array
tau_time_ss = (unix_date_ss/3600 - 131496)
tau_time_ss = tau2yymmdd(tau_time_ss, /GEOS1)

;Convert unix time to decimal year
time_ss = (unix_date_ss/3600 - 131496)/8760 + 1985

;=====Finish reading in Simpson et al. data files=====
;============================================================================================================

;Analyzing the Simpson et al. data, the data is available from 1996 to 2009. Will be picking out the spcific latitudes in the lat_var_GT30 array.


lat_var_GT30 = [ 71.3000, 64.4900, 60.7500, 57.8300, 57.8000, 42.8300, 30.4000, 29.9400, 23.4400, 22.8700, 20.2600, 15.2000, $
	13.3600, 9.52000, 7.05000, 6.92000, 1.42000, -8.52000, -14.2300, -21.2300, -29.0200, -34.9000, -36.8200, -42.7200, -43.8800 ]
	
bin01 = fltarr( 3 , n_elements( lat_var_GT30 ) )
bin02 = fltarr( 3 , n_elements( lat_var_GT30 ) )
bin03 = fltarr( 3 , n_elements( lat_var_GT30 ) )

for i = 0, n_elements( lat_var_GT30 ) - 1 do begin
	a = where( lat1_ss[ * ] EQ lat_var_GT30[ i ] )
	temp_lat = lat1_ss[ a ]
	temp_lon = lon1_ss[ a ]
	temp_avg = avg_ss[ a ] 
	temp_bin01 = where( tau_time_ss.year[ a ] EQ 1996 OR $
		tau_time_ss.year[ a ] EQ 1997 OR $
		tau_time_ss.year[ a ] EQ 1998 OR $
		tau_time_ss.year[ a ] EQ 1999 OR $
		tau_time_ss.year[ a ] EQ 2000 , count )
	if ( count GT 3 ) then begin 
		bin01[ 2 , i ] = mean( temp_avg[ temp_bin01 ], /NAN )	
		bin01[ 0 , i ] = temp_lat[ 0 ]
		bin01[ 1 , i ] = temp_lon[ 0 ]
		endif else begin	
			bin01[ 2 , i ] = !Values.F_NAN
			bin01[ 0 , i ] = temp_lat[ 0 ]
			bin01[ 1 , i ] = temp_lon[ 0 ]
			endelse

	temp_bin02 = where( tau_time_ss.year[ a ] EQ 2001 OR $
		tau_time_ss.year[ a ] EQ 2002 OR $
		tau_time_ss.year[ a ] EQ 2003 OR $
		tau_time_ss.year[ a ] EQ 2004 OR $
		tau_time_ss.year[ a ] EQ 2005 , count )
	if ( count GT 3 ) then begin 
		bin02[ 2 , i ] = mean( temp_avg[ temp_bin02 ], /NAN )	
		bin02[ 0 , i ] = temp_lat[ 0 ]
		bin02[ 1 , i ] = temp_lon[ 0 ]
		endif else begin	
			bin02[ 2 , i ] = !Values.F_NAN
			bin02[ 0 , i ] = temp_lat[ 0 ]
			bin02[ 1 , i ] = temp_lon[ 0 ]
			endelse
			
	temp_bin03 = where( tau_time_ss.year[ a ] EQ 2006 OR $
		tau_time_ss.year[ a ] EQ 2007 OR $
		tau_time_ss.year[ a ] EQ 2008 OR $
		tau_time_ss.year[ a ] EQ 2009 , count )
	if ( count GT 3 ) then begin 
		bin03[ 2 , i ] = mean( temp_avg[ temp_bin03 ], /NAN )	
		bin03[ 0 , i ] = temp_lat[ 0 ]
		bin03[ 1 , i ] = temp_lon[ 0 ]
		endif else begin	
			bin03[ 2 , i ] = !Values.F_NAN
			bin03[ 0 , i ] = temp_lat[ 0 ]
			bin03[ 1 , i ] = temp_lon[ 0 ]
			endelse
			
endfor 
			
;finish calculating the averages for the Simpson et al. data

;bringing in the deseasonal data to plot uncertainty for the data
;Note: for the Simpson et al. data, data is only available in March, June, September and December
;The process is the same as sub part three-one, just different files names and variables

;Create link to the data file
temp_deseason = "/home/jakechung/IDL/myIDL/temp_folder/temp_residual_Simpson.dat"
openr, deseason, temp_deseason, /get_lun

deseasonal_lines = file_lines(temp_deseason)
deseasonal_arr_ss = fltarr(5, deseasonal_lines)
for u = 0, deseasonal_lines - 1 do begin
	readf, deseason, deseasonal_arr0, lat0, lon0, year0, month0
	deseasonal_arr_ss[ 0 , u ] = deseasonal_arr0
	deseasonal_arr_ss[ 1 , u ] = lat0
	deseasonal_arr_ss[ 2 , u ] = lon0
	deseasonal_arr_ss[ 3 , u ] = year0
	deseasonal_arr_ss[ 4 , u ] = month0
endfor

free_lun, deseason

;again, using the same method to seperate the months out from the deseasonal temp array from the deseason_all_noaa program
infile_bin01_de = "/home/jakechung/IDL/myIDL/temp_folder/temp_bin01_SS_deseason.dat"
infile_bin02_de = "/home/jakechung/IDL/myIDL/temp_folder/temp_bin02_SS_deseason.dat"
infile_bin03_de = "/home/jakechung/IDL/myIDL/temp_folder/temp_bin03_SS_deseason.dat"

;open to write in the temp files
openw , lun1 , infile_bin01_de , /get_lun
openw , lun2 , infile_bin02_de , /get_lun
openw , lun3 , infile_bin03_de , /get_lun

;loop to pull out the months, sort by latitude
for j = 0, deseasonal_lines - 1 do begin
	if (j NE deseasonal_lines - 1) then b = j + 1 else b = j
		if ( deseasonal_arr_ss[ 1 , j ] EQ deseasonal_arr_ss[ 1 , b ] ) then begin
			a = where( deseasonal_arr_ss[ 1 , * ] EQ deseasonal_arr_ss[ 1 , j ] )
			temp_de_year = deseasonal_arr_ss[ 3 , a ]
			temp_de_month = deseasonal_arr_ss[ 4 , a ]
			temp_de_var = deseasonal_arr_ss[ 0 , a ]
			bin01_de_index = where( temp_de_year EQ 1996 OR $
				temp_de_year EQ 1997 OR $
				temp_de_year EQ 1998 OR $
				temp_de_year EQ 1999 OR $
				temp_de_year EQ 2000 ) 
			bin02_de_index = where( temp_de_year EQ 2001 OR $
				temp_de_year EQ 2002 OR $
				temp_de_year EQ 2003 OR $
				temp_de_year EQ 2004 OR $
				temp_de_year EQ 2005 ) 
			bin03_de_index = where( temp_de_year EQ 2006 OR $
				temp_de_year EQ 2007 OR $
				temp_de_year EQ 2008 OR $
				temp_de_year EQ 2009 ) 
			;---------------------------------------------------------------------------------
			bin01_de_var = stddev( temp_de_var [ bin01_de_index ] ) / sqrt( n_elements( bin01_de_index ) )
			bin02_de_var = stddev( temp_de_var [ bin02_de_index ] ) / sqrt( n_elements( bin02_de_index ) )
			bin03_de_var = stddev( temp_de_var [ bin03_de_index ] ) / sqrt( n_elements( bin03_de_index ) )
			;----------------------------------------------------------------------------------
			printf, lun1, deseasonal_arr_ss[ 1 , j ], bin01_de_var
			printf, lun2, deseasonal_arr_ss[ 1 , j ], bin02_de_var
			printf, lun3, deseasonal_arr_ss[ 1 , j ], bin03_de_var
			;------------------------------------------------------------------------------------
			j = a[n_elements(a) - 1]					
			a = 0
		endif else if ( deseasonal_arr_ss[ 1 , j ] NE deseasonal_arr_ss[ 1 , b ] ) then j = j + 1
endfor

free_lun, lun1, lun2, lun3

;refresh the temp files to write the data into arrays for plotting
openr , lun1 , infile_bin01_de , /get_lun
openr , lun2 , infile_bin02_de , /get_lun
openr , lun3 , infile_bin03_de , /get_lun

bin01_de_size = file_lines( infile_bin01_de )
bin02_de_size = file_lines( infile_bin02_de )
bin03_de_size = file_lines( infile_bin03_de )

de_obs_bin01_ss = fltarr( 2 , bin01_de_size)
de_obs_bin02_ss = fltarr( 2 , bin02_de_size)
de_obs_bin03_ss = fltarr( 2 , bin03_de_size)

for t = 0, bin01_de_size - 1 do begin
	readf, lun1, lat_var_de, stderr
	de_obs_bin01_ss[ 0 , t ] = lat_var_de
	de_obs_bin01_ss[ 1 , t ] = stderr
endfor
for t = 0, bin02_de_size - 1 do begin
	readf, lun2, lat_var_de, stderr
	de_obs_bin02_ss[ 0 , t ] = lat_var_de
	de_obs_bin02_ss[ 1 , t ] = stderr
endfor
for t = 0, bin03_de_size - 1 do begin
	readf, lun3, lat_var_de, stderr
	de_obs_bin03_ss[ 0 , t ] = lat_var_de
	de_obs_bin03_ss[ 1 , t ] = stderr
endfor

free_lun, lun1, lun2, lun3

;remove the NAN values so array is less messy and NAN cannot stay for further calculations
i_nan = where( ~finite( bin01[ 2 , * ] ), count, /null )
if ( count GT 0 ) then begin
	bin01 = RemoveRows( bin01, i_nan )
	de_obs_bin01_ss = RemoveRows( de_obs_bin01_ss, i_nan )
endif
i_nan = where( ~finite( bin02[ 2 , * ] ), count, /null )
if ( count GT 0 ) then begin
	bin02 = RemoveRows( bin02, i_nan )
	de_obs_bin02_ss = RemoveRows( de_obs_bin02_ss, i_nan )
endif
i_nan = where( ~finite( bin03[ 2 , * ] ), count, /null )
if ( count GT 0 ) then begin
	bin03 = RemoveRows( bin03, i_nan )
	de_obs_bin03_ss = RemoveRows( de_obs_bin03_ss, i_nan )
endif
;=====================================
;===Start extracting simulated data===
;=====================================
; Define the filename

;sim_title = 'Aydin emissions scaled to Xiao et al over 1996-2003'
;filename = "/home/excluded-from-backup/data/C2H6/trac_avg.AydinSF_1981_2015.bpch"

;sim_title = 'Unscaled (absolute) Aydin emissions'
;filename = "/home/excluded-from-backup/data/C2H6/trac_avg.AydinAbsSF_1981_2015.bpch"

;sim_title = 'PSU emissions scaled to Xiao et al over 1996-2003'
;filename = "/home/excluded-from-backup/data/C2H6/trac_avg.PSUSF_1981_2015.bpch"

;sim_title = 'Unscaled PSU emissions using methane-to-ethane (MER) ratio of 5 for biomass burning, and 18 for fossil fuel'
;filename = "/home/excluded-from-backup/data/C2H6/trac_avg.PSU_MER3BB_MER18FF_1981_2015.bpch"

;sim_title = 'Unscaled emissions from Simpson et al.'
;filename = "/home/excluded-from-backup/data/C2H6/trac_avg.SimpsonSF_1981_2015.bpch"

;sim_title = 'Constant default base emissions '
;filename = "/home/excluded-from-backup/data/C2H6/trac_avg.spinup_GFED4_MAVG_1981_2015.198101010000"

;read in the file
;ctm_get_data, datainfo, 'IJ-AVG-$', tracer=1, filename=filename, tau0=tau0

; Get MODELINFO and GRIDINFO structures, xmid, ymid hold lon/lat centers 
;getmodelandgridinfo, datainfo[0], modelinfo, gridinfo, lon=xmid, lat=ymid
;if (n_elements(datainfo) eq 0) then print, 'No data records found.'

;nt is the number of data points in the simulation.
;nt = n_elements(datainfo)
;globalavg= fltarr(nt)
;print, 'Simulation ran for ', nt, ' months'

;convert tau0 values to year month day array and store in sim_yymmdd
;sim_yymmdd stores the entire time element of the input simulation
;sim_yymmdd= tau2yymmdd(tau0, /GEOS1)

;bin01 - - 2005 - 2006
;bin02 - - 2007 - 2008
;bin03 - - 2009 and 2010
;bin04 - - 2011 and 2012
;bin05 - - 2013 and 2014
;locate indexes of all months of the simulated data 
;sim_bin01_idx = where( sim_yymmdd.year EQ 1996 OR sim_yymmdd.year EQ 1997 )
;sim_bin02_idx = where( sim_yymmdd.year EQ 1998 OR sim_yymmdd.year EQ 1999 )
;sim_bin03_idx = where( sim_yymmdd.year EQ 2000 OR sim_yymmdd.year EQ 2001 )
;sim_bin04_idx = where( sim_yymmdd.year EQ 2002 OR sim_yymmdd.year EQ 2003 )
;sim_bin05_idx = where( sim_yymmdd.year EQ 2004 OR sim_yymmdd.year EQ 2005 )
;sim_bin06_idx = where( sim_yymmdd.year EQ 2006 OR sim_yymmdd.year EQ 2007 )
;sim_bin07_idx = where( sim_yymmdd.year EQ 2008 OR sim_yymmdd.year EQ 2009 )

;extracting the simulated data for the corresponding months and years and locations
;sim_bin01_var = fltarr(n_elements(bin01[0 , *]))
;for i = 0, n_elements(bin01[0 , *]) - 1 do begin
;	for t = 0, nt - 1 do begin
;		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
;		gridinfo= gridinfo, lat= [bin01[0, i] - 0.1, bin01[0, i] + 0.1], $
;		lon= [bin01[1, i] - 0.1, bin01[1, i] + 0.1], alrange = [1, 3], average=7)
;		globalavg[ t ]= data/2 * 1000
;	endfor
;	sim_bin01_var[ i ] = mean(globalavg[ sim_bin01_idx ])
;endfor

;sim_bin02_var = fltarr(n_elements(bin02[0 , *]))
;for i = 0, n_elements(bin02[0 , *]) - 1 do begin
;	for t = 0, nt - 1 do begin
;		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
;		gridinfo= gridinfo, lat= [bin02[0, i] - 0.1, bin02[0, i] + 0.1], $
;		lon= [bin02[1, i] - 0.1, bin02[1, i] + 0.1], alrange = [1, 3], average=7)
;		globalavg[ t ]= data/2 * 1000
;	endfor
;	sim_bin02_var[ i ] = mean(globalavg[ sim_bin02_idx ])
;endfor

;sim_bin03_var = fltarr(n_elements(bin03[0 , *]))
;for i = 0, n_elements(bin03[0 , *]) - 1 do begin
;	for t = 0, nt - 1 do begin
;		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
;		gridinfo= gridinfo, lat= [bin03[0, i] - 0.1, bin03[0, i] + 0.1], $
;		lon= [bin03[1, i] - 0.1, bin03[1, i] + 0.1], alrange = [1, 3], average=7)
;		globalavg[ t ]= data/2 * 1000
;	endfor
;	sim_bin03_var[ i ] = mean(globalavg[ sim_bin03_idx ])
;endfor

;sim_bin04_var = fltarr(n_elements(bin04[0 , *]))
;for i = 0, n_elements(bin04[0 , *]) - 1 do begin
;	for t = 0, nt - 1 do begin
;		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
;		gridinfo= gridinfo, lat= [bin04[0, i] - 0.1, bin04[0, i] + 0.1], $
;		lon= [bin04[1, i] - 0.1, bin04[1, i] + 0.1], alrange = [1, 3], average=7)
;		globalavg[ t ]= data/2 * 1000
;	endfor
;	sim_bin04_var[ i ] = mean(globalavg[ sim_bin04_idx ])
;endfor

;sim_bin05_var = fltarr(n_elements(bin05[0 , *]))
;for i = 0, n_elements(bin05[0 , *]) - 1 do begin
;	for t = 0, nt - 1 do begin
;		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
;		gridinfo= gridinfo, lat= [bin05[0, i] - 0.1, bin05[0, i] + 0.1], $
;		lon= [bin05[1, i] - 0.1, bin05[1, i] + 0.1], alrange = [1, 3], average=7)
;		globalavg[ t ]= data/2 * 1000
;	endfor
;	sim_bin05_var[ i ] = mean(globalavg[ sim_bin05_idx ])
;endfor

;sim_bin06_var = fltarr(n_elements(bin06[0 , *]))
;for i = 0, n_elements(bin06[0 , *]) - 1 do begin
;	for t = 0, nt - 1 do begin
;		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
;		gridinfo= gridinfo, lat= [bin06[0, i] - 0.1, bin06[0, i] + 0.1], $
;		lon= [bin06[1, i] - 0.1, bin06[1, i] + 0.1], alrange = [1, 3], average=7)
;		globalavg[ t ]= data/2 * 1000
;	endfor
;	sim_bin06_var[ i ] = mean(globalavg[ sim_bin06_idx ])
;endfor

;sim_bin07_var = fltarr(n_elements(bin07[0 , *]))
;for i = 0, n_elements(bin07[0 , *]) - 1 do begin
;	for t = 0, nt - 1 do begin
;		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
;		gridinfo= gridinfo, lat= [bin07[0, i] - 0.1, bin07[0, i] + 0.1], $
;		lon= [bin07[1, i] - 0.1, bin07[1, i] + 0.1], alrange = [1, 3], average=7)
;		globalavg[ t ]= data/2 * 1000
;	endfor
;	sim_bin07_var[ i ] = mean(globalavg[ sim_bin07_idx ])
;endfor

;Start plotting procudures
!P.Multi = 0
cgDisplay, 1100 , 970


	
;start plotting the data
;cgPlot, bin01[ 0 , * ], bin01[ 2 , * ], /overplot, psym = 7, err_yhigh = de_obs_bin01_ss[ 1 , * ], $
;	err_ylow = de_obs_bin01_ss[ 1 , * ], color = 'red'
;cgPlot, bin02[ 0 , * ], bin02[ 2 , * ], /overplot, psym = 7, err_yhigh = de_obs_bin02_ss[ 1 , * ], $
;	err_ylow = de_obs_bin02_ss[ 1 , * ], color = 'forest green'
;cgPlot, bin03[ 0 , * ], bin03[ 2 , * ], /overplot, psym = 6, err_yhigh = de_obs_bin03_ss[ 1 , * ], $
;	err_ylow = de_obs_bin03_ss[ 1 , * ], color = 'red'
;cgPlot, bin04[ 0 , * ], bin04[ 2 , * ], /overplot, psym = 6, err_yhigh = de_obs_bin04_ss[ 1 , * ], $
;	err_ylow = de_obs_bin04_ss[ 1 , * ], color = 'forest green'
;cgPlot, bin05[ 0 , * ], bin05[ 2 , * ], /overplot, psym = 5, err_yhigh = de_obs_bin05_ss[ 1 , * ], $
;	err_ylow = de_obs_bin05_ss[ 1 , * ], color = 'red'
;cgPlot, bin06[ 0 , * ], bin06[ 2 , * ], /overplot, psym = 5, err_yhigh = de_obs_bin06_ss[ 1 , * ], $
;	err_ylow = de_obs_bin06_ss[ 1 , * ], color = 'forest green'
;cgPlot, bin07[ 0 , * ], bin07[ 2 , * ], /overplot, psym = 9, err_yhigh = de_obs_bin07_ss[ 1 , * ], $
;	err_ylow = de_obs_bin07_ss[ 1 , * ], color = 'blue', symsize = 2

;Legend for the data points
;AL_LEGEND, [ '1996-1997' , '1998-1999' , '2000-2001' , '2002-2003' , '2004-2005' , '2006-2007' , '2008-2009'], $
;	color = ['red' , 'forest green' , 'red' , 'forest green' , 'red' , 'forest green' , 'blue'], $
;	psym = [ 7 , 7 , 6 , 6 , 5 , 5 , 9 ], $
;	position = [ -50 ,  3000 ], charsize = 2.0
	

;AL_LEGEND, [ '1996-1997' , '1998-1999' , '2000-2001' , '2002-2003' , '2004-2005' , '2006-2007' , '2008-2009'], $
;	color = ['red' , 'forest green' , 'pink' , 'green' , 'black' , 'maroon' , 'blue'], linestyle = [ 0 , 0 , 0 , 0 , 0 , 0 , 0 ]


;poly_fit the shit out of the data
measure_errors = de_obs_bin01_ss[ 1 , * ]
p_fit01 = poly_fit( bin01[ 0 , * ], bin01[ 2 , * ], 6, measure_errors = measure_errors, $
	sigma = sigma01, yfit = yfit01, /double )

measure_errors = de_obs_bin02_ss[ 1 , * ]
p_fit02 = poly_fit( bin02[ 0 , * ], bin02[ 2 , * ], 6, measure_errors = measure_errors, $
	sigma = sigma02, yfit = yfit02, /double )

measure_errors = de_obs_bin03_ss[ 1 , * ]
p_fit03 = poly_fit( bin03[ 0 , * ], bin03[ 2 , * ], 6, measure_errors = measure_errors, $
	sigma = sigma03, yfit = yfit03, /double )

;Set up the empty plot first
cgPlot, bin03[ 2 , * ], /nodata, ytitle = 'UCI Ethane mixing ratio (ppptv)', xtitle = 'Latitudes', $
	xticklen = 1, yticklen = 1, xGridStyle = 1, yGridStyle = 1, xrange = [ -90, 90], $
	yrange = [0, 2000]
cgPlot, bin01[ 0 , * ], bin01[ 2 , * ], color = 'blue', psym = 16, err_yhigh = de_obs_bin01_ss, $
	err_ylow = de_obs_bin01_ss, /overplot
cgPlot, bin02[ 0 , * ], bin02[ 2 , * ], color = 'forest green', psym = 16, err_yhigh = de_obs_bin02_ss, $
	err_ylow = de_obs_bin02_ss, /overplot
cgPlot, bin03[ 0 , * ], bin03[ 2 , * ], color = 'red', psym = 16, err_yhigh = de_obs_bin03_ss, $
	err_ylow = de_obs_bin03_ss, /overplot
cgPlot, bin01[ 0 , * ], yfit01, color = 'blue', /overplot
cgPlot, bin02[ 0 , * ], yfit02, color = 'forest green', /overplot
cgPlot, bin03[ 0 , * ], yfit03, color = 'red', /overplot

AL_LEGEND, ['1996-2000', '2001-2005', '2006-2009', '1996-2000 fit', '2001-2005 fit', '2006-2009 fit'], $
	psym = [ 16 , 16 , 16 , -3 , -3 , -3 ], $
	linestyle = [ 0,0,0,0,0,0 ], $
	position = [ -75, 1800], color = ['blue', 'forest green', 'red', 'blue', 'forest green', 'red']
;bin01_residual 
;poly fit the sim data
;p_fit01_sim = poly_fit( bin01[ 0 , * ], sim_bin01_var[ * ], 6, sigma = sigma01_sim, $
;	yband = yband01_sim, yfit = yfit01_sim, /double )

;p_fit02_sim = poly_fit( bin02[ 0 , * ], sim_bin02_var[ * ], 6, sigma = sigma02_sim, $
;	yband = yband02_sim, yfit = yfit02_sim, /double )
	
;p_fit03_sim = poly_fit( bin03[ 0 , * ], sim_bin03_var[ * ], 6, sigma = sigma03_sim, $
;	yband = yband03_sim, yfit = yfit03_sim, /double )
	
;make a super plot
;cgPlot, bin01[ 0 , * ], yfit01, yrange = [ 0 , 2000 ], xrange = [ -100 , 100 ], $
;	ytitle = 'UCI Ethane ratio (pptv)', yticklen = 1, xticklen = 1, $
;	yGridStyle = 1, xGridStyle = 1, Position=[0.15, 0.72, 0.90, 0.90], color = 'red', $
;	psym = -3, /NoErase, XTickformat = '(A1)'
;cgPlot, bin02[ 0 , * ], yfit02, color = 'blue', /overplot, psym = -3
;cgPlot, bin03[ 0 , * ], yfit03, color = 'forest green', /overplot, psym = -3
;cgPlot, bin01[ 0 , * ], bin01[ 2 , * ], psym = 3, err_yhigh = de_obs_bin01_ss[ 1 , * ], $
;	err_ylow = de_obs_bin01_ss[ 1 , * ], color = 'red', /overplot
;cgPlot, bin02[ 0 , * ], bin02[ 2 , * ], psym = 3, err_yhigh = de_obs_bin02_ss[ 1 , * ], $
;	err_ylow = de_obs_bin02_ss[ 1 , * ], color = 'blue', /overplot
;cgPlot, bin03[ 0 , * ], bin03[ 2 , * ], psym = 3, err_yhigh = de_obs_bin03_ss[ 1 , * ], $
;	err_ylow = de_obs_bin03_ss[ 1 , * ], color = 'forest green', /overplot
;cgLEGEND, title = ['1996-1997', '1998-1999', '2000-2001'], $
;	color = ['red', 'blue', 'forest green'], $
;;	linestyle = [ 0, 0, 0 ], $
;	psym = [ -16, -16, -16 ], $
;	location = [0.2, 0.87], $
;	/Box, /Background, $
;	BG_Color = 'rose', $
;	/Center_Sym

;cgPlot, bin01[ 0 , * ], yfit01_sim, XTickformat = '(A1)', yrange = [ 0 , 2000 ], $
;	position = [0.15, 0.52, 0.90, 0.70], yTicklen = 1, xTicklen = 1, yGridStyle = 1, $
;	xGridStyle = 1, xrange = [ -100, 100 ], ytitle = 'Sim. Ethane ratio (pptv)', $
;	color = 'red', linestyle = 2, err_yhigh = yband01_sim, err_ylow = yband01_sim, $
;	thick = 1.5, xStyle = 8, /NoErase
;cgPlot, bin02[ 0 , * ], yfit02_sim, linestyle = 2, color = 'blue', err_yhigh = yband02_sim, $
;	err_ylow = yband02_sim, /overplot, 	thick = 1.5
;cgPlot, bin03[ 0 , * ], yfit03_sim, linestyle = 2, color = 'forest green', err_yhigh = yband03_sim, $
;	err_ylow = yband03_sim, /overplot, 	thick = 1.5


;cgLEGEND, title =  ['1996-1997', '1998-1999', '2000-2001'], $
;	color = ['red', 'blue', 'forest green'], $
;	linestyle = [ 2, 2, 2 ], $
;	location = [0.2, 0.67], $
;	/Box, /Background, $
;	BG_Color = 'rose'
;Plots, [0.15, 0.15], [0.70, 0.72], /Normal ; Fix left axis.
;Plots, [0.90, 0.90], [0.70, 0.72], /Normal ; Fix right axis.

;make a super plot
;cgPlot, bin04[ 0 , * ], yfit04, yrange = [ 0 , 2000 ], xrange = [ -100 , 100 ], $
;	ytitle = 'UCI Ethane ratio (pptv)', yticklen = 1, xticklen = 1, $
;	yGridStyle = 1, xGridStyle = 1, Position=[0.15, 0.31, 0.90, 0.50], color = 'red', $
;	psym = -3, /NoErase, xStyle = 8, XTickformat = '(A1)'
;cgPlot, bin05[ 0 , * ], yfit05, color = 'blue', /overplot, psym = -3
;cgPlot, bin06[ 0 , * ], yfit06, color = 'forest green', /overplot, psym = -3
;cgPlot, bin07[ 0 , * ], yfit07, color = 'maroon', /overplot, psym = -3
;cgPlot, bin04[ 0 , * ], bin04[ 2 , * ], psym = 3, err_yhigh = de_obs_bin04_ss[ 1 , * ], $
;	err_ylow = de_obs_bin04_ss[ 1 , * ], color = 'red', /overplot
;cgPlot, bin05[ 0 , * ], bin05[ 2 , * ], psym = 3, err_yhigh = de_obs_bin05_ss[ 1 , * ], $
;	err_ylow = de_obs_bin05_ss[ 1 , * ], color = 'blue', /overplot
;cgPlot, bin06[ 0 , * ], bin06[ 2 , * ], psym = 3, err_yhigh = de_obs_bin06_ss[ 1 , * ], $
;	err_ylow = de_obs_bin06_ss[ 1 , * ], color = 'forest green', /overplot
;cgPlot, bin07[ 0 , * ], bin07[ 2 , * ], psym = 3, err_yhigh = de_obs_bin07_ss[ 1 , * ], $
;	err_ylow = de_obs_bin07_ss[ 1 , * ], color = 'maroon', /overplot
;cgLEGEND, title =  ['2002-2003', '2004-2005', '2006-2007', '2008-2009'], $
;	color = ['red', 'blue', 'forest green', 'maroon'], $
;	linestyle = [ 0, 0, 0, 0 ], $
;	psym = [ -16, -16, -16, -16 ], $
;	location = [0.2, 0.47], $
;	/Box, /Background, $
;	BG_Color = 'rose', $
;	/Center_Sym
	
;Plots, [0.15, 0.15], [0.28, 0.31], /Normal ; Fix left axis.
;Plots, [0.90, 0.90], [0.28, 0.31], /Normal ; Fix right axis.
;
;cgPlot, bin04[ 0 , * ], yfit04_sim, yrange = [ 0 , 2000 ], $
;	position = [0.15, 0.10, 0.90, 0.29], yTicklen = 1, xTicklen = 1, yGridStyle = 1, $
;	xGridStyle = 1, xrange = [ -100, 100 ], ytitle = 'Sim. Ethane ratio (pptv)', $
;	color = 'red', linestyle = 2, err_yhigh = yband04_sim, err_ylow = yband04_sim, $
;	thick = 1.5,  xtitle = 'Latitudes', /NoErase, xStyle = 8
;cgPlot, bin05[ 0 , * ], yfit05_sim, linestyle = 2, color = 'blue', err_yhigh = yband05_sim, $
;	err_ylow = yband05_sim, /overplot, 	thick = 1.5
;cgPlot, bin06[ 0 , * ], yfit06_sim, linestyle = 2, color = 'forest green', err_yhigh = yband06_sim, $
;	err_ylow = yband06_sim, /overplot, 	thick = 1.5
;cgPlot, bin07[ 0 , * ], yfit07_sim, linestyle = 2, color = 'maroon', err_yhigh = yband07_sim, $
;	err_ylow = yband07_sim, /overplot, 	thick = 1.5

;cgLEGEND, title =  ['2002-2003', '2004-2005', '2006-2007', '2008-2009'], $
;	color = ['red', 'blue', 'forest green', 'maroon'], $
;	linestyle = [ 2, 2, 2, 2 ], $
;	location = [0.2, 0.25], $
;	/Box, /Background, $
;	BG_Color = 'rose'
;Plots, [0.15, 0.15], [0.50, 0.52], /Normal ; Fix left axis.
;Plots, [0.90, 0.90], [0.50, 0.52], /Normal ; Fix right axis.


;cgText, 0.2, 0.3, sim_title, ALIGNMENT= 0.5


;Make a JPEG for the graphs
READ, option, PROMPT= 'Do I need to save this graph as a jpeg so I can print it out later? Enter "1" for HELLYEAH or "0" for NOPE       '
name = strarr(1)
IF (option EQ 1) THEN BEGIN
	READ, name, PROMPT= 'Name of the JPEG file? (No need to include the jpeg extension, just name is fine)   '
	SCREEN2JPG, name
	print, 'JPEG has been saved to current directory'
ENDIF

;plot the poly fitted line on the graph
;cgPlot, bin01[ 0 , * ], yfit01, color = 'red', /overplot ;1996-1997
;cgPlot, bin02[ 0 , * ], yfit02, color = 'maroon', /overplot ;1998-1999
;cgPlot, bin03[ 0 , * ], yfit03, color = 'yellow', /overplot ;2000-2001
;cgPlot, bin04[ 0 , * ], yfit04, color = 'forest green', /overplot ;2002-2003
;cgPlot, bin05[ 0 , * ], yfit05, color = 'green', /overplot ;2004-2005
;cgPlot, bin06[ 0 , * ], yfit06, color = 'blue', /overplot ;2006-2007
;cgPlot, bin07[ 0 , * ], yfit07, color = 'violet', /overplot ;2008-2009

end
