pro annual_avg_lat_dis_noaa

;This program calculate the annual means of the NOAA data set and plot all of it on a graph of latitude vs concentration

;There are four major parts
;<<PART ONE>> 
;	Extract the NOAA data from the ASCII files
;<<PART TWO>>
;	Classify the data into bins to calculate the means 
;<<PART THREE>>
;	Read in the deseasonal data to calculate the uncertainty of the graph
;<<PART FOUR>>
;	Deal with simulated data
;<<PART FIVE>>
;	Plotting routines, note: might want to hold off on plotting the error bars or else it will be very cluttering 


compile_opt IDL2				;set compile options


;============================================================================================================
;PART ONE===================================================================================================
;=====Read in NOAA data=====
;get the files
infile1= "/home/jakechung/IDL/myIDL/NOAA_data.dat"
infile2= "/home/jakechung/IDL/myIDL/NOAA_locations.dat"

;get number of lines in the data file
n1_noaa = file_lines(infile1)
n2_noaa = file_lines(infile2)

print, 'The NOAA data has ', n1_noaa, '   data lines with', n2_noaa, '  locations'

;create an empty array to store the data that was read in lon, lat, avg, date and locations
lat1_noaa= fltarr(n1_noaa)
lon1_noaa= fltarr(n1_noaa)
avg_noaa= fltarr(n1_noaa)
location_noaa= strarr(n2_noaa)
year_noaa= fltarr(n1_noaa)
month_noaa= fltarr(n1_noaa)
day_noaa= fltarr(n1_noaa)
hour_noaa= fltarr(n1_noaa)
alt_noaa= fltarr(n1_noaa)
time_noaa= fltarr(n1_noaa)

;FLOAT input variable
sample_site_code= ' '
sample_latitude= 0.0
sample_longitude= 0.0
analysis_value= 0.0
sample_year= 0.0
sample_month= 0.0
sample_day= 0.0
sample_hour= 0.0
sample_altitude= 0.0

;open input files
openr, iunit1, infile1, /get_lun
openr, iunit2, infile2, /get_lun

;read float data of the data file
for i = 0, n1_noaa-1 do begin
	;read one line from the file
	readf, iunit1, sample_year, sample_month, sample_day, sample_hour, analysis_value,	$
				sample_latitude, sample_longitude, sample_altitude
	;store the values
	year_noaa[i] = sample_year
	month_noaa[i] = sample_month
	day_noaa[i] = sample_day
	hour_noaa[i] = sample_hour
	avg_noaa[i] = analysis_value
	lat1_noaa[i] = sample_latitude
	lon1_noaa[i] = sample_longitude
	alt_noaa[i] = sample_altitude
endfor

;read float data of the location file
for i = 0, n2_noaa-1 do begin
	readf, iunit2, sample_site_code
	location_noaa[i] = sample_site_code
endfor

;close input file
free_lun, iunit1, iunit2


;combine the year, month, day, hour into a decimal year values
;in this script this line is not necesary but keep it here anyway
for i = 0, n1_noaa-1 do begin	
	time_noaa[i] = year_noaa[i] + hour_noaa[i]/8544 + day_noaa[i]/365 + month_noaa[i]/12
endfor

;Remove the -999 values in the C2H6 arrays
neg999 = where(avg_noaa LT 0, count)
location_noaa = RemoveRows(rotate(location_noaa, 1), neg999)
time_noaa = RemoveRows(rotate(time_noaa, 1), neg999)
lat1_noaa = RemoveRows(rotate(lat1_noaa, 1), neg999)
lon1_noaa = RemoveRows(rotate(lon1_noaa, 1), neg999)
avg_noaa = RemoveRows(rotate(avg_noaa, 1), neg999)
year_noaa = RemoveRows(rotate(year_noaa, 1), neg999)
month_noaa = RemoveRows(rotate(month_noaa, 1), neg999)
day_noaa = RemoveRows(rotate(day_noaa, 1), neg999)

;check the new size of the two arrays, n1 and n2 should be the same.
n2_noaa = n_elements(location_noaa)
n1_noaa = n_elements(avg_noaa)
print, n1_noaa, n2_noaa

;FINISH PART ONE======================================================================================


;START PART TWO=======================================================================================

;link to file
noaa_ids_file = "/home/jakechung/IDL/myIDL/NOAA_stations_IDs.dat"

;check how many lines in the file 
noaa_ids_size = file_lines( noaa_ids_file )

;create string array to store those stations ids
noaa_ids = strarr( noaa_ids_size )
ids = ' ' 

;open to read
openr, ids_lun, noaa_ids_file, /get_lun

;read the data into the noaa_ids array
for i = 0, noaa_ids_size - 1 do begin
	readf, ids_lun, ids
	noaa_ids [ i ] = ids 
endfor 

free_lun, ids_lun
;bin01 2005 2006
;bin02 2007 2008
;bin03 2009 2010
;bin04 2011 2012
;bin05 2013 2014
;bin06 2015 
;setting up the arrays that will hold the annual values
bin01 = fltarr( 3 , noaa_ids_size )
bin02 = fltarr( 3 , noaa_ids_size )
bin03 = fltarr( 3 , noaa_ids_size )
bin04 = fltarr( 3 , noaa_ids_size )
bin05 = fltarr( 3 , noaa_ids_size )
bin06 = fltarr( 3 , noaa_ids_size )

;loop to find the station and extract the year.
;version 2.0: change the variable naming scheme to bin numbers instead of year numbers so if needed be, I can do a 2 - year or 5 - year mean
for i = 0, noaa_ids_size - 1 do begin
	a = where( strmatch ( location_noaa[ * ], noaa_ids[ i ], /FOLD_CASE) EQ 1 )
	temp_lat = lat1_noaa[ a ]
	temp_lon = lon1_noaa[ a ]
	temp_avg = avg_noaa[ a ]
	temp_bin01 = where( year_noaa[ a ] EQ 2005 OR year_noaa[ a ] EQ 2006 , count )
	;having a if-then filter here to make sure that when a station does not have a certain year data, 
	;it will not send the script into frenzy and writing false values
	if ( count GT 0 ) then begin
		bin01[ 2 , i ] = mean( temp_avg[ temp_bin01 ], /NAN )	
		bin01[ 0 , i ] = temp_lat[ 0 ]
		bin01[ 1 , i ] = temp_lon[ 0 ]
		endif else begin	
			bin01[ 2 , i ] = !Values.F_NAN
			bin01[ 0 , i ] = temp_lat[ 0 ]
			bin01[ 1 , i ] = temp_lon[ 0 ]
			endelse
	
	temp_bin02 = where( year_noaa[ a ] EQ 2007 OR year_noaa[ a ] EQ 2008 , count )
	if ( count GT 0 ) then begin
		bin02[ 2 , i ] = mean( temp_avg[ temp_bin02 ], /NAN )	
		bin02[ 0 , i ] = temp_lat[ 0 ]
		bin02[ 1 , i ] = temp_lon[ 0 ]
		endif else begin	
			bin02[ 2 , i ] = !Values.F_NAN
			bin02[ 0 , i ] = temp_lat[ 0 ]
			bin02[ 1 , i ] = temp_lon[ 0 ]
			endelse
		
	temp_bin03 = where( year_noaa[ a ] EQ 2009 OR year_noaa[ a ] EQ 2010 , count )
	if ( count GT 0 ) then begin
		bin03[ 2 , i ] = mean( temp_avg[ temp_bin03 ], /NAN )	
		bin03[ 0 , i ] = temp_lat[ 0 ]
		bin03[ 1 , i ] = temp_lon[ 0 ]
		endif else begin	
			bin03[ 2 , i ] = !Values.F_NAN
			bin03[ 0 , i ] = temp_lat[ 0 ]
			bin03[ 1 , i ] = temp_lon[ 0 ]
			endelse

	temp_bin04 = where( year_noaa[ a ] EQ 2011 OR year_noaa[ a ] EQ 2012 , count )
	if ( count GT 0 ) then begin
		bin04[ 2 , i ] = mean( temp_avg[ temp_bin04 ], /NAN )	
		bin04[ 0 , i ] = temp_lat[ 0 ]
		bin04[ 1 , i ] = temp_lon[ 0 ]
		endif else begin	
			bin04[ 2 , i ] = !Values.F_NAN
			bin04[ 0 , i ] = temp_lat[ 0 ]
			bin04[ 1 , i ] = temp_lon[ 0 ]
			endelse

	temp_bin05 = where( year_noaa[ a ] EQ 2013 OR year_noaa[ a ] EQ 2014 , count )
	if ( count GT 0 ) then begin
		bin05[ 2 , i ] = mean( temp_avg[ temp_bin05 ] )	
		bin05[ 0 , i ] = temp_lat[ 0 ]
		bin05[ 1 , i ] = temp_lon[ 0 ]
		endif else begin	
			bin05[ 2 , i ] = !Values.F_NAN
			bin05[ 0 , i ] = temp_lat[ 0 ]
			bin05[ 1 , i ] = temp_lon[ 0 ]
			endelse
			
	temp_bin06 = where( year_noaa[ a ] EQ 2015 , count )
	if ( count GT 0 ) then begin
		bin06[ 2 , i ] = mean( temp_avg[ temp_bin06 ] )	
		bin06[ 0 , i ] = temp_lat[ 0 ]
		bin06[ 1 , i ] = temp_lon[ 0 ]
		endif else begin	
			bin06[ 2 , i ] = !Values.F_NAN
			bin06[ 0 , i ] = temp_lat[ 0 ]
			bin06[ 1 , i ] = temp_lon[ 0 ]
			endelse
		
endfor 

;finish PART TWO=========================================================================================


;start PART THREE========================================================================================

;This next section will read in the deseasonal values from the deseason_all_noaa_v2 program ASCII file
;Create link to the data file
temp_deseason = "/home/jakechung/IDL/myIDL/temp_folder/temp_residual_noaa.dat"
openr, deseason, temp_deseason, /get_lun

deseasonal_lines = file_lines(temp_deseason)
deseasonal_arr_noaa = fltarr(5, deseasonal_lines)
for u = 0, deseasonal_lines - 1 do begin
	readf, deseason, deseasonal_arr0, lat0, lon0, year0, month0
	deseasonal_arr_noaa[ 0 , u ] = deseasonal_arr0
	deseasonal_arr_noaa[ 1 , u ] = lat0
	deseasonal_arr_noaa[ 2 , u ] = lon0
	deseasonal_arr_noaa[ 3 , u ] = year0
	deseasonal_arr_noaa[ 4 , u ] = month0
endfor

free_lun, deseason

;again, using the same method to seperate the months out from the deseasonal temp array from the deseason_all_noaa program
infile_bin01_de = "/home/jakechung/IDL/myIDL/temp_folder/temp_bin01_NOAA_deseason.dat"
infile_bin02_de = "/home/jakechung/IDL/myIDL/temp_folder/temp_bin02_NOAA_deseason.dat"
infile_bin03_de = "/home/jakechung/IDL/myIDL/temp_folder/temp_bin03_NOAA_deseason.dat"
infile_bin04_de = "/home/jakechung/IDL/myIDL/temp_folder/temp_bin04_NOAA_deseason.dat"
infile_bin05_de = "/home/jakechung/IDL/myIDL/temp_folder/temp_bin05_NOAA_deseason.dat"
infile_bin06_de = "/home/jakechung/IDL/myIDL/temp_folder/temp_bin06_NOAA_deseason.dat"

;open to write in the temp files
openw , lun1 , infile_bin01_de , /get_lun
openw , lun2 , infile_bin02_de , /get_lun
openw , lun3 , infile_bin03_de , /get_lun
openw , lun4 , infile_bin04_de , /get_lun
openw , lun5 , infile_bin05_de , /get_lun
openw , lun6 , infile_bin06_de , /get_lun

;loop to pull out the months, sort by latitude
for j = 0, deseasonal_lines - 1 do begin
	if (j NE deseasonal_lines - 1) then b = j + 1 else b = j
		if ( deseasonal_arr_noaa[ 1 , j ] EQ deseasonal_arr_noaa[ 1 , b ] ) then begin
			a = where( deseasonal_arr_noaa[ 1 , * ] EQ deseasonal_arr_noaa[ 1 , j ] )
			temp_de_year = deseasonal_arr_noaa[ 3 , a ]
			temp_de_month = deseasonal_arr_noaa[ 4 , a ]
			temp_de_var = deseasonal_arr_noaa[ 0 , a ]
			bin01_de_idx = where( temp_de_year EQ 2005 OR temp_de_year EQ 2006)
			bin02_de_idx = where( temp_de_year EQ 2007 OR temp_de_year EQ 2008)
			bin03_de_idx = where( temp_de_year EQ 2009 OR temp_de_year EQ 2010) 
			bin04_de_idx = where( temp_de_year EQ 2011 OR temp_de_year EQ 2012 ) 
			bin05_de_idx = where( temp_de_year EQ 2013 OR temp_de_year EQ 2014 ) 
			bin06_de_idx = where( temp_de_year EQ 2015 ) 
			;---------------------------------------------------------------------------------
			bin01_de_var = stddev( temp_de_var [ bin01_de_idx ] ) / sqrt( n_elements( bin01_de_idx ) )
			bin02_de_var = stddev( temp_de_var [ bin02_de_idx ] ) / sqrt( n_elements( bin02_de_idx ) )
			bin03_de_var = stddev( temp_de_var [ bin03_de_idx ] ) / sqrt( n_elements( bin03_de_idx ) )
			bin04_de_var = stddev( temp_de_var [ bin04_de_idx ] ) / sqrt( n_elements( bin04_de_idx ) )
			bin05_de_var = stddev( temp_de_var [ bin05_de_idx ] ) / sqrt( n_elements( bin05_de_idx ) )
			bin06_de_var = stddev( temp_de_var [ bin06_de_idx ] ) / sqrt( n_elements( bin06_de_idx ) )
			;----------------------------------------------------------------------------------
			printf, lun1, deseasonal_arr_noaa[ 1 , j ], deseasonal_arr_noaa[ 3 , j ], bin01_de_var
			printf, lun2, deseasonal_arr_noaa[ 1 , j ], deseasonal_arr_noaa[ 3 , j ], bin02_de_var
			printf, lun3, deseasonal_arr_noaa[ 1 , j ], deseasonal_arr_noaa[ 3 , j ], bin03_de_var
			printf, lun4, deseasonal_arr_noaa[ 1 , j ], deseasonal_arr_noaa[ 3 , j ], bin04_de_var
			printf, lun5, deseasonal_arr_noaa[ 1 , j ], deseasonal_arr_noaa[ 3 , j ], bin05_de_var
			printf, lun6, deseasonal_arr_noaa[ 1 , j ], deseasonal_arr_noaa[ 3 , j ], bin06_de_var
			;------------------------------------------------------------------------------------
			j = a[n_elements(a) - 1]					
			a = 0
		endif else if ( deseasonal_arr_noaa[ 1 , j ] NE deseasonal_arr_noaa[ 1 , b ] ) then j = j + 1
endfor

free_lun, lun1, lun2, lun3, lun4, lun5, lun6

;refresh the temp files to write the data into arrays for plotting
openr , lun1 , infile_bin01_de , /get_lun
openr , lun2 , infile_bin02_de , /get_lun
openr , lun3 , infile_bin03_de , /get_lun
openr , lun4 , infile_bin04_de , /get_lun
openr , lun5 , infile_bin05_de , /get_lun
openr , lun6 , infile_bin06_de , /get_lun

bin01_de_size = file_lines( infile_bin01_de )
bin02_de_size = file_lines( infile_bin02_de )
bin03_de_size = file_lines( infile_bin03_de )
bin04_de_size = file_lines( infile_bin04_de )
bin05_de_size = file_lines( infile_bin05_de )
bin06_de_size = file_lines( infile_bin06_de )

de_obs_bin01_noaa = fltarr( 3 , bin01_de_size)
de_obs_bin02_noaa = fltarr( 3 , bin02_de_size)
de_obs_bin03_noaa = fltarr( 3 , bin03_de_size)
de_obs_bin04_noaa = fltarr( 3 , bin04_de_size)
de_obs_bin05_noaa = fltarr( 3 , bin05_de_size)
de_obs_bin06_noaa = fltarr( 3 , bin06_de_size)

for t = 0, bin01_de_size - 1 do begin
	readf, lun1, lat_var_de, year_var, stderr
	de_obs_bin01_noaa[ 0 , t ] = lat_var_de
	de_obs_bin01_noaa[ 1 , t ] = year_var
	de_obs_bin01_noaa[ 2 , t ] = stderr
endfor
for t = 0, bin02_de_size - 1 do begin
	readf, lun2, lat_var_de, year_var, stderr
	de_obs_bin02_noaa[ 0 , t ] = lat_var_de
	de_obs_bin02_noaa[ 1 , t ] = year_var
	de_obs_bin02_noaa[ 2 , t ] = stderr
endfor
for t = 0, bin03_de_size - 1 do begin
	readf, lun3, lat_var_de, year_var, stderr
	de_obs_bin03_noaa[ 0 , t ] = lat_var_de
	de_obs_bin03_noaa[ 1 , t ] = year_var
	de_obs_bin03_noaa[ 2 , t ] = stderr
endfor
for t = 0, bin04_de_size - 1 do begin
	readf, lun4, lat_var_de, year_var, stderr
	de_obs_bin04_noaa[ 0 , t ] = lat_var_de
	de_obs_bin04_noaa[ 1 , t ] = year_var
	de_obs_bin04_noaa[ 2 , t ] = stderr
endfor
for t = 0, bin05_de_size - 1 do begin
	readf, lun5, lat_var_de, year_var, stderr
	de_obs_bin05_noaa[ 0 , t ] = lat_var_de
	de_obs_bin05_noaa[ 1 , t ] = year_var
	de_obs_bin05_noaa[ 2 , t ] = stderr
endfor
for t = 0, bin06_de_size - 1 do begin
	readf, lun6, lat_var_de, year_var, stderr
	de_obs_bin06_noaa[ 0 , t ] = lat_var_de
	de_obs_bin06_noaa[ 1 , t ] = year_var
	de_obs_bin06_noaa[ 2 , t ] = stderr
endfor

free_lun, lun1, lun2, lun3, lun4, lun5, lun6


;start PART FOUR=========================================================================================
;dealing with simulated data. sim_data will be extracted out and average yearly or every 2 years, depends on the bins

;remove the NAN values so the poly_fit function would work and to pull out the correct coordinates for the simulated data 
i_nan = where( ~finite( bin01[ 2 , * ] ), /null )
bin01 = RemoveRows( bin01, i_nan )
de_obs_bin01_noaa = RemoveRows( de_obs_bin01_noaa, i_nan )
i_nan = where( ~finite( bin02[ 2 , * ] ), /null )
bin02 = RemoveRows( bin02, i_nan )
de_obs_bin02_noaa = RemoveRows( de_obs_bin02_noaa, i_nan )
i_nan = where( ~finite( bin03[ 2 , * ] ), /null )
bin03 = RemoveRows( bin03, i_nan )
de_obs_bin03_noaa = RemoveRows( de_obs_bin03_noaa, i_nan )
i_nan = where( ~finite( bin04[ 2 , * ] ), /null )
bin04 = RemoveRows( bin04, i_nan )
de_obs_bin04_noaa = RemoveRows( de_obs_bin04_noaa, i_nan )
i_nan = where( ~finite( bin05[ 2 , * ] ), /null )
bin05 = RemoveRows( bin05, i_nan )
de_obs_bin05_noaa = RemoveRows( de_obs_bin05_noaa, i_nan )
i_nan = where( ~finite( bin06[ 2 , * ] ), /null )
bin06 = RemoveRows( bin06, i_nan )
de_obs_bin06_noaa = RemoveRows( de_obs_bin06_noaa, i_nan )

;PART FOUR==================================================================================================================

;=====================================
;===Start extracting simulated data===
;=====================================
; Define the filename

;Aydin emissions scaled to Xiao et al over 1996-2003
;filename = "/home/excluded-from-backup/data/C2H6/trac_avg.AydinSF_1981_2015.bpch"

;Unscaled (absolute) Aydin emissions
;ilename = "/home/excluded-from-backup/data/C2H6/trac_avg.AydinAbsSF_1981_2015.bpch"

; PSU emissions scaled to Xiao et al over 1996-2003
;filename = "/home/excluded-from-backup/data/C2H6/trac_avg.PSUSF_1981_2015.bpch"

; Unscaled PSU emissions using methane-to-ethane (MER) ratio of 5 for biomass burning, and 18 for fossil fuel
;filename = "/home/excluded-from-backup/data/C2H6/trac_avg.PSU_MER3BB_MER18FF_1981_2015.bpch"

; Unscaled emissions from Simpson et al.
;filename = "/home/excluded-from-backup/data/C2H6/trac_avg.SimpsonSF_1981_2015.bpch"

;Constant default base emissions 
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
;sim_bin01_idx = where( sim_yymmdd.year EQ 2005 )
;sim_bin02_idx = where( sim_yymmdd.year EQ 2006 )
;sim_bin03_idx = where( sim_yymmdd.year EQ 2007 )
;sim_bin04_idx = where( sim_yymmdd.year EQ 2008 )
;sim_bin05_idx = where( sim_yymmdd.year EQ 2009 )
;sim_bin06_idx = where( sim_yymmdd.year EQ 2010 )
;sim_bin07_idx = where( sim_yymmdd.year EQ 2011 )
;sim_bin08_idx = where( sim_yymmdd.year EQ 2012 )
;sim_bin09_idx = where( sim_yymmdd.year EQ 2013 )
;sim_bin10_idx = where( sim_yymmdd.year EQ 2014 )
;sim_bin11_idx = where( sim_yymmdd.year EQ 2015 )

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
;
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
;
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
;
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
;
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
;
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
;
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
;
;sim_bin08_var = fltarr(n_elements(bin08[0 , *]))
;for i = 0, n_elements(bin08[0 , *]) - 1 do begin
;	for t = 0, nt - 1 do begin
;		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
;		gridinfo= gridinfo, lat= [bin08[0, i] - 0.1, bin08[0, i] + 0.1], $
;		lon= [bin08[1, i] - 0.1, bin08[1, i] + 0.1], alrange = [1, 3], average=7)
;		globalavg[ t ]= data/2 * 1000
;	endfor
;	sim_bin08_var[ i ] = mean(globalavg[ sim_bin08_idx ])
;endfor
;
;sim_bin09_var = fltarr(n_elements(bin09[0 , *]))
;for i = 0, n_elements(bin09[0 , *]) - 1 do begin
;	for t = 0, nt - 1 do begin
;		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
;		gridinfo= gridinfo, lat= [bin09[0, i] - 0.1, bin09[0, i] + 0.1], $
;		lon= [bin09[1, i] - 0.1, bin09[1, i] + 0.1], alrange = [1, 3], average=7)
;		globalavg[ t ]= data/2 * 1000
;	endfor
;	sim_bin09_var[ i ] = mean(globalavg[ sim_bin09_idx ])
;endfor
;
;sim_bin10_var = fltarr(n_elements(bin10[0 , *]))
;for i = 0, n_elements(bin10[0 , *]) - 1 do begin
;	for t = 0, nt - 1 do begin
;		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
;		gridinfo= gridinfo, lat= [bin10[0, i] - 0.1, bin10[0, i] + 0.1], $
;		lon= [bin10[1, i] - 0.1, bin10[1, i] + 0.1], alrange = [1, 3], average=7)
;		globalavg[ t ]= data/2 * 1000
;	endfor
;	sim_bin10_var[ i ] = mean(globalavg[ sim_bin10_idx ])
;endfor
;
;sim_bin11_var = fltarr(n_elements(bin11[0 , *]))
;for i = 0, n_elements(bin11[0 , *]) - 1 do begin
;	for t = 0, nt - 1 do begin
;		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
;		gridinfo= gridinfo, lat= [bin11[0, i] - 0.1, bin11[0, i] + 0.1], $
;		lon= [bin11[1, i] - 0.1, bin11[1, i] + 0.1], alrange = [1, 3], average=7)
;		globalavg[ t ]= data/2 * 1000
;	endfor
;	sim_bin11_var[ i ] = mean(globalavg[ sim_bin11_idx ])
;endfor

;Finish PART FOUR=================================================================================================
;preliminary data processing before plotting routines

;There is an outlier in bin05, need to remove it manually
a = where( bin05[ 0 , * ] EQ 54.95)
bin05 = RemoveRows( bin05, a )
de_obs_bin05_noaa = RemoveRows( de_obs_bin05_noaa, a )

!P.Multi = [ 0 , 1 , 2 , 0 , 0 ]
cgDisplay, 1200 , 1000
title = 'Lat. Distribution of NOAA data averaging every two years'
;plotting annually
;setting up the plot
cgPlot, bin01[ 2 , * ], xtitle = 'Latitudes (degrees)', ytitle = 'Ethane mixing ratio (pptv)', title = title, $
	yrange = [ 0 , 4000 ], xrange = [ -90, 90], xticklen = 1.0, yticklen = 1.0, xGridstyle = 1, yGridstyle = 1, /nodata
	;================================================================================================
cgPlot, bin01[ 0 , * ], bin01[ 2 , * ], psym = 16, color = 'blue', /overplot, err_yhigh = de_obs_bin01_noaa[ 2 , * ], $
	err_ylow = de_obs_bin01_noaa[ 2 , * ]
cgPlot, bin02[ 0 , * ], bin02[ 2 , * ], psym = 16, color = 'forest green', /overplot, err_yhigh = de_obs_bin02_noaa[ 2 , * ], $
	err_ylow = de_obs_bin02_noaa[ 2 , * ]
cgPlot, bin03[ 0 , * ], bin03[ 2 , * ], psym = 16, color = 'red', /overplot, err_yhigh = de_obs_bin03_noaa[ 2 , * ], $
	err_ylow = de_obs_bin03_noaa[ 2 , * ]
;poly_fit the shit out of the data, yeah!!!! all the data, ALL OF THEM!!!!!


measure_errors = de_obs_bin01_noaa[ 2 , * ]
p_fit01 = poly_fit( bin01[ 0 , * ], bin01[ 2 , * ], 6, measure_errors = measure_errors, sigma = sigma01, yfit = yfit01, /double )
measure_errors = de_obs_bin02_noaa[ 2 , * ]
p_fit02 = poly_fit( bin02[ 0 , * ], bin02[ 2 , * ], 6, measure_errors = measure_errors, sigma = sigma02, yfit = yfit02, /double )
measure_errors = de_obs_bin03_noaa[ 2 , * ]
p_fit03 = poly_fit( bin03[ 0 , * ], bin03[ 2 , * ], 6, measure_errors = measure_errors, sigma = sigma03, yfit = yfit03, /double )
;poly_fit the simulated data
;p_fit01_sim = poly_fit( bin01[ 0 , * ], sim_bin01_var[ * ], 6, sigma = sigma01_sim, $
;	yfit = yfit01_sim, /double, yband = yband_sim01 )
;p_fit02_sim = poly_fit( bin02[ 0 , * ], sim_bin02_var[ * ], 6, sigma = sigma02_sim, $
;	yfit = yfit02_sim, /double, yband = yband_sim02 )
;p_fit03_sim = poly_fit( bin03[ 0 , * ], sim_bin03_var[ * ], 6, sigma = sigma03_sim, $
;	yfit = yfit03_sim, /double, yband = yband_sim03 )

;calculate the residual of the measured data and its poly fit to estimate the uncertainty
bin01_residual = yfit01 - bin01[ 2 , * ]	
bin02_residual = yfit02 - bin02[ 2 , * ]	
bin03_residual = yfit03 - bin03[ 2 , * ]	
	
;plot the poly fit
cgPlot, bin01[ 0 , sort( bin01[ 0 , * ] ) ], yfit01[ sort( bin01[ 0 , * ] ) ], color = 'blue', /overplot
;	err_yhigh = abs( bin01_residual[ sort( bin01[ 0 , * ] ) ] ), err_ylow = abs( bin01_residual[ sort( bin01[ 0 , * ] ) ] ), /err_clip
cgPlot, bin02[ 0 , sort( bin02[ 0 , * ] ) ], yfit02[ sort( bin02[ 0 , * ] ) ], color = 'forest green', /overplot
;	err_yhigh = abs( bin02_residual[ sort( bin02[ 0 , * ] ) ] ), err_ylow = abs( bin02_residual[ sort( bin02[ 0 , * ] ) ] ), /err_clip
cgPlot, bin03[ 0 , sort( bin03[ 0 , * ] ) ], yfit03[ sort( bin03[ 0 , * ] ) ], color = 'red', /overplot
;	err_yhigh = abs( bin03_residual[ sort( bin03[ 0 , * ] ) ] ), err_ylow = abs( bin03_residual[ sort( bin03[ 0 , * ] ) ] ), /err_clip
;plot the simulated data
;cgPlot, bin01[ 0 , sort( bin01[ 0 , * ] ) ], yfit01_sim[ sort( bin01[ 0 , * ] ) ], color = 'blue', linestyle = 2, /overplot, $
;	err_yhigh = yband_sim01[ sort( bin01[ 0 , * ] ) ], err_ylow = yband_sim01[ sort( bin01[ 0 , * ] ) ], thick = 1.5
;cgPlot, bin02[ 0 , sort( bin02[ 0 , * ] ) ], yfit02_sim[ sort( bin02[ 0 , * ] ) ], color = 'forest green', linestyle = 2, /overplot, $
;	err_yhigh = yband_sim02[ sort( bin02[ 0 , * ] ) ], err_ylow = yband_sim02[ sort( bin02[ 0 , * ] ) ], thick = 1.5
;cgPlot, bin03[ 0 , sort( bin03[ 0 , * ] ) ], yfit03_sim[ sort( bin03[ 0 , * ] ) ], color = 'red', linestyle = 2, /overplot, $
;	err_yhigh = yband_sim03[ sort( bin03[ 0 , * ] ) ], err_ylow = yband_sim03[ sort( bin03[ 0 , * ] ) ], thick = 1.5
;plot the simulated data
;cgPlot, bin01[ 0 , sort( bin01[ 0 , * ] ) ], sim_bin01_var[ sort( bin01[ 0 , * ] ) ], color = 'blue', psym = 1, /overplot
;cgPlot, bin02[ 0 , sort( bin02[ 0 , * ] ) ], sim_bin02_var[ sort( bin02[ 0 , * ] ) ], color = 'forest green', psym = 1, /overplot
;cgPlot, bin03[ 0 , sort( bin03[ 0 , * ] ) ], sim_bin03_var[ sort( bin03[ 0 , * ] ) ], color = 'red', psym = 1, /overplot
;legend
AL_LEGEND, ['2005-2006', '2007-2008', '2009-2010', '2005-2006 fit ', '2007-2008 fit ', '2009-2010 fit'], psym = [ 16 , 16 , 16 , -3 , -3 , -3], $
	color = ['blue', 'forest green', 'red', 'blue', 'forest green', 'red'], linestyle = [ 0 , 0 , 0 , 0 , 0 , 0 ], position = [ - 80, 2500 ], $
	/box, background_color = ['rose']
;=====================================================================================
	
;setting up the plot
cgPlot, bin04[ 2 , * ], xtitle = 'Latitudes (degrees)', ytitle = 'Ethane mixing ratio (pptv)', title = title, $
	yrange = [ 0 , 4000 ], xrange = [ -90, 90], xticklen = 1.0, yticklen = 1.0, xGridstyle = 1, yGridstyle = 1, /nodata
	;================================================================================================
cgPlot, bin04[ 0 , * ], bin04[ 2 , * ], psym = 16, color = 'blue', /overplot, err_yhigh = de_obs_bin04_noaa[ 2 , * ], $
	err_ylow = de_obs_bin04_noaa[ 2 , * ]
cgPlot, bin05[ 0 , * ], bin05[ 2 , * ], /overplot, psym = 16, color = 'forest green', err_yhigh = de_obs_bin05_noaa[ 2 , * ], $
	err_ylow = de_obs_bin05_noaa[ 2 , * ]
cgPlot, bin06[ 0 , * ], bin06[ 2 , * ], /overplot, psym = 16, color = 'red', err_yhigh = de_obs_bin06_noaa[ 2 , * ], $
	err_ylow = de_obs_bin06_noaa[ 2 , * ]
;poly_fit the shit out of the data
measure_errors = de_obs_bin04_noaa[ 2 , * ]
p_fit04 = poly_fit( bin04[ 0 , * ], bin04[ 2 , * ], 6, measure_errors = measure_errors, sigma = sigma04, yfit = yfit04, /double )
measure_errors = de_obs_bin05_noaa[ 2 , * ]
p_fit05 = poly_fit( bin05[ 0 , * ], bin05[ 2 , * ], 6, measure_errors = measure_errors, sigma = sigma05, yfit = yfit05, /double )
measure_errors = de_obs_bin06_noaa[ 2 , * ]
p_fit06 = poly_fit( bin06[ 0 , * ], bin06[ 2 , * ], 6, measure_errors = measure_errors, sigma = sigma06, yfit = yfit06, /double )
;poly_fit the simulated data
;p_fit04_sim = poly_fit( bin04[ 0 , * ], sim_bin04_var[ * ], 6, sigma = sigma04_sim, $
;	yfit = yfit04_sim, /double, yband = yband_sim04 )
;p_fit05_sim = poly_fit( bin05[ 0 , * ], sim_bin05_var[ * ], 6, sigma = sigma05_sim, $
;	yfit = yfit05_sim, /double, yband = yband_sim05 )
;p_fit06_sim = poly_fit( bin06[ 0 , * ], sim_bin06_var[ * ], 6, sigma = sigma06_sim, $
;	yfit = yfit06_sim, /double, yband = yband_sim06 )

;calculate the residual of the measured data and its poly fit to estimate the uncertainty
bin04_residual = yfit04 - bin04[ 2 , * ]	
bin05_residual = yfit05 - bin05[ 2 , * ]	
bin06_residual = yfit06 - bin06[ 2 , * ]	

;plot the poly fit
cgPlot, bin04[ 0 , sort( bin04[ 0 , * ] ) ], yfit04[ sort( bin04[ 0 , * ] ) ], color = 'blue', /overplot
;	err_yhigh = abs( bin04_residual[ sort( bin04[ 0 , * ] ) ] ), err_ylow = abs( bin04_residual[ sort( bin04[ 0 , * ] ) ] ), /err_clip
cgPlot, bin05[ 0 , sort( bin05[ 0 , * ] ) ], yfit05[ sort( bin05[ 0 , * ] ) ], color = 'forest green', /overplot
;	err_yhigh = abs( bin05_residual[ sort( bin05[ 0 , * ] ) ] ), err_ylow = abs( bin05_residual[ sort( bin05[ 0 , * ] ) ] ), /err_clip
cgPlot, bin06[ 0 , sort( bin06[ 0 , * ] ) ], yfit06[ sort( bin06[ 0 , * ] ) ], color = 'red', /overplot
;	err_yhigh = abs( bin06_residual[ sort( bin06[ 0 , * ] ) ] ), err_ylow = abs( bin06_residual[ sort( bin06[ 0 , * ] ) ] ), /err_clip
;plot the simulated data
;cgPlot, bin04[ 0 , sort( bin04[ 0 , * ] ) ], yfit04_sim[ sort( bin04[ 0 , * ] ) ], color = 'blue', linestyle = 2, /overplot, $
;	err_yhigh = yband_sim04[ sort( bin04[ 0 , * ] ) ], err_ylow = yband_sim04[ sort( bin04[ 0 , * ] ) ], thick = 1.5
;cgPlot, bin05[ 0 , sort( bin05[ 0 , * ] ) ], yfit05_sim[ sort( bin05[ 0 , * ] ) ], color = 'forest green', linestyle = 2, /overplot, $
;	err_yhigh = yband_sim05[ sort( bin05[ 0 , * ] ) ], err_ylow = yband_sim05[ sort( bin05[ 0 , * ] ) ], thick = 1.5
;cgPlot, bin06[ 0 , sort( bin06[ 0 , * ] ) ], yfit06_sim[ sort( bin06[ 0 , * ] ) ], color = 'red', linestyle = 2, /overplot, $
;	err_yhigh = yband_sim06[ sort( bin06[ 0 , * ] ) ], err_ylow = yband_sim06[ sort( bin06[ 0 , * ] ) ], thick = 1.5
;plot the simulated data
;cgPlot, bin04[ 0 , sort( bin04[ 0 , * ] ) ], sim_bin04_var[ sort( bin04[ 0 , * ] ) ], color = 'blue', psym = 1, /overplot
;cgPlot, bin05[ 0 , sort( bin05[ 0 , * ] ) ], sim_bin05_var[ sort( bin05[ 0 , * ] ) ], color = 'forest green', psym = 1, /overplot
;cgPlot, bin06[ 0 , sort( bin06[ 0 , * ] ) ], sim_bin06_var[ sort( bin06[ 0 , * ] ) ], color = 'red', psym = 1, /overplot
;legend
AL_LEGEND, ['2011-2012', '2013-2014', '2015', '2011-2012 fit ', '2013-2014 fit ', '2015 fit'], psym = [ 16 , 16 , 16 , -3 , -3 , -3], $
	color = ['blue', 'forest green', 'red', 'blue', 'forest green', 'red'], linestyle = [ 0 , 0 , 0 , 0 , 0 , 0 ], position = [ - 80, 2500 ], $
	/box, background_color = ['rose']
;=====================================================================================
	
	

	
	
;Make a JPEG for the graphs
READ, option, PROMPT= 'Do I need to save this graph as a jpeg so I can print it out later? Enter "1" for HELLYEAH or "0" for NOPE       '
name = strarr(1)
IF (option EQ 1) THEN BEGIN
	READ, name, PROMPT= 'Name of the JPEG file? (No need to include the jpeg extension, just name is fine)   '
	SCREEN2JPG, name
	print, 'JPEG has been saved to current directory'
ENDIF





end