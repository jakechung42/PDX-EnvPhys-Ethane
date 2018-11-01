PRO TEMPORAL_IHD_V4

;This v4 is built upon the resulted algorithm from the test version.
;It uses a new algorithm to do v1, v2 and v3.

compile_opt idl2

;=====Read in NOAA data=====
;get the files
infile1= "/home/jakechung/IDL/myIDL/NOAA_data.dat"
infile2= "/home/jakechung/IDL/myIDL/NOAA_locations.dat"

;get number of lines in the data file
n1_noaa = file_lines(infile1)
n2_noaa = file_lines(infile2)

print, 'The NOAA data has ', n1_noaa, '   data lines with', n2_noaa, '  locations'

;create arrays to store the data that was read in from the NOAA data file.
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



;=====Start reading in Simpson et al. data files=====
;get the files
infile1= "/home/jakechung/IDL/myIDL/data_set.dat"
infile2= "/home/jakechung/IDL/myIDL/location_set.dat"

;get number of lines in the data file
n1_ss = file_lines(infile1)
n2_ss = file_lines(infile2)

print, 'The input has ', n1_ss, '   data lines with', n2_ss, '  locations'

;create an empty array to store the data that was read in from the data files
;since the time data in the UCI data set is in a different format, so the time data when read into
;IDL will be handled differently compared to NOAA data.
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



;=====Start reading in OGI data files=====
;get the files
infile1= "/home/jakechung/IDL/myIDL/KhalilEtAl_data_set.dat"
infile2= "/home/jakechung/IDL/myIDL/KhalilEtAl_locations.dat"

;get number of lines in the data file
n1_ogi = file_lines(infile1)
n2_ogi = file_lines(infile2)

print, 'The input has ', n1_ogi, '   data lines with', n2_ogi, '  locations'

;create an empty array to store the data that was read in from the data files
;the time data is handled the same as UCI data.
unix_date_ogi = fltarr(n1_ogi)
avg_ogi = fltarr(n1_ogi)
location_ogi = strarr(n2_ogi)
lat1_ogi = fltarr(n1_ogi)
lon1_ogi = fltarr(n1_ogi)

;FLOAT input variable
location0= ' '
date0= 0.0
avg0= 0.0
lat0= 0.0
lon0=0.0

;open input files
openr, iunit1, infile1, /get_lun
openr, iunit2, infile2, /get_lun

;read float data of the data file
for i = 0, n1_ogi-1 do begin
	;read one line from the file
	readf, iunit1, date0, avg0, lat0, lon0
	;store the values
	unix_date_ogi[ i ] = date0
	avg_ogi[ i ] = avg0
	lat1_ogi[ i ] = lat0
	lon1_ogi[ i ] = lon0
endfor

;read float data of the location file
for i = 0, n2_ogi-1 do begin
	readf, iunit2, location0
	location_ogi[ i ] = location0
endfor

;close input file
free_lun, iunit1
free_lun, iunit2

;Because time in the OGI data set denoted time as Unix time, these lines convert Unix time to Geos_chem time and stored in tau_time array
tau_time_ogi = (unix_date_ogi/3600 - 131496)
tau_time_ogi = tau2yymmdd(tau_time_ogi, /GEOS1)



;--------------------------Calculating the NOAA data-----------------------
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
;bin01: 2006 + 2007
;bin02: 2008 + 2009
;bin03: 2010 + 2011
;bin04: 2012 + 2013
;bin05: 2014 

;setting up the arrays that will hold the averages
bin01_noaa = fltarr( 3 , noaa_ids_size )
bin02_noaa = fltarr( 3 , noaa_ids_size )
bin03_noaa = fltarr( 3 , noaa_ids_size )
bin04_noaa = fltarr( 3 , noaa_ids_size )
bin05_noaa = fltarr( 3 , noaa_ids_size )

;loop to find the station and extract the year.
;change the variable naming scheme to bin numbers instead of year numbers so if needed be, I can do a 2 - year or 5 - year mean
for i = 0, noaa_ids_size - 1 do begin
	a = where( strmatch ( location_noaa[ * ], noaa_ids[ i ], /FOLD_CASE) EQ 1 )
	temp_lat = lat1_noaa[ a ]
	temp_lon = lon1_noaa[ a ]
	temp_avg = avg_noaa[ a ]
	temp_bin01 = where( year_noaa[ a ] EQ 2006 OR year_noaa[ a ] EQ 2007 , count )
	;having a if-then filter here to make sure that when a station does not have a certain year data, 
	;it will not send the script into frenzy and writing false values
	if ( count GT 0 ) then begin
		bin01_noaa[ 2 , i ] = mean( temp_avg[ temp_bin01 ], /NAN )	
		bin01_noaa[ 0 , i ] = temp_lat[ 0 ]
		bin01_noaa[ 1 , i ] = temp_lon[ 0 ]
	endif 

	temp_bin02 = where( year_noaa[ a ] EQ 2008 OR year_noaa[ a ] EQ 2009 , count )
	if ( count GT 0 ) then begin
		bin02_noaa[ 2 , i ] = mean( temp_avg[ temp_bin02 ], /NAN )	
		bin02_noaa[ 0 , i ] = temp_lat[ 0 ]
		bin02_noaa[ 1 , i ] = temp_lon[ 0 ]
	endif 

	temp_bin03 = where( year_noaa[ a ] EQ 2010 OR year_noaa[ a ] EQ 2011 , count )
	if ( count GT 0 ) then begin
		bin03_noaa[ 2 , i ] = mean( temp_avg[ temp_bin03 ], /NAN )	
		bin03_noaa[ 0 , i ] = temp_lat[ 0 ]
		bin03_noaa[ 1 , i ] = temp_lon[ 0 ]
	endif 

	temp_bin04 = where( year_noaa[ a ] EQ 2012 OR year_noaa[ a ] EQ 2013 , count )
	if ( count GT 0 ) then begin
		bin04_noaa[ 2 , i ] = mean( temp_avg[ temp_bin04 ], /NAN )	
		bin04_noaa[ 0 , i ] = temp_lat[ 0 ]
		bin04_noaa[ 1 , i ] = temp_lon[ 0 ]
	endif

	temp_bin05 = where( year_noaa[ a ] EQ 2014 , count )
	if ( count GT 0 ) then begin
		bin05_noaa[ 2 , i ] = mean( temp_avg[ temp_bin05 ] )	
		bin05_noaa[ 0 , i ] = temp_lat[ 0 ]
		bin05_noaa[ 1 , i ] = temp_lon[ 0 ]
	endif
endfor 


;----------Calculations for the UCI data--------------

;array contains the lat values of the stations with more than 30 data points, there are 25 stations so 25 latitude values
lat_var_GT30 = [ 71.3000, 64.4900, 60.7500, 57.8300, 57.8000, 42.8300, 30.4000, 29.9400, 23.4400, 22.8700, 20.2600, 15.2000, $
	13.3600, 9.52000, 7.05000, 6.92000, 1.42000, -8.52000, -14.2300, -21.2300, -29.0200, -34.9000, -36.8200, -42.7200, -43.8800 ]
	
;bin01: 1996 - 2000
;bin02: 2001 - 2005
;bin03: 2006 - 2009
bin01_ss = fltarr( 3 , n_elements( lat_var_GT30 ) )
bin02_ss = fltarr( 3 , n_elements( lat_var_GT30 ) )
bin03_ss = fltarr( 3 , n_elements( lat_var_GT30 ) )

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
		bin01_ss[ 2 , i ] = mean( temp_avg[ temp_bin01 ], /NAN )	
		bin01_ss[ 0 , i ] = temp_lat[ 0 ]
		bin01_ss[ 1 , i ] = temp_lon[ 0 ]
		endif

	temp_bin02 = where( tau_time_ss.year[ a ] EQ 2001 OR $
		tau_time_ss.year[ a ] EQ 2002 OR $
		tau_time_ss.year[ a ] EQ 2003 OR $
		tau_time_ss.year[ a ] EQ 2004 OR $
		tau_time_ss.year[ a ] EQ 2005 , count )
	if ( count GT 3 ) then begin 
		bin02_ss[ 2 , i ] = mean( temp_avg[ temp_bin02 ], /NAN )	
		bin02_ss[ 0 , i ] = temp_lat[ 0 ]
		bin02_ss[ 1 , i ] = temp_lon[ 0 ]
		endif

	temp_bin03 = where( tau_time_ss.year[ a ] EQ 2006 OR $
		tau_time_ss.year[ a ] EQ 2007 OR $
		tau_time_ss.year[ a ] EQ 2008 OR $
		tau_time_ss.year[ a ] EQ 2009 , count )
	if ( count GT 3 ) then begin 
		bin03_ss[ 2 , i ] = mean( temp_avg[ temp_bin03 ], /NAN )	
		bin03_ss[ 0 , i ] = temp_lat[ 0 ]
		bin03_ss[ 1 , i ] = temp_lon[ 0 ]
		endif 

endfor 


;run calculation for the OGI data
ogi_lat_var = [ 71.16, 45.50, 21.08, -14.10, -42.0, -90 ]
bin01_ogi = fltarr( 3 , 6 )
for i = 0, 5 do begin
	a = where( lat1_ogi EQ ogi_lat_var[ i ] )
	bin01_ogi[ 1 , i ] = mean( avg_ogi[ a ] )
	bin01_ogi[ 2 , i ] = lon1_ogi[ a[ 0 ] ]
	bin01_ogi[ 0 , i ] = lat1_ogi[ a[ 0 ] ] 
endfor 




;This next section will read in the deseasonal values of the NOAA data from the deseason_all_noaa_v2 program ASCII file
;Create link to the data file
temp_deseason = "/home/jakechung/IDL/myIDL/temp_folder/temp_residual_noaa.dat"
openr, deseason, temp_deseason, /get_lun

;there should be 39 lines since the NOAA data has 39 stations
deseasonal_lines = file_lines(temp_deseason)
deseasonal_arr_noaa = fltarr(5, deseasonal_lines)
for u = 0, deseasonal_lines - 1 do begin
	readf, deseason, deseasonal_var0, lat0, lon0, year0, month0
	deseasonal_arr_noaa[ 0 , u ] = deseasonal_var0
	deseasonal_arr_noaa[ 1 , u ] = lat0
	deseasonal_arr_noaa[ 2 , u ] = lon0
	deseasonal_arr_noaa[ 3 , u ] = year0
	deseasonal_arr_noaa[ 4 , u ] = month0
endfor

free_lun, deseason

;the NOAA data and UCI data will be handled differently compared to the OGI data because the OGI data only has 1 bin.
;again, using the same method to seperate the years out from the deseasonal temp array from the deseason_all_noaa program
infile_bin01_de = "/home/jakechung/IDL/myIDL/temp_folder/temp_bin01_NOAA_deseason.dat"
infile_bin02_de = "/home/jakechung/IDL/myIDL/temp_folder/temp_bin02_NOAA_deseason.dat"
infile_bin03_de = "/home/jakechung/IDL/myIDL/temp_folder/temp_bin03_NOAA_deseason.dat"
infile_bin04_de = "/home/jakechung/IDL/myIDL/temp_folder/temp_bin04_NOAA_deseason.dat"
infile_bin05_de = "/home/jakechung/IDL/myIDL/temp_folder/temp_bin05_NOAA_deseason.dat"

;open to write in the temp files
openw , lun1 , infile_bin01_de , /get_lun
openw , lun2 , infile_bin02_de , /get_lun
openw , lun3 , infile_bin03_de , /get_lun
openw , lun4 , infile_bin04_de , /get_lun
openw , lun5 , infile_bin05_de , /get_lun

;loop to pull out the year 
for j = 0, deseasonal_lines - 1 do begin
	if (j NE deseasonal_lines - 1) then b = j + 1 else b = j
		if ( deseasonal_arr_noaa[ 1 , j ] EQ deseasonal_arr_noaa[ 1 , b ] ) then begin
			a = where( deseasonal_arr_noaa[ 1 , * ] EQ deseasonal_arr_noaa[ 1 , j ] )
			temp_de_year = deseasonal_arr_noaa[ 3 , a ]
			temp_de_month = deseasonal_arr_noaa[ 4 , a ]
			temp_de_var = deseasonal_arr_noaa[ 0 , a ]
			bin01_de_idx = where( temp_de_year EQ 2006 OR temp_de_year EQ 2007)
			bin02_de_idx = where( temp_de_year EQ 2008 OR temp_de_year EQ 2009)
			bin03_de_idx = where( temp_de_year EQ 2010 OR temp_de_year EQ 2011) 
			bin04_de_idx = where( temp_de_year EQ 2012 OR temp_de_year EQ 2013) 
			bin05_de_idx = where( temp_de_year EQ 2014 ) 
			;---------------------------------------------------------------------------------
			bin01_de_var = stddev( temp_de_var [ bin01_de_idx ] ) / sqrt( n_elements( bin01_de_idx ) )
			bin02_de_var = stddev( temp_de_var [ bin02_de_idx ] ) / sqrt( n_elements( bin02_de_idx ) )
			bin03_de_var = stddev( temp_de_var [ bin03_de_idx ] ) / sqrt( n_elements( bin03_de_idx ) )
			bin04_de_var = stddev( temp_de_var [ bin04_de_idx ] ) / sqrt( n_elements( bin04_de_idx ) )
			bin05_de_var = stddev( temp_de_var [ bin05_de_idx ] ) / sqrt( n_elements( bin05_de_idx ) )
			;----------------------------------------------------------------------------------
			printf, lun1, deseasonal_arr_noaa[ 1 , j ], deseasonal_arr_noaa[ 3 , j ], bin01_de_var
			printf, lun2, deseasonal_arr_noaa[ 1 , j ], deseasonal_arr_noaa[ 3 , j ], bin02_de_var
			printf, lun3, deseasonal_arr_noaa[ 1 , j ], deseasonal_arr_noaa[ 3 , j ], bin03_de_var
			printf, lun4, deseasonal_arr_noaa[ 1 , j ], deseasonal_arr_noaa[ 3 , j ], bin04_de_var
			printf, lun5, deseasonal_arr_noaa[ 1 , j ], deseasonal_arr_noaa[ 3 , j ], bin05_de_var
			;------------------------------------------------------------------------------------
			j = a[n_elements(a) - 1]					
			a = 0
		endif else if ( deseasonal_arr_noaa[ 1 , j ] NE deseasonal_arr_noaa[ 1 , b ] ) then j = j + 1
endfor

free_lun, lun1, lun2, lun3, lun4, lun5

;refresh the temp files to write the data into arrays for plotting
openr , lun1 , infile_bin01_de , /get_lun
openr , lun2 , infile_bin02_de , /get_lun
openr , lun3 , infile_bin03_de , /get_lun
openr , lun4 , infile_bin04_de , /get_lun
openr , lun5 , infile_bin05_de , /get_lun

bin01_de_size = file_lines( infile_bin01_de )
bin02_de_size = file_lines( infile_bin02_de )
bin03_de_size = file_lines( infile_bin03_de )
bin04_de_size = file_lines( infile_bin04_de )
bin05_de_size = file_lines( infile_bin05_de )

de_obs_bin01_noaa = fltarr( 3 , bin01_de_size)
de_obs_bin02_noaa = fltarr( 3 , bin02_de_size)
de_obs_bin03_noaa = fltarr( 3 , bin03_de_size)
de_obs_bin04_noaa = fltarr( 3 , bin04_de_size)
de_obs_bin05_noaa = fltarr( 3 , bin05_de_size)

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

free_lun, lun1, lun2, lun3, lun4, lun5

;remove 0 values in the deseasonal and the obs NOAA data 
de_obs_bin01_noaa = RemoveRows( de_obs_bin01_noaa, where( bin01_noaa[ 0 , * ] EQ 0 ) )
bin01_noaa = RemoveRows( bin01_noaa, where( bin01_noaa[ 0 , * ] EQ 0 ) )
de_obs_bin02_noaa = RemoveRows( de_obs_bin02_noaa, where( bin02_noaa[ 0 , * ] EQ 0 ) )
bin02_noaa = RemoveRows( bin02_noaa, where( bin02_noaa[ 0 , * ] EQ 0 ) )
de_obs_bin03_noaa = RemoveRows( de_obs_bin03_noaa, where( bin03_noaa[ 0 , * ] EQ 0 ) )
bin03_noaa = RemoveRows( bin03_noaa, where( bin03_noaa[ 0 , * ] EQ 0 ) ) 
de_obs_bin04_noaa = RemoveRows( de_obs_bin04_noaa, where( bin04_noaa[ 0 , * ] EQ 0 ) )
bin04_noaa = RemoveRows( bin04_noaa, where( bin04_noaa[ 0 , * ] EQ 0 ) )
de_obs_bin05_noaa = RemoveRows( de_obs_bin05_noaa, where( bin05_noaa[ 0 , * ] EQ 0 ) )
bin05_noaa = RemoveRows( bin05_noaa, where( bin05_noaa[ 0 , * ] EQ 0 ) )



;Dealing with  the deseasonal OGI data
;Create link to the data file
temp_deseason = "/home/jakechung/IDL/myIDL/temp_folder/temp_residual_OGI.dat"
openr, deseason, temp_deseason, /get_lun

deseasonal_lines = file_lines(temp_deseason)
deseasonal_arr_ogi = fltarr(5, deseasonal_lines)
for u = 0, deseasonal_lines - 1 do begin
	readf, deseason, deseasonal_var0, lat0, lon0, year0, month0
	deseasonal_arr_ogi[ 0 , u ] = deseasonal_var0
	deseasonal_arr_ogi[ 1 , u ] = lat0
	deseasonal_arr_ogi[ 2 , u ] = lon0
	deseasonal_arr_ogi[ 3 , u ] = year0
	deseasonal_arr_ogi[ 4 , u ] = month0
endfor

free_lun, deseason

;de_obs_bin01_ogi array contains the deseasonal values of each station of the OGI data
de_obs_bin01_ogi = fltarr( 2, 6 )
for i = 0, 5 do begin
	a = where( deseasonal_arr_ogi[ 1 , * ] EQ ogi_lat_var[ i ] )
	de_obs_bin01_ogi[ 0 , i ] = deseasonal_arr_ogi[ 1 , a[ 0 ] ] 
	de_obs_bin01_ogi[ 1 , i ] = stddev( deseasonal_arr_ogi[ 0 , a ] ) / sqrt( n_elements( a ) )
endfor 


;this section deal with the UCI deseasoned data
;Create link to the data file
temp_deseason = "/home/jakechung/IDL/myIDL/temp_folder/temp_residual_Simpson.dat"
openr, deseason, temp_deseason, /get_lun

deseasonal_lines = file_lines(temp_deseason)
deseasonal_arr_ss = fltarr(5, deseasonal_lines)
for u = 0, deseasonal_lines - 1 do begin
	readf, deseason, deseasonal_var0, lat0, lon0, year0, month0
	deseasonal_arr_ss[ 0 , u ] = deseasonal_var0
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

;loop to pull out the years
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

;refresh the temp files to write the data into arrays
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

de_obs_bin01_ss = RemoveRows( de_obs_bin01_ss, where( bin01_ss[ 0 , * ] EQ 0 ) )
bin01_ss = RemoveRows( bin01_ss, where( bin01_ss[ 0 , * ] EQ 0 ) )
de_obs_bin02_ss = RemoveRows( de_obs_bin02_ss, where( bin02_ss[ 0 , * ] EQ 0 ) )
bin02_ss = RemoveRows( bin02_ss, where( bin02_ss[ 0 , * ] EQ 0 ) )
de_obs_bin03_ss = RemoveRows( de_obs_bin03_ss, where( bin03_ss[ 0 , * ] EQ 0 ) )
bin03_ss = RemoveRows( bin03_ss, where( bin03_ss[ 0 , * ] EQ 0 ) ) 
;preliminary data processing for all data before calculations
;sort the NOAA data by latitudes
a = sort( bin01_noaa[ 0 , * ] )
bin01_noaa[ 0 , * ] = bin01_noaa[ 0 , a ]
bin01_noaa[ 1 , * ] = bin01_noaa[ 1 , a ]
bin01_noaa[ 2 , * ] = bin01_noaa[ 2 , a ]
de_obs_bin01_noaa[ 0 , * ] = de_obs_bin01_noaa[ 0 , a ]
de_obs_bin01_noaa[ 2 , * ] = de_obs_bin01_noaa[ 2 , a ]

a = sort( bin02_noaa[ 0 , * ] )
bin02_noaa[ 0 , * ] = bin02_noaa[ 0 , a ]
bin02_noaa[ 1 , * ] = bin02_noaa[ 1 , a ]
bin02_noaa[ 2 , * ] = bin02_noaa[ 2 , a ]
de_obs_bin02_noaa[ 0 , * ] = de_obs_bin02_noaa[ 0 , a ]
de_obs_bin02_noaa[ 2 , * ] = de_obs_bin02_noaa[ 2 , a ]

a = sort( bin03_noaa[ 0 , * ] )
bin03_noaa[ 0 , * ] = bin03_noaa[ 0 , a ]
bin03_noaa[ 1 , * ] = bin03_noaa[ 1 , a ]
bin03_noaa[ 2 , * ] = bin03_noaa[ 2 , a ]
de_obs_bin03_noaa[ 0 , * ] = de_obs_bin03_noaa[ 0 , a ]
de_obs_bin03_noaa[ 2 , * ] = de_obs_bin03_noaa[ 2 , a ]

a = sort( bin04_noaa[ 0 , * ] )
bin04_noaa[ 0 , * ] = bin04_noaa[ 0 , a ]
bin04_noaa[ 1 , * ] = bin04_noaa[ 1 , a ]
bin04_noaa[ 2 , * ] = bin04_noaa[ 2 , a ]
de_obs_bin04_noaa[ 0 , * ] = de_obs_bin04_noaa[ 0 , a ]
de_obs_bin04_noaa[ 2 , * ] = de_obs_bin04_noaa[ 2 , a ]

a = sort( bin05_noaa[ 0 , * ] )
bin05_noaa[ 0 , * ] = bin05_noaa[ 0 , a ]
bin05_noaa[ 1 , * ] = bin05_noaa[ 1 , a ]
bin05_noaa[ 2 , * ] = bin05_noaa[ 2 , a ]
de_obs_bin05_noaa[ 0 , * ] = de_obs_bin05_noaa[ 0 , a ]
de_obs_bin05_noaa[ 2 , * ] = de_obs_bin05_noaa[ 2 , a ]
;sort the OGI and UCI data
a = sort( bin01_ss[ 0 , * ] )
bin01_ss[ 0 , * ] = bin01_ss[ 0 , a ]
bin01_ss[ 1 , * ] = bin01_ss[ 1 , a ]
bin01_ss[ 2 , * ] = bin01_ss[ 2 , a ]
de_obs_bin01_ss[ 0 , * ] = de_obs_bin01_ss[ 0 , a ]
de_obs_bin01_ss[ 1 , * ] = de_obs_bin01_ss[ 1 , a ]

a = sort( bin02_ss[ 0 , * ] )
bin02_ss[ 0 , * ] = bin02_ss[ 0 , a ]
bin02_ss[ 1 , * ] = bin02_ss[ 1 , a ]
bin02_ss[ 2 , * ] = bin02_ss[ 2 , a ]
de_obs_bin02_ss[ 0 , * ] = de_obs_bin02_ss[ 0 , a ]
de_obs_bin02_ss[ 1 , * ] = de_obs_bin02_ss[ 1 , a ]

a = sort( bin03_ss[ 0 , * ] )
bin03_ss[ 0 , * ] = bin03_ss[ 0 , a ]
bin03_ss[ 1 , * ] = bin03_ss[ 1 , a ]
bin03_ss[ 2 , * ] = bin03_ss[ 2 , a ]
de_obs_bin03_ss[ 0 , * ] = de_obs_bin03_ss[ 0 , a ]
de_obs_bin03_ss[ 1 , * ] = de_obs_bin03_ss[ 1 , a ]

;sort the OGI data
a = sort( bin01_ogi[ 0 , * ] )
bin01_ogi[ 0 , * ] = bin01_ogi[ 0 , a ]
bin01_ogi[ 1 , * ] = bin01_ogi[ 1 , a ]
bin01_ogi[ 2 , * ] = bin01_ogi[ 2 , a ]
de_obs_bin01_ogi[ 0 , * ] = de_obs_bin01_ogi[ 0 , a ]
de_obs_bin01_ogi[ 1 , * ] = de_obs_bin01_ogi[ 1 , a ]
;-----------This new section using a new algorithm to call the sim data-----------------

;sim_title = 'PSU emissions scaled to Xiao et al over 1996-2003'
;filename = "/home/excluded-from-backup/data/C2H6/trac_avg.PSUSF_1981_2015.bpch"

;sim_title = 'Unscaled emissions from Simpson et al.'
;filename = "/home/excluded-from-backup/data/C2H6/trac_avg.SimpsonSF_1981_2015.bpch"

;sim_title = 'Constant default base emissions' 
;filename = "/home/excluded-from-backup/data/C2H6/trac_avg.spinup_GFED4_MAVG_1981_2015.198101010000"

;sim_title = 'Aydin et al. no normalization to Xiao et al.'
;filename = "/home/excluded-from-backup/data/C2H6/trac_avg.AydinAbsSF_1981_2015.bpch"

;sim_title = 'Simpson et al. no normalization to Xiao et al.'
;filename = "/home/excluded-from-backup/data/C2H6/trac_avg.SimpsonSF_1981_2015.bpch"

sim_title = 'Unscaled PSU emission, MER3BB, MER18FF'
filename = "/home/excluded-from-backup/data/C2H6/trac_avg.PSU_MER3BB_MER18FF_1981_2015.bpch"

;read in the file
ctm_get_data, datainfo, 'IJ-AVG-$', tracer=1, filename=filename, tau0=tau0

;Get MODELINFO and GRIDINFO structures, xmid, ymid hold lon/lat centers
getmodelandgridinfo, datainfo[0], modelinfo, gridinfo

;nt is the number of data points in the simulation.
nt = n_elements(datainfo)
globalavg= fltarr(nt)

;sim_yymmdd stores the entire time element of the input simulation
sim_yymmdd= tau2yymmdd(tau0, /GEOS1)

;locate indexes of all months of the simulated data for noaa
sim_bin01_idx_noaa = where( sim_yymmdd.year EQ 2006 OR sim_yymmdd.year EQ 2007 )
sim_bin02_idx_noaa = where( sim_yymmdd.year EQ 2008 OR sim_yymmdd.year EQ 2009 )
sim_bin03_idx_noaa = where( sim_yymmdd.year EQ 2010 OR sim_yymmdd.year EQ 2011 )
sim_bin04_idx_noaa = where( sim_yymmdd.year EQ 2012 OR sim_yymmdd.year EQ 2013 )
sim_bin05_idx_noaa = where( sim_yymmdd.year EQ 2014 )

sim_bin01_idx_ss = where( sim_yymmdd.year EQ 1996 OR sim_yymmdd.year EQ 1997	$
	OR sim_yymmdd.year EQ 1998 OR sim_yymmdd.year EQ 1999 OR sim_yymmdd.year EQ 2000 )
sim_bin02_idx_ss = where( sim_yymmdd.year EQ 2001 OR sim_yymmdd.year EQ 2002	$
	OR sim_yymmdd.year EQ 2003 OR sim_yymmdd.year EQ 2004 OR sim_yymmdd.year EQ 2005 )
sim_bin03_idx_ss = where( sim_yymmdd.year EQ 2006 OR sim_yymmdd.year EQ 2007	$
	OR sim_yymmdd.year EQ 2008 OR sim_yymmdd.year EQ 2009 )

sim_bin01_idx_ogi = where( sim_yymmdd.year EQ 1983 $
	OR sim_yymmdd.year EQ 1984 $
	OR sim_yymmdd.year EQ 1985 $
	OR sim_yymmdd.year EQ 1986 $
	OR sim_yymmdd.year EQ 1987 )

;this section call out the data blocks of the simulated data and store it in simArr
;simArr is a 3-D array with the following attributes:
;simArr[ *, 0, 0 ]: the time dimension
;simArr[ 0, *, 0 ]: longitudes
;simArr[ 0, 0, * ]: latitudes
simArr = fltarr( nt, 144, 91 )
help, simArr
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

;this changes the lat and lon values to index values to call the data
;out from the sim array 
lat_idx_noaa_bin01 = round( ( bin01_noaa[ 0 , * ] - (-90) ) / 2 )
lat_idx_noaa_bin02 = round( ( bin02_noaa[ 0 , * ] - (-90) ) / 2 )
lat_idx_noaa_bin03 = round( ( bin03_noaa[ 0 , * ] - (-90) ) / 2 )
lat_idx_noaa_bin04 = round( ( bin04_noaa[ 0 , * ] - (-90) ) / 2 )
lat_idx_noaa_bin05 = round( ( bin05_noaa[ 0 , * ] - (-90) ) / 2 )

lon_idx_noaa_bin01 = round( ( bin01_noaa[ 1 , * ] - (-180) ) / 2.5 )
lon_idx_noaa_bin02 = round( ( bin02_noaa[ 1 , * ] - (-180) ) / 2.5 )
lon_idx_noaa_bin03 = round( ( bin03_noaa[ 1 , * ] - (-180) ) / 2.5 )
lon_idx_noaa_bin04 = round( ( bin04_noaa[ 1 , * ] - (-180) ) / 2.5 )
lon_idx_noaa_bin05 = round( ( bin05_noaa[ 1 , * ] - (-180) ) / 2.5 )

lat_idx_ss_bin01 = round( ( bin01_ss[ 0 , * ] - (-90) ) / 2 )
lat_idx_ss_bin02 = round( ( bin02_ss[ 0 , * ] - (-90) ) / 2 )
lat_idx_ss_bin03 = round( ( bin03_ss[ 0 , * ] - (-90) ) / 2 )

lon_idx_ss_bin01 = round( ( bin01_ss[ 1 , * ] - (-180) ) / 2.5 )
lon_idx_ss_bin02 = round( ( bin02_ss[ 1 , * ] - (-180) ) / 2.5 )
lon_idx_ss_bin03 = round( ( bin03_ss[ 1 , * ] - (-180) ) / 2.5 )

lat_idx_ogi_bin01 = round( ( bin01_ogi[ 0 , * ] - (-90) ) / 2 )

lon_idx_ogi_bin01 = round( ( bin01_ogi[ 1 , * ] - (-180) ) / 2.5 )

;call the sim data out from the simArr, based on coordinates and time
;First, call out the correct time period.
bin01_sim_noaa = mean( simArr[ sim_bin01_idx_noaa , * , * ] , 1 )
bin02_sim_noaa = mean( simArr[ sim_bin02_idx_noaa , * , * ] , 1 )
bin03_sim_noaa = mean( simArr[ sim_bin03_idx_noaa , * , * ] , 1 )
bin04_sim_noaa = mean( simArr[ sim_bin04_idx_noaa , * , * ] , 1 )
bin05_sim_noaa = mean( simArr[ sim_bin05_idx_noaa , * , * ] , 1 )

bin01_sim_ss = mean( simArr[ sim_bin01_idx_ss , * , * ] , 1 )
bin02_sim_ss = mean( simArr[ sim_bin02_idx_ss , * , * ] , 1 )
bin03_sim_ss = mean( simArr[ sim_bin03_idx_ss , * , * ] , 1 )

bin01_sim_ogi = mean( simArr[ sim_bin01_idx_ogi , * , * ] , 1 )

;create some arrays to store the mixing ratio values from the simulated data
bin01_sim_noaa_var = fltarr( 2 , n_elements( lat_idx_noaa_bin01 ) )
bin02_sim_noaa_var = fltarr( 2 , n_elements( lat_idx_noaa_bin02 ) )
bin03_sim_noaa_var = fltarr( 2 , n_elements( lat_idx_noaa_bin03 ) )
bin04_sim_noaa_var = fltarr( 2 , n_elements( lat_idx_noaa_bin04 ) )
bin05_sim_noaa_var = fltarr( 2 , n_elements( lat_idx_noaa_bin05 ) )

bin01_sim_ss_var = fltarr( 2 , n_elements( lat_idx_ss_bin01 ) )
bin02_sim_ss_var = fltarr( 2 , n_elements( lat_idx_ss_bin02 ) )
bin03_sim_ss_var = fltarr( 2 , n_elements( lat_idx_ss_bin03 ) )

bin01_sim_ogi_var = fltarr( 2 , n_elements( lat_idx_ogi_bin01 ) )

;Second, call out the correct coordinate of the simulated data
;bin01_sim_noaa_var[ 0 , * ] stores the latitudes 
;bin01_sim_noaa_var[ 1 , * ] stores the mixing ratio values
bin01_sim_noaa_var[ 0 , * ] = bin01_noaa[ 0 , * ]
bin01_sim_noaa_var[ 1 , * ] = bin01_sim_noaa[ lon_idx_noaa_bin01 , lat_idx_noaa_bin01 ] / 2 * 1000
bin02_sim_noaa_var[ 0 , * ] = bin02_noaa[ 0 , * ]
bin02_sim_noaa_var[ 1 , * ] = bin02_sim_noaa[ lon_idx_noaa_bin02 , lat_idx_noaa_bin02 ] / 2 * 1000
bin03_sim_noaa_var[ 0 , * ] = bin03_noaa[ 0 , * ]
bin03_sim_noaa_var[ 1 , * ] = bin03_sim_noaa[ lon_idx_noaa_bin03 , lat_idx_noaa_bin03 ] / 2 * 1000
bin04_sim_noaa_var[ 0 , * ] = bin04_noaa[ 0 , * ]
bin04_sim_noaa_var[ 1 , * ] = bin04_sim_noaa[ lon_idx_noaa_bin04 , lat_idx_noaa_bin04 ] / 2 * 1000
bin05_sim_noaa_var[ 0 , * ] = bin05_noaa[ 0 , * ]
bin05_sim_noaa_var[ 1 , * ] = bin05_sim_noaa[ lon_idx_noaa_bin05 , lat_idx_noaa_bin05 ] / 2 * 1000

bin01_sim_ss_var[ 0 , * ] = bin01_ss[ 0 , * ]
bin01_sim_ss_var[ 1 , * ] = bin01_sim_ss[ lon_idx_ss_bin01 , lat_idx_ss_bin01 ] / 2 * 1000
bin02_sim_ss_var[ 0 , * ] = bin02_ss[ 0 , * ]
bin02_sim_ss_var[ 1 , * ] = bin02_sim_ss[ lon_idx_ss_bin02 , lat_idx_ss_bin02 ] / 2 * 1000
bin03_sim_ss_var[ 0 , * ] = bin03_ss[ 0 , * ]
bin03_sim_ss_var[ 1 , * ] = bin03_sim_ss[ lon_idx_ss_bin03 , lat_idx_ss_bin03 ] / 2 * 1000

bin01_sim_ogi_var[ 0 , * ] = bin01_ogi[ 0 , * ]
bin01_sim_ogi_var[ 1 , * ] = bin01_sim_ogi[ lon_idx_ogi_bin01 , lat_idx_ogi_bin01 ] / 2 * 1000

;convert the latitude values of the obs to sine values in order to do the poly fit, same thing for 
;sim data
bin01_noaa[ 0 , * ] = sin( bin01_noaa[ 0 , * ] * !pi / 180 )
bin02_noaa[ 0 , * ] = sin( bin02_noaa[ 0 , * ] * !pi / 180 )
bin03_noaa[ 0 , * ] = sin( bin03_noaa[ 0 , * ] * !pi / 180 )
bin04_noaa[ 0 , * ] = sin( bin04_noaa[ 0 , * ] * !pi / 180 )
bin05_noaa[ 0 , * ] = sin( bin05_noaa[ 0 , * ] * !pi / 180 )
bin01_sim_noaa_var[ 0 , * ] = sin( bin01_sim_noaa_var[ 0 , * ] * !pi / 180 )
bin02_sim_noaa_var[ 0 , * ] = sin( bin02_sim_noaa_var[ 0 , * ] * !pi / 180 )
bin03_sim_noaa_var[ 0 , * ] = sin( bin03_sim_noaa_var[ 0 , * ] * !pi / 180 )
bin04_sim_noaa_var[ 0 , * ] = sin( bin04_sim_noaa_var[ 0 , * ] * !pi / 180 )
bin05_sim_noaa_var[ 0 , * ] = sin( bin05_sim_noaa_var[ 0 , * ] * !pi / 180 )

bin01_ss[ 0 , * ] = sin( bin01_ss[ 0 , * ] * !pi / 180 )
bin02_ss[ 0 , * ] = sin( bin02_ss[ 0 , * ] * !pi / 180 )
bin03_ss[ 0 , * ] = sin( bin03_ss[ 0 , * ] * !pi / 180 )
bin01_sim_ss_var[ 0 , * ] = sin( bin01_sim_ss_var[ 0 , * ] * !pi / 180 )
bin02_sim_ss_var[ 0 , * ] = sin( bin02_sim_ss_var[ 0 , * ] * !pi / 180 )
bin03_sim_ss_var[ 0 , * ] = sin( bin03_sim_ss_var[ 0 , * ] * !pi / 180 )

bin01_ogi[ 0 , * ] = sin( bin01_ogi[ 0 , * ] * !pi / 180 )
bin01_sim_ogi_var[ 0 , * ] = sin( bin01_sim_ogi_var[ 0 , * ] * !pi / 180 )

;poly fit the obs data using the deseasonal value to weight the fit.
measure_errors = de_obs_bin01_noaa[ 2 , sort( bin01_noaa[ 0 , * ] ) ]
bin01_noaa_pfit = poly_fit( bin01_noaa [ 0 , * ] , $
	bin01_noaa [ 2 , * ], 4, measure_errors = measure_errors, $
	/double, yfit = yfit_bin01_noaa )
measure_errors = de_obs_bin02_noaa[ 2 , sort( bin02_noaa[ 0 , * ] ) ]
bin02_noaa_pfit =  poly_fit( bin02_noaa [ 0 , * ] , $
	bin02_noaa [ 2 , * ], 4, measure_errors = measure_errors, $
	/double, yfit = yfit_bin02_noaa )
measure_errors = de_obs_bin03_noaa[ 2 , sort( bin03_noaa[ 0 , * ] ) ]
bin03_noaa_pfit =  poly_fit( bin03_noaa [ 0 , * ] , $
	bin03_noaa [ 2 , * ], 4, measure_errors = measure_errors, $
	/double, yfit = yfit_bin03_noaa )
measure_errors = de_obs_bin04_noaa[ 2 , sort( bin04_noaa[ 0 , * ] ) ]
bin04_noaa_pfit =  poly_fit( bin04_noaa [ 0 , * ] , $
	bin04_noaa [ 2 , * ], 4, measure_errors = measure_errors, $
	/double, yfit = yfit_bin04_noaa )
measure_errors = de_obs_bin05_noaa[ 2 , sort( bin05_noaa[ 0 , * ] ) ]
bin05_noaa_pfit =  poly_fit( bin05_noaa [ 0 , * ] , $
	bin05_noaa [ 2 , * ], 4, measure_errors = measure_errors, $
	/double, yfit = yfit_bin05_noaa )

measure_errors = de_obs_bin01_ss[ 1 , * ]
bin01_ss_pfit = poly_fit( bin01_ss [ 0 , * ] , $
	bin01_ss [ 2 , * ], 4, measure_errors = measure_errors, $
	/double, yfit = yfit_bin01_ss )
measure_errors = de_obs_bin02_ss[ 1 , * ]
bin02_ss_pfit =  poly_fit( bin02_ss [ 0 , * ] , $
	bin02_ss [ 2 , * ], 4, measure_errors = measure_errors, $
	/double, yfit = yfit_bin02_ss )
measure_errors = de_obs_bin03_ss[ 1 , * ]
bin03_ss_pfit =  poly_fit( bin03_ss [ 0 , * ] , $
	bin03_ss [ 2 , * ], 4, measure_errors = measure_errors, $
	/double, yfit = yfit_bin03_ss )
	
measure_errors = de_obs_bin01_ogi[ 1 , * ]
bin01_ogi_pfit = poly_fit( bin01_ogi [ 0 , * ] , $
	bin01_ogi [ 1 , * ], 4, measure_errors = measure_errors, $
	/double, yfit = yfit_bin01_ogi )
;poly fit the sim data.
bin01_sim_noaa_pfit = poly_fit( bin01_sim_noaa_var[ 0 , * ], $
	bin01_sim_noaa_var[ 1 , * ], 4, yfit = yfit_bin01_sim_noaa, /double )
bin02_sim_noaa_pfit = poly_fit( bin02_sim_noaa_var[ 0 , * ], $
	bin02_sim_noaa_var[ 1 , * ], 4, yfit = yfit_bin02_sim_noaa, /double )
bin03_sim_noaa_pfit = poly_fit( bin03_sim_noaa_var[ 0 , * ], $
	bin03_sim_noaa_var[ 1 , * ], 4, yfit = yfit_bin03_sim_noaa, /double )
bin04_sim_noaa_pfit = poly_fit( bin04_sim_noaa_var[ 0 , * ], $
	bin04_sim_noaa_var[ 1 , * ], 4, yfit = yfit_bin04_sim_noaa, /double )
bin05_sim_noaa_pfit = poly_fit( bin05_sim_noaa_var[ 0 , * ], $
	bin05_sim_noaa_var[ 1 , * ], 4, yfit = yfit_bin05_sim_noaa, /double )
	
bin01_sim_ss_pfit = poly_fit( bin01_sim_ss_var[ 0 , * ], $
	bin01_sim_ss_var[ 1 , * ], 4, yfit = yfit_bin01_sim_ss, /double )
bin02_sim_ss_pfit = poly_fit( bin02_sim_ss_var[ 0 , * ], $
	bin02_sim_ss_var[ 1 , * ], 4, yfit = yfit_bin02_sim_ss, /double )
bin03_sim_ss_pfit = poly_fit( bin03_sim_ss_var[ 0 , * ], $
	bin03_sim_ss_var[ 1 , * ], 4, yfit = yfit_bin03_sim_ss, /double )
	
bin01_sim_ogi_pfit = poly_fit( bin01_sim_ogi_var[ 0 , * ], $
	bin01_sim_ogi_var[ 1 , * ], 4, yfit = yfit_bin01_sim_ogi, /double )

;calculating the area under the curve of the fit for obs data using int_tabulated routines
;the where function is to narrow down to the latitudinal restriction when calculating 
;the area under the curve
rad = !pi / 180
nor_noaa = fltarr( 5 )
;a = where( bin01_noaa[ 0 , * ] gt 0 and bin01_noaa[ 0 , * ] lt sin( 70*rad ) )
;nor_noaa[ 0 ] = int_tabulated( bin01_noaa[ 0 , a ], yfit_bin01_noaa[ a ], /double )
;a = where( bin02_noaa[ 0 , * ] gt 0 and bin02_noaa[ 0 , * ] lt sin( 70*rad ) )
;nor_noaa[ 1 ] = int_tabulated( bin02_noaa[ 0 , a ], yfit_bin02_noaa[ a ], /double )
;a = where( bin03_noaa[ 0 , * ] gt 0 and bin03_noaa[ 0 , * ] lt sin( 70*rad ) )
;nor_noaa[ 2 ] = int_tabulated( bin03_noaa[ 0 , a ], yfit_bin03_noaa[ a ], /double )
;a = where( bin04_noaa[ 0 , * ] gt 0 and bin04_noaa[ 0 , * ] lt sin( 70*rad ) )
;nor_noaa[ 3 ] = int_tabulated( bin04_noaa[ 0 , a ], yfit_bin04_noaa[ a ], /double )
;a = where( bin05_noaa[ 0 , * ] gt 0 and bin05_noaa[ 0 , * ] lt sin( 70*rad ) )
;nor_noaa[ 4 ] = int_tabulated( bin05_noaa[ 0 , a ], yfit_bin05_noaa[ a ], /double )
;change the methond of calculating the area under the curve to doing a algebraic 
;analysis rather than using the int__tabulated algorithm
nor_noaa[ 0 ] = def_integrate_poly4( bin01_noaa_pfit[0], bin01_noaa_pfit[1], bin01_noaa_pfit[2], $
	bin01_noaa_pfit[3], bin01_noaa_pfit[4], sin( 0 * rad ), sin( 70 * rad ))
nor_noaa[ 1 ] = def_integrate_poly4( bin02_noaa_pfit[0], bin02_noaa_pfit[1], bin02_noaa_pfit[2], $
	bin02_noaa_pfit[3], bin02_noaa_pfit[4], sin( 0 * rad ), sin( 70 * rad ))
nor_noaa[ 2 ] = def_integrate_poly4( bin03_noaa_pfit[0], bin03_noaa_pfit[1], bin03_noaa_pfit[2], $
	bin03_noaa_pfit[3], bin03_noaa_pfit[4], sin( 0 * rad ), sin( 70 * rad ))
nor_noaa[ 3 ] = def_integrate_poly4( bin04_noaa_pfit[0], bin04_noaa_pfit[1], bin04_noaa_pfit[2], $
	bin04_noaa_pfit[3], bin04_noaa_pfit[4], sin( 0 * rad ), sin( 70 * rad ))
nor_noaa[ 4 ] = def_integrate_poly4( bin05_noaa_pfit[0], bin05_noaa_pfit[1], bin05_noaa_pfit[2], $
	bin05_noaa_pfit[3], bin05_noaa_pfit[4], sin( 0 * rad ), sin( 70 * rad ))

sou_noaa = fltarr( 5 )
;a = where( bin01_noaa[ 0 , * ] lt 0 and bin01_noaa[ 0 , * ] gt sin( -30*rad ) )
;sou_noaa[ 0 ] = int_tabulated( bin01_noaa[ 0 , a ], yfit_bin01_noaa[ a ], /double )
;a = where( bin02_noaa[ 0 , * ] lt 0 and bin02_noaa[ 0 , * ] gt sin( -30*rad ) )
;sou_noaa[ 1 ] = int_tabulated( bin02_noaa[ 0 , a ], yfit_bin02_noaa[ a ], /double )
;a = where( bin03_noaa[ 0 , * ] lt 0 and bin03_noaa[ 0 , * ] gt sin( -30*rad ) )
;sou_noaa[ 2 ] = int_tabulated( bin03_noaa[ 0 , a ], yfit_bin03_noaa[ a ], /double )
;a = where( bin04_noaa[ 0 , * ] lt 0 and bin04_noaa[ 0 , * ] gt sin( -30*rad ) )
;sou_noaa[ 3 ] = int_tabulated( bin04_noaa[ 0 , a ], yfit_bin04_noaa[ a ], /double )
;a = where( bin05_noaa[ 0 , * ] lt 0 and bin05_noaa[ 0 , * ] gt sin( -30*rad ) )
;sou_noaa[ 4 ] = int_tabulated( bin05_noaa[ 0 , a ], yfit_bin05_noaa[ a ], /double )
sou_noaa[ 0 ] = def_integrate_poly4( bin01_noaa_pfit[0], bin01_noaa_pfit[1], bin01_noaa_pfit[2], $
	bin01_noaa_pfit[3], bin01_noaa_pfit[4], sin( (-30) * rad ), sin( 0 * rad ))
sou_noaa[ 1 ] = def_integrate_poly4( bin02_noaa_pfit[0], bin02_noaa_pfit[1], bin02_noaa_pfit[2], $
	bin02_noaa_pfit[3], bin02_noaa_pfit[4], sin( (-30) * rad ), sin( 0 * rad ))
sou_noaa[ 2 ] = def_integrate_poly4( bin03_noaa_pfit[0], bin03_noaa_pfit[1], bin03_noaa_pfit[2], $
	bin03_noaa_pfit[3], bin03_noaa_pfit[4], sin( (-30) * rad ), sin( 0 * rad ))
sou_noaa[ 3 ] = def_integrate_poly4( bin04_noaa_pfit[0], bin04_noaa_pfit[1], bin04_noaa_pfit[2], $
	bin04_noaa_pfit[3], bin04_noaa_pfit[4], sin( (-30) * rad ), sin( 0 * rad ))
sou_noaa[ 4 ] = def_integrate_poly4( bin05_noaa_pfit[0], bin05_noaa_pfit[1], bin05_noaa_pfit[2], $
	bin05_noaa_pfit[3], bin05_noaa_pfit[4], sin( (-30) * rad ), sin( 0 * rad ))
	
nor_ss = fltarr( 3 )
;a = where( bin01_ss[ 0 , * ] gt 0 and bin01_ss[ 0 , * ] lt sin( 70*rad ) )
;nor_ss[ 0 ] = int_tabulated( bin01_ss[ 0 , a ], yfit_bin01_ss[ a ], /double )
;a = where( bin02_ss[ 0 , * ] gt 0 and bin02_ss[ 0 , * ] lt sin( 70*rad ) )
;nor_ss[ 1 ] = int_tabulated( bin02_ss[ 0 , a ], yfit_bin02_ss[ a ], /double )
;a = where( bin03_ss[ 0 , * ] gt 0 and bin03_ss[ 0 , * ] lt sin( 70*rad ) )
;nor_ss[ 2 ] = int_tabulated( bin03_ss[ 0 , a ], yfit_bin03_ss[ a ], /double )
nor_ss[ 0 ] = def_integrate_poly4( bin01_ss_pfit[0], bin01_ss_pfit[1], bin01_ss_pfit[2], $
	bin01_ss_pfit[3], bin01_ss_pfit[4], sin( 0 * rad ), sin( 70 * rad ))
nor_ss[ 1 ] = def_integrate_poly4( bin02_ss_pfit[0], bin02_ss_pfit[1], bin02_ss_pfit[2], $
	bin02_ss_pfit[3], bin02_ss_pfit[4], sin( 0 * rad ), sin( 70 * rad ))
nor_ss[ 2 ] = def_integrate_poly4( bin03_ss_pfit[0], bin03_ss_pfit[1], bin03_ss_pfit[2], $
	bin03_ss_pfit[3], bin03_ss_pfit[4], sin( 0 * rad ), sin( 70 * rad ))

sou_ss = fltarr( 3 )
;a = where( bin01_ss[ 0 , * ] lt 0 and bin01_ss[ 0 , * ] gt sin( -30*rad ) )
;sou_ss[ 0 ] = int_tabulated( bin01_ss[ 0 , a ], yfit_bin01_ss[ a ], /double )
;a = where( bin02_ss[ 0 , * ] lt 0 and bin02_ss[ 0 , * ] gt sin( -30*rad ) )
;sou_ss[ 1 ] = int_tabulated( bin02_ss[ 0 , a ], yfit_bin02_ss[ a ], /double )
;a = where( bin03_ss[ 0 , * ] lt 0 and bin03_ss[ 0 , * ] gt sin( -30*rad ) )
;sou_ss[ 2 ] = int_tabulated( bin03_ss[ 0 , a ], yfit_bin03_ss[ a ], /double )
sou_ss[ 0 ] = def_integrate_poly4( bin01_ss_pfit[0], bin01_ss_pfit[1], bin01_ss_pfit[2], $
	bin01_ss_pfit[3], bin01_ss_pfit[4], sin( (-30) * rad ), sin( 0 * rad ))
sou_ss[ 1 ] = def_integrate_poly4( bin02_ss_pfit[0], bin02_ss_pfit[1], bin02_ss_pfit[2], $
	bin02_ss_pfit[3], bin02_ss_pfit[4], sin( (-30) * rad ), sin( 0 * rad ))
sou_ss[ 2 ] = def_integrate_poly4( bin03_ss_pfit[0], bin03_ss_pfit[1], bin03_ss_pfit[2], $
	bin03_ss_pfit[3], bin03_ss_pfit[4], sin( (-30) * rad ), sin( 0 * rad ))

;a = where( bin01_ogi[ 0 , * ] gt 0 and bin01_ogi[ 0 , * ] lt sin( 70*rad ) )
;nor_ogi = int_tabulated( bin01_ogi[ 0 , * ], yfit_bin01_ogi[ * ], /double )
nor_ogi = def_integrate_poly4( bin01_ogi_pfit[0], bin01_ogi_pfit[1], bin01_ogi_pfit[2], $
	bin01_ogi_pfit[3], bin01_ogi_pfit[4], sin( 0 * rad ), sin( 70 * rad ))

;a = where( bin01_ogi[ 0 , * ] lt 0 and bin01_ogi[ 0 , * ] gt sin( -30*rad ) )
;sou_ogi = int_tabulated( bin01_ogi[ 0 , * ], yfit_bin01_ogi[ * ], /double )
sou_ogi = def_integrate_poly4( bin01_ogi_pfit[0], bin01_ogi_pfit[1], bin01_ogi_pfit[2], $
	bin01_ogi_pfit[3], bin01_ogi_pfit[4], sin( (-30) * rad ), sin( 0 * rad ))
;calculating the area under the curve of the fit for the sim data
nor_sim_noaa = fltarr( 5 )
;a = where( bin01_noaa[ 0 , * ] gt 0 and bin01_noaa[ 0 , * ] lt sin( 70*rad ) )
;nor_sim_noaa[ 0 ] = int_tabulated( bin01_noaa[ 0 , a ], yfit_bin01_sim_noaa[ a ], /double )
;a = where( bin02_noaa[ 0 , * ] gt 0 and bin02_noaa[ 0 , * ] lt sin( 70*rad ) )
;nor_sim_noaa[ 1 ] = int_tabulated( bin02_noaa[ 0 , a ], yfit_bin02_sim_noaa[ a ], /double )
;a = where( bin03_noaa[ 0 , * ] gt 0 and bin03_noaa[ 0 , * ] lt sin( 70*rad ) )
;nor_sim_noaa[ 2 ] = int_tabulated( bin03_noaa[ 0 , a ], yfit_bin03_sim_noaa[ a ], /double )
;a = where( bin04_noaa[ 0 , * ] gt 0 and bin04_noaa[ 0 , * ] lt sin( 70*rad ) )
;nor_sim_noaa[ 3 ] = int_tabulated( bin04_noaa[ 0 , a ], yfit_bin04_sim_noaa[ a ], /double )
;a = where( bin05_noaa[ 0 , * ] gt 0 and bin05_noaa[ 0 , * ] lt sin( 70*rad ) )
;nor_sim_noaa[ 4 ] = int_tabulated( bin05_noaa[ 0 , a ], yfit_bin05_sim_noaa[ a ], /double )
nor_sim_noaa[ 0 ] = def_integrate_poly4( bin01_sim_noaa_pfit[0], bin01_sim_noaa_pfit[1], bin01_sim_noaa_pfit[2], $
	bin01_sim_noaa_pfit[3], bin01_sim_noaa_pfit[4], sin( 0 * rad ), sin( 70 * rad ))
nor_sim_noaa[ 1 ] = def_integrate_poly4( bin02_sim_noaa_pfit[0], bin02_sim_noaa_pfit[1], bin02_sim_noaa_pfit[2], $
	bin02_sim_noaa_pfit[3], bin02_sim_noaa_pfit[4], sin( 0 * rad ), sin( 70 * rad ))
nor_sim_noaa[ 2 ] = def_integrate_poly4( bin03_sim_noaa_pfit[0], bin03_sim_noaa_pfit[1], bin03_sim_noaa_pfit[2], $
	bin03_sim_noaa_pfit[3], bin03_sim_noaa_pfit[4], sin( 0 * rad ), sin( 70 * rad ))
nor_sim_noaa[ 3 ] = def_integrate_poly4( bin04_sim_noaa_pfit[0], bin04_sim_noaa_pfit[1], bin04_sim_noaa_pfit[2], $
	bin04_sim_noaa_pfit[3], bin04_sim_noaa_pfit[4], sin( 0 * rad ), sin( 70 * rad ))
nor_sim_noaa[ 4 ] = def_integrate_poly4( bin05_sim_noaa_pfit[0], bin05_sim_noaa_pfit[1], bin05_sim_noaa_pfit[2], $
	bin05_sim_noaa_pfit[3], bin05_sim_noaa_pfit[4], sin( 0 * rad ), sin( 70 * rad ))
	
sou_sim_noaa = fltarr( 5 )
;a = where( bin01_noaa[ 0 , * ] lt 0 and bin01_noaa[ 0 , * ] gt sin( -30*rad ) )
;sou_sim_noaa[ 0 ] = int_tabulated( bin01_noaa[ 0 , a ], yfit_bin01_sim_noaa[ a ], /double )
;a = where( bin02_noaa[ 0 , * ] lt 0 and bin02_noaa[ 0 , * ] gt sin( -30*rad ) )
;sou_sim_noaa[ 1 ] = int_tabulated( bin02_noaa[ 0 , a ], yfit_bin02_sim_noaa[ a ], /double )
;a = where( bin03_noaa[ 0 , * ] lt 0 and bin03_noaa[ 0 , * ] gt sin( -30*rad ) )
;sou_sim_noaa[ 2 ] = int_tabulated( bin03_noaa[ 0 , a ], yfit_bin03_sim_noaa[ a ], /double )
;a = where( bin04_noaa[ 0 , * ] lt 0 and bin04_noaa[ 0 , * ] gt sin( -30*rad ) )
;sou_sim_noaa[ 3 ] = int_tabulated( bin04_noaa[ 0 , a ], yfit_bin04_sim_noaa[ a ], /double )
;a = where( bin05_noaa[ 0 , * ] lt 0 and bin05_noaa[ 0 , * ] gt sin( -30*rad ) )
;sou_sim_noaa[ 4 ] = int_tabulated( bin05_noaa[ 0 , a ], yfit_bin05_sim_noaa[ a ], /double )
sou_sim_noaa[ 0 ] = def_integrate_poly4( bin01_sim_noaa_pfit[0], bin01_sim_noaa_pfit[1], bin01_sim_noaa_pfit[2], $
	bin01_sim_noaa_pfit[3], bin01_sim_noaa_pfit[4], sin( (-30) * rad ), sin( 0 * rad ))
sou_sim_noaa[ 1 ] = def_integrate_poly4( bin02_sim_noaa_pfit[0], bin02_sim_noaa_pfit[1], bin02_sim_noaa_pfit[2], $
	bin02_sim_noaa_pfit[3], bin02_sim_noaa_pfit[4], sin( (-30) * rad ), sin( 0 * rad ))
sou_sim_noaa[ 2 ] = def_integrate_poly4( bin03_sim_noaa_pfit[0], bin03_sim_noaa_pfit[1], bin03_sim_noaa_pfit[2], $
	bin03_sim_noaa_pfit[3], bin03_sim_noaa_pfit[4], sin( (-30) * rad ), sin( 0 * rad ))
sou_sim_noaa[ 3 ] = def_integrate_poly4( bin04_sim_noaa_pfit[0], bin04_sim_noaa_pfit[1], bin04_sim_noaa_pfit[2], $
	bin04_sim_noaa_pfit[3], bin04_sim_noaa_pfit[4], sin( (-30) * rad ), sin( 0 * rad ))
sou_sim_noaa[ 4 ] = def_integrate_poly4( bin05_sim_noaa_pfit[0], bin05_sim_noaa_pfit[1], bin05_sim_noaa_pfit[2], $
	bin05_sim_noaa_pfit[3], bin05_sim_noaa_pfit[4], sin( (-30) * rad ), sin( 0 * rad ))

nor_sim_ss = fltarr( 3 )
;a = where( bin01_ss[ 0 , * ] gt 0 and bin01_ss[ 0 , * ] lt sin( 70*rad ) )
;nor_sim_ss[ 0 ] = int_tabulated( bin01_ss[ 0 , a ], yfit_bin01_sim_ss[ a ], /double )
;a = where( bin02_ss[ 0 , * ] gt 0 and bin02_ss[ 0 , * ] lt sin( 70*rad ) )
;nor_sim_ss[ 1 ] = int_tabulated( bin02_ss[ 0 , a ], yfit_bin02_sim_ss[ a ], /double )
;a = where( bin03_ss[ 0 , * ] gt 0 and bin03_ss[ 0 , * ] lt sin( 70*rad ) )
;nor_sim_ss[ 2 ] = int_tabulated( bin03_ss[ 0 , a ], yfit_bin03_sim_ss[ a ], /double )
nor_sim_ss[ 0 ] = def_integrate_poly4( bin01_sim_ss_pfit[0], bin01_sim_ss_pfit[1], bin01_sim_ss_pfit[2], $
	bin01_sim_ss_pfit[3], bin01_sim_ss_pfit[4], sin( 0 * rad ), sin( 70 * rad ))
nor_sim_ss[ 1 ] = def_integrate_poly4( bin02_sim_ss_pfit[0], bin02_sim_ss_pfit[1], bin02_sim_ss_pfit[2], $
	bin02_sim_ss_pfit[3], bin02_sim_ss_pfit[4], sin( 0 * rad ), sin( 70 * rad ))
nor_sim_ss[ 2 ] = def_integrate_poly4( bin03_sim_ss_pfit[0], bin03_sim_ss_pfit[1], bin03_sim_ss_pfit[2], $
	bin03_sim_ss_pfit[3], bin03_sim_ss_pfit[4], sin( 0 * rad ), sin( 70 * rad ))

sou_sim_ss = fltarr( 5 )
;a = where( bin01_ss[ 0 , * ] lt 0 and bin01_ss[ 0 , * ] gt sin( -30*rad ) )
;sou_sim_ss[ 0 ] = int_tabulated( bin01_ss[ 0 , a ], yfit_bin01_sim_ss[ a ], /double )
;a = where( bin02_ss[ 0 , * ] lt 0 and bin02_ss[ 0 , * ] gt sin( -30*rad ) )
;sou_sim_ss[ 1 ] = int_tabulated( bin02_ss[ 0 , a ], yfit_bin02_sim_ss[ a ], /double )
;a = where( bin03_ss[ 0 , * ] lt 0 and bin03_ss[ 0 , * ] gt sin( -30*rad ) )
;sou_sim_ss[ 2 ] = int_tabulated( bin03_ss[ 0 , a ], yfit_bin03_sim_ss[ a ], /double )
sou_sim_ss[ 0 ] = def_integrate_poly4( bin01_sim_ss_pfit[0], bin01_sim_ss_pfit[1], bin01_sim_ss_pfit[2], $
	bin01_sim_ss_pfit[3], bin01_sim_ss_pfit[4], sin( (-30) * rad ), sin( 0 * rad ))
sou_sim_ss[ 1 ] = def_integrate_poly4( bin02_sim_ss_pfit[0], bin02_sim_ss_pfit[1], bin02_sim_ss_pfit[2], $
	bin02_sim_ss_pfit[3], bin02_sim_ss_pfit[4], sin( (-30) * rad ), sin( 0 * rad ))
sou_sim_ss[ 2 ] = def_integrate_poly4( bin03_sim_ss_pfit[0], bin03_sim_ss_pfit[1], bin03_sim_ss_pfit[2], $
	bin03_sim_ss_pfit[3], bin03_sim_ss_pfit[4], sin( (-30) * rad ), sin( 0 * rad ))

;a = where( bin01_ogi[ 0 , * ] gt 0 and bin01_ogi[ 0 , * ] lt sin( 70*rad ) )
;nor_sim_ogi = int_tabulated( bin01_ogi[ 0 , * ], yfit_bin01_sim_ogi[ * ], /double )
nor_sim_ogi = def_integrate_poly4( bin01_sim_ogi_pfit[0], bin01_sim_ogi_pfit[1], bin01_sim_ogi_pfit[2], $
	bin01_sim_ogi_pfit[3], bin01_sim_ogi_pfit[4], sin( 0 * rad ), sin( 70 * rad ))

;a = where( bin01_ogi[ 0 , * ] lt 0 and bin01_ogi[ 0 , * ] gt sin( -30*rad ) )
;sou_sim_ogi = int_tabulated( bin01_ogi[ 0 , * ], yfit_bin01_sim_ogi[ * ], /double )
sou_sim_ogi = def_integrate_poly4( bin01_sim_ogi_pfit[0], bin01_sim_ogi_pfit[1], bin01_sim_ogi_pfit[2], $
	bin01_sim_ogi_pfit[3], bin01_sim_ogi_pfit[4], sin( (-30) * rad ), sin( 0 * rad ))
	
;calculating the error of the NHM and SHM
nor_noaa_err = fltarr( 5 )
a = where( bin01_noaa[ 0 , * ] gt 0 and bin01_noaa[ 0 , * ] lt sin( 70*rad ) )
nor_noaa_err[ 0 ] = stddev( bin01_noaa[ 2 , a ] - yfit_bin01_noaa[ a ] ) / $
	sqrt( n_elements( a ) ) 
a = where( bin02_noaa[ 0 , * ] gt 0 and bin02_noaa[ 0 , * ] lt sin( 70*rad ) )
nor_noaa_err[ 1 ] = stddev( bin02_noaa[ 2 , a ] - yfit_bin02_noaa[ a ] ) / $
	sqrt( n_elements( a ) ) 
a = where( bin03_noaa[ 0 , * ] gt 0 and bin03_noaa[ 0 , * ] lt sin( 70*rad ) )
nor_noaa_err[ 2 ] = stddev( bin03_noaa[ 2 , a ] - yfit_bin03_noaa[ a ] ) / $
	sqrt( n_elements( a ) ) 
a = where( bin04_noaa[ 0 , * ] gt 0 and bin04_noaa[ 0 , * ] lt sin( 70*rad ) )
nor_noaa_err[ 3 ] = stddev( bin04_noaa[ 2 , a ] - yfit_bin04_noaa[ a ] ) / $
	sqrt( n_elements( a ) ) 
a = where( bin05_noaa[ 0 , * ] gt 0 and bin05_noaa[ 0 , * ] lt sin( 70*rad ) )
nor_noaa_err[ 4 ] = stddev( bin05_noaa[ 2 , a ] - yfit_bin05_noaa[ a ] ) / $
	sqrt( n_elements( a ) )
	
sou_noaa_err = fltarr( 5 )
a = where( bin01_noaa[ 0 , * ] lt 0 and bin01_noaa[ 0 , * ] gt sin( -30*rad ) )
sou_noaa_err[ 0 ] = stddev( bin01_noaa[ 2 , a ] - yfit_bin01_noaa[ a ] ) / $
	sqrt( n_elements( a ) ) 
a = where( bin02_noaa[ 0 , * ] lt 0 and bin02_noaa[ 0 , * ] gt sin( -30*rad ) )
sou_noaa_err[ 1 ] = stddev( bin02_noaa[ 2 , a ] - yfit_bin02_noaa[ a ] ) / $
	sqrt( n_elements( a ) ) 
a = where( bin03_noaa[ 0 , * ] lt 0 and bin03_noaa[ 0 , * ] gt sin( -30*rad ) )
sou_noaa_err[ 2 ] = stddev( bin03_noaa[ 2 , a ] - yfit_bin03_noaa[ a ] ) / $
	sqrt( n_elements( a ) ) 
a = where( bin04_noaa[ 0 , * ] lt 0 and bin04_noaa[ 0 , * ] gt sin( -30*rad ) )
sou_noaa_err[ 3 ] = stddev( bin04_noaa[ 2 , a ] - yfit_bin04_noaa[ a ] ) / $
	sqrt( n_elements( a ) ) 
a = where( bin05_noaa[ 0 , * ] lt 0 and bin05_noaa[ 0 , * ] gt sin( -30*rad ) )
sou_noaa_err[ 4 ] = stddev( bin05_noaa[ 2 , a ] - yfit_bin05_noaa[ a ] ) / $
	sqrt( n_elements( a ) ) 
	
nor_ss_err = fltarr( 3 )
a = where( bin01_ss[ 0 , * ] gt 0 and bin01_ss[ 0 , * ] lt sin( 70*rad ) )
nor_ss_err[ 0 ] = stddev( bin01_ss[ 2 , a ] - yfit_bin01_ss[ a ] ) / $
	sqrt( n_elements( a ) ) 
a = where( bin02_ss[ 0 , * ] gt 0 and bin02_ss[ 0 , * ] lt sin( 70*rad ) )
nor_ss_err[ 1 ] = stddev( bin02_ss[ 2 , a ] - yfit_bin02_ss[ a ] ) / $
	sqrt( n_elements( a ) ) 
a = where( bin03_ss[ 0 , * ] gt 0 and bin03_ss[ 0 , * ] lt sin( 70*rad ) )
nor_ss_err[ 2 ] = stddev( bin03_ss[ 2 , a ] - yfit_bin03_ss[ a ] ) / $
	sqrt( n_elements( a ) ) 

sou_ss_err = fltarr( 3 )
a = where( bin01_ss[ 0 , * ] lt 0 and bin01_ss[ 0 , * ] gt sin( -30*rad ) )
sou_ss_err[ 0 ] = stddev( bin01_ss[ 2 , a ] - yfit_bin01_ss[ a ] ) / $
	sqrt( n_elements( a ) ) 
a = where( bin02_ss[ 0 , * ] lt 0 and bin02_ss[ 0 , * ] gt sin( -30*rad ) )
sou_ss_err[ 1 ] = stddev( bin02_ss[ 2 , a ] - yfit_bin02_ss[ a ] ) / $
	sqrt( n_elements( a ) ) 
a = where( bin03_ss[ 0 , * ] lt 0 and bin03_ss[ 0 , * ] gt sin( -30*rad ) )
sou_ss_err[ 2 ] = stddev( bin03_ss[ 2 , a ] - yfit_bin03_ss[ a ] ) / $
	sqrt( n_elements( a ) ) 
	
a = where( bin01_ogi[ 0 , * ] gt 0 and bin01_ogi[ 0 , * ] lt sin( 70*rad ) )
nor_ogi_err = stddev( bin01_ogi[ 2 , a ] - yfit_bin01_ogi[ a ] ) / $
	sqrt( n_elements( a ) ) 

a = where( bin01_ogi[ 0 , * ] lt 0 and bin01_ogi[ 0 , * ] gt sin( -30*rad ) )
sou_ogi_err = stddev( bin01_ogi[ 2 , a ] - yfit_bin01_ogi[ a ] ) / $
	sqrt( n_elements( a ) ) 	
;calculating the interhemispheric differences (IHD) from the NHM and SHM
noaa_ihd = nor_noaa - sou_noaa
noaa_sim_ihd = nor_sim_noaa - sou_sim_noaa

ss_ihd = nor_ss - sou_ss
ss_sim_ihd = nor_sim_ss - sou_sim_ss

ogi_ihd = nor_ogi - sou_ogi
ogi_sim_ihd = nor_sim_ogi - sou_sim_ogi

;--------------------------------------------------------------------------------------
;This part is to retrieve all the coordinates of the NOAA stations
;the array noaa_coord contains the lat and lon values of all NOAA stations
;noaa_coord[ 0 , * ] is latitude
;noaa_coord[ 1 , * ] is longitude
noaa_coord = fltarr( 2, n_elements( noaa_ids ) )
for i = 0, n_elements( noaa_ids ) - 1 do begin 
	a = where( strmatch( location_noaa[ * ], noaa_ids[ i ], /FOLD_CASE ) EQ 1 )
	noaa_coord[ 0 , i ] = lat1_noaa[ a [ 0 ] ]
	noaa_coord[ 1 , i ] = lon1_noaa[ a [ 0 ] ]
endfor 

;converting all the NOAA locations to indecies in order to input to simArr
noaa_coord_idx = fltarr( 2 , n_elements( noaa_ids) )
noaa_coord_idx[ 0 , * ] = round( ( noaa_coord[ 0 , * ] - (-90) ) / 2 ) 
noaa_coord_idx[ 1 , * ] = round( ( noaa_coord[ 1 , * ] - (-180) ) / 2.5 ) 

a = sort( noaa_coord[ 0 , * ] )
noaa_coord[ 0 , * ] = noaa_coord[ 0 , a ]
noaa_coord[ 1 , * ] = noaa_coord[ 1 , a ]
noaa_coord_idx[ 0 , * ] = noaa_coord_idx[ 0 , a ]
noaa_coord_idx[ 1 , * ] = noaa_coord_idx[ 1 , a ]

;--------------------------------------------------------------------------------------
;This section extract and convert coordinates of UCI data to indecies to use
;in the 3-D array that stores the simulated data
;ss_coord stores the coordinates of the UCI data
;ss_coord[ 0 , * ] stores latitudes
;ss_coord[ 1 , * ] stores longitudes
ss_coord = fltarr( 2 , n_elements( lat_var_GT30 ) )
for i = 0, n_elements( lat_var_GT30 ) - 1 do begin
	a = where( lat1_ss[ * ] eq lat_var_GT30[ i ] )
	ss_coord[ 0 , i ] = lat1_ss[ a[ 0 ] ]
	ss_coord[ 1 , i ] = lon1_ss[ a[ 0 ] ]
endfor

ss_coord_idx = fltarr( 2 , n_elements( lat_var_GT30 ) )
ss_coord_idx[ 0 , * ] = round( ( ss_coord[ 0 , * ] - (-90 ) ) / 2 ) 
ss_coord_idx[ 1 , * ] = round( ( ss_coord[ 1 , * ] - (-180) ) / 2.5 ) 

a = sort( ss_coord[ 0 , * ] )
ss_coord[ 0 , * ] = ss_coord[ 0 , a ]
ss_coord[ 1 , * ] = ss_coord[ 1 , a ]
ss_coord_idx[ 0 , * ] = ss_coord_idx[ 0 , a ]
ss_coord_idx[ 1 , * ] = ss_coord_idx[ 1 , a ]

;fix the idx where it is out of range
ss_coord_idx[ 1 , 7 ] = 143
;------------------------------------------------------------------------------------------

;this following is a special loop to average the simArr annually
;creating an array to store the annual values
annualSimArr = fltarr( 34, 144 , 91 )
for i = 0, 33 do begin	;there are 34 years in the simulation, so from 0 to 33
	j = ( i ) * 12 
	j0 = findgen( 12 ) + j
	annualSimArr[ i , * , * ] = mean( simArr[ j0 , * , * ], 1 ) / 2 * 1000
endfor


;-------------------------------------------------------------------------------------------
;this next section pull out the simulated values using the coordinates indecies of the 
;UCI data
;annual_sim_avg[ * , 1 , * ] is the averages
;annual_sim_avg[ * , 0 , * ] is the latitudes
annual_ss_sim_avg = fltarr( 34 , 2 , n_elements( ss_coord_idx[ 1 , * ] ) )
for i = 0, 33 do begin
	for j = 0, n_elements( lat_var_GT30 ) - 1 do begin
		annual_ss_sim_avg[ i , 1 , j ] = annualSimArr[ i , ss_coord_idx[ 1 , j ] , ss_coord_idx[ 0 , j ] ]
		annual_ss_sim_avg[ i , 0 , j ] = ss_coord[ 0 , j ]
	endfor
endfor

;polynomial fit the simulated data
annual_ss_sim_avg[ * , 0 , * ] = sin( annual_ss_sim_avg[ * , 0 , * ] * !pi / 180 )
nor_ss_sim_avg_var = fltarr( 34 )
sou_ss_sim_avg_var = fltarr( 34 )
sim_fit_coeff = fltarr( 34 , 5 , n_elements( lat_var_GT30) ) ;i comp is years, j comp is number of 
;coeffcients, k comp is the number of stations of the UCI data.
for i = 0, 33 do begin
	a = poly_fit( annual_ss_sim_avg[ i , 0 , * ], annual_ss_sim_avg[ i , 1 , * ], 4, /double )
	;this section is for calculating the area under the curve of the fit using algebraic analysis
	temp = def_integrate_poly4( a[0], a[1], a[2], a[3], a[4], sin( 0 * rad ), sin( 70 * rad ) )
	nor_ss_sim_avg_var[ i ] = temp
	temp = def_integrate_poly4( a[0], a[1], a[2], a[3], a[4], sin( -30 * rad ), sin( 0 * rad ) )
	sou_ss_sim_avg_var[ i ] = temp
endfor


;calculating the IHD for all the sim data that is plotted outside the available observed data
sim_ss_allyear_ihd = fltarr( 34 )
sim_ss_allyear_ihd = nor_ss_sim_avg_var - sou_ss_sim_avg_var


;-------------------------------------------------------------------------------------------
;this next section pull out the simulated values using the coordinates indecies of the 
;NOAA
;annual_sim_avg[ * , 1 , * ] is the averages
;annual_sim_avg[ * , 0 , * ] is the latitudes
annual_noaa_sim_avg = fltarr( 34 , 2 , n_elements( noaa_coord_idx[ 1 , * ] ) )
for i = 0, 33 do begin
	for j = 0, n_elements( noaa_ids ) - 1 do begin
		annual_noaa_sim_avg[ i , 1 , j ] = annualSimArr[ i , noaa_coord_idx[ 1 , j ] , noaa_coord_idx[ 0 , j ] ]
		annual_noaa_sim_avg[ i , 0 , j ] = noaa_coord[ 0 , j ]
	endfor
endfor

;polynomial fit the simulated data
annual_noaa_sim_avg[ * , 0 , * ] = sin( annual_noaa_sim_avg[ * , 0 , * ] * !pi / 180 )
nor_noaa_sim_avg_var = fltarr( 34 )
sou_noaa_sim_avg_var = fltarr( 34 )
sim_fit_coeff = fltarr( 34 , 5 , n_elements( noaa_ids) ) ;i comp is years, j comp is number of 
;coeffcients, k comp is the number of stations of the UCI data.
for i = 0, 33 do begin
	a = poly_fit( annual_noaa_sim_avg[ i , 0 , * ], annual_noaa_sim_avg[ i , 1 , * ], 4, /double )
	;this section is for calculating the area under the curve of the fit using algebraic analysis
	temp = def_integrate_poly4( a[0], a[1], a[2], a[3], a[4], sin( 0 * rad ), sin( 70 * rad ) )
	nor_noaa_sim_avg_var[ i ] = temp
	temp = def_integrate_poly4( a[0], a[1], a[2], a[3], a[4], sin( -30 * rad ), sin( 0 * rad ) )
	sou_noaa_sim_avg_var[ i ] = temp
endfor

;calculating the IHD for all the sim data that is plotted outside the available observed data
sim_noaa_allyear_ihd = fltarr( 34 )
sim_noaa_allyear_ihd = nor_noaa_sim_avg_var - sou_noaa_sim_avg_var

;calculating the IHD uncertainties from NHM and SHM
noaa_ihd_err = sqrt( nor_noaa_err ^ 2 + sou_noaa_err ^ 2 )

ss_ihd_err = sqrt( nor_ss_err ^ 2 + sou_ss_err ^ 2 )

ogi_ihd_err = sqrt( nor_ogi_err ^ 2 + sou_ogi_err ^ 2 )

ogi_time = [ 1985 ]
ss_time = [ 1998 , 2003 , 2008 ]
noaa_time = [2007 , 2009 , 2011 , 2013 , 2014]
all_time = findgen( 34 ) + 1981

!P.Multi = [ 0 , 1 , 2 , 0 , 0 ]
cgplot, noaa_time, nor_noaa, psym = 2, color = 'red', yrange = [ 0 , 1200 ], $
	xrange = [ 1980 , 2015 ], err_yhigh = nor_noaa_err, err_ylow = nor_noaa_err
cgplot, noaa_time, sou_noaa, psym = 1, color = 'red', /overplot, $
	err_yhigh = sou_noaa_err, err_ylow = sou_noaa_err
cgplot, ss_time, nor_ss, psym = 2, color = 'blue', /overplot, $
	err_yhigh = nor_ss_err, err_ylow = nor_ss_err
cgplot, ss_time, sou_ss, psym = 1, color = 'blue', /overplot, $
	err_yhigh = sou_ss_err, err_ylow = sou_ss_err
cgplot, ogi_time, nor_ogi, psym = 2, color = 'violet', /overplot, $
	err_yhigh = nor_ogi_err, err_ylow = nor_ogi_err
cgplot, ogi_time, sou_ogi, psym = 1, colot = 'violet', /overplot, $
	err_yhigh = sou_ogi_err, err_ylow = sou_ogi_err
cgplot, noaa_time, nor_sim_noaa, color = 'red', /overplot
cgplot, noaa_time, sou_sim_noaa, color = 'red', /overplot 
cgplot, ss_time, nor_sim_ss, color = 'blue', /overplot
cgplot, ss_time, sou_sim_ss, color = 'blue', /overplot
cgplot, all_time, nor_ss_sim_avg_var, /overplot, color = 'brown'
cgplot, all_time, sou_ss_sim_avg_var, /overplot, color = 'brown'
cgplot, all_time, nor_noaa_sim_avg_var, /overplot, color = 'forest green'
cgplot, all_time, sou_noaa_sim_avg_var, /overplot, color = 'forest green'

cgplot, noaa_time, noaa_ihd, psym = 4, yrange = [ 500, 1100 ], xrange = [ 1980 , 2015 ], $
	err_yhigh = noaa_ihd_err, err_ylow = noaa_ihd_err, color = 'red', title = sim_title
cgplot, ogi_time, ogi_ihd, psym = 4, color = 'violet', /overplot
cgplot, noaa_time, noaa_sim_ihd, /overplot, color = 'red'
cgplot, all_time, sim_ss_allyear_ihd, /overplot, color = 'brown'
cgplot, all_time, sim_noaa_allyear_ihd, /overplot, color = 'forest green'
cgplot, ss_time, ss_ihd, /overplot, psym = 4, color = 'blue', $
	err_yhigh = ss_ihd_err, err_ylow = ss_ihd_err
cgplot, ss_time, ss_sim_ihd, /overplot, color = 'blue'

;Make a JPEG for the graphs
READ, option, PROMPT= 'Do I need to save this graph as a jpeg so I can print it out later? Enter "1" for HELLYEAH or "0" for NOPE       '
name = strarr(1)
IF (option EQ 1) THEN BEGIN
	READ, name, PROMPT= 'Name of the JPEG file? (No need to include the jpeg extension, just name is fine)   '
	SCREEN2JPG, name
	print, 'JPEG has been saved to current directory'
ENDIF
;making some plots to verify the fits of the obs and sim data
;!P.Multi = [ 0 , 2 , 3 , 0 , 0 ]
;cgplot, bin01_sim_noaa_var[ 0 , * ], yfit_bin01_sim_noaa, color = 'blue'
;cgplot, bin01_sim_noaa_var[ 0 , * ], bin01_sim_noaa_var[ 1 , * ], color = 'blue', /overplot, psym = 2
;cgplot, bin01_noaa[ 0 , * ], yfit_bin01_noaa, color = 'red', /overplot
;cgplot, bin01_noaa[ 0 , * ], bin01_noaa[ 2 , * ], color = 'red', /overplot, psym = 2

;cgplot, bin02_sim_noaa_var[ 0 , * ], yfit_bin02_sim_noaa, color = 'blue'
;cgplot, bin02_sim_noaa_var[ 0 , * ], bin02_sim_noaa_var[ 1 , * ], color = 'blue', /overplot, psym = 2
;cgplot, bin02_noaa[ 0 , * ], yfit_bin02_noaa, color = 'red', /overplot
;cgplot, bin02_noaa[ 0 , * ], bin02_noaa[ 2 , * ], color = 'red', /overplot, psym = 2

;cgplot, bin03_sim_noaa_var[ 0 , * ], yfit_bin03_sim_noaa, color = 'blue'
;cgplot, bin03_sim_noaa_var[ 0 , * ], bin03_sim_noaa_var[ 1 , * ], color = 'blue', /overplot, psym = 2
;cgplot, bin03_noaa[ 0 , * ], yfit_bin03_noaa, color = 'red', /overplot
;cgplot, bin03_noaa[ 0 , * ], bin03_noaa[ 2 , * ], color = 'red', /overplot, psym = 2

;cgplot, bin04_sim_noaa_var[ 0 , * ], yfit_bin04_sim_noaa, color = 'blue'
;cgplot, bin04_sim_noaa_var[ 0 , * ], bin04_sim_noaa_var[ 1 , * ], color = 'blue', /overplot, psym = 2
;cgplot, bin04_noaa[ 0 , * ], yfit_bin04_noaa, color = 'red', /overplot
;cgplot, bin04_noaa[ 0 , * ], bin04_noaa[ 2 , * ], color = 'red', /overplot, psym = 2

;cgplot, bin05_sim_noaa_var[ 0 , * ], yfit_bin05_sim_noaa, color = 'blue'
;cgplot, bin05_sim_noaa_var[ 0 , * ], bin05_sim_noaa_var[ 1 , * ], color = 'blue', /overplot, psym = 2
;cgplot, bin05_noaa[ 0 , * ], yfit_bin05_noaa, color = 'red'
;cgplot, bin05_noaa[ 0 , * ], bin05_noaa[ 2 , * ], color = 'red', /overplot, psym = 2


END

