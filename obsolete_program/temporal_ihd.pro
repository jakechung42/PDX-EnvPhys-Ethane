PRO TEMPORAL_IHD

;this program plot the temporal trend of the UCI, NOAA, OGI data
;the averages are calculated by the sine weighted integral of the poly fit.

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
print, n1_noaa, n2_noaa
;=====Finish reading in NOAA data files=====
;============================================================================================================
;============================================================================================================
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

;this line is not needed for this program, but keep it here anyway
;Convert unix time to decimal year
time_ogi = (unix_date_ogi/3600 - 131496)/8760 + 1985
;=====Finish reading in OGI data files=====

;run calculation for the NOAA file
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

;NOAA data is averaged for every two years, therefore there will be 6 bins for the NOAA data
;bin01: 2005 + 2006
;bin02: 2007 + 2008
;bin03: 2009 + 2010
;bin04: 2011 + 2012
;bin05: 2013 + 2014
;bin06: 2015
;setting up the arrays that will hold the averages
bin01_noaa = fltarr( 3 , noaa_ids_size )
bin02_noaa = fltarr( 3 , noaa_ids_size )
bin03_noaa = fltarr( 3 , noaa_ids_size )
bin04_noaa = fltarr( 3 , noaa_ids_size )
bin05_noaa = fltarr( 3 , noaa_ids_size )
bin06_noaa = fltarr( 3 , noaa_ids_size )

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
		bin01_noaa[ 2 , i ] = mean( temp_avg[ temp_bin01 ], /NAN )	
		bin01_noaa[ 0 , i ] = temp_lat[ 0 ]
		bin01_noaa[ 1 , i ] = temp_lon[ 0 ]
	endif 
	
	temp_bin02 = where( year_noaa[ a ] EQ 2007 OR year_noaa[ a ] EQ 2008 , count )
	if ( count GT 0 ) then begin
		bin02_noaa[ 2 , i ] = mean( temp_avg[ temp_bin02 ], /NAN )	
		bin02_noaa[ 0 , i ] = temp_lat[ 0 ]
		bin02_noaa[ 1 , i ] = temp_lon[ 0 ]
	endif 
		
	temp_bin03 = where( year_noaa[ a ] EQ 2009 OR year_noaa[ a ] EQ 2010 , count )
	if ( count GT 0 ) then begin
		bin03_noaa[ 2 , i ] = mean( temp_avg[ temp_bin03 ], /NAN )	
		bin03_noaa[ 0 , i ] = temp_lat[ 0 ]
		bin03_noaa[ 1 , i ] = temp_lon[ 0 ]
	endif 

	temp_bin04 = where( year_noaa[ a ] EQ 2011 OR year_noaa[ a ] EQ 2012 , count )
	if ( count GT 0 ) then begin
		bin04_noaa[ 2 , i ] = mean( temp_avg[ temp_bin04 ], /NAN )	
		bin04_noaa[ 0 , i ] = temp_lat[ 0 ]
		bin04_noaa[ 1 , i ] = temp_lon[ 0 ]
	endif 

	temp_bin05 = where( year_noaa[ a ] EQ 2013 OR year_noaa[ a ] EQ 2014 , count )
	if ( count GT 0 ) then begin
		bin05_noaa[ 2 , i ] = mean( temp_avg[ temp_bin05 ] )	
		bin05_noaa[ 0 , i ] = temp_lat[ 0 ]
		bin05_noaa[ 1 , i ] = temp_lon[ 0 ]
	endif 
			
	temp_bin06 = where( year_noaa[ a ] EQ 2015 , count )
	if ( count GT 0 ) then begin
		bin06_noaa[ 2 , i ] = mean( temp_avg[ temp_bin06 ] )	
		bin06_noaa[ 0 , i ] = temp_lat[ 0 ]
		bin06_noaa[ 1 , i ] = temp_lon[ 0 ]
	endif 
		
endfor 


;=========================================================================================


;the NOAA data runs from 2005 to 2015
;the SS data runs from 1996 - 2009

;Dealing with the SS data 
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
infile_bin06_de = "/home/jakechung/IDL/myIDL/temp_folder/temp_bin06_NOAA_deseason.dat"

;open to write in the temp files
openw , lun1 , infile_bin01_de , /get_lun
openw , lun2 , infile_bin02_de , /get_lun
openw , lun3 , infile_bin03_de , /get_lun
openw , lun4 , infile_bin04_de , /get_lun
openw , lun5 , infile_bin05_de , /get_lun
openw , lun6 , infile_bin06_de , /get_lun

;loop to pull out the year 
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

print, 'bin01_noaa', bin01_noaa
print, 'bin02_noaa', bin02_noaa
print, 'de_obs_bin01_noaa', de_obs_bin01_noaa
print, 'de_obs_bin02_noaa', de_obs_bin02_noaa
help, bin01_noaa, bin02_noaa
help, de_obs_bin01_noaa, de_obs_bin02_noaa

;this was where I made the mistake that lead to the wrong conclusion on S6 ems scenario being the best fit.
;I did not remove the 0 values  for the deseasonal data thus causing some station values to have incorrect deseasonal
;data which influenced the poly fit thus causes the averages to be calculated incorrectly.
;Now, bug has been fixed.
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
de_obs_bin06_noaa = RemoveRows( de_obs_bin06_noaa, where( bin06_noaa[ 0 , * ] EQ 0 ) )
bin06_noaa = RemoveRows( bin06_noaa, where( bin06_noaa[ 0 , * ] EQ 0 ) )


de_obs_bin01_ss = RemoveRows( de_obs_bin01_ss, where( bin01_ss[ 0 , * ] EQ 0 ) )
bin01_ss = RemoveRows( bin01_ss, where( bin01_ss[ 0 , * ] EQ 0 ) )
de_obs_bin02_ss = RemoveRows( de_obs_bin02_ss, where( bin02_ss[ 0 , * ] EQ 0 ) )
bin02_ss = RemoveRows( bin02_ss, where( bin02_ss[ 0 , * ] EQ 0 ) )
de_obs_bin03_ss = RemoveRows( de_obs_bin03_ss, where( bin03_ss[ 0 , * ] EQ 0 ) )
bin03_ss = RemoveRows( bin03_ss, where( bin03_ss[ 0 , * ] EQ 0 ) ) 




;Deal with sim data, this is going to be longggggggg!!!
;=====================================
;===Start extracting simulated data===
;=====================================
;Define the filename

;sim_title = 'Aydin emissions scaled to Xiao et al over 1996-2003'
;filename = "/home/excluded-from-backup/data/C2H6/trac_avg.AydinSF_1981_2015.bpch"

;sim_title = 'Unscaled (absolute) Aydin emissions'
;filename = "/home/excluded-from-backup/data/C2H6/trac_avg.AydinAbsSF_1981_2015.bpch"

;sim_title = 'PSU emissions scaled to Xiao et al over 1996-2003'
filename = "/home/excluded-from-backup/data/C2H6/trac_avg.PSUSF_1981_2015.bpch"

;sim_title = 'Unscaled PSU emissions using methane-to-ethane (MER) ratio of 5 for biomass burning, and 18 for fossil fuel'
;filename = "/home/excluded-from-backup/data/C2H6/trac_avg.PSU_MER3BB_MER18FF_1981_2015.bpch"

;sim_title = 'Unscaled emissions from Simpson et al.'
;filename = "/home/excluded-from-backup/data/C2H6/trac_avg.SimpsonSF_1981_2015.bpch"

;sim_title = 'Constant default base emissions '
;filename = "/home/excluded-from-backup/data/C2H6/trac_avg.spinup_GFED4_MAVG_1981_2015.198101010000"

;read in the file
ctm_get_data, datainfo, 'IJ-AVG-$', tracer=1, filename=filename, tau0=tau0

;Get MODELINFO and GRIDINFO structures, xmid, ymid hold lon/lat centers 
getmodelandgridinfo, datainfo[0], modelinfo, gridinfo, lon=xmid, lat=ymid
if (n_elements(datainfo) eq 0) then print, 'No data records found.'

;nt is the number of data points in the simulation.
nt = n_elements(datainfo)
globalavg= fltarr(nt)

;convert tau0 values to year month day array and store in sim_yymmdd
;sim_yymmdd stores the entire time element of the input simulation
sim_yymmdd= tau2yymmdd(tau0, /GEOS1)

;---------------------------------------------------------------------------------------------------------------------------------------------
;locate indexes of all months of the simulated data for noaa
sim_bin01_idx_noaa = where( sim_yymmdd.year EQ 2005 OR sim_yymmdd.year EQ 2006 )
sim_bin02_idx_noaa = where( sim_yymmdd.year EQ 2007 OR sim_yymmdd.year EQ 2008 )
sim_bin03_idx_noaa = where( sim_yymmdd.year EQ 2009 OR sim_yymmdd.year EQ 2010 )
sim_bin04_idx_noaa = where( sim_yymmdd.year EQ 2011 OR sim_yymmdd.year EQ 2012 )
sim_bin05_idx_noaa = where( sim_yymmdd.year EQ 2013 OR sim_yymmdd.year EQ 2014 )
sim_bin06_idx_noaa = where( sim_yymmdd.year EQ 2015 )

;extracting the simulated data for the corresponding months and years and locations
sim_bin01_var_noaa = fltarr(n_elements(bin01_noaa[0 , *]))
for i = 0, n_elements(bin01_noaa[0 , *]) - 1 do begin
	for t = 0, nt - 1 do begin
		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
		gridinfo= gridinfo, lat= [bin01_noaa[0, i] - 0.1, bin01_noaa[0, i] + 0.1], $
		lon= [bin01_noaa[1, i] - 0.1, bin01_noaa[1, i] + 0.1], alrange = [1, 3], average=7)
		globalavg[ t ]= data/2 * 1000
	endfor
	sim_bin01_var_noaa[ i ] = mean(globalavg[ sim_bin01_idx_noaa ])
endfor

sim_bin02_var_noaa = fltarr(n_elements(bin02_noaa[0 , *]))
for i = 0, n_elements(bin02_noaa[0 , *]) - 1 do begin
	for t = 0, nt - 1 do begin
		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
		gridinfo= gridinfo, lat= [bin02_noaa[0, i] - 0.1, bin02_noaa[0, i] + 0.1], $
		lon= [bin02_noaa[1, i] - 0.1, bin02_noaa[1, i] + 0.1], alrange = [1, 3], average=7)
		globalavg[ t ]= data/2 * 1000
	endfor
	sim_bin02_var_noaa[ i ] = mean(globalavg[ sim_bin02_idx_noaa ])
endfor

sim_bin03_var_noaa = fltarr(n_elements(bin03_noaa[0 , *]))
for i = 0, n_elements(bin03_noaa[0 , *]) - 1 do begin
	for t = 0, nt - 1 do begin
		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
		gridinfo= gridinfo, lat= [bin03_noaa[0, i] - 0.1, bin03_noaa[0, i] + 0.1], $
		lon= [bin03_noaa[1, i] - 0.1, bin03_noaa[1, i] + 0.1], alrange = [1, 3], average=7)
		globalavg[ t ]= data/2 * 1000
	endfor
	sim_bin03_var_noaa[ i ] = mean(globalavg[ sim_bin03_idx_noaa ])
endfor

sim_bin04_var_noaa = fltarr(n_elements(bin04_noaa[0 , *]))
for i = 0, n_elements(bin04_noaa[0 , *]) - 1 do begin
	for t = 0, nt - 1 do begin
		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
		gridinfo= gridinfo, lat= [bin04_noaa[0, i] - 0.1, bin04_noaa[0, i] + 0.1], $
		lon= [bin04_noaa[1, i] - 0.1, bin04_noaa[1, i] + 0.1], alrange = [1, 3], average=7)
		globalavg[ t ]= data/2 * 1000
	endfor
	sim_bin04_var_noaa[ i ] = mean(globalavg[ sim_bin04_idx_noaa ])
endfor

sim_bin05_var_noaa = fltarr(n_elements(bin05_noaa[0 , *]))
for i = 0, n_elements(bin05_noaa[0 , *]) - 1 do begin
	for t = 0, nt - 1 do begin
		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
		gridinfo= gridinfo, lat= [bin05_noaa[0, i] - 0.1, bin05_noaa[0, i] + 0.1], $
		lon= [bin05_noaa[1, i] - 0.1, bin05_noaa[1, i] + 0.1], alrange = [1, 3], average=7)
		globalavg[ t ]= data/2 * 1000
	endfor
	sim_bin05_var_noaa[ i ] = mean(globalavg[ sim_bin05_idx_noaa ])
endfor

sim_bin06_var_noaa = fltarr(n_elements(bin06_noaa[0 , *]))
for i = 0, n_elements(bin06_noaa[0 , *]) - 1 do begin
	for t = 0, nt - 1 do begin
		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
		gridinfo= gridinfo, lat= [bin06_noaa[0, i] - 0.1, bin06_noaa[0, i] + 0.1], $
		lon= [bin06_noaa[1, i] - 0.1, bin06_noaa[1, i] + 0.1], alrange = [1, 3], average=7)
		globalavg[ t ]= data/2 * 1000
	endfor
	sim_bin06_var_noaa[ i ] = mean(globalavg[ sim_bin06_idx_noaa ])
endfor







;---------------------------------------------------------------------------------------------------------------------------------------------
;locate indexes of all months of the simulated data for UCI
sim_bin01_idx_ss = where( sim_yymmdd.year EQ 1996 $
	OR sim_yymmdd.year EQ 1997 $
	OR sim_yymmdd.year EQ 1998 $
	OR sim_yymmdd.year EQ 1999 $
	OR sim_yymmdd.year EQ 2000 )
sim_bin02_idx_ss = where( sim_yymmdd.year EQ 2001 $
	OR sim_yymmdd.year EQ 2002 $
	OR sim_yymmdd.year EQ 2003 $
	OR sim_yymmdd.year EQ 2004 $
	OR sim_yymmdd.year EQ 2005 )
sim_bin03_idx_ss = where( sim_yymmdd.year EQ 2006 $
	OR sim_yymmdd.year EQ 2007 $
	OR sim_yymmdd.year EQ 2008 $
	OR sim_yymmdd.year EQ 2009 )

;extracting the simulated data for the corresponding months and years and locations
sim_bin01_var_ss = fltarr(n_elements(bin01_ss[0 , *]))
for i = 0, n_elements(bin01_ss[0 , *]) - 1 do begin
	for t = 0, nt - 1 do begin
		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
		gridinfo= gridinfo, lat= [bin01_ss[0, i] - 0.1, bin01_ss[0, i] + 0.1], $
		lon= [bin01_ss[1, i] - 0.1, bin01_ss[1, i] + 0.1], alrange = [1, 3], average=7)
		globalavg[ t ]= data/2 * 1000
	endfor
	sim_bin01_var_ss[ i ] = mean(globalavg[ sim_bin01_idx_ss ])
endfor

sim_bin02_var_ss = fltarr(n_elements(bin02_ss[0 , *]))
for i = 0, n_elements(bin02_ss[0 , *]) - 1 do begin
	for t = 0, nt - 1 do begin
		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
		gridinfo= gridinfo, lat= [bin02_ss[0, i] - 0.1, bin02_ss[0, i] + 0.1], $
		lon= [bin02_ss[1, i] - 0.1, bin02_ss[1, i] + 0.1], alrange = [1, 3], average=7)
		globalavg[ t ]= data/2 * 1000
	endfor
	sim_bin02_var_ss[ i ] = mean(globalavg[ sim_bin02_idx_ss ])
endfor

sim_bin03_var_ss = fltarr(n_elements(bin03_ss[0 , *]))
for i = 0, n_elements(bin03_ss[0 , *]) - 1 do begin
	for t = 0, nt - 1 do begin
		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
		gridinfo= gridinfo, lat= [bin03_ss[0, i] - 0.1, bin03_ss[0, i] + 0.1], $
		lon= [bin03_ss[1, i] - 0.1, bin03_ss[1, i] + 0.1], alrange = [1, 3], average=7)
		globalavg[ t ]= data/2 * 1000
	endfor
	sim_bin03_var_ss[ i ] = mean(globalavg[ sim_bin03_idx_ss ])
endfor





;---------------------------------------------------------------------------------------------------------------------------------------------
;locate indexes of all months of the simulated data for OGI
sim_bin01_idx_ogi = where( sim_yymmdd.year EQ 1983 $
	OR sim_yymmdd.year EQ 1984 $
	OR sim_yymmdd.year EQ 1985 $
	OR sim_yymmdd.year EQ 1986 $
	OR sim_yymmdd.year EQ 1987 )


;extracting the simulated data for the corresponding months and years and locations
sim_bin01_var_ogi = fltarr(n_elements(bin01_ogi[0 , *]))
for i = 0, n_elements(bin01_ogi[0 , *]) - 1 do begin
	for t = 0, nt - 1 do begin
		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
		gridinfo= gridinfo, lat= [bin01_ogi[0, i] - 0.1, bin01_ogi[0, i] + 0.1], $
		lon= [bin01_ogi[2, i] - 0.1, bin01_ogi[2, i] + 0.1], alrange = [1, 3], average=7)
		globalavg[ t ]= data/2 * 1000
	endfor
	sim_bin01_var_ogi[ i ] = mean(globalavg[ sim_bin01_idx_ogi ])
endfor




;sim_title = 'Aydin emissions scaled to Xiao et al over 1996-2003'
filename = "/home/excluded-from-backup/data/C2H6/trac_avg.AydinSF_1981_2015.bpch"
;read in the file
ctm_get_data, datainfo, 'IJ-AVG-$', tracer=1, filename=filename, tau0=tau0

;Get MODELINFO and GRIDINFO structures, xmid, ymid hold lon/lat centers 
getmodelandgridinfo, datainfo[0], modelinfo, gridinfo, lon=xmid, lat=ymid
if (n_elements(datainfo) eq 0) then print, 'No data records found.'

;extracting the simulated data for the corresponding months and years and locations
sim_bin01_var_noaa_01 = fltarr(n_elements(bin01_noaa[0 , *]))
for i = 0, n_elements(bin01_noaa[0 , *]) - 1 do begin
	for t = 0, nt - 1 do begin
		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
		gridinfo= gridinfo, lat= [bin01_noaa[0, i] - 0.1, bin01_noaa[0, i] + 0.1], $
		lon= [bin01_noaa[1, i] - 0.1, bin01_noaa[1, i] + 0.1], alrange = [1, 3], average=7)
		globalavg[ t ]= data/2 * 1000
	endfor
	sim_bin01_var_noaa_01[ i ] = mean(globalavg[ sim_bin01_idx_noaa ])
endfor

sim_bin02_var_noaa_01 = fltarr(n_elements(bin02_noaa[0 , *]))
for i = 0, n_elements(bin02_noaa[0 , *]) - 1 do begin
	for t = 0, nt - 1 do begin
		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
		gridinfo= gridinfo, lat= [bin02_noaa[0, i] - 0.1, bin02_noaa[0, i] + 0.1], $
		lon= [bin02_noaa[1, i] - 0.1, bin02_noaa[1, i] + 0.1], alrange = [1, 3], average=7)
		globalavg[ t ]= data/2 * 1000
	endfor
	sim_bin02_var_noaa_01[ i ] = mean(globalavg[ sim_bin02_idx_noaa ])
endfor

sim_bin03_var_noaa_01 = fltarr(n_elements(bin03_noaa[0 , *]))
for i = 0, n_elements(bin03_noaa[0 , *]) - 1 do begin
	for t = 0, nt - 1 do begin
		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
		gridinfo= gridinfo, lat= [bin03_noaa[0, i] - 0.1, bin03_noaa[0, i] + 0.1], $
		lon= [bin03_noaa[1, i] - 0.1, bin03_noaa[1, i] + 0.1], alrange = [1, 3], average=7)
		globalavg[ t ]= data/2 * 1000
	endfor
	sim_bin03_var_noaa_01[ i ] = mean(globalavg[ sim_bin03_idx_noaa ])
endfor

sim_bin04_var_noaa_01 = fltarr(n_elements(bin04_noaa[0 , *]))
for i = 0, n_elements(bin04_noaa[0 , *]) - 1 do begin
	for t = 0, nt - 1 do begin
		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
		gridinfo= gridinfo, lat= [bin04_noaa[0, i] - 0.1, bin04_noaa[0, i] + 0.1], $
		lon= [bin04_noaa[1, i] - 0.1, bin04_noaa[1, i] + 0.1], alrange = [1, 3], average=7)
		globalavg[ t ]= data/2 * 1000
	endfor
	sim_bin04_var_noaa_01[ i ] = mean(globalavg[ sim_bin04_idx_noaa ])
endfor

sim_bin05_var_noaa_01 = fltarr(n_elements(bin05_noaa[0 , *]))
for i = 0, n_elements(bin05_noaa[0 , *]) - 1 do begin
	for t = 0, nt - 1 do begin
		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
		gridinfo= gridinfo, lat= [bin05_noaa[0, i] - 0.1, bin05_noaa[0, i] + 0.1], $
		lon= [bin05_noaa[1, i] - 0.1, bin05_noaa[1, i] + 0.1], alrange = [1, 3], average=7)
		globalavg[ t ]= data/2 * 1000
	endfor
	sim_bin05_var_noaa_01[ i ] = mean(globalavg[ sim_bin05_idx_noaa ])
endfor

sim_bin06_var_noaa_01 = fltarr(n_elements(bin06_noaa[0 , *]))
for i = 0, n_elements(bin06_noaa[0 , *]) - 1 do begin
	for t = 0, nt - 1 do begin
		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
		gridinfo= gridinfo, lat= [bin06_noaa[0, i] - 0.1, bin06_noaa[0, i] + 0.1], $
		lon= [bin06_noaa[1, i] - 0.1, bin06_noaa[1, i] + 0.1], alrange = [1, 3], average=7)
		globalavg[ t ]= data/2 * 1000
	endfor
	sim_bin06_var_noaa_01[ i ] = mean(globalavg[ sim_bin06_idx_noaa ])
endfor





;extracting the simulated data for the corresponding months and years and locations
sim_bin01_var_ss_01 = fltarr(n_elements(bin01_ss[0 , *]))
for i = 0, n_elements(bin01_ss[0 , *]) - 1 do begin
	for t = 0, nt - 1 do begin
		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
		gridinfo= gridinfo, lat= [bin01_ss[0, i] - 0.1, bin01_ss[0, i] + 0.1], $
		lon= [bin01_ss[1, i] - 0.1, bin01_ss[1, i] + 0.1], alrange = [1, 3], average=7)
		globalavg[ t ]= data/2 * 1000
	endfor
	sim_bin01_var_ss_01[ i ] = mean(globalavg[ sim_bin01_idx_ss ])
endfor

sim_bin02_var_ss_01 = fltarr(n_elements(bin02_ss[0 , *]))
for i = 0, n_elements(bin02_ss[0 , *]) - 1 do begin
	for t = 0, nt - 1 do begin
		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
		gridinfo= gridinfo, lat= [bin02_ss[0, i] - 0.1, bin02_ss[0, i] + 0.1], $
		lon= [bin02_ss[1, i] - 0.1, bin02_ss[1, i] + 0.1], alrange = [1, 3], average=7)
		globalavg[ t ]= data/2 * 1000
	endfor
	sim_bin02_var_ss_01[ i ] = mean(globalavg[ sim_bin02_idx_ss ])
endfor

sim_bin03_var_ss_01 = fltarr(n_elements(bin03_ss[0 , *]))
for i = 0, n_elements(bin03_ss[0 , *]) - 1 do begin
	for t = 0, nt - 1 do begin
		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
		gridinfo= gridinfo, lat= [bin03_ss[0, i] - 0.1, bin03_ss[0, i] + 0.1], $
		lon= [bin03_ss[1, i] - 0.1, bin03_ss[1, i] + 0.1], alrange = [1, 3], average=7)
		globalavg[ t ]= data/2 * 1000
	endfor
	sim_bin03_var_ss_01[ i ] = mean(globalavg[ sim_bin03_idx_ss ])
endfor





;extracting the simulated data for the corresponding months and years and locations
sim_bin01_var_ogi_01 = fltarr(n_elements(bin01_ogi[0 , *]))
for i = 0, n_elements(bin01_ogi[0 , *]) - 1 do begin
	for t = 0, nt - 1 do begin
		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
		gridinfo= gridinfo, lat= [bin01_ogi[0, i] - 0.1, bin01_ogi[0, i] + 0.1], $
		lon= [bin01_ogi[2, i] - 0.1, bin01_ogi[2, i] + 0.1], alrange = [1, 3], average=7)
		globalavg[ t ]= data/2 * 1000
	endfor
	sim_bin01_var_ogi_01[ i ] = mean(globalavg[ sim_bin01_idx_ogi ])
endfor






;sim_title = 'Constant default base emissions '
filename = "/home/excluded-from-backup/data/C2H6/trac_avg.spinup_GFED4_MAVG_1981_2015.198101010000"

;read in the file
ctm_get_data, datainfo, 'IJ-AVG-$', tracer=1, filename=filename, tau0=tau0

;Get MODELINFO and GRIDINFO structures, xmid, ymid hold lon/lat centers 
getmodelandgridinfo, datainfo[0], modelinfo, gridinfo, lon=xmid, lat=ymid
if (n_elements(datainfo) eq 0) then print, 'No data records found.'

;extracting the simulated data for the corresponding months and years and locations
sim_bin01_var_noaa_02 = fltarr(n_elements(bin01_noaa[0 , *]))
for i = 0, n_elements(bin01_noaa[0 , *]) - 1 do begin
	for t = 0, nt - 1 do begin
		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
		gridinfo= gridinfo, lat= [bin01_noaa[0, i] - 0.1, bin01_noaa[0, i] + 0.1], $
		lon= [bin01_noaa[1, i] - 0.1, bin01_noaa[1, i] + 0.1], alrange = [1, 3], average=7)
		globalavg[ t ]= data/2 * 1000
	endfor
	sim_bin01_var_noaa_02[ i ] = mean(globalavg[ sim_bin01_idx_noaa ])
endfor

sim_bin02_var_noaa_02 = fltarr(n_elements(bin02_noaa[0 , *]))
for i = 0, n_elements(bin02_noaa[0 , *]) - 1 do begin
	for t = 0, nt - 1 do begin
		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
		gridinfo= gridinfo, lat= [bin02_noaa[0, i] - 0.1, bin02_noaa[0, i] + 0.1], $
		lon= [bin02_noaa[1, i] - 0.1, bin02_noaa[1, i] + 0.1], alrange = [1, 3], average=7)
		globalavg[ t ]= data/2 * 1000
	endfor
	sim_bin02_var_noaa_02[ i ] = mean(globalavg[ sim_bin02_idx_noaa ])
endfor

sim_bin03_var_noaa_02 = fltarr(n_elements(bin03_noaa[0 , *]))
for i = 0, n_elements(bin03_noaa[0 , *]) - 1 do begin
	for t = 0, nt - 1 do begin
		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
		gridinfo= gridinfo, lat= [bin03_noaa[0, i] - 0.1, bin03_noaa[0, i] + 0.1], $
		lon= [bin03_noaa[1, i] - 0.1, bin03_noaa[1, i] + 0.1], alrange = [1, 3], average=7)
		globalavg[ t ]= data/2 * 1000
	endfor
	sim_bin03_var_noaa_02[ i ] = mean(globalavg[ sim_bin03_idx_noaa ])
endfor

sim_bin04_var_noaa_02 = fltarr(n_elements(bin04_noaa[0 , *]))
for i = 0, n_elements(bin04_noaa[0 , *]) - 1 do begin
	for t = 0, nt - 1 do begin
		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
		gridinfo= gridinfo, lat= [bin04_noaa[0, i] - 0.1, bin04_noaa[0, i] + 0.1], $
		lon= [bin04_noaa[1, i] - 0.1, bin04_noaa[1, i] + 0.1], alrange = [1, 3], average=7)
		globalavg[ t ]= data/2 * 1000
	endfor
	sim_bin04_var_noaa_02[ i ] = mean(globalavg[ sim_bin04_idx_noaa ])
endfor

sim_bin05_var_noaa_02 = fltarr(n_elements(bin05_noaa[0 , *]))
for i = 0, n_elements(bin05_noaa[0 , *]) - 1 do begin
	for t = 0, nt - 1 do begin
		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
		gridinfo= gridinfo, lat= [bin05_noaa[0, i] - 0.1, bin05_noaa[0, i] + 0.1], $
		lon= [bin05_noaa[1, i] - 0.1, bin05_noaa[1, i] + 0.1], alrange = [1, 3], average=7)
		globalavg[ t ]= data/2 * 1000
	endfor
	sim_bin05_var_noaa_02[ i ] = mean(globalavg[ sim_bin05_idx_noaa ])
endfor

sim_bin06_var_noaa_02 = fltarr(n_elements(bin06_noaa[0 , *]))
for i = 0, n_elements(bin06_noaa[0 , *]) - 1 do begin
	for t = 0, nt - 1 do begin
		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
		gridinfo= gridinfo, lat= [bin06_noaa[0, i] - 0.1, bin06_noaa[0, i] + 0.1], $
		lon= [bin06_noaa[1, i] - 0.1, bin06_noaa[1, i] + 0.1], alrange = [1, 3], average=7)
		globalavg[ t ]= data/2 * 1000
	endfor
	sim_bin06_var_noaa_02[ i ] = mean(globalavg[ sim_bin06_idx_noaa ])
endfor





;extracting the simulated data for the corresponding months and years and locations
sim_bin01_var_ss_02 = fltarr(n_elements(bin01_ss[0 , *]))
for i = 0, n_elements(bin01_ss[0 , *]) - 1 do begin
	for t = 0, nt - 1 do begin
		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
		gridinfo= gridinfo, lat= [bin01_ss[0, i] - 0.1, bin01_ss[0, i] + 0.1], $
		lon= [bin01_ss[1, i] - 0.1, bin01_ss[1, i] + 0.1], alrange = [1, 3], average=7)
		globalavg[ t ]= data/2 * 1000
	endfor
	sim_bin01_var_ss_02[ i ] = mean(globalavg[ sim_bin01_idx_ss ])
endfor

sim_bin02_var_ss_02 = fltarr(n_elements(bin02_ss[0 , *]))
for i = 0, n_elements(bin02_ss[0 , *]) - 1 do begin
	for t = 0, nt - 1 do begin
		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
		gridinfo= gridinfo, lat= [bin02_ss[0, i] - 0.1, bin02_ss[0, i] + 0.1], $
		lon= [bin02_ss[1, i] - 0.1, bin02_ss[1, i] + 0.1], alrange = [1, 3], average=7)
		globalavg[ t ]= data/2 * 1000
	endfor
	sim_bin02_var_ss_02[ i ] = mean(globalavg[ sim_bin02_idx_ss ])
endfor

sim_bin03_var_ss_02 = fltarr(n_elements(bin03_ss[0 , *]))
for i = 0, n_elements(bin03_ss[0 , *]) - 1 do begin
	for t = 0, nt - 1 do begin
		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
		gridinfo= gridinfo, lat= [bin03_ss[0, i] - 0.1, bin03_ss[0, i] + 0.1], $
		lon= [bin03_ss[1, i] - 0.1, bin03_ss[1, i] + 0.1], alrange = [1, 3], average=7)
		globalavg[ t ]= data/2 * 1000
	endfor
	sim_bin03_var_ss_02[ i ] = mean(globalavg[ sim_bin03_idx_ss ])
endfor





;extracting the simulated data for the corresponding months and years and locations
sim_bin01_var_ogi_02 = fltarr(n_elements(bin01_ogi[0 , *]))
for i = 0, n_elements(bin01_ogi[0 , *]) - 1 do begin
	for t = 0, nt - 1 do begin
		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
		gridinfo= gridinfo, lat= [bin01_ogi[0, i] - 0.1, bin01_ogi[0, i] + 0.1], $
		lon= [bin01_ogi[2, i] - 0.1, bin01_ogi[2, i] + 0.1], alrange = [1, 3], average=7)
		globalavg[ t ]= data/2 * 1000
	endfor
	sim_bin01_var_ogi_02[ i ] = mean(globalavg[ sim_bin01_idx_ogi ])
endfor








;sim_title = ''Aydin et al. no normalization to Xiao et al.
filename = "/home/excluded-from-backup/data/C2H6/trac_avg.AydinAbsSF_1981_2015.bpch"

;read in the file
ctm_get_data, datainfo, 'IJ-AVG-$', tracer=1, filename=filename, tau0=tau0

;Get MODELINFO and GRIDINFO structures, xmid, ymid hold lon/lat centers 
getmodelandgridinfo, datainfo[0], modelinfo, gridinfo, lon=xmid, lat=ymid
if (n_elements(datainfo) eq 0) then print, 'No data records found.'

;extracting the simulated data for the corresponding months and years and locations
sim_bin01_var_noaa_03 = fltarr(n_elements(bin01_noaa[0 , *]))
for i = 0, n_elements(bin01_noaa[0 , *]) - 1 do begin
	for t = 0, nt - 1 do begin
		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
		gridinfo= gridinfo, lat= [bin01_noaa[0, i] - 0.1, bin01_noaa[0, i] + 0.1], $
		lon= [bin01_noaa[1, i] - 0.1, bin01_noaa[1, i] + 0.1], alrange = [1, 3], average=7)
		globalavg[ t ]= data/2 * 1000
	endfor
	sim_bin01_var_noaa_03[ i ] = mean(globalavg[ sim_bin01_idx_noaa ])
endfor

sim_bin02_var_noaa_03 = fltarr(n_elements(bin02_noaa[0 , *]))
for i = 0, n_elements(bin02_noaa[0 , *]) - 1 do begin
	for t = 0, nt - 1 do begin
		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
		gridinfo= gridinfo, lat= [bin02_noaa[0, i] - 0.1, bin02_noaa[0, i] + 0.1], $
		lon= [bin02_noaa[1, i] - 0.1, bin02_noaa[1, i] + 0.1], alrange = [1, 3], average=7)
		globalavg[ t ]= data/2 * 1000
	endfor
	sim_bin02_var_noaa_03[ i ] = mean(globalavg[ sim_bin02_idx_noaa ])
endfor

sim_bin03_var_noaa_03 = fltarr(n_elements(bin03_noaa[0 , *]))
for i = 0, n_elements(bin03_noaa[0 , *]) - 1 do begin
	for t = 0, nt - 1 do begin
		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
		gridinfo= gridinfo, lat= [bin03_noaa[0, i] - 0.1, bin03_noaa[0, i] + 0.1], $
		lon= [bin03_noaa[1, i] - 0.1, bin03_noaa[1, i] + 0.1], alrange = [1, 3], average=7)
		globalavg[ t ]= data/2 * 1000
	endfor
	sim_bin03_var_noaa_03[ i ] = mean(globalavg[ sim_bin03_idx_noaa ])
endfor

sim_bin04_var_noaa_03 = fltarr(n_elements(bin04_noaa[0 , *]))
for i = 0, n_elements(bin04_noaa[0 , *]) - 1 do begin
	for t = 0, nt - 1 do begin
		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
		gridinfo= gridinfo, lat= [bin04_noaa[0, i] - 0.1, bin04_noaa[0, i] + 0.1], $
		lon= [bin04_noaa[1, i] - 0.1, bin04_noaa[1, i] + 0.1], alrange = [1, 3], average=7)
		globalavg[ t ]= data/2 * 1000
	endfor
	sim_bin04_var_noaa_03[ i ] = mean(globalavg[ sim_bin04_idx_noaa ])
endfor

sim_bin05_var_noaa_03 = fltarr(n_elements(bin05_noaa[0 , *]))
for i = 0, n_elements(bin05_noaa[0 , *]) - 1 do begin
	for t = 0, nt - 1 do begin
		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
		gridinfo= gridinfo, lat= [bin05_noaa[0, i] - 0.1, bin05_noaa[0, i] + 0.1], $
		lon= [bin05_noaa[1, i] - 0.1, bin05_noaa[1, i] + 0.1], alrange = [1, 3], average=7)
		globalavg[ t ]= data/2 * 1000
	endfor
	sim_bin05_var_noaa_03[ i ] = mean(globalavg[ sim_bin05_idx_noaa ])
endfor

sim_bin06_var_noaa_03 = fltarr(n_elements(bin06_noaa[0 , *]))
for i = 0, n_elements(bin06_noaa[0 , *]) - 1 do begin
	for t = 0, nt - 1 do begin
		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
		gridinfo= gridinfo, lat= [bin06_noaa[0, i] - 0.1, bin06_noaa[0, i] + 0.1], $
		lon= [bin06_noaa[1, i] - 0.1, bin06_noaa[1, i] + 0.1], alrange = [1, 3], average=7)
		globalavg[ t ]= data/2 * 1000
	endfor
	sim_bin06_var_noaa_03[ i ] = mean(globalavg[ sim_bin06_idx_noaa ])
endfor





;extracting the simulated data for the corresponding months and years and locations
sim_bin01_var_ss_03 = fltarr(n_elements(bin01_ss[0 , *]))
for i = 0, n_elements(bin01_ss[0 , *]) - 1 do begin
	for t = 0, nt - 1 do begin
		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
		gridinfo= gridinfo, lat= [bin01_ss[0, i] - 0.1, bin01_ss[0, i] + 0.1], $
		lon= [bin01_ss[1, i] - 0.1, bin01_ss[1, i] + 0.1], alrange = [1, 3], average=7)
		globalavg[ t ]= data/2 * 1000
	endfor
	sim_bin01_var_ss_03[ i ] = mean(globalavg[ sim_bin01_idx_ss ])
endfor

sim_bin02_var_ss_03 = fltarr(n_elements(bin02_ss[0 , *]))
for i = 0, n_elements(bin02_ss[0 , *]) - 1 do begin
	for t = 0, nt - 1 do begin
		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
		gridinfo= gridinfo, lat= [bin02_ss[0, i] - 0.1, bin02_ss[0, i] + 0.1], $
		lon= [bin02_ss[1, i] - 0.1, bin02_ss[1, i] + 0.1], alrange = [1, 3], average=7)
		globalavg[ t ]= data/2 * 1000
	endfor
	sim_bin02_var_ss_03[ i ] = mean(globalavg[ sim_bin02_idx_ss ])
endfor

sim_bin03_var_ss_03 = fltarr(n_elements(bin03_ss[0 , *]))
for i = 0, n_elements(bin03_ss[0 , *]) - 1 do begin
	for t = 0, nt - 1 do begin
		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
		gridinfo= gridinfo, lat= [bin03_ss[0, i] - 0.1, bin03_ss[0, i] + 0.1], $
		lon= [bin03_ss[1, i] - 0.1, bin03_ss[1, i] + 0.1], alrange = [1, 3], average=7)
		globalavg[ t ]= data/2 * 1000
	endfor
	sim_bin03_var_ss_03[ i ] = mean(globalavg[ sim_bin03_idx_ss ])
endfor





;extracting the simulated data for the corresponding months and years and locations
sim_bin01_var_ogi_03 = fltarr(n_elements(bin01_ogi[0 , *]))
for i = 0, n_elements(bin01_ogi[0 , *]) - 1 do begin
	for t = 0, nt - 1 do begin
		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
		gridinfo= gridinfo, lat= [bin01_ogi[0, i] - 0.1, bin01_ogi[0, i] + 0.1], $
		lon= [bin01_ogi[2, i] - 0.1, bin01_ogi[2, i] + 0.1], alrange = [1, 3], average=7)
		globalavg[ t ]= data/2 * 1000
	endfor
	sim_bin01_var_ogi_03[ i ] = mean(globalavg[ sim_bin01_idx_ogi ])
endfor





;sim_title = ''Aydin et al. no normalization to Xiao et al.
filename = "/home/excluded-from-backup/data/C2H6/trac_avg.SimpsonSF_1981_2015.bpch"

;read in the file
ctm_get_data, datainfo, 'IJ-AVG-$', tracer=1, filename=filename, tau0=tau0

;Get MODELINFO and GRIDINFO structures, xmid, ymid hold lon/lat centers 
getmodelandgridinfo, datainfo[0], modelinfo, gridinfo, lon=xmid, lat=ymid
if (n_elements(datainfo) eq 0) then print, 'No data records found.'

;extracting the simulated data for the corresponding months and years and locations
sim_bin01_var_noaa_04 = fltarr(n_elements(bin01_noaa[0 , *]))
for i = 0, n_elements(bin01_noaa[0 , *]) - 1 do begin
	for t = 0, nt - 1 do begin
		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
		gridinfo= gridinfo, lat= [bin01_noaa[0, i] - 0.1, bin01_noaa[0, i] + 0.1], $
		lon= [bin01_noaa[1, i] - 0.1, bin01_noaa[1, i] + 0.1], alrange = [1, 3], average=7)
		globalavg[ t ]= data/2 * 1000
	endfor
	sim_bin01_var_noaa_04[ i ] = mean(globalavg[ sim_bin01_idx_noaa ])
endfor

sim_bin02_var_noaa_04 = fltarr(n_elements(bin02_noaa[0 , *]))
for i = 0, n_elements(bin02_noaa[0 , *]) - 1 do begin
	for t = 0, nt - 1 do begin
		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
		gridinfo= gridinfo, lat= [bin02_noaa[0, i] - 0.1, bin02_noaa[0, i] + 0.1], $
		lon= [bin02_noaa[1, i] - 0.1, bin02_noaa[1, i] + 0.1], alrange = [1, 3], average=7)
		globalavg[ t ]= data/2 * 1000
	endfor
	sim_bin02_var_noaa_04[ i ] = mean(globalavg[ sim_bin02_idx_noaa ])
endfor

sim_bin03_var_noaa_04 = fltarr(n_elements(bin03_noaa[0 , *]))
for i = 0, n_elements(bin03_noaa[0 , *]) - 1 do begin
	for t = 0, nt - 1 do begin
		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
		gridinfo= gridinfo, lat= [bin03_noaa[0, i] - 0.1, bin03_noaa[0, i] + 0.1], $
		lon= [bin03_noaa[1, i] - 0.1, bin03_noaa[1, i] + 0.1], alrange = [1, 3], average=7)
		globalavg[ t ]= data/2 * 1000
	endfor
	sim_bin03_var_noaa_04[ i ] = mean(globalavg[ sim_bin03_idx_noaa ])
endfor

sim_bin04_var_noaa_04 = fltarr(n_elements(bin04_noaa[0 , *]))
for i = 0, n_elements(bin04_noaa[0 , *]) - 1 do begin
	for t = 0, nt - 1 do begin
		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
		gridinfo= gridinfo, lat= [bin04_noaa[0, i] - 0.1, bin04_noaa[0, i] + 0.1], $
		lon= [bin04_noaa[1, i] - 0.1, bin04_noaa[1, i] + 0.1], alrange = [1, 3], average=7)
		globalavg[ t ]= data/2 * 1000
	endfor
	sim_bin04_var_noaa_04[ i ] = mean(globalavg[ sim_bin04_idx_noaa ])
endfor

sim_bin05_var_noaa_04 = fltarr(n_elements(bin05_noaa[0 , *]))
for i = 0, n_elements(bin05_noaa[0 , *]) - 1 do begin
	for t = 0, nt - 1 do begin
		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
		gridinfo= gridinfo, lat= [bin05_noaa[0, i] - 0.1, bin05_noaa[0, i] + 0.1], $
		lon= [bin05_noaa[1, i] - 0.1, bin05_noaa[1, i] + 0.1], alrange = [1, 3], average=7)
		globalavg[ t ]= data/2 * 1000
	endfor
	sim_bin05_var_noaa_04[ i ] = mean(globalavg[ sim_bin05_idx_noaa ])
endfor

sim_bin06_var_noaa_04 = fltarr(n_elements(bin06_noaa[0 , *]))
for i = 0, n_elements(bin06_noaa[0 , *]) - 1 do begin
	for t = 0, nt - 1 do begin
		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
		gridinfo= gridinfo, lat= [bin06_noaa[0, i] - 0.1, bin06_noaa[0, i] + 0.1], $
		lon= [bin06_noaa[1, i] - 0.1, bin06_noaa[1, i] + 0.1], alrange = [1, 3], average=7)
		globalavg[ t ]= data/2 * 1000
	endfor
	sim_bin06_var_noaa_04[ i ] = mean(globalavg[ sim_bin06_idx_noaa ])
endfor





;extracting the simulated data for the corresponding months and years and locations
sim_bin01_var_ss_04 = fltarr(n_elements(bin01_ss[0 , *]))
for i = 0, n_elements(bin01_ss[0 , *]) - 1 do begin
	for t = 0, nt - 1 do begin
		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
		gridinfo= gridinfo, lat= [bin01_ss[0, i] - 0.1, bin01_ss[0, i] + 0.1], $
		lon= [bin01_ss[1, i] - 0.1, bin01_ss[1, i] + 0.1], alrange = [1, 3], average=7)
		globalavg[ t ]= data/2 * 1000
	endfor
	sim_bin01_var_ss_04[ i ] = mean(globalavg[ sim_bin01_idx_ss ])
endfor

sim_bin02_var_ss_04 = fltarr(n_elements(bin02_ss[0 , *]))
for i = 0, n_elements(bin02_ss[0 , *]) - 1 do begin
	for t = 0, nt - 1 do begin
		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
		gridinfo= gridinfo, lat= [bin02_ss[0, i] - 0.1, bin02_ss[0, i] + 0.1], $
		lon= [bin02_ss[1, i] - 0.1, bin02_ss[1, i] + 0.1], alrange = [1, 3], average=7)
		globalavg[ t ]= data/2 * 1000
	endfor
	sim_bin02_var_ss_04[ i ] = mean(globalavg[ sim_bin02_idx_ss ])
endfor

sim_bin03_var_ss_04 = fltarr(n_elements(bin03_ss[0 , *]))
for i = 0, n_elements(bin03_ss[0 , *]) - 1 do begin
	for t = 0, nt - 1 do begin
		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
		gridinfo= gridinfo, lat= [bin03_ss[0, i] - 0.1, bin03_ss[0, i] + 0.1], $
		lon= [bin03_ss[1, i] - 0.1, bin03_ss[1, i] + 0.1], alrange = [1, 3], average=7)
		globalavg[ t ]= data/2 * 1000
	endfor
	sim_bin03_var_ss_04[ i ] = mean(globalavg[ sim_bin03_idx_ss ])
endfor





;extracting the simulated data for the corresponding months and years and locations
sim_bin01_var_ogi_04 = fltarr(n_elements(bin01_ogi[0 , *]))
for i = 0, n_elements(bin01_ogi[0 , *]) - 1 do begin
	for t = 0, nt - 1 do begin
		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
		gridinfo= gridinfo, lat= [bin01_ogi[0, i] - 0.1, bin01_ogi[0, i] + 0.1], $
		lon= [bin01_ogi[2, i] - 0.1, bin01_ogi[2, i] + 0.1], alrange = [1, 3], average=7)
		globalavg[ t ]= data/2 * 1000
	endfor
	sim_bin01_var_ogi_04[ i ] = mean(globalavg[ sim_bin01_idx_ogi ])
endfor






;sim_title = ''Aydin et al. no normalization to Xiao et al.
filename = "/home/excluded-from-backup/data/C2H6/trac_avg.PSU_MER3BB_MER18FF_1981_2015.bpch"

;read in the file
ctm_get_data, datainfo, 'IJ-AVG-$', tracer=1, filename=filename, tau0=tau0

;Get MODELINFO and GRIDINFO structures, xmid, ymid hold lon/lat centers 
getmodelandgridinfo, datainfo[0], modelinfo, gridinfo, lon=xmid, lat=ymid
if (n_elements(datainfo) eq 0) then print, 'No data records found.'

;extracting the simulated data for the corresponding months and years and locations
sim_bin01_var_noaa_05 = fltarr(n_elements(bin01_noaa[0 , *]))
for i = 0, n_elements(bin01_noaa[0 , *]) - 1 do begin
	for t = 0, nt - 1 do begin
		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
		gridinfo= gridinfo, lat= [bin01_noaa[0, i] - 0.1, bin01_noaa[0, i] + 0.1], $
		lon= [bin01_noaa[1, i] - 0.1, bin01_noaa[1, i] + 0.1], alrange = [1, 3], average=7)
		globalavg[ t ]= data/2 * 1000
	endfor
	sim_bin01_var_noaa_05[ i ] = mean(globalavg[ sim_bin01_idx_noaa ])
endfor

sim_bin02_var_noaa_05 = fltarr(n_elements(bin02_noaa[0 , *]))
for i = 0, n_elements(bin02_noaa[0 , *]) - 1 do begin
	for t = 0, nt - 1 do begin
		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
		gridinfo= gridinfo, lat= [bin02_noaa[0, i] - 0.1, bin02_noaa[0, i] + 0.1], $
		lon= [bin02_noaa[1, i] - 0.1, bin02_noaa[1, i] + 0.1], alrange = [1, 3], average=7)
		globalavg[ t ]= data/2 * 1000
	endfor
	sim_bin02_var_noaa_05[ i ] = mean(globalavg[ sim_bin02_idx_noaa ])
endfor

sim_bin03_var_noaa_05 = fltarr(n_elements(bin03_noaa[0 , *]))
for i = 0, n_elements(bin03_noaa[0 , *]) - 1 do begin
	for t = 0, nt - 1 do begin
		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
		gridinfo= gridinfo, lat= [bin03_noaa[0, i] - 0.1, bin03_noaa[0, i] + 0.1], $
		lon= [bin03_noaa[1, i] - 0.1, bin03_noaa[1, i] + 0.1], alrange = [1, 3], average=7)
		globalavg[ t ]= data/2 * 1000
	endfor
	sim_bin03_var_noaa_05[ i ] = mean(globalavg[ sim_bin03_idx_noaa ])
endfor

sim_bin04_var_noaa_05 = fltarr(n_elements(bin04_noaa[0 , *]))
for i = 0, n_elements(bin04_noaa[0 , *]) - 1 do begin
	for t = 0, nt - 1 do begin
		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
		gridinfo= gridinfo, lat= [bin04_noaa[0, i] - 0.1, bin04_noaa[0, i] + 0.1], $
		lon= [bin04_noaa[1, i] - 0.1, bin04_noaa[1, i] + 0.1], alrange = [1, 3], average=7)
		globalavg[ t ]= data/2 * 1000
	endfor
	sim_bin04_var_noaa_05[ i ] = mean(globalavg[ sim_bin04_idx_noaa ])
endfor

sim_bin05_var_noaa_05 = fltarr(n_elements(bin05_noaa[0 , *]))
for i = 0, n_elements(bin05_noaa[0 , *]) - 1 do begin
	for t = 0, nt - 1 do begin
		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
		gridinfo= gridinfo, lat= [bin05_noaa[0, i] - 0.1, bin05_noaa[0, i] + 0.1], $
		lon= [bin05_noaa[1, i] - 0.1, bin05_noaa[1, i] + 0.1], alrange = [1, 3], average=7)
		globalavg[ t ]= data/2 * 1000
	endfor
	sim_bin05_var_noaa_05[ i ] = mean(globalavg[ sim_bin05_idx_noaa ])
endfor

sim_bin06_var_noaa_05 = fltarr(n_elements(bin06_noaa[0 , *]))
for i = 0, n_elements(bin06_noaa[0 , *]) - 1 do begin
	for t = 0, nt - 1 do begin
		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
		gridinfo= gridinfo, lat= [bin06_noaa[0, i] - 0.1, bin06_noaa[0, i] + 0.1], $
		lon= [bin06_noaa[1, i] - 0.1, bin06_noaa[1, i] + 0.1], alrange = [1, 3], average=7)
		globalavg[ t ]= data/2 * 1000
	endfor
	sim_bin06_var_noaa_05[ i ] = mean(globalavg[ sim_bin06_idx_noaa ])
endfor





;extracting the simulated data for the corresponding months and years and locations
sim_bin01_var_ss_05 = fltarr(n_elements(bin01_ss[0 , *]))
for i = 0, n_elements(bin01_ss[0 , *]) - 1 do begin
	for t = 0, nt - 1 do begin
		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
		gridinfo= gridinfo, lat= [bin01_ss[0, i] - 0.1, bin01_ss[0, i] + 0.1], $
		lon= [bin01_ss[1, i] - 0.1, bin01_ss[1, i] + 0.1], alrange = [1, 3], average=7)
		globalavg[ t ]= data/2 * 1000
	endfor
	sim_bin01_var_ss_05[ i ] = mean(globalavg[ sim_bin01_idx_ss ])
endfor

sim_bin02_var_ss_05 = fltarr(n_elements(bin02_ss[0 , *]))
for i = 0, n_elements(bin02_ss[0 , *]) - 1 do begin
	for t = 0, nt - 1 do begin
		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
		gridinfo= gridinfo, lat= [bin02_ss[0, i] - 0.1, bin02_ss[0, i] + 0.1], $
		lon= [bin02_ss[1, i] - 0.1, bin02_ss[1, i] + 0.1], alrange = [1, 3], average=7)
		globalavg[ t ]= data/2 * 1000
	endfor
	sim_bin02_var_ss_05[ i ] = mean(globalavg[ sim_bin02_idx_ss ])
endfor

sim_bin03_var_ss_05 = fltarr(n_elements(bin03_ss[0 , *]))
for i = 0, n_elements(bin03_ss[0 , *]) - 1 do begin
	for t = 0, nt - 1 do begin
		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
		gridinfo= gridinfo, lat= [bin03_ss[0, i] - 0.1, bin03_ss[0, i] + 0.1], $
		lon= [bin03_ss[1, i] - 0.1, bin03_ss[1, i] + 0.1], alrange = [1, 3], average=7)
		globalavg[ t ]= data/2 * 1000
	endfor
	sim_bin03_var_ss_05[ i ] = mean(globalavg[ sim_bin03_idx_ss ])
endfor





;extracting the simulated data for the corresponding months and years and locations
sim_bin01_var_ogi_05 = fltarr(n_elements(bin01_ogi[0 , *]))
for i = 0, n_elements(bin01_ogi[0 , *]) - 1 do begin
	for t = 0, nt - 1 do begin
		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
		gridinfo= gridinfo, lat= [bin01_ogi[0, i] - 0.1, bin01_ogi[0, i] + 0.1], $
		lon= [bin01_ogi[2, i] - 0.1, bin01_ogi[2, i] + 0.1], alrange = [1, 3], average=7)
		globalavg[ t ]= data/2 * 1000
	endfor
	sim_bin01_var_ogi_05[ i ] = mean(globalavg[ sim_bin01_idx_ogi ])
endfor





;This section turns the latitudes to sine values

bin01_ss[ 0 , * ] = sin( bin01_ss[ 0 , * ] * !pi / 180 )
bin02_ss[ 0 , * ] = sin( bin02_ss[ 0 , * ] * !pi / 180 )
bin03_ss[ 0 , * ] = sin( bin03_ss[ 0 , * ] * !pi / 180 )

bin01_noaa[ 0 , * ] = sin( bin01_noaa[ 0 , * ] * !pi / 180 )
bin02_noaa[ 0 , * ] = sin( bin02_noaa[ 0 , * ] * !pi / 180 )
bin03_noaa[ 0 , * ] = sin( bin03_noaa[ 0 , * ] * !pi / 180 )
bin04_noaa[ 0 , * ] = sin( bin04_noaa[ 0 , * ] * !pi / 180 )
bin05_noaa[ 0 , * ] = sin( bin05_noaa[ 0 , * ] * !pi / 180 )
bin06_noaa[ 0 , * ] = sin( bin06_noaa[ 0 , * ] * !pi / 180 )

bin01_ogi[ 0 , * ] = sin( bin01_ogi[ 0 , * ] * !pi / 180 )

;begin doing poly fit for the latitudinal distribution
measure_errors = de_obs_bin01_ogi[ 1 , * ]
bin01_ogi_pfit = poly_fit( bin01_ogi [ 0 , * ] , bin01_ogi [ 1 , * ], 3, measure_errors = measure_errors, $
	/double, yfit = yfit_bin01_ogi )

measure_errors = de_obs_bin01_ss[ 1 , * ]
bin01_ss_pfit = poly_fit( bin01_ss [ 0 , * ] , bin01_ss [ 2 , * ], 6, measure_errors = measure_errors, $
	/double, yfit = yfit_bin01_ss )
measure_errors = de_obs_bin02_ss[ 1 , * ]
bin02_ss_pfit = poly_fit( bin02_ss [ 0 , * ] , bin02_ss [ 2 , * ], 6, measure_errors = measure_errors, $
	/double, yfit = yfit_bin02_ss ) 
measure_errors = de_obs_bin03_ss[ 1 , * ]
bin03_ss_pfit = poly_fit( bin03_ss [ 0 , * ] , bin03_ss [ 2 , * ], 6, measure_errors = measure_errors, $
	/double, yfit = yfit_bin03_ss ) 


measure_errors = de_obs_bin01_noaa[ 2 , sort( bin01_noaa[ 0 , * ] ) ]
bin01_noaa_pfit = poly_fit( bin01_noaa [ 0 , sort( bin01_noaa[ 0 , * ] ) ] , $
	bin01_noaa [ 2 , sort( bin01_noaa[ 0 , * ] ) ], 6, measure_errors = measure_errors, $
	/double, yfit = yfit_bin01_noaa )
measure_errors = de_obs_bin02_noaa[ 2 , sort( bin02_noaa[ 0 , * ] ) ]
bin02_noaa_pfit = poly_fit( bin02_noaa [ 0 , sort( bin02_noaa[ 0 , * ] ) ] , $
	bin02_noaa [ 2 , sort( bin02_noaa[ 0 , * ] ) ], 6, measure_errors = measure_errors, $
	/double, yfit = yfit_bin02_noaa ) 
measure_errors = de_obs_bin03_noaa[ 2 , sort( bin03_noaa[ 0 , * ] ) ]
bin03_noaa_pfit = poly_fit( bin03_noaa [ 0 , sort( bin03_noaa[ 0 , * ] ) ] , $
	bin03_noaa [ 2 , sort( bin03_noaa[ 0 , * ] ) ], 6, measure_errors = measure_errors, $
	/double, yfit = yfit_bin03_noaa ) 
measure_errors = de_obs_bin04_noaa[ 2 , sort( bin04_noaa[ 0 , * ] ) ]
bin04_noaa_pfit = poly_fit( bin04_noaa [ 0 , sort( bin04_noaa[ 0 , * ] ) ] , $
	bin04_noaa [ 2 , sort( bin04_noaa[ 0 , * ] ) ], 6, measure_errors = measure_errors, $
	/double, yfit = yfit_bin04_noaa ) 
measure_errors = de_obs_bin05_noaa[ 2 , sort( bin05_noaa[ 0 , * ] ) ]
bin05_noaa_pfit = poly_fit( bin05_noaa [ 0 , sort( bin05_noaa[ 0 , * ] ) ] , $
	bin05_noaa [ 2 , sort( bin05_noaa[ 0 , * ] ) ], 6, measure_errors = measure_errors, $
	/double, yfit = yfit_bin05_noaa ) 
measure_errors = de_obs_bin06_noaa[ 2 , sort( bin06_noaa[ 0 , * ] ) ]
bin06_noaa_pfit = poly_fit( bin06_noaa [ 0 , sort( bin06_noaa[ 0 , * ] ) ] , $
	bin06_noaa [ 2 , sort( bin06_noaa[ 0 , * ] ) ], 6, measure_errors = measure_errors, $
	/double, yfit = yfit_bin06_noaa ) 

;this section does poly fit for the simulated data
bin01_noaa_sim_pfit = poly_fit( bin01_noaa[ 0 , sort( bin01_noaa[ 0 , * ] ) ], $
	sim_bin01_var_noaa[ sort( bin01_noaa[ 0 , * ] ) ], 6, yband = yband_sim01_noaa, yfit = yfit_sim01_noaa )
bin02_noaa_sim_pfit = poly_fit( bin02_noaa[ 0 , sort( bin02_noaa[ 0 , * ] ) ], $
	sim_bin02_var_noaa[ sort( bin02_noaa[ 0 , * ] ) ], 6, yband = yband_sim02_noaa, yfit = yfit_sim02_noaa )
bin03_noaa_sim_pfit = poly_fit( bin03_noaa[ 0 , sort( bin03_noaa[ 0 , * ] ) ], $
	sim_bin03_var_noaa[ sort( bin03_noaa[ 0 , * ] ) ], 6, yband = yband_sim03_noaa, yfit = yfit_sim03_noaa )
bin04_noaa_sim_pfit = poly_fit( bin04_noaa[ 0 , sort( bin04_noaa[ 0 , * ] ) ], $
	sim_bin04_var_noaa[ sort( bin04_noaa[ 0 , * ] ) ], 6, yband = yband_sim04_noaa, yfit = yfit_sim04_noaa )
bin05_noaa_sim_pfit = poly_fit( bin05_noaa[ 0 , sort( bin05_noaa[ 0 , * ] ) ], $
	sim_bin05_var_noaa[ sort( bin05_noaa[ 0 , * ] ) ], 6, yband = yband_sim05_noaa, yfit = yfit_sim05_noaa )
bin06_noaa_sim_pfit = poly_fit( bin06_noaa[ 0 , sort( bin06_noaa[ 0 , * ] ) ], $
	sim_bin06_var_noaa[ sort( bin06_noaa[ 0 , * ] ) ], 6, yband = yband_sim06_noaa, yfit = yfit_sim06_noaa )
	
	
bin01_noaa_sim_pfit_01 = poly_fit( bin01_noaa[ 0 , sort( bin01_noaa[ 0 , * ] ) ], $
	sim_bin01_var_noaa_01[ sort( bin01_noaa[ 0 , * ] ) ], 6, yband = yband_sim01_noaa_01, yfit = yfit_sim01_noaa_01 )
bin02_noaa_sim_pfit_01 = poly_fit( bin02_noaa[ 0 , sort( bin02_noaa[ 0 , * ] ) ], $
	sim_bin02_var_noaa_01[ sort( bin02_noaa[ 0 , * ] ) ], 6, yband = yband_sim02_noaa_01, yfit = yfit_sim02_noaa_01 )
bin03_noaa_sim_pfit_01 = poly_fit( bin03_noaa[ 0 , sort( bin03_noaa[ 0 , * ] ) ], $
	sim_bin03_var_noaa_01[ sort( bin03_noaa[ 0 , * ] ) ], 6, yband = yband_sim03_noaa_01, yfit = yfit_sim03_noaa_01 )
bin04_noaa_sim_pfit_01 = poly_fit( bin04_noaa[ 0 , sort( bin04_noaa[ 0 , * ] ) ], $
	sim_bin04_var_noaa_01[ sort( bin04_noaa[ 0 , * ] ) ], 6, yband = yband_sim04_noaa_01, yfit = yfit_sim04_noaa_01 )
bin05_noaa_sim_pfit_01 = poly_fit( bin05_noaa[ 0 , sort( bin05_noaa[ 0 , * ] ) ], $
	sim_bin05_var_noaa_01[ sort( bin05_noaa[ 0 , * ] ) ], 6, yband = yband_sim05_noaa_01, yfit = yfit_sim05_noaa_01 )
bin06_noaa_sim_pfit_01 = poly_fit( bin06_noaa[ 0 , sort( bin06_noaa[ 0 , * ] ) ], $
	sim_bin06_var_noaa_01[ sort( bin06_noaa[ 0 , * ] ) ], 6, yband = yband_sim06_noaa_01, yfit = yfit_sim06_noaa_01 )
	
	
bin01_noaa_sim_pfit_02 = poly_fit( bin01_noaa[ 0 , sort( bin01_noaa[ 0 , * ] ) ], $
	sim_bin01_var_noaa_02[ sort( bin01_noaa[ 0 , * ] ) ], 6, yband = yband_sim01_noaa_02, yfit = yfit_sim01_noaa_02 )
bin02_noaa_sim_pfit_02 = poly_fit( bin02_noaa[ 0 , sort( bin02_noaa[ 0 , * ] ) ], $
	sim_bin02_var_noaa_02[ sort( bin02_noaa[ 0 , * ] ) ], 6, yband = yband_sim02_noaa_02, yfit = yfit_sim02_noaa_02 )
bin03_noaa_sim_pfit_02 = poly_fit( bin03_noaa[ 0 , sort( bin03_noaa[ 0 , * ] ) ], $
	sim_bin03_var_noaa_02[ sort( bin03_noaa[ 0 , * ] ) ], 6, yband = yband_sim03_noaa_02, yfit = yfit_sim03_noaa_02 )
bin04_noaa_sim_pfit_02 = poly_fit( bin04_noaa[ 0 , sort( bin04_noaa[ 0 , * ] ) ], $
	sim_bin04_var_noaa_02[ sort( bin04_noaa[ 0 , * ] ) ], 6, yband = yband_sim04_noaa_02, yfit = yfit_sim04_noaa_02 )
bin05_noaa_sim_pfit_02 = poly_fit( bin05_noaa[ 0 , sort( bin05_noaa[ 0 , * ] ) ], $
	sim_bin05_var_noaa_02[ sort( bin05_noaa[ 0 , * ] ) ], 6, yband = yband_sim05_noaa_02, yfit = yfit_sim05_noaa_02 )
bin06_noaa_sim_pfit_02 = poly_fit( bin06_noaa[ 0 , sort( bin06_noaa[ 0 , * ] ) ], $
	sim_bin06_var_noaa_02[ sort( bin06_noaa[ 0 , * ] ) ], 6, yband = yband_sim06_noaa_02, yfit = yfit_sim06_noaa_02 )
	

bin01_noaa_sim_pfit_03 = poly_fit( bin01_noaa[ 0 , sort( bin01_noaa[ 0 , * ] ) ], $
	sim_bin01_var_noaa_03[ sort( bin01_noaa[ 0 , * ] ) ], 6, yband = yband_sim01_noaa_03, yfit = yfit_sim01_noaa_03 )
bin02_noaa_sim_pfit_03 = poly_fit( bin02_noaa[ 0 , sort( bin02_noaa[ 0 , * ] ) ], $
	sim_bin02_var_noaa_03[ sort( bin02_noaa[ 0 , * ] ) ], 6, yband = yband_sim02_noaa_03, yfit = yfit_sim02_noaa_03 )
bin03_noaa_sim_pfit_03 = poly_fit( bin03_noaa[ 0 , sort( bin03_noaa[ 0 , * ] ) ], $
	sim_bin03_var_noaa_03[ sort( bin03_noaa[ 0 , * ] ) ], 6, yband = yband_sim03_noaa_03, yfit = yfit_sim03_noaa_03 )
bin04_noaa_sim_pfit_03 = poly_fit( bin04_noaa[ 0 , sort( bin04_noaa[ 0 , * ] ) ], $
	sim_bin04_var_noaa_03[ sort( bin04_noaa[ 0 , * ] ) ], 6, yband = yband_sim04_noaa_03, yfit = yfit_sim04_noaa_03 )
bin05_noaa_sim_pfit_03 = poly_fit( bin05_noaa[ 0 , sort( bin05_noaa[ 0 , * ] ) ], $
	sim_bin05_var_noaa_03[ sort( bin05_noaa[ 0 , * ] ) ], 6, yband = yband_sim05_noaa_03, yfit = yfit_sim05_noaa_03 )
bin06_noaa_sim_pfit_03 = poly_fit( bin06_noaa[ 0 , sort( bin06_noaa[ 0 , * ] ) ], $
	sim_bin06_var_noaa_03[ sort( bin06_noaa[ 0 , * ] ) ], 6, yband = yband_sim06_noaa_03, yfit = yfit_sim06_noaa_03 )
	
	
bin01_noaa_sim_pfit_04 = poly_fit( bin01_noaa[ 0 , sort( bin01_noaa[ 0 , * ] ) ], $
	sim_bin01_var_noaa_04[ sort( bin01_noaa[ 0 , * ] ) ], 6, yband = yband_sim01_noaa_04, yfit = yfit_sim01_noaa_04 )
bin02_noaa_sim_pfit_04 = poly_fit( bin02_noaa[ 0 , sort( bin02_noaa[ 0 , * ] ) ], $
	sim_bin02_var_noaa_04[ sort( bin02_noaa[ 0 , * ] ) ], 6, yband = yband_sim02_noaa_04, yfit = yfit_sim02_noaa_04 )
bin03_noaa_sim_pfit_04 = poly_fit( bin03_noaa[ 0 , sort( bin03_noaa[ 0 , * ] ) ], $
	sim_bin03_var_noaa_04[ sort( bin03_noaa[ 0 , * ] ) ], 6, yband = yband_sim03_noaa_04, yfit = yfit_sim03_noaa_04 )
bin04_noaa_sim_pfit_04 = poly_fit( bin04_noaa[ 0 , sort( bin04_noaa[ 0 , * ] ) ], $
	sim_bin04_var_noaa_04[ sort( bin04_noaa[ 0 , * ] ) ], 6, yband = yband_sim04_noaa_04, yfit = yfit_sim04_noaa_04 )
bin05_noaa_sim_pfit_04 = poly_fit( bin05_noaa[ 0 , sort( bin05_noaa[ 0 , * ] ) ], $
	sim_bin05_var_noaa_04[ sort( bin05_noaa[ 0 , * ] ) ], 6, yband = yband_sim05_noaa_04, yfit = yfit_sim05_noaa_04 )
bin06_noaa_sim_pfit_04 = poly_fit( bin06_noaa[ 0 , sort( bin06_noaa[ 0 , * ] ) ], $
	sim_bin06_var_noaa_04[ sort( bin06_noaa[ 0 , * ] ) ], 6, yband = yband_sim06_noaa_04, yfit = yfit_sim06_noaa_04 )
	
	
bin01_noaa_sim_pfit_05 = poly_fit( bin01_noaa[ 0 , sort( bin01_noaa[ 0 , * ] ) ], $
	sim_bin01_var_noaa_05[ sort( bin01_noaa[ 0 , * ] ) ], 6, yband = yband_sim01_noaa_05, yfit = yfit_sim01_noaa_05 )
bin02_noaa_sim_pfit_05 = poly_fit( bin02_noaa[ 0 , sort( bin02_noaa[ 0 , * ] ) ], $
	sim_bin02_var_noaa_05[ sort( bin02_noaa[ 0 , * ] ) ], 6, yband = yband_sim02_noaa_05, yfit = yfit_sim02_noaa_05 )
bin03_noaa_sim_pfit_05 = poly_fit( bin03_noaa[ 0 , sort( bin03_noaa[ 0 , * ] ) ], $
	sim_bin03_var_noaa_05[ sort( bin03_noaa[ 0 , * ] ) ], 6, yband = yband_sim03_noaa_05, yfit = yfit_sim03_noaa_05 )
bin04_noaa_sim_pfit_05 = poly_fit( bin04_noaa[ 0 , sort( bin04_noaa[ 0 , * ] ) ], $
	sim_bin04_var_noaa_05[ sort( bin04_noaa[ 0 , * ] ) ], 6, yband = yband_sim04_noaa_05, yfit = yfit_sim04_noaa_05 )
bin05_noaa_sim_pfit_05 = poly_fit( bin05_noaa[ 0 , sort( bin05_noaa[ 0 , * ] ) ], $
	sim_bin05_var_noaa_05[ sort( bin05_noaa[ 0 , * ] ) ], 6, yband = yband_sim05_noaa_05, yfit = yfit_sim05_noaa_05 )
bin06_noaa_sim_pfit_05 = poly_fit( bin06_noaa[ 0 , sort( bin06_noaa[ 0 , * ] ) ], $
	sim_bin06_var_noaa_05[ sort( bin06_noaa[ 0 , * ] ) ], 6, yband = yband_sim06_noaa_05, yfit = yfit_sim06_noaa_05 )
	
	
	
bin01_ss_sim_pfit = poly_fit( bin01_ss[ 0 , * ], sim_bin01_var_ss, 6, yband = yband_sim01_ss, $
	yfit = yfit_sim01_ss )
bin02_ss_sim_pfit = poly_fit( bin02_ss[ 0 , * ], sim_bin02_var_ss, 6, yband = yband_sim02_ss, $
	yfit = yfit_sim02_ss )
bin03_ss_sim_pfit = poly_fit( bin03_ss[ 0 , * ], sim_bin03_var_ss, 6, yband = yband_sim03_ss, $
	yfit = yfit_sim03_ss )
	
	
bin01_ss_sim_pfit_01 = poly_fit( bin01_ss[ 0 , * ], sim_bin01_var_ss_01, 6, yband = yband_sim01_ss_01, $
	yfit = yfit_sim01_ss_01 )
bin02_ss_sim_pfit_01 = poly_fit( bin02_ss[ 0 , * ], sim_bin02_var_ss_01, 6, yband = yband_sim02_ss_01, $
	yfit = yfit_sim02_ss_01 )
bin03_ss_sim_pfit_01 = poly_fit( bin03_ss[ 0 , * ], sim_bin03_var_ss_01, 6, yband = yband_sim03_ss_01, $
	yfit = yfit_sim03_ss_01 )
	

bin01_ss_sim_pfit_02 = poly_fit( bin01_ss[ 0 , * ], sim_bin01_var_ss_02, 6, yband = yband_sim01_ss_02, $
	yfit = yfit_sim01_ss_02 )
bin02_ss_sim_pfit_02 = poly_fit( bin02_ss[ 0 , * ], sim_bin02_var_ss_02, 6, yband = yband_sim02_ss_02, $
	yfit = yfit_sim02_ss_02 )
bin03_ss_sim_pfit_02 = poly_fit( bin03_ss[ 0 , * ], sim_bin03_var_ss_02, 6, yband = yband_sim03_ss_02, $
	yfit = yfit_sim03_ss_02 )
	

bin01_ss_sim_pfit_03 = poly_fit( bin01_ss[ 0 , * ], sim_bin01_var_ss_03, 6, yband = yband_sim01_ss_03, $
	yfit = yfit_sim01_ss_03 )
bin02_ss_sim_pfit_03 = poly_fit( bin02_ss[ 0 , * ], sim_bin02_var_ss_03, 6, yband = yband_sim02_ss_03, $
	yfit = yfit_sim02_ss_03 )
bin03_ss_sim_pfit_03 = poly_fit( bin03_ss[ 0 , * ], sim_bin03_var_ss_03, 6, yband = yband_sim03_ss_03, $
	yfit = yfit_sim03_ss_03 )
	
	
bin01_ss_sim_pfit_04 = poly_fit( bin01_ss[ 0 , * ], sim_bin01_var_ss_04, 6, yband = yband_sim01_ss_04, $
	yfit = yfit_sim01_ss_04 )
bin02_ss_sim_pfit_04 = poly_fit( bin02_ss[ 0 , * ], sim_bin02_var_ss_04, 6, yband = yband_sim02_ss_04, $
	yfit = yfit_sim02_ss_04 )
bin03_ss_sim_pfit_04 = poly_fit( bin03_ss[ 0 , * ], sim_bin03_var_ss_04, 6, yband = yband_sim03_ss_04, $
	yfit = yfit_sim03_ss_04 )
	
	
bin01_ss_sim_pfit_05 = poly_fit( bin01_ss[ 0 , * ], sim_bin01_var_ss_05, 6, yband = yband_sim01_ss_05, $
	yfit = yfit_sim01_ss_05 )
bin02_ss_sim_pfit_05 = poly_fit( bin02_ss[ 0 , * ], sim_bin02_var_ss_05, 6, yband = yband_sim02_ss_05, $
	yfit = yfit_sim02_ss_05 )
bin03_ss_sim_pfit_05 = poly_fit( bin03_ss[ 0 , * ], sim_bin03_var_ss_05, 6, yband = yband_sim03_ss_05, $
	yfit = yfit_sim03_ss_05 )
	
	
	
	
	
bin01_ogi_sim_pfit = poly_fit( bin01_ogi[ 0 , * ], sim_bin01_var_ogi, 3, yband = yband_sim01_ogi, $
	yfit = yfit_sim01_ogi )
	
bin01_ogi_sim_pfit_01 = poly_fit( bin01_ogi[ 0 , * ], sim_bin01_var_ogi_01, 3, yband = yband_sim01_ogi_01, $
	yfit = yfit_sim01_ogi_01 )
	
bin01_ogi_sim_pfit_02 = poly_fit( bin01_ogi[ 0 , * ], sim_bin01_var_ogi_02, 3, yband = yband_sim01_ogi_02, $
	yfit = yfit_sim01_ogi_02 )
	
bin01_ogi_sim_pfit_03 = poly_fit( bin01_ogi[ 0 , * ], sim_bin01_var_ogi_03, 3, yband = yband_sim01_ogi_03, $
	yfit = yfit_sim01_ogi_03 )
	
bin01_ogi_sim_pfit_04 = poly_fit( bin01_ogi[ 0 , * ], sim_bin01_var_ogi_04, 3, yband = yband_sim01_ogi_04, $
	yfit = yfit_sim01_ogi_04 )
	
bin01_ogi_sim_pfit_05 = poly_fit( bin01_ogi[ 0 , * ], sim_bin01_var_ogi_05, 3, yband = yband_sim01_ogi_05, $
	yfit = yfit_sim01_ogi_05 )
;print, yfit_bin01_noaa, yfit_bin02_noaa, yfit_bin03_noaa, yfit_bin04_noaa, yfit_bin05_noaa, yfit_bin06_noaa, yfit_bin07_noaa
;alright, I'm going to use brute force for this, this equation will be very long, and complex, here ya go, get ready
;this equation will assign the weights, calculate the weighted mean of the integral of the poly fit function

;calculate the uncertainty of the residual of the yfit and the original values
nor_err_ogi = stddev( bin01_ogi[ 1 , where( bin01_ogi[ 0 , * ] GT 0) ] - yfit_bin01_ogi[ where( bin01_ogi[ 0 , * ] GT 0 ) ] ) / $
	sqrt( n_elements ( where( bin01_ogi[ 0 , * ] GT 0 ) ) )

sou_err_ogi = stddev( bin01_ogi[ 1 , where( bin01_ogi[ 0 , * ] LT 0) ] - yfit_bin01_ogi[ where( bin01_ogi[ 0 , * ] LT 0 ) ] ) / $
	sqrt( n_elements ( where( bin01_ogi[ 0 , * ] LT 0 ) ) )


nor_err_ss = fltarr( 3 )
sou_err_ss = fltarr( 3 )
nor_err_ss[ 0 ] = stddev( bin01_ss[ 2 , where( bin01_ss[ 0 , * ] GT 0) ] - yfit_bin01_ss[ where( bin01_ss[ 0 , * ] GT 0 ) ] ) / $
	sqrt( n_elements ( where( bin01_ss[ 0 , * ] GT 0 ) ) )
nor_err_ss[ 1 ] = stddev( bin02_ss[ 2 , where( bin02_ss[ 0 , * ] GT 0) ] - yfit_bin02_ss[ where( bin02_ss[ 0 , * ] GT 0 ) ] ) / $
	sqrt( n_elements ( where( bin02_ss[ 0 , * ] GT 0 ) ) )
nor_err_ss[ 2 ] = stddev( bin03_ss[ 2 , where( bin03_ss[ 0 , * ] GT 0) ] - yfit_bin03_ss[ where( bin03_ss[ 0 , * ] GT 0 ) ] ) / $
	sqrt( n_elements ( where( bin03_ss[ 0 , * ] GT 0 ) ) )


	
sou_err_ss[ 0 ] = stddev( bin01_ss[ 2 , where( bin01_ss[ 0 , * ] LT 0 ) ] - yfit_bin01_ss[ where( bin01_ss[ 0 , * ] LT 0 ) ] ) / $
	sqrt( n_elements ( where( bin01_ss[ 0 , * ] LT 0 ) ) )
sou_err_ss[ 1 ] = stddev( bin02_ss[ 2 , where( bin02_ss[ 0 , * ] LT 0) ] - yfit_bin02_ss[ where( bin02_ss[ 0 , * ] LT 0 ) ] ) / $
	sqrt( n_elements ( where( bin02_ss[ 0 , * ] LT 0 ) ) )
sou_err_ss[ 2 ] = stddev( bin03_ss[ 2 , where( bin03_ss[ 0 , * ] LT 0) ] - yfit_bin03_ss[ where( bin03_ss[ 0 , * ] LT 0 ) ] ) / $
	sqrt( n_elements ( where( bin03_ss[ 0 , * ] LT 0 ) ) )



temp_idx = sort( bin01_noaa[ 0 , * ] )
bin01_noaa[ 0 , * ] = bin01_noaa[ 0 , temp_idx ]
bin01_noaa[ 1 , * ] = bin01_noaa[ 1 , temp_idx ]
bin01_noaa[ 2 , * ] = bin01_noaa[ 2 , temp_idx ]

temp_idx = sort( bin02_noaa[ 0 , * ] )
bin02_noaa[ 0 , * ] = bin02_noaa[ 0 , temp_idx ]
bin02_noaa[ 1 , * ] = bin02_noaa[ 1 , temp_idx ]
bin02_noaa[ 2 , * ] = bin02_noaa[ 2 , temp_idx ]

temp_idx = sort( bin03_noaa[ 0 , * ] )
bin03_noaa[ 0 , * ] = bin03_noaa[ 0 , temp_idx ]
bin03_noaa[ 1 , * ] = bin03_noaa[ 1 , temp_idx ]
bin03_noaa[ 2 , * ] = bin03_noaa[ 2 , temp_idx ]

temp_idx = sort( bin04_noaa[ 0 , * ] )
bin04_noaa[ 0 , * ] = bin04_noaa[ 0 , temp_idx ]
bin04_noaa[ 1 , * ] = bin04_noaa[ 1 , temp_idx ]
bin04_noaa[ 2 , * ] = bin04_noaa[ 2 , temp_idx ]

temp_idx = sort( bin05_noaa[ 0 , * ] )
bin05_noaa[ 0 , * ] = bin05_noaa[ 0 , temp_idx ]
bin05_noaa[ 1 , * ] = bin05_noaa[ 1 , temp_idx ]
bin05_noaa[ 2 , * ] = bin05_noaa[ 2 , temp_idx ]

temp_idx = sort( bin06_noaa[ 0 , * ] )
bin06_noaa[ 0 , * ] = bin06_noaa[ 0 , temp_idx ]
bin06_noaa[ 1 , * ] = bin06_noaa[ 1 , temp_idx ]
bin06_noaa[ 2 , * ] = bin06_noaa[ 2 , temp_idx ]



nor_err_noaa = fltarr( 6 )
sou_err_noaa = fltarr( 6 )

nor_err_noaa[ 0 ] = stddev( bin01_noaa[ 2 , * ] - $
	yfit_bin01_noaa[ where( bin01_noaa[ 0 , * ] GT 0 ) ] ) / $
	sqrt( n_elements ( where( bin01_noaa[ 0 , * ] GT 0 ) ) )
	
nor_err_noaa[ 1 ] = stddev( bin02_noaa[ 2 , * ] - $
	yfit_bin02_noaa[ where( bin02_noaa[ 0 , * ] GT 0 ) ] ) / $
	sqrt( n_elements ( where( bin02_noaa[ 0 , * ] GT 0 ) ) )
	
nor_err_noaa[ 2 ] = stddev( bin03_noaa[ 2 , * ] - $
	yfit_bin03_noaa[ where( bin03_noaa[ 0 , * ] GT 0 ) ] ) / $
	sqrt( n_elements ( where( bin03_noaa[ 0 , * ] GT 0 ) ) )
	
nor_err_noaa[ 3 ] = stddev( bin04_noaa[ 2 , * ] - $
	yfit_bin04_noaa[ where( bin04_noaa[ 0 , * ] GT 0 ) ] ) / $
	sqrt( n_elements ( where( bin04_noaa[ 0 , * ] GT 0 ) ) )
	
nor_err_noaa[ 4 ] = stddev( bin05_noaa[ 2 , * ] - $
	yfit_bin05_noaa[ where( bin05_noaa[ 0 , * ] GT 0 ) ] ) / $
	sqrt( n_elements ( where( bin05_noaa[ 0 , * ] GT 0 ) ) )
	
nor_err_noaa[ 5 ] = stddev( bin06_noaa[ 2 , * ] - $
	yfit_bin06_noaa[ where( bin06_noaa[ 0 , * ] GT 0 ) ] ) / $
	sqrt( n_elements ( where( bin06_noaa[ 0 , * ] GT 0 ) ) )
	


sou_err_noaa[ 0 ] = stddev( bin01_noaa[ 2 , * ] - $
	yfit_bin01_noaa[ where( bin01_noaa[ 0 , * ] LT 0 ) ] ) / $
	sqrt( n_elements ( where( bin01_noaa[ 0 , * ] LT 0 ) ) )
	
sou_err_noaa[ 1 ] = stddev( bin02_noaa[ 2 , * ] - $
	yfit_bin02_noaa[ where( bin02_noaa[ 0 , * ] LT 0 ) ] ) / $
	sqrt( n_elements ( where( bin02_noaa[ 0 , * ] LT 0 ) ) )
	
sou_err_noaa[ 2 ] = stddev( bin03_noaa[ 2 , * ] - $
	yfit_bin03_noaa[ where( bin03_noaa[ 0 , * ] LT 0 ) ] ) / $
	sqrt( n_elements ( where( bin03_noaa[ 0 , * ] LT 0 ) ) )
	
sou_err_noaa[ 3 ] = stddev( bin04_noaa[ 2 , * ] - $
	yfit_bin04_noaa[ where( bin04_noaa[ 0 , * ] LT 0 ) ] ) / $
	sqrt( n_elements ( where( bin04_noaa[ 0 , * ] LT 0 ) ) )
	
sou_err_noaa[ 4 ] = stddev( bin05_noaa[ 2 , * ] - $
	yfit_bin05_noaa[ where( bin05_noaa[ 0 , * ] LT 0 ) ] ) / $
	sqrt( n_elements ( where( bin05_noaa[ 0 , * ] LT 0 ) ) )
	
sou_err_noaa[ 5 ] = stddev( bin06_noaa[ 2 , * ] - $
	yfit_bin06_noaa[ where( bin06_noaa[ 0 , * ] LT 0 ) ] ) / $
	sqrt( n_elements ( where( bin06_noaa[ 0 , * ] LT 0 ) ) )
	

nor_ss_sim = fltarr( 3 )
sou_ss_sim = fltarr( 3 )
nor_noaa_sim = fltarr( 6 )
sou_noaa_sim = fltarr( 6 )
sou_ss = fltarr( 3 )
nor_ss = fltarr( 3 )
sou_noaa = fltarr( 6 )
nor_noaa = fltarr( 6 )

nor_ss_sim_01 = fltarr( 3 )
sou_ss_sim_01 = fltarr( 3 )
nor_noaa_sim_01 = fltarr( 6 )
sou_noaa_sim_01 = fltarr( 6 )

nor_ss_sim_02 = fltarr( 3 )
sou_ss_sim_02 = fltarr( 3 )
nor_noaa_sim_02 = fltarr( 6 )
sou_noaa_sim_02 = fltarr( 6 )

nor_ss_sim_03 = fltarr( 3 )
sou_ss_sim_03 = fltarr( 3 )
nor_noaa_sim_03 = fltarr( 6 )
sou_noaa_sim_03 = fltarr( 6 )

nor_ss_sim_04 = fltarr( 3 )
sou_ss_sim_04 = fltarr( 3 )
nor_noaa_sim_04 = fltarr( 6 )
sou_noaa_sim_04 = fltarr( 6 )

nor_ss_sim_05 = fltarr( 3 )
sou_ss_sim_05 = fltarr( 3 )
nor_noaa_sim_05 = fltarr( 6 )
sou_noaa_sim_05 = fltarr( 6 )

;x = findgen(400) / 50 - 5
;y = bin01_ogi_pfit[ 0 ] + bin01_ogi_pfit[ 1 ] * x + bin01_ogi_pfit[ 2 ] * x ^ 2 + bin01_ogi_pfit[ 3 ] * x ^ 3 + bin01_ogi_pfit[ 3 ] * x ^ 4 + $
;	bin01_ogi_pfit[ 4 ] * x ^ 5 + bin01_ogi_pfit[ 5 ] * x ^ 6


rad = !pi / 180
;I don't have a script to calculate the integral of a parabola so here I'm going to estimate it using int_intabulated, not the best, but it'll do
;nor_ogi = def_integrate_poly6( bin01_ogi_pfit[ 0 ], bin01_ogi_pfit[ 1 ], bin01_ogi_pfit[ 2 ], bin01_ogi_pfit[ 3 ], $
;	bin01_ogi_pfit[ 4 ], bin01_ogi_pfit[ 5 ], bin01_ogi_pfit[ 6 ], sin( 71*rad ) , sin( 0*rad) )

nor_ogi = int_tabulated( bin01_ogi[ 0 , 0 : 2], yfit_bin01_ogi[ 0 : 2 ] )
sou_ogi = int_tabulated( bin01_ogi[ 0 , 3 : 5], yfit_bin01_ogi[ 3 : 5 ] ) 
nor_ogi_sim = int_tabulated( bin01_ogi[ 0 , 0 : 2], yfit_sim01_ogi[ 0 : 2 ] )
sou_ogi_sim = int_tabulated( bin01_ogi[ 0 , 3 : 5], yfit_sim01_ogi[ 3 : 5 ] )
nor_ogi_sim_01 = int_tabulated( bin01_ogi[ 0 , 0 : 2], yfit_sim01_ogi_01[ 0 : 2 ] )
sou_ogi_sim_01 = int_tabulated( bin01_ogi[ 0 , 3 : 5], yfit_sim01_ogi_01[ 3 : 5 ] )
nor_ogi_sim_02 = int_tabulated( bin01_ogi[ 0 , 0 : 2], yfit_sim01_ogi_02[ 0 : 2 ] )
sou_ogi_sim_02 = int_tabulated( bin01_ogi[ 0 , 3 : 5], yfit_sim01_ogi_02[ 3 : 5 ] )
nor_ogi_sim_03 = int_tabulated( bin01_ogi[ 0 , 0 : 2], yfit_sim01_ogi_03[ 0 : 2 ] )
sou_ogi_sim_03 = int_tabulated( bin01_ogi[ 0 , 3 : 5], yfit_sim01_ogi_03[ 3 : 5 ] )
nor_ogi_sim_04 = int_tabulated( bin01_ogi[ 0 , 0 : 2], yfit_sim01_ogi_04[ 0 : 2 ] )
sou_ogi_sim_04 = int_tabulated( bin01_ogi[ 0 , 3 : 5], yfit_sim01_ogi_04[ 3 : 5 ] )
nor_ogi_sim_05 = int_tabulated( bin01_ogi[ 0 , 0 : 2], yfit_sim01_ogi_05[ 0 : 2 ] )
sou_ogi_sim_05 = int_tabulated( bin01_ogi[ 0 , 3 : 5], yfit_sim01_ogi_05[ 3 : 5 ] )
;sou_ogi = def_integrate_poly6( bin01_ogi_pfit[ 0 ], bin01_ogi_pfit[ 1 ], bin01_ogi_pfit[ 2 ], bin01_ogi_pfit[ 3 ], $
;	bin01_ogi_pfit[ 4 ], bin01_ogi_pfit[ 5 ], bin01_ogi_pfit[ 6 ], sin( 0*rad ) , sin( -90*rad) )
	
	
	
nor_ss[ 0 ] = ( $
	def_integrate_poly6( bin01_ss_pfit[ 0 ], bin01_ss_pfit[ 1 ], bin01_ss_pfit[ 2 ], bin01_ss_pfit[ 3 ], $
	bin01_ss_pfit[ 4 ], bin01_ss_pfit[ 5 ], bin01_ss_pfit[ 6 ], sin( 70*rad ) , sin( 0*rad) ) )
nor_ss[ 1 ] = ( $
	def_integrate_poly6( bin02_ss_pfit[ 0 ], bin02_ss_pfit[ 1 ], bin02_ss_pfit[ 2 ], bin02_ss_pfit[ 3 ], $
	bin02_ss_pfit[ 4 ], bin02_ss_pfit[ 5 ], bin02_ss_pfit[ 6 ], sin( 70*rad ) , sin( 0*rad) ) ) 
nor_ss[ 2 ] = ( $
	def_integrate_poly6( bin03_ss_pfit[ 0 ], bin03_ss_pfit[ 1 ], bin03_ss_pfit[ 2 ], bin03_ss_pfit[ 3 ], $
	bin03_ss_pfit[ 4 ], bin03_ss_pfit[ 5 ], bin03_ss_pfit[ 6 ], sin( 70*rad ) , sin( 0*rad) ) ) 
			
			
nor_ss_sim[ 0 ] = def_integrate_poly6( bin01_ss_sim_pfit[ 0 ], bin01_ss_sim_pfit[ 1 ], bin01_ss_sim_pfit[ 2 ], $
	bin01_ss_sim_pfit[ 3 ], bin01_ss_sim_pfit[ 4 ], bin01_ss_sim_pfit[ 5 ], bin01_ss_sim_pfit[ 6 ], $
	sin( 70*rad ), sin( 0*rad ) )
nor_ss_sim[ 1 ] = def_integrate_poly6( bin02_ss_sim_pfit[ 0 ], bin02_ss_sim_pfit[ 1 ], bin02_ss_sim_pfit[ 2 ], $
	bin02_ss_sim_pfit[ 3 ], bin02_ss_sim_pfit[ 4 ], bin02_ss_sim_pfit[ 5 ], bin02_ss_sim_pfit[ 6 ], $
	sin( 70*rad ), sin( 0*rad ) )
nor_ss_sim[ 2 ] = def_integrate_poly6( bin03_ss_sim_pfit[ 0 ], bin03_ss_sim_pfit[ 1 ], bin03_ss_sim_pfit[ 2 ], $
	bin03_ss_sim_pfit[ 3 ], bin03_ss_sim_pfit[ 4 ], bin03_ss_sim_pfit[ 5 ], bin03_ss_sim_pfit[ 6 ], $
	sin( 70*rad ), sin( 0*rad ) )
	
	
	
nor_ss_sim_01[ 0 ] = def_integrate_poly6( bin01_ss_sim_pfit_01[ 0 ], bin01_ss_sim_pfit_01[ 1 ], bin01_ss_sim_pfit_01[ 2 ], $
	bin01_ss_sim_pfit_01[ 3 ], bin01_ss_sim_pfit_01[ 4 ], bin01_ss_sim_pfit_01[ 5 ], bin01_ss_sim_pfit_01[ 6 ], $
	sin( 70*rad ), sin( 0*rad ) )
nor_ss_sim_01[ 1 ] = def_integrate_poly6( bin02_ss_sim_pfit_01[ 0 ], bin02_ss_sim_pfit_01[ 1 ], bin02_ss_sim_pfit_01[ 2 ], $
	bin02_ss_sim_pfit_01[ 3 ], bin02_ss_sim_pfit_01[ 4 ], bin02_ss_sim_pfit_01[ 5 ], bin02_ss_sim_pfit_01[ 6 ], $
	sin( 70*rad ), sin( 0*rad ) )
nor_ss_sim_01[ 2 ] = def_integrate_poly6( bin03_ss_sim_pfit_01[ 0 ], bin03_ss_sim_pfit_01[ 1 ], bin03_ss_sim_pfit_01[ 2 ], $
	bin03_ss_sim_pfit_01[ 3 ], bin03_ss_sim_pfit_01[ 4 ], bin03_ss_sim_pfit_01[ 5 ], bin03_ss_sim_pfit_01[ 6 ], $
	sin( 70*rad ), sin( 0*rad ) )
	
	
	
nor_ss_sim_02[ 0 ] = def_integrate_poly6( bin01_ss_sim_pfit_02[ 0 ], bin01_ss_sim_pfit_02[ 1 ], bin01_ss_sim_pfit_02[ 2 ], $
	bin01_ss_sim_pfit_02[ 3 ], bin01_ss_sim_pfit_02[ 4 ], bin01_ss_sim_pfit_02[ 5 ], bin01_ss_sim_pfit_02[ 6 ], $
	sin( 70*rad ), sin( 0*rad ) )
nor_ss_sim_02[ 1 ] = def_integrate_poly6( bin02_ss_sim_pfit_02[ 0 ], bin02_ss_sim_pfit_02[ 1 ], bin02_ss_sim_pfit_02[ 2 ], $
	bin02_ss_sim_pfit_02[ 3 ], bin02_ss_sim_pfit_02[ 4 ], bin02_ss_sim_pfit_02[ 5 ], bin02_ss_sim_pfit_02[ 6 ], $
	sin( 70*rad ), sin( 0*rad ) )
nor_ss_sim_02[ 2 ] = def_integrate_poly6( bin03_ss_sim_pfit_02[ 0 ], bin03_ss_sim_pfit_02[ 1 ], bin03_ss_sim_pfit_02[ 2 ], $
	bin03_ss_sim_pfit_02[ 3 ], bin03_ss_sim_pfit_02[ 4 ], bin03_ss_sim_pfit_02[ 5 ], bin03_ss_sim_pfit_02[ 6 ], $
	sin( 70*rad ), sin( 0*rad ) )
	

nor_ss_sim_03[ 0 ] = def_integrate_poly6( bin01_ss_sim_pfit_03[ 0 ], bin01_ss_sim_pfit_03[ 1 ], bin01_ss_sim_pfit_03[ 2 ], $
	bin01_ss_sim_pfit_03[ 3 ], bin01_ss_sim_pfit_03[ 4 ], bin01_ss_sim_pfit_03[ 5 ], bin01_ss_sim_pfit_03[ 6 ], $
	sin( 70*rad ), sin( 0*rad ) )
nor_ss_sim_03[ 1 ] = def_integrate_poly6( bin02_ss_sim_pfit_03[ 0 ], bin02_ss_sim_pfit_03[ 1 ], bin02_ss_sim_pfit_03[ 2 ], $
	bin02_ss_sim_pfit_03[ 3 ], bin02_ss_sim_pfit_03[ 4 ], bin02_ss_sim_pfit_03[ 5 ], bin02_ss_sim_pfit_03[ 6 ], $
	sin( 70*rad ), sin( 0*rad ) )
nor_ss_sim_03[ 2 ] = def_integrate_poly6( bin03_ss_sim_pfit_03[ 0 ], bin03_ss_sim_pfit_03[ 1 ], bin03_ss_sim_pfit_03[ 2 ], $
	bin03_ss_sim_pfit_03[ 3 ], bin03_ss_sim_pfit_03[ 4 ], bin03_ss_sim_pfit_03[ 5 ], bin03_ss_sim_pfit_03[ 6 ], $
	sin( 70*rad ), sin( 0*rad ) )
	
	
nor_ss_sim_04[ 0 ] = def_integrate_poly6( bin01_ss_sim_pfit_04[ 0 ], bin01_ss_sim_pfit_04[ 1 ], bin01_ss_sim_pfit_04[ 2 ], $
	bin01_ss_sim_pfit_04[ 3 ], bin01_ss_sim_pfit_04[ 4 ], bin01_ss_sim_pfit_04[ 5 ], bin01_ss_sim_pfit_04[ 6 ], $
	sin( 70*rad ), sin( 0*rad ) )
nor_ss_sim_04[ 1 ] = def_integrate_poly6( bin02_ss_sim_pfit_04[ 0 ], bin02_ss_sim_pfit_04[ 1 ], bin02_ss_sim_pfit_04[ 2 ], $
	bin02_ss_sim_pfit_04[ 3 ], bin02_ss_sim_pfit_04[ 4 ], bin02_ss_sim_pfit_04[ 5 ], bin02_ss_sim_pfit_04[ 6 ], $
	sin( 70*rad ), sin( 0*rad ) )
nor_ss_sim_04[ 2 ] = def_integrate_poly6( bin03_ss_sim_pfit_04[ 0 ], bin03_ss_sim_pfit_04[ 1 ], bin03_ss_sim_pfit_04[ 2 ], $
	bin03_ss_sim_pfit_04[ 3 ], bin03_ss_sim_pfit_04[ 4 ], bin03_ss_sim_pfit_04[ 5 ], bin03_ss_sim_pfit_04[ 6 ], $
	sin( 70*rad ), sin( 0*rad ) )
	
	
nor_ss_sim_05[ 0 ] = def_integrate_poly6( bin01_ss_sim_pfit_05[ 0 ], bin01_ss_sim_pfit_05[ 1 ], bin01_ss_sim_pfit_05[ 2 ], $
	bin01_ss_sim_pfit_05[ 3 ], bin01_ss_sim_pfit_05[ 4 ], bin01_ss_sim_pfit_05[ 5 ], bin01_ss_sim_pfit_05[ 6 ], $
	sin( 70*rad ), sin( 0*rad ) )
nor_ss_sim_05[ 1 ] = def_integrate_poly6( bin02_ss_sim_pfit_05[ 0 ], bin02_ss_sim_pfit_05[ 1 ], bin02_ss_sim_pfit_05[ 2 ], $
	bin02_ss_sim_pfit_05[ 3 ], bin02_ss_sim_pfit_05[ 4 ], bin02_ss_sim_pfit_05[ 5 ], bin02_ss_sim_pfit_05[ 6 ], $
	sin( 70*rad ), sin( 0*rad ) )
nor_ss_sim_05[ 2 ] = def_integrate_poly6( bin03_ss_sim_pfit_05[ 0 ], bin03_ss_sim_pfit_05[ 1 ], bin03_ss_sim_pfit_05[ 2 ], $
	bin03_ss_sim_pfit_05[ 3 ], bin03_ss_sim_pfit_05[ 4 ], bin03_ss_sim_pfit_05[ 5 ], bin03_ss_sim_pfit_05[ 6 ], $
	sin( 70*rad ), sin( 0*rad ) )
	
	
	
sou_ss[ 0 ] = ( $
	def_integrate_poly6( bin01_ss_pfit[ 0 ], bin01_ss_pfit[ 1 ], bin01_ss_pfit[ 2 ], bin01_ss_pfit[ 3 ], $
	bin01_ss_pfit[ 4 ], bin01_ss_pfit[ 5 ], bin01_ss_pfit[ 6 ], sin( 0*rad ), sin( -30*rad ) ) ) 
sou_ss[ 1 ] = ( $
	def_integrate_poly6( bin02_ss_pfit[ 0 ], bin02_ss_pfit[ 1 ], bin02_ss_pfit[ 2 ], bin02_ss_pfit[ 3 ], $
	bin02_ss_pfit[ 4 ], bin02_ss_pfit[ 5 ], bin02_ss_pfit[ 6 ], sin( 0*rad ), sin( -30*rad ) ) ) 
sou_ss[ 2 ] = ( $
	def_integrate_poly6( bin03_ss_pfit[ 0 ], bin03_ss_pfit[ 1 ], bin03_ss_pfit[ 2 ], bin03_ss_pfit[ 3 ], $
	bin03_ss_pfit[ 4 ], bin03_ss_pfit[ 5 ], bin03_ss_pfit[ 6 ], sin( 0*rad ), sin( -30*rad ) ) ) 

	
	
sou_ss_sim[ 0 ] = def_integrate_poly6( bin01_ss_sim_pfit[ 0 ], bin01_ss_sim_pfit[ 1 ], bin01_ss_sim_pfit[ 2 ], $
	bin01_ss_sim_pfit[ 3 ], bin01_ss_sim_pfit[ 4 ], bin01_ss_sim_pfit[ 5 ], bin01_ss_sim_pfit[ 6 ], $
	sin( 0*rad ), sin( -30*rad ) )
sou_ss_sim[ 1 ] = def_integrate_poly6( bin02_ss_sim_pfit[ 0 ], bin02_ss_sim_pfit[ 1 ], bin02_ss_sim_pfit[ 2 ], $
	bin02_ss_sim_pfit[ 3 ], bin02_ss_sim_pfit[ 4 ], bin02_ss_sim_pfit[ 5 ], bin02_ss_sim_pfit[ 6 ], $
	sin( 0*rad ), sin( -30*rad ) )
sou_ss_sim[ 2 ] = def_integrate_poly6( bin03_ss_sim_pfit[ 0 ], bin03_ss_sim_pfit[ 1 ], bin03_ss_sim_pfit[ 2 ], $
	bin03_ss_sim_pfit[ 3 ], bin03_ss_sim_pfit[ 4 ], bin03_ss_sim_pfit[ 5 ], bin03_ss_sim_pfit[ 6 ], $
	sin( 0*rad ), sin( -30*rad ) )
	
	
	
sou_ss_sim_01[ 0 ] = def_integrate_poly6( bin01_ss_sim_pfit_01[ 0 ], bin01_ss_sim_pfit_01[ 1 ], bin01_ss_sim_pfit_01[ 2 ], $
	bin01_ss_sim_pfit_01[ 3 ], bin01_ss_sim_pfit_01[ 4 ], bin01_ss_sim_pfit_01[ 5 ], bin01_ss_sim_pfit_01[ 6 ], $
	sin( 0*rad ), sin( -30*rad ) )
sou_ss_sim_01[ 1 ] = def_integrate_poly6( bin02_ss_sim_pfit_01[ 0 ], bin02_ss_sim_pfit_01[ 1 ], bin02_ss_sim_pfit_01[ 2 ], $
	bin02_ss_sim_pfit_01[ 3 ], bin02_ss_sim_pfit_01[ 4 ], bin02_ss_sim_pfit_01[ 5 ], bin02_ss_sim_pfit_01[ 6 ], $
	sin( 0*rad ), sin( -30*rad ) )
sou_ss_sim_01[ 2 ] = def_integrate_poly6( bin03_ss_sim_pfit_01[ 0 ], bin03_ss_sim_pfit_01[ 1 ], bin03_ss_sim_pfit_01[ 2 ], $
	bin03_ss_sim_pfit_01[ 3 ], bin03_ss_sim_pfit_01[ 4 ], bin03_ss_sim_pfit_01[ 5 ], bin03_ss_sim_pfit_01[ 6 ], $
	sin( 0*rad ), sin( -30*rad ) )
	
	
sou_ss_sim_02[ 0 ] = def_integrate_poly6( bin01_ss_sim_pfit_02[ 0 ], bin01_ss_sim_pfit_02[ 1 ], bin01_ss_sim_pfit_02[ 2 ], $
	bin01_ss_sim_pfit_02[ 3 ], bin01_ss_sim_pfit_02[ 4 ], bin01_ss_sim_pfit_02[ 5 ], bin01_ss_sim_pfit_02[ 6 ], $
	sin( 0*rad ), sin( -30*rad ) )
sou_ss_sim_02[ 1 ] = def_integrate_poly6( bin02_ss_sim_pfit_02[ 0 ], bin02_ss_sim_pfit_02[ 1 ], bin02_ss_sim_pfit_02[ 2 ], $
	bin02_ss_sim_pfit_02[ 3 ], bin02_ss_sim_pfit_02[ 4 ], bin02_ss_sim_pfit_02[ 5 ], bin02_ss_sim_pfit_02[ 6 ], $
	sin( 0*rad ), sin( -30*rad ) )
sou_ss_sim_02[ 2 ] = def_integrate_poly6( bin03_ss_sim_pfit_02[ 0 ], bin03_ss_sim_pfit_02[ 1 ], bin03_ss_sim_pfit_02[ 2 ], $
	bin03_ss_sim_pfit_02[ 3 ], bin03_ss_sim_pfit_02[ 4 ], bin03_ss_sim_pfit_02[ 5 ], bin03_ss_sim_pfit_02[ 6 ], $
	sin( 0*rad ), sin( -30*rad ) )
	
	
sou_ss_sim_03[ 0 ] = def_integrate_poly6( bin01_ss_sim_pfit_03[ 0 ], bin01_ss_sim_pfit_03[ 1 ], bin01_ss_sim_pfit_03[ 2 ], $
	bin01_ss_sim_pfit_03[ 3 ], bin01_ss_sim_pfit_03[ 4 ], bin01_ss_sim_pfit_03[ 5 ], bin01_ss_sim_pfit_03[ 6 ], $
	sin( 0*rad ), sin( -30*rad ) )
sou_ss_sim_03[ 1 ] = def_integrate_poly6( bin02_ss_sim_pfit_03[ 0 ], bin02_ss_sim_pfit_03[ 1 ], bin02_ss_sim_pfit_03[ 2 ], $
	bin02_ss_sim_pfit_03[ 3 ], bin02_ss_sim_pfit_03[ 4 ], bin02_ss_sim_pfit_03[ 5 ], bin02_ss_sim_pfit_03[ 6 ], $
	sin( 0*rad ), sin( -30*rad ) )
sou_ss_sim_03[ 2 ] = def_integrate_poly6( bin03_ss_sim_pfit_03[ 0 ], bin03_ss_sim_pfit_03[ 1 ], bin03_ss_sim_pfit_03[ 2 ], $
	bin03_ss_sim_pfit_03[ 3 ], bin03_ss_sim_pfit_03[ 4 ], bin03_ss_sim_pfit_03[ 5 ], bin03_ss_sim_pfit_03[ 6 ], $
	sin( 0*rad ), sin( -30*rad ) )
	
	
sou_ss_sim_04[ 0 ] = def_integrate_poly6( bin01_ss_sim_pfit_04[ 0 ], bin01_ss_sim_pfit_04[ 1 ], bin01_ss_sim_pfit_04[ 2 ], $
	bin01_ss_sim_pfit_04[ 3 ], bin01_ss_sim_pfit_04[ 4 ], bin01_ss_sim_pfit_04[ 5 ], bin01_ss_sim_pfit_04[ 6 ], $
	sin( 0*rad ), sin( -30*rad ) )
sou_ss_sim_04[ 1 ] = def_integrate_poly6( bin02_ss_sim_pfit_04[ 0 ], bin02_ss_sim_pfit_04[ 1 ], bin02_ss_sim_pfit_04[ 2 ], $
	bin02_ss_sim_pfit_04[ 3 ], bin02_ss_sim_pfit_04[ 4 ], bin02_ss_sim_pfit_04[ 5 ], bin02_ss_sim_pfit_04[ 6 ], $
	sin( 0*rad ), sin( -30*rad ) )
sou_ss_sim_04[ 2 ] = def_integrate_poly6( bin03_ss_sim_pfit_04[ 0 ], bin03_ss_sim_pfit_04[ 1 ], bin03_ss_sim_pfit_04[ 2 ], $
	bin03_ss_sim_pfit_04[ 3 ], bin03_ss_sim_pfit_04[ 4 ], bin03_ss_sim_pfit_04[ 5 ], bin03_ss_sim_pfit_04[ 6 ], $
	sin( 0*rad ), sin( -30*rad ) )
	
	
sou_ss_sim_05[ 0 ] = def_integrate_poly6( bin01_ss_sim_pfit_05[ 0 ], bin01_ss_sim_pfit_05[ 1 ], bin01_ss_sim_pfit_05[ 2 ], $
	bin01_ss_sim_pfit_05[ 3 ], bin01_ss_sim_pfit_05[ 4 ], bin01_ss_sim_pfit_05[ 5 ], bin01_ss_sim_pfit_05[ 6 ], $
	sin( 0*rad ), sin( -30*rad ) )
sou_ss_sim_05[ 1 ] = def_integrate_poly6( bin02_ss_sim_pfit_05[ 0 ], bin02_ss_sim_pfit_05[ 1 ], bin02_ss_sim_pfit_05[ 2 ], $
	bin02_ss_sim_pfit_05[ 3 ], bin02_ss_sim_pfit_05[ 4 ], bin02_ss_sim_pfit_05[ 5 ], bin02_ss_sim_pfit_05[ 6 ], $
	sin( 0*rad ), sin( -30*rad ) )
sou_ss_sim_05[ 2 ] = def_integrate_poly6( bin03_ss_sim_pfit_05[ 0 ], bin03_ss_sim_pfit_05[ 1 ], bin03_ss_sim_pfit_05[ 2 ], $
	bin03_ss_sim_pfit_05[ 3 ], bin03_ss_sim_pfit_05[ 4 ], bin03_ss_sim_pfit_05[ 5 ], bin03_ss_sim_pfit_05[ 6 ], $
	sin( 0*rad ), sin( -30*rad ) )
			
			
			

nor_noaa[ 0 ] = ( $
	def_integrate_poly6( bin01_noaa_pfit[ 0 ], bin01_noaa_pfit[ 1 ], bin01_noaa_pfit[ 2 ], bin01_noaa_pfit[ 3 ], $
	bin01_noaa_pfit[ 4 ], bin01_noaa_pfit[ 5 ], bin01_noaa_pfit[ 6 ], sin( 70*rad ) , sin( 0*rad) ) ) 
			
nor_noaa[ 1 ] = ( $
	def_integrate_poly6( bin02_noaa_pfit[ 0 ], bin02_noaa_pfit[ 1 ], bin02_noaa_pfit[ 2 ], bin02_noaa_pfit[ 3 ], $
	bin02_noaa_pfit[ 4 ], bin02_noaa_pfit[ 5 ], bin02_noaa_pfit[ 6 ], sin( 70*rad ) , sin( 0*rad) ) ) 
			
nor_noaa[ 2 ] = ( $
	def_integrate_poly6( bin03_noaa_pfit[ 0 ], bin03_noaa_pfit[ 1 ], bin03_noaa_pfit[ 2 ], bin03_noaa_pfit[ 3 ], $
	bin03_noaa_pfit[ 4 ], bin03_noaa_pfit[ 5 ], bin03_noaa_pfit[ 6 ], sin( 70*rad ) , sin( 0*rad) ) ) 
			
nor_noaa[ 3 ] = ( $
	def_integrate_poly6( bin04_noaa_pfit[ 0 ], bin04_noaa_pfit[ 1 ], bin04_noaa_pfit[ 2 ], bin04_noaa_pfit[ 3 ], $
	bin04_noaa_pfit[ 4 ], bin04_noaa_pfit[ 5 ], bin04_noaa_pfit[ 6 ], sin( 70*rad ) , sin( 0*rad) ) ) 
			
nor_noaa[ 4 ] = ( $
	def_integrate_poly6( bin05_noaa_pfit[ 0 ], bin05_noaa_pfit[ 1 ], bin05_noaa_pfit[ 2 ], bin05_noaa_pfit[ 3 ], $
	bin05_noaa_pfit[ 4 ], bin05_noaa_pfit[ 5 ], bin05_noaa_pfit[ 6 ], sin( 70*rad ) , sin( 0*rad) ) ) 
			
nor_noaa[ 5 ] = ( $
	def_integrate_poly6( bin06_noaa_pfit[ 0 ], bin06_noaa_pfit[ 1 ], bin06_noaa_pfit[ 2 ], bin06_noaa_pfit[ 3 ], $
	bin06_noaa_pfit[ 4 ], bin06_noaa_pfit[ 5 ], bin06_noaa_pfit[ 6 ], sin( 70*rad ) , sin( 0*rad) ) )
			
			
			
nor_noaa_sim[ 0 ] = def_integrate_poly6( bin01_noaa_sim_pfit[ 0 ], bin01_noaa_sim_pfit[ 1 ], bin01_noaa_sim_pfit[ 2 ], $
	bin01_noaa_sim_pfit[ 3 ], bin01_noaa_sim_pfit[ 4 ], bin01_noaa_sim_pfit[ 5 ], bin01_noaa_sim_pfit[ 6 ], $
	sin( 70*rad ), sin( 0*rad ) )
nor_noaa_sim[ 1 ] = def_integrate_poly6( bin02_noaa_sim_pfit[ 0 ], bin02_noaa_sim_pfit[ 1 ], bin02_noaa_sim_pfit[ 2 ], $
	bin02_noaa_sim_pfit[ 3 ], bin02_noaa_sim_pfit[ 4 ], bin02_noaa_sim_pfit[ 5 ], bin02_noaa_sim_pfit[ 6 ], $
	sin( 70*rad ), sin( 0*rad ) )
nor_noaa_sim[ 2 ] = def_integrate_poly6( bin03_noaa_sim_pfit[ 0 ], bin03_noaa_sim_pfit[ 1 ], bin03_noaa_sim_pfit[ 2 ], $
	bin03_noaa_sim_pfit[ 3 ], bin03_noaa_sim_pfit[ 4 ], bin03_noaa_sim_pfit[ 5 ], bin03_noaa_sim_pfit[ 6 ], $
	sin( 70*rad ), sin( 0*rad ) )
nor_noaa_sim[ 3 ] = def_integrate_poly6( bin04_noaa_sim_pfit[ 0 ], bin04_noaa_sim_pfit[ 1 ], bin04_noaa_sim_pfit[ 2 ], $
	bin04_noaa_sim_pfit[ 3 ], bin04_noaa_sim_pfit[ 4 ], bin04_noaa_sim_pfit[ 5 ], bin04_noaa_sim_pfit[ 6 ], $
	sin( 70*rad ), sin( 0*rad ) )
nor_noaa_sim[ 4 ] = def_integrate_poly6( bin05_noaa_sim_pfit[ 0 ], bin05_noaa_sim_pfit[ 1 ], bin05_noaa_sim_pfit[ 2 ], $
	bin05_noaa_sim_pfit[ 3 ], bin05_noaa_sim_pfit[ 4 ], bin05_noaa_sim_pfit[ 5 ], bin05_noaa_sim_pfit[ 6 ], $
	sin( 70*rad ), sin( 0*rad ) )
nor_noaa_sim[ 5 ] = def_integrate_poly6( bin06_noaa_sim_pfit[ 0 ], bin06_noaa_sim_pfit[ 1 ], bin06_noaa_sim_pfit[ 2 ], $
	bin06_noaa_sim_pfit[ 3 ], bin06_noaa_sim_pfit[ 4 ], bin06_noaa_sim_pfit[ 5 ], bin06_noaa_sim_pfit[ 6 ], $
	sin( 70*rad ), sin( 0*rad ) )
			
sou_noaa_sim[ 0 ] = def_integrate_poly6( bin01_noaa_sim_pfit[ 0 ], bin01_noaa_sim_pfit[ 1 ], bin01_noaa_sim_pfit[ 2 ], $
	bin01_noaa_sim_pfit[ 3 ], bin01_noaa_sim_pfit[ 4 ], bin01_noaa_sim_pfit[ 5 ], bin01_noaa_sim_pfit[ 6 ], $
	sin( 0*rad ), sin( -30*rad ) )
sou_noaa_sim[ 1 ] = def_integrate_poly6( bin02_noaa_sim_pfit[ 0 ], bin02_noaa_sim_pfit[ 1 ], bin02_noaa_sim_pfit[ 2 ], $
	bin02_noaa_sim_pfit[ 3 ], bin02_noaa_sim_pfit[ 4 ], bin02_noaa_sim_pfit[ 5 ], bin02_noaa_sim_pfit[ 6 ], $
	sin( 0*rad ), sin( -30*rad ) )
sou_noaa_sim[ 2 ] = def_integrate_poly6( bin03_noaa_sim_pfit[ 0 ], bin03_noaa_sim_pfit[ 1 ], bin03_noaa_sim_pfit[ 2 ], $
	bin03_noaa_sim_pfit[ 3 ], bin03_noaa_sim_pfit[ 4 ], bin03_noaa_sim_pfit[ 5 ], bin03_noaa_sim_pfit[ 6 ], $
	sin( 0*rad ), sin( -30*rad ) )
sou_noaa_sim[ 3 ] = def_integrate_poly6( bin04_noaa_sim_pfit[ 0 ], bin04_noaa_sim_pfit[ 1 ], bin04_noaa_sim_pfit[ 2 ], $
	bin04_noaa_sim_pfit[ 3 ], bin04_noaa_sim_pfit[ 4 ], bin04_noaa_sim_pfit[ 5 ], bin04_noaa_sim_pfit[ 6 ], $
	sin( 0*rad ), sin( -30*rad ) )
sou_noaa_sim[ 4 ] = def_integrate_poly6( bin05_noaa_sim_pfit[ 0 ], bin05_noaa_sim_pfit[ 1 ], bin05_noaa_sim_pfit[ 2 ], $
	bin05_noaa_sim_pfit[ 3 ], bin05_noaa_sim_pfit[ 4 ], bin05_noaa_sim_pfit[ 5 ], bin05_noaa_sim_pfit[ 6 ], $
	sin( 0*rad ), sin( -30*rad ) )
sou_noaa_sim[ 5 ] = def_integrate_poly6( bin06_noaa_sim_pfit[ 0 ], bin06_noaa_sim_pfit[ 1 ], bin06_noaa_sim_pfit[ 2 ], $
	bin06_noaa_sim_pfit[ 3 ], bin06_noaa_sim_pfit[ 4 ], bin06_noaa_sim_pfit[ 5 ], bin06_noaa_sim_pfit[ 6 ], $
	sin( 0*rad ), sin( -30*rad ) )
	
	

			
nor_noaa_sim_01[ 0 ] = def_integrate_poly6( bin01_noaa_sim_pfit_01[ 0 ], bin01_noaa_sim_pfit_01[ 1 ], bin01_noaa_sim_pfit_01[ 2 ], $
	bin01_noaa_sim_pfit_01[ 3 ], bin01_noaa_sim_pfit_01[ 4 ], bin01_noaa_sim_pfit_01[ 5 ], bin01_noaa_sim_pfit_01[ 6 ], $
	sin( 70*rad ), sin( 0*rad ) )
nor_noaa_sim_01[ 1 ] = def_integrate_poly6( bin02_noaa_sim_pfit_01[ 0 ], bin02_noaa_sim_pfit_01[ 1 ], bin02_noaa_sim_pfit_01[ 2 ], $
	bin02_noaa_sim_pfit_01[ 3 ], bin02_noaa_sim_pfit_01[ 4 ], bin02_noaa_sim_pfit_01[ 5 ], bin02_noaa_sim_pfit_01[ 6 ], $
	sin( 70*rad ), sin( 0*rad ) )
nor_noaa_sim_01[ 2 ] = def_integrate_poly6( bin03_noaa_sim_pfit_01[ 0 ], bin03_noaa_sim_pfit_01[ 1 ], bin03_noaa_sim_pfit_01[ 2 ], $
	bin03_noaa_sim_pfit_01[ 3 ], bin03_noaa_sim_pfit_01[ 4 ], bin03_noaa_sim_pfit_01[ 5 ], bin03_noaa_sim_pfit_01[ 6 ], $
	sin( 70*rad ), sin( 0*rad ) )
nor_noaa_sim_01[ 3 ] = def_integrate_poly6( bin04_noaa_sim_pfit_01[ 0 ], bin04_noaa_sim_pfit_01[ 1 ], bin04_noaa_sim_pfit_01[ 2 ], $
	bin04_noaa_sim_pfit_01[ 3 ], bin04_noaa_sim_pfit_01[ 4 ], bin04_noaa_sim_pfit_01[ 5 ], bin04_noaa_sim_pfit_01[ 6 ], $
	sin( 70*rad ), sin( 0*rad ) )
nor_noaa_sim_01[ 4 ] = def_integrate_poly6( bin05_noaa_sim_pfit_01[ 0 ], bin05_noaa_sim_pfit_01[ 1 ], bin05_noaa_sim_pfit_01[ 2 ], $
	bin05_noaa_sim_pfit_01[ 3 ], bin05_noaa_sim_pfit_01[ 4 ], bin05_noaa_sim_pfit_01[ 5 ], bin05_noaa_sim_pfit_01[ 6 ], $
	sin( 70*rad ), sin( 0*rad ) )
nor_noaa_sim_01[ 5 ] = def_integrate_poly6( bin06_noaa_sim_pfit_01[ 0 ], bin06_noaa_sim_pfit_01[ 1 ], bin06_noaa_sim_pfit_01[ 2 ], $
	bin06_noaa_sim_pfit_01[ 3 ], bin06_noaa_sim_pfit_01[ 4 ], bin06_noaa_sim_pfit_01[ 5 ], bin06_noaa_sim_pfit_01[ 6 ], $
	sin( 70*rad ), sin( 0*rad ) )
			
			
			
sou_noaa_sim_01[ 0 ] = def_integrate_poly6( bin01_noaa_sim_pfit_01[ 0 ], bin01_noaa_sim_pfit_01[ 1 ], bin01_noaa_sim_pfit_01[ 2 ], $
	bin01_noaa_sim_pfit_01[ 3 ], bin01_noaa_sim_pfit_01[ 4 ], bin01_noaa_sim_pfit_01[ 5 ], bin01_noaa_sim_pfit_01[ 6 ], $
	sin( 0*rad ), sin( -30*rad ) )
sou_noaa_sim_01[ 1 ] = def_integrate_poly6( bin02_noaa_sim_pfit_01[ 0 ], bin02_noaa_sim_pfit_01[ 1 ], bin02_noaa_sim_pfit_01[ 2 ], $
	bin02_noaa_sim_pfit_01[ 3 ], bin02_noaa_sim_pfit_01[ 4 ], bin02_noaa_sim_pfit_01[ 5 ], bin02_noaa_sim_pfit_01[ 6 ], $
	sin( 0*rad ), sin( -30*rad ) )
sou_noaa_sim_01[ 2 ] = def_integrate_poly6( bin03_noaa_sim_pfit_01[ 0 ], bin03_noaa_sim_pfit_01[ 1 ], bin03_noaa_sim_pfit_01[ 2 ], $
	bin03_noaa_sim_pfit_01[ 3 ], bin03_noaa_sim_pfit_01[ 4 ], bin03_noaa_sim_pfit_01[ 5 ], bin03_noaa_sim_pfit_01[ 6 ], $
	sin( 0*rad ), sin( -30*rad ) )
sou_noaa_sim_01[ 3 ] = def_integrate_poly6( bin04_noaa_sim_pfit_01[ 0 ], bin04_noaa_sim_pfit_01[ 1 ], bin04_noaa_sim_pfit_01[ 2 ], $
	bin04_noaa_sim_pfit_01[ 3 ], bin04_noaa_sim_pfit_01[ 4 ], bin04_noaa_sim_pfit_01[ 5 ], bin04_noaa_sim_pfit_01[ 6 ], $
	sin( 0*rad ), sin( -30*rad ) )
sou_noaa_sim_01[ 4 ] = def_integrate_poly6( bin05_noaa_sim_pfit_01[ 0 ], bin05_noaa_sim_pfit_01[ 1 ], bin05_noaa_sim_pfit_01[ 2 ], $
	bin05_noaa_sim_pfit_01[ 3 ], bin05_noaa_sim_pfit_01[ 4 ], bin05_noaa_sim_pfit_01[ 5 ], bin05_noaa_sim_pfit_01[ 6 ], $
	sin( 0*rad ), sin( -30*rad ) )
sou_noaa_sim_01[ 5 ] = def_integrate_poly6( bin06_noaa_sim_pfit_01[ 0 ], bin06_noaa_sim_pfit_01[ 1 ], bin06_noaa_sim_pfit_01[ 2 ], $
	bin06_noaa_sim_pfit_01[ 3 ], bin06_noaa_sim_pfit_01[ 4 ], bin06_noaa_sim_pfit_01[ 5 ], bin06_noaa_sim_pfit_01[ 6 ], $
	sin( 0*rad ), sin( -30*rad ) )

	

nor_noaa_sim_02[ 0 ] = def_integrate_poly6( bin01_noaa_sim_pfit_02[ 0 ], bin01_noaa_sim_pfit_02[ 1 ], bin01_noaa_sim_pfit_02[ 2 ], $
	bin01_noaa_sim_pfit_02[ 3 ], bin01_noaa_sim_pfit_02[ 4 ], bin01_noaa_sim_pfit_02[ 5 ], bin01_noaa_sim_pfit_02[ 6 ], $
	sin( 70*rad ), sin( 0*rad ) )
nor_noaa_sim_02[ 1 ] = def_integrate_poly6( bin02_noaa_sim_pfit_02[ 0 ], bin02_noaa_sim_pfit_02[ 1 ], bin02_noaa_sim_pfit_02[ 2 ], $
	bin02_noaa_sim_pfit_02[ 3 ], bin02_noaa_sim_pfit_02[ 4 ], bin02_noaa_sim_pfit_02[ 5 ], bin02_noaa_sim_pfit_02[ 6 ], $
	sin( 70*rad ), sin( 0*rad ) )
nor_noaa_sim_02[ 2 ] = def_integrate_poly6( bin03_noaa_sim_pfit_02[ 0 ], bin03_noaa_sim_pfit_02[ 1 ], bin03_noaa_sim_pfit_02[ 2 ], $
	bin03_noaa_sim_pfit_02[ 3 ], bin03_noaa_sim_pfit_02[ 4 ], bin03_noaa_sim_pfit_02[ 5 ], bin03_noaa_sim_pfit_02[ 6 ], $
	sin( 70*rad ), sin( 0*rad ) )
nor_noaa_sim_02[ 3 ] = def_integrate_poly6( bin04_noaa_sim_pfit_02[ 0 ], bin04_noaa_sim_pfit_02[ 1 ], bin04_noaa_sim_pfit_02[ 2 ], $
	bin04_noaa_sim_pfit_02[ 3 ], bin04_noaa_sim_pfit_02[ 4 ], bin04_noaa_sim_pfit_02[ 5 ], bin04_noaa_sim_pfit_02[ 6 ], $
	sin( 70*rad ), sin( 0*rad ) )
nor_noaa_sim_02[ 4 ] = def_integrate_poly6( bin05_noaa_sim_pfit_02[ 0 ], bin05_noaa_sim_pfit_02[ 1 ], bin05_noaa_sim_pfit_02[ 2 ], $
	bin05_noaa_sim_pfit_02[ 3 ], bin05_noaa_sim_pfit_02[ 4 ], bin05_noaa_sim_pfit_02[ 5 ], bin05_noaa_sim_pfit_02[ 6 ], $
	sin( 70*rad ), sin( 0*rad ) )
nor_noaa_sim_02[ 5 ] = def_integrate_poly6( bin06_noaa_sim_pfit_02[ 0 ], bin06_noaa_sim_pfit_02[ 1 ], bin06_noaa_sim_pfit_02[ 2 ], $
	bin06_noaa_sim_pfit_02[ 3 ], bin06_noaa_sim_pfit_02[ 4 ], bin06_noaa_sim_pfit_02[ 5 ], bin06_noaa_sim_pfit_02[ 6 ], $
	sin( 70*rad ), sin( 0*rad ) )
			
sou_noaa_sim_02[ 0 ] = def_integrate_poly6( bin01_noaa_sim_pfit_02[ 0 ], bin01_noaa_sim_pfit_02[ 1 ], bin01_noaa_sim_pfit_02[ 2 ], $
	bin01_noaa_sim_pfit_02[ 3 ], bin01_noaa_sim_pfit_02[ 4 ], bin01_noaa_sim_pfit_02[ 5 ], bin01_noaa_sim_pfit_02[ 6 ], $
	sin( 0*rad ), sin( -30*rad ) )
sou_noaa_sim_02[ 1 ] = def_integrate_poly6( bin02_noaa_sim_pfit_02[ 0 ], bin02_noaa_sim_pfit_02[ 1 ], bin02_noaa_sim_pfit_02[ 2 ], $
	bin02_noaa_sim_pfit_02[ 3 ], bin02_noaa_sim_pfit_02[ 4 ], bin02_noaa_sim_pfit_02[ 5 ], bin02_noaa_sim_pfit_02[ 6 ], $
	sin( 0*rad ), sin( -30*rad ) )
sou_noaa_sim_02[ 2 ] = def_integrate_poly6( bin03_noaa_sim_pfit_02[ 0 ], bin03_noaa_sim_pfit_02[ 1 ], bin03_noaa_sim_pfit_02[ 2 ], $
	bin03_noaa_sim_pfit_02[ 3 ], bin03_noaa_sim_pfit_02[ 4 ], bin03_noaa_sim_pfit_02[ 5 ], bin03_noaa_sim_pfit_02[ 6 ], $
	sin( 0*rad ), sin( -30*rad ) )
sou_noaa_sim_02[ 3 ] = def_integrate_poly6( bin04_noaa_sim_pfit_02[ 0 ], bin04_noaa_sim_pfit_02[ 1 ], bin04_noaa_sim_pfit_02[ 2 ], $
	bin04_noaa_sim_pfit_02[ 3 ], bin04_noaa_sim_pfit_02[ 4 ], bin04_noaa_sim_pfit_02[ 5 ], bin04_noaa_sim_pfit_02[ 6 ], $
	sin( 0*rad ), sin( -30*rad ) )
sou_noaa_sim_02[ 4 ] = def_integrate_poly6( bin05_noaa_sim_pfit_02[ 0 ], bin05_noaa_sim_pfit_02[ 1 ], bin05_noaa_sim_pfit_02[ 2 ], $
	bin05_noaa_sim_pfit_02[ 3 ], bin05_noaa_sim_pfit_02[ 4 ], bin05_noaa_sim_pfit_02[ 5 ], bin05_noaa_sim_pfit_02[ 6 ], $
	sin( 0*rad ), sin( -30*rad ) )
sou_noaa_sim_02[ 5 ] = def_integrate_poly6( bin06_noaa_sim_pfit_02[ 0 ], bin06_noaa_sim_pfit_02[ 1 ], bin06_noaa_sim_pfit_02[ 2 ], $
	bin06_noaa_sim_pfit_02[ 3 ], bin06_noaa_sim_pfit_02[ 4 ], bin06_noaa_sim_pfit_02[ 5 ], bin06_noaa_sim_pfit_02[ 6 ], $
	sin( 0*rad ), sin( -30*rad ) )
	
	
	
nor_noaa_sim_03[ 0 ] = def_integrate_poly6( bin01_noaa_sim_pfit_03[ 0 ], bin01_noaa_sim_pfit_03[ 1 ], bin01_noaa_sim_pfit_03[ 2 ], $
	bin01_noaa_sim_pfit_03[ 3 ], bin01_noaa_sim_pfit_03[ 4 ], bin01_noaa_sim_pfit_03[ 5 ], bin01_noaa_sim_pfit_03[ 6 ], $
	sin( 70*rad ), sin( 0*rad ) )
nor_noaa_sim_03[ 1 ] = def_integrate_poly6( bin02_noaa_sim_pfit_03[ 0 ], bin02_noaa_sim_pfit_03[ 1 ], bin02_noaa_sim_pfit_03[ 2 ], $
	bin02_noaa_sim_pfit_03[ 3 ], bin02_noaa_sim_pfit_03[ 4 ], bin02_noaa_sim_pfit_03[ 5 ], bin02_noaa_sim_pfit_03[ 6 ], $
	sin( 70*rad ), sin( 0*rad ) )
nor_noaa_sim_03[ 2 ] = def_integrate_poly6( bin03_noaa_sim_pfit_03[ 0 ], bin03_noaa_sim_pfit_03[ 1 ], bin03_noaa_sim_pfit_03[ 2 ], $
	bin03_noaa_sim_pfit_03[ 3 ], bin03_noaa_sim_pfit_03[ 4 ], bin03_noaa_sim_pfit_03[ 5 ], bin03_noaa_sim_pfit_03[ 6 ], $
	sin( 70*rad ), sin( 0*rad ) )
nor_noaa_sim_03[ 3 ] = def_integrate_poly6( bin04_noaa_sim_pfit_03[ 0 ], bin04_noaa_sim_pfit_03[ 1 ], bin04_noaa_sim_pfit_03[ 2 ], $
	bin04_noaa_sim_pfit_03[ 3 ], bin04_noaa_sim_pfit_03[ 4 ], bin04_noaa_sim_pfit_03[ 5 ], bin04_noaa_sim_pfit_03[ 6 ], $
	sin( 70*rad ), sin( 0*rad ) )
nor_noaa_sim_03[ 4 ] = def_integrate_poly6( bin05_noaa_sim_pfit_03[ 0 ], bin05_noaa_sim_pfit_03[ 1 ], bin05_noaa_sim_pfit_03[ 2 ], $
	bin05_noaa_sim_pfit_03[ 3 ], bin05_noaa_sim_pfit_03[ 4 ], bin05_noaa_sim_pfit_03[ 5 ], bin05_noaa_sim_pfit_03[ 6 ], $
	sin( 70*rad ), sin( 0*rad ) )
nor_noaa_sim_03[ 5 ] = def_integrate_poly6( bin06_noaa_sim_pfit_03[ 0 ], bin06_noaa_sim_pfit_03[ 1 ], bin06_noaa_sim_pfit_03[ 2 ], $
	bin06_noaa_sim_pfit_03[ 3 ], bin06_noaa_sim_pfit_03[ 4 ], bin06_noaa_sim_pfit_03[ 5 ], bin06_noaa_sim_pfit_03[ 6 ], $
	sin( 70*rad ), sin( 0*rad ) )
			
sou_noaa_sim_03[ 0 ] = def_integrate_poly6( bin01_noaa_sim_pfit_03[ 0 ], bin01_noaa_sim_pfit_03[ 1 ], bin01_noaa_sim_pfit_03[ 2 ], $
	bin01_noaa_sim_pfit_03[ 3 ], bin01_noaa_sim_pfit_03[ 4 ], bin01_noaa_sim_pfit_03[ 5 ], bin01_noaa_sim_pfit_03[ 6 ], $
	sin( 0*rad ), sin( -30*rad ) )
sou_noaa_sim_03[ 1 ] = def_integrate_poly6( bin02_noaa_sim_pfit_03[ 0 ], bin02_noaa_sim_pfit_03[ 1 ], bin02_noaa_sim_pfit_03[ 2 ], $
	bin02_noaa_sim_pfit_03[ 3 ], bin02_noaa_sim_pfit_03[ 4 ], bin02_noaa_sim_pfit_03[ 5 ], bin02_noaa_sim_pfit_03[ 6 ], $
	sin( 0*rad ), sin( -30*rad ) )
sou_noaa_sim_03[ 2 ] = def_integrate_poly6( bin03_noaa_sim_pfit_03[ 0 ], bin03_noaa_sim_pfit_03[ 1 ], bin03_noaa_sim_pfit_03[ 2 ], $
	bin03_noaa_sim_pfit_03[ 3 ], bin03_noaa_sim_pfit_03[ 4 ], bin03_noaa_sim_pfit_03[ 5 ], bin03_noaa_sim_pfit_03[ 6 ], $
	sin( 0*rad ), sin( -30*rad ) )
sou_noaa_sim_03[ 3 ] = def_integrate_poly6( bin04_noaa_sim_pfit_03[ 0 ], bin04_noaa_sim_pfit_03[ 1 ], bin04_noaa_sim_pfit_03[ 2 ], $
	bin04_noaa_sim_pfit_03[ 3 ], bin04_noaa_sim_pfit_03[ 4 ], bin04_noaa_sim_pfit_03[ 5 ], bin04_noaa_sim_pfit_03[ 6 ], $
	sin( 0*rad ), sin( -30*rad ) )
sou_noaa_sim_03[ 4 ] = def_integrate_poly6( bin05_noaa_sim_pfit_03[ 0 ], bin05_noaa_sim_pfit_03[ 1 ], bin05_noaa_sim_pfit_03[ 2 ], $
	bin05_noaa_sim_pfit_03[ 3 ], bin05_noaa_sim_pfit_03[ 4 ], bin05_noaa_sim_pfit_03[ 5 ], bin05_noaa_sim_pfit_03[ 6 ], $
	sin( 0*rad ), sin( -30*rad ) )
sou_noaa_sim_03[ 5 ] = def_integrate_poly6( bin06_noaa_sim_pfit_03[ 0 ], bin06_noaa_sim_pfit_03[ 1 ], bin06_noaa_sim_pfit_03[ 2 ], $
	bin06_noaa_sim_pfit_03[ 3 ], bin06_noaa_sim_pfit_03[ 4 ], bin06_noaa_sim_pfit_03[ 5 ], bin06_noaa_sim_pfit_03[ 6 ], $
	sin( 0*rad ), sin( -30*rad ) )
	
	
	
nor_noaa_sim_04[ 0 ] = def_integrate_poly6( bin01_noaa_sim_pfit_04[ 0 ], bin01_noaa_sim_pfit_04[ 1 ], bin01_noaa_sim_pfit_04[ 2 ], $
	bin01_noaa_sim_pfit_04[ 3 ], bin01_noaa_sim_pfit_04[ 4 ], bin01_noaa_sim_pfit_04[ 5 ], bin01_noaa_sim_pfit_04[ 6 ], $
	sin( 70*rad ), sin( 0*rad ) )
nor_noaa_sim_04[ 1 ] = def_integrate_poly6( bin02_noaa_sim_pfit_04[ 0 ], bin02_noaa_sim_pfit_04[ 1 ], bin02_noaa_sim_pfit_04[ 2 ], $
	bin02_noaa_sim_pfit_04[ 3 ], bin02_noaa_sim_pfit_04[ 4 ], bin02_noaa_sim_pfit_04[ 5 ], bin02_noaa_sim_pfit_04[ 6 ], $
	sin( 70*rad ), sin( 0*rad ) )
nor_noaa_sim_04[ 2 ] = def_integrate_poly6( bin03_noaa_sim_pfit_04[ 0 ], bin03_noaa_sim_pfit_04[ 1 ], bin03_noaa_sim_pfit_04[ 2 ], $
	bin03_noaa_sim_pfit_04[ 3 ], bin03_noaa_sim_pfit_04[ 4 ], bin03_noaa_sim_pfit_04[ 5 ], bin03_noaa_sim_pfit_04[ 6 ], $
	sin( 70*rad ), sin( 0*rad ) )
nor_noaa_sim_04[ 3 ] = def_integrate_poly6( bin04_noaa_sim_pfit_04[ 0 ], bin04_noaa_sim_pfit_04[ 1 ], bin04_noaa_sim_pfit_04[ 2 ], $
	bin04_noaa_sim_pfit_04[ 3 ], bin04_noaa_sim_pfit_04[ 4 ], bin04_noaa_sim_pfit_04[ 5 ], bin04_noaa_sim_pfit_04[ 6 ], $
	sin( 70*rad ), sin( 0*rad ) )
nor_noaa_sim_04[ 4 ] = def_integrate_poly6( bin05_noaa_sim_pfit_04[ 0 ], bin05_noaa_sim_pfit_04[ 1 ], bin05_noaa_sim_pfit_04[ 2 ], $
	bin05_noaa_sim_pfit_04[ 3 ], bin05_noaa_sim_pfit_04[ 4 ], bin05_noaa_sim_pfit_04[ 5 ], bin05_noaa_sim_pfit_04[ 6 ], $
	sin( 70*rad ), sin( 0*rad ) )
nor_noaa_sim_04[ 5 ] = def_integrate_poly6( bin06_noaa_sim_pfit_04[ 0 ], bin06_noaa_sim_pfit_04[ 1 ], bin06_noaa_sim_pfit_04[ 2 ], $
	bin06_noaa_sim_pfit_04[ 3 ], bin06_noaa_sim_pfit_04[ 4 ], bin06_noaa_sim_pfit_04[ 5 ], bin06_noaa_sim_pfit_04[ 6 ], $
	sin( 70*rad ), sin( 0*rad ) )
			
sou_noaa_sim_04[ 0 ] = def_integrate_poly6( bin01_noaa_sim_pfit_04[ 0 ], bin01_noaa_sim_pfit_04[ 1 ], bin01_noaa_sim_pfit_04[ 2 ], $
	bin01_noaa_sim_pfit_04[ 3 ], bin01_noaa_sim_pfit_04[ 4 ], bin01_noaa_sim_pfit_04[ 5 ], bin01_noaa_sim_pfit_04[ 6 ], $
	sin( 0*rad ), sin( -30*rad ) )
sou_noaa_sim_04[ 1 ] = def_integrate_poly6( bin02_noaa_sim_pfit_04[ 0 ], bin02_noaa_sim_pfit_04[ 1 ], bin02_noaa_sim_pfit_04[ 2 ], $
	bin02_noaa_sim_pfit_04[ 3 ], bin02_noaa_sim_pfit_04[ 4 ], bin02_noaa_sim_pfit_04[ 5 ], bin02_noaa_sim_pfit_04[ 6 ], $
	sin( 0*rad ), sin( -30*rad ) )
sou_noaa_sim_04[ 2 ] = def_integrate_poly6( bin03_noaa_sim_pfit_04[ 0 ], bin03_noaa_sim_pfit_04[ 1 ], bin03_noaa_sim_pfit_04[ 2 ], $
	bin03_noaa_sim_pfit_04[ 3 ], bin03_noaa_sim_pfit_04[ 4 ], bin03_noaa_sim_pfit_04[ 5 ], bin03_noaa_sim_pfit_04[ 6 ], $
	sin( 0*rad ), sin( -30*rad ) )
sou_noaa_sim_04[ 3 ] = def_integrate_poly6( bin04_noaa_sim_pfit_04[ 0 ], bin04_noaa_sim_pfit_04[ 1 ], bin04_noaa_sim_pfit_04[ 2 ], $
	bin04_noaa_sim_pfit_04[ 3 ], bin04_noaa_sim_pfit_04[ 4 ], bin04_noaa_sim_pfit_04[ 5 ], bin04_noaa_sim_pfit_04[ 6 ], $
	sin( 0*rad ), sin( -30*rad ) )
sou_noaa_sim_04[ 4 ] = def_integrate_poly6( bin05_noaa_sim_pfit_04[ 0 ], bin05_noaa_sim_pfit_04[ 1 ], bin05_noaa_sim_pfit_04[ 2 ], $
	bin05_noaa_sim_pfit_04[ 3 ], bin05_noaa_sim_pfit_04[ 4 ], bin05_noaa_sim_pfit_04[ 5 ], bin05_noaa_sim_pfit_04[ 6 ], $
	sin( 0*rad ), sin( -30*rad ) )
sou_noaa_sim_04[ 5 ] = def_integrate_poly6( bin06_noaa_sim_pfit_04[ 0 ], bin06_noaa_sim_pfit_04[ 1 ], bin06_noaa_sim_pfit_04[ 2 ], $
	bin06_noaa_sim_pfit_04[ 3 ], bin06_noaa_sim_pfit_04[ 4 ], bin06_noaa_sim_pfit_04[ 5 ], bin06_noaa_sim_pfit_04[ 6 ], $
	sin( 0*rad ), sin( -30*rad ) )
	
	
	
nor_noaa_sim_05[ 0 ] = def_integrate_poly6( bin01_noaa_sim_pfit_05[ 0 ], bin01_noaa_sim_pfit_05[ 1 ], bin01_noaa_sim_pfit_05[ 2 ], $
	bin01_noaa_sim_pfit_05[ 3 ], bin01_noaa_sim_pfit_05[ 4 ], bin01_noaa_sim_pfit_05[ 5 ], bin01_noaa_sim_pfit_05[ 6 ], $
	sin( 70*rad ), sin( 0*rad ) )
nor_noaa_sim_05[ 1 ] = def_integrate_poly6( bin02_noaa_sim_pfit_05[ 0 ], bin02_noaa_sim_pfit_05[ 1 ], bin02_noaa_sim_pfit_05[ 2 ], $
	bin02_noaa_sim_pfit_05[ 3 ], bin02_noaa_sim_pfit_05[ 4 ], bin02_noaa_sim_pfit_05[ 5 ], bin02_noaa_sim_pfit_05[ 6 ], $
	sin( 70*rad ), sin( 0*rad ) )
nor_noaa_sim_05[ 2 ] = def_integrate_poly6( bin03_noaa_sim_pfit_05[ 0 ], bin03_noaa_sim_pfit_05[ 1 ], bin03_noaa_sim_pfit_05[ 2 ], $
	bin03_noaa_sim_pfit_05[ 3 ], bin03_noaa_sim_pfit_05[ 4 ], bin03_noaa_sim_pfit_05[ 5 ], bin03_noaa_sim_pfit_05[ 6 ], $
	sin( 70*rad ), sin( 0*rad ) )
nor_noaa_sim_05[ 3 ] = def_integrate_poly6( bin04_noaa_sim_pfit_05[ 0 ], bin04_noaa_sim_pfit_05[ 1 ], bin04_noaa_sim_pfit_05[ 2 ], $
	bin04_noaa_sim_pfit_05[ 3 ], bin04_noaa_sim_pfit_05[ 4 ], bin04_noaa_sim_pfit_05[ 5 ], bin04_noaa_sim_pfit_05[ 6 ], $
	sin( 70*rad ), sin( 0*rad ) )
nor_noaa_sim_05[ 4 ] = def_integrate_poly6( bin05_noaa_sim_pfit_05[ 0 ], bin05_noaa_sim_pfit_05[ 1 ], bin05_noaa_sim_pfit_05[ 2 ], $
	bin05_noaa_sim_pfit_05[ 3 ], bin05_noaa_sim_pfit_05[ 4 ], bin05_noaa_sim_pfit_05[ 5 ], bin05_noaa_sim_pfit_05[ 6 ], $
	sin( 70*rad ), sin( 0*rad ) )
nor_noaa_sim_05[ 5 ] = def_integrate_poly6( bin06_noaa_sim_pfit_05[ 0 ], bin06_noaa_sim_pfit_05[ 1 ], bin06_noaa_sim_pfit_05[ 2 ], $
	bin06_noaa_sim_pfit_05[ 3 ], bin06_noaa_sim_pfit_05[ 4 ], bin06_noaa_sim_pfit_05[ 5 ], bin06_noaa_sim_pfit_05[ 6 ], $
	sin( 70*rad ), sin( 0*rad ) )
			
sou_noaa_sim_05[ 0 ] = def_integrate_poly6( bin01_noaa_sim_pfit_05[ 0 ], bin01_noaa_sim_pfit_05[ 1 ], bin01_noaa_sim_pfit_05[ 2 ], $
	bin01_noaa_sim_pfit_05[ 3 ], bin01_noaa_sim_pfit_05[ 4 ], bin01_noaa_sim_pfit_05[ 5 ], bin01_noaa_sim_pfit_05[ 6 ], $
	sin( 0*rad ), sin( -30*rad ) )
sou_noaa_sim_05[ 1 ] = def_integrate_poly6( bin02_noaa_sim_pfit_05[ 0 ], bin02_noaa_sim_pfit_05[ 1 ], bin02_noaa_sim_pfit_05[ 2 ], $
	bin02_noaa_sim_pfit_05[ 3 ], bin02_noaa_sim_pfit_05[ 4 ], bin02_noaa_sim_pfit_05[ 5 ], bin02_noaa_sim_pfit_05[ 6 ], $
	sin( 0*rad ), sin( -30*rad ) )
sou_noaa_sim_05[ 2 ] = def_integrate_poly6( bin03_noaa_sim_pfit_05[ 0 ], bin03_noaa_sim_pfit_05[ 1 ], bin03_noaa_sim_pfit_05[ 2 ], $
	bin03_noaa_sim_pfit_05[ 3 ], bin03_noaa_sim_pfit_05[ 4 ], bin03_noaa_sim_pfit_05[ 5 ], bin03_noaa_sim_pfit_05[ 6 ], $
	sin( 0*rad ), sin( -30*rad ) )
sou_noaa_sim_05[ 3 ] = def_integrate_poly6( bin04_noaa_sim_pfit_05[ 0 ], bin04_noaa_sim_pfit_05[ 1 ], bin04_noaa_sim_pfit_05[ 2 ], $
	bin04_noaa_sim_pfit_05[ 3 ], bin04_noaa_sim_pfit_05[ 4 ], bin04_noaa_sim_pfit_05[ 5 ], bin04_noaa_sim_pfit_05[ 6 ], $
	sin( 0*rad ), sin( -30*rad ) )
sou_noaa_sim_05[ 4 ] = def_integrate_poly6( bin05_noaa_sim_pfit_05[ 0 ], bin05_noaa_sim_pfit_05[ 1 ], bin05_noaa_sim_pfit_05[ 2 ], $
	bin05_noaa_sim_pfit_05[ 3 ], bin05_noaa_sim_pfit_05[ 4 ], bin05_noaa_sim_pfit_05[ 5 ], bin05_noaa_sim_pfit_05[ 6 ], $
	sin( 0*rad ), sin( -30*rad ) )
sou_noaa_sim_05[ 5 ] = def_integrate_poly6( bin06_noaa_sim_pfit_05[ 0 ], bin06_noaa_sim_pfit_05[ 1 ], bin06_noaa_sim_pfit_05[ 2 ], $
	bin06_noaa_sim_pfit_05[ 3 ], bin06_noaa_sim_pfit_05[ 4 ], bin06_noaa_sim_pfit_05[ 5 ], bin06_noaa_sim_pfit_05[ 6 ], $
	sin( 0*rad ), sin( -30*rad ) )
	
	
	
sou_noaa[ 0 ] = ( $
	def_integrate_poly6( bin01_noaa_pfit[ 0 ], bin01_noaa_pfit[ 1 ], bin01_noaa_pfit[ 2 ], bin01_noaa_pfit[ 3 ], $
	bin01_noaa_pfit[ 4 ], bin01_noaa_pfit[ 5 ], bin01_noaa_pfit[ 6 ], sin( 0*rad ), sin( -30*rad ) ) ) 
	
sou_noaa[ 1 ] = ( $
	def_integrate_poly6( bin02_noaa_pfit[ 0 ], bin02_noaa_pfit[ 1 ], bin02_noaa_pfit[ 2 ], bin02_noaa_pfit[ 3 ], $
	bin02_noaa_pfit[ 4 ], bin02_noaa_pfit[ 5 ], bin02_noaa_pfit[ 6 ], sin( 0*rad ), sin( -30*rad ) ) )
			
sou_noaa[ 2 ] = ( $
	def_integrate_poly6( bin03_noaa_pfit[ 0 ], bin03_noaa_pfit[ 1 ], bin03_noaa_pfit[ 2 ], bin03_noaa_pfit[ 3 ], $
	bin03_noaa_pfit[ 4 ], bin03_noaa_pfit[ 5 ], bin03_noaa_pfit[ 6 ], sin( 0*rad ), sin( -30*rad ) ) )
		
sou_noaa[ 3 ] = ( $
	def_integrate_poly6( bin04_noaa_pfit[ 0 ], bin04_noaa_pfit[ 1 ], bin04_noaa_pfit[ 2 ], bin04_noaa_pfit[ 3 ], $
	bin04_noaa_pfit[ 4 ], bin04_noaa_pfit[ 5 ], bin04_noaa_pfit[ 6 ], sin( 0*rad ), sin( -30*rad ) ) ) 
			
sou_noaa[ 4 ] = ( $
	def_integrate_poly6( bin05_noaa_pfit[ 0 ], bin05_noaa_pfit[ 1 ], bin05_noaa_pfit[ 2 ], bin05_noaa_pfit[ 3 ], $
	bin05_noaa_pfit[ 4 ], bin05_noaa_pfit[ 5 ], bin05_noaa_pfit[ 6 ], sin( 0*rad ), sin( -30*rad ) ) ) 
			
sou_noaa[ 5 ] = ( $
	def_integrate_poly6( bin06_noaa_pfit[ 0 ], bin06_noaa_pfit[ 1 ], bin06_noaa_pfit[ 2 ], bin06_noaa_pfit[ 3 ], $
	bin06_noaa_pfit[ 4 ], bin06_noaa_pfit[ 5 ], bin06_noaa_pfit[ 6 ], sin( 0*rad ), sin( -30*rad ) ) ) 
			

year_series_ogi = 1985
year_series_ss = [ 1998.5 , 2003.5 , 2008 ]
year_series_noaa = [ 2006 , 2008 , 2010 , 2012 , 2014 , 2015 ]
;because we don't have 2015 sim data, so remove year 2015 from the plot
year_series_noaa_sim = [ 2006, 2008, 2010, 2012, 2014 ]
			
;the fomular for calculating the uncertainty sigma squared for the ratio=
ratio_err_noaa = sqrt( ( ( nor_err_noaa ) ^ 2 ) / ( ( sou_noaa ) ^ 2 ) + $
	( ( nor_noaa ^ 2 ) * ( sou_err_noaa ^ 2 ) ) / ( sou_noaa ^ 4 ) )

ratio_err_ogi = sqrt( ( ( nor_err_ogi ) ^ 2 ) / ( ( sou_ogi ) ^ 2 ) + $
	( ( nor_ogi ^ 2 ) * ( sou_err_ogi ^ 2 ) ) / ( sou_ogi ^ 4 ) )

ratio_err_ss = sqrt( ( ( nor_err_ss ) ^ 2 ) / ( ( sou_ss ) ^ 2 ) + $
	( ( nor_ss ^ 2 ) * ( sou_err_ss ^ 2 ) ) / ( sou_ss ^ 4 ) )
print, 'nor_err_noaa:', nor_err_noaa
print, 'sou_noaa:', sou_noaa
print, 'nor_noaa:', nor_noaa
print, 'sou_err_noaa:', sou_err_noaa
print, 'Error: ', ratio_err_ss, ratio_err_noaa, ratio_err_ogi
	
;remove the last element of the sim arrays because we don't have 2015 data, only for the NOAA data set
sou_noaa_sim = RemoveRows( rotate( sou_noaa_sim, 1 ), n_elements( sou_noaa_sim ) - 1 )
nor_noaa_sim = RemoveRows( rotate( nor_noaa_sim, 1 ), n_elements( nor_noaa_sim ) - 1 )
sou_noaa_sim_01 = RemoveRows( rotate( sou_noaa_sim_01, 1 ), n_elements( sou_noaa_sim_01 ) - 1 )
nor_noaa_sim_01 = RemoveRows( rotate( nor_noaa_sim_01, 1 ), n_elements( nor_noaa_sim_01 ) - 1 )
sou_noaa_sim_02 = RemoveRows( rotate( sou_noaa_sim_02, 1 ), n_elements( sou_noaa_sim_02 ) - 1 )
nor_noaa_sim_02 = RemoveRows( rotate( nor_noaa_sim_02, 1 ), n_elements( nor_noaa_sim_02 ) - 1 )
sou_noaa_sim_03 = RemoveRows( rotate( sou_noaa_sim_03, 1 ), n_elements( sou_noaa_sim_03 ) - 1 )
nor_noaa_sim_03 = RemoveRows( rotate( nor_noaa_sim_03, 1 ), n_elements( nor_noaa_sim_03 ) - 1 )
sou_noaa_sim_04 = RemoveRows( rotate( sou_noaa_sim_04, 1 ), n_elements( sou_noaa_sim_04 ) - 1 )
nor_noaa_sim_04 = RemoveRows( rotate( nor_noaa_sim_04, 1 ), n_elements( nor_noaa_sim_04 ) - 1 )
sou_noaa_sim_05 = RemoveRows( rotate( sou_noaa_sim_05, 1 ), n_elements( sou_noaa_sim_05 ) - 1 )
nor_noaa_sim_05 = RemoveRows( rotate( nor_noaa_sim_05, 1 ), n_elements( nor_noaa_sim_05 ) - 1 )


cgDisplay, 1200, 850
;set up the plot
title = 'Temporal trend of Ethane from 1997 to 2015 from NOAA and UCI'
ytitle = 'Weighted mean of the integral (pptv)'
xtitle = 'Years'
cgPlot, nor_ss, /nodata, xrange = [ 1980 , 2016 ], xtitle = 'Years', yrange = [ 0 , 1400 ], $
	xStyle = 8, position = [0.15, 0.15, 0.9, 0.50], charsize = 2, $
	ytitle = 'Ethane mixing ratio (pptv)'

;plot nor hemisphere NOAA sim
cgPlot, year_series_noaa_sim, nor_noaa_sim, color = 'violet', /overplot, /NoErase, linestyle = 2, thick = 2.0
;plot sou hemisphere NOAA sim
cgPlot, year_series_noaa_sim, sou_noaa_sim, color = 'violet', /overplot, /NoErase, thick = 2.0
;plot nor hemisphere SS sim
cgPlot, year_series_ss, nor_ss_sim, color = 'violet', /overplot, /NoErase, linestyle = 2, thick = 2.0
;plot sou hemisphere SS sim
cgPlot, year_series_ss, sou_ss_sim, color = 'violet', /overplot, /NoErase, thick = 2.0
;plot nor hemisphere OGI sim
cgPlot, year_series_ogi, nor_ogi_sim, color = 'violet', /overplot, /NoErase, linestyle = 2, thick = 2.0
;plot sou hemisphere OGI sim
cgPlot, year_series_ogi, sou_ogi_sim, color = 'violet', /overplot, /NoErase, thick = 2.0
;plot nor hemisphere NOAA sim01
cgPlot, year_series_noaa_sim, nor_noaa_sim_01, color = 'forest green', /overplot, /NoErase, linestyle = 2, thick = 2.0
;plot sou hemisphere NOAA sim01
cgPlot, year_series_noaa_sim, sou_noaa_sim_01, color = 'forest green', /overplot, /NoErase, thick = 2.0
;plot nor hemisphere OGI sim01
cgPlot, year_series_ogi, nor_ogi_sim_01, color = 'forest green', /overplot, /NoErase, linestyle = 2, thick = 1.5, psym = 2
;plot sou hemisphere OGI sim01
cgPlot, year_series_ogi, sou_ogi_sim_01, color = 'forest green', /overplot, /NoErase, psym = 2
;plot nor hemisphere SS sim01
cgPlot, year_series_ss, nor_ss_sim_01, color = 'forest green', /overplot, /NoErase, linestyle = 2, thick = 2.0
;plot  sou hemisphere SS sim01
cgPlot, year_series_ss, sou_ss_sim_01, color = 'forest green', /overplot, /NoErase , thick = 2.0
;plot nor hemisphere NOAA sim02
cgPlot, year_series_noaa_sim, nor_noaa_sim_02, color = 'tomato', /overplot, /NoErase, linestyle = 2, thick = 2.0
;plot sou hemisphere NOAA sim02
cgPlot, year_series_noaa_sim, sou_noaa_sim_02, color = 'tomato', /overplot, /NoErase, thick = 2.0
;plot nor hemisphere OGI sim02
cgPlot, year_series_ogi, nor_ogi_sim_02, color = 'tomato', /overplot, /NoErase, linestyle = 2, thick = 1.5, psym = 2
;plot sou hemisphere OGI sim02
cgPlot, year_series_ogi, sou_ogi_sim_02, color = 'tomato', /overplot, /NoErase, psym = 2
;plot nor hemisphere SS sim02
cgPlot, year_series_ss, nor_ss_sim_02, color = 'tomato', /overplot, /NoErase, linestyle = 2, thick = 2.0
;plot  sou hemisphere SS sim02
cgPlot, year_series_ss, sou_ss_sim_02, color = 'tomato', /overplot, /NoErase , thick = 2.0
;plot nor hemisphere NOAA sim03
cgPlot, year_series_noaa_sim, nor_noaa_sim_03, color = 'brown', /overplot, /NoErase, linestyle = 2, thick = 2.0
;plot sou hemisphere NOAA sim03
cgPlot, year_series_noaa_sim, sou_noaa_sim_03, color = 'brown', /overplot, /NoErase, thick = 2.0
;plot nor hemisphere OGI sim03
cgPlot, year_series_ogi, nor_ogi_sim_03, color = 'brown', /overplot, /NoErase, linestyle = 2, thick = 1.5, psym = 2
;plot sou hemisphere OGI sim03
cgPlot, year_series_ogi, sou_ogi_sim_03, color = 'brown', /overplot, /NoErase, psym = 2
;plot nor hemisphere SS sim03
cgPlot, year_series_ss, nor_ss_sim_03, color = 'brown', /overplot, /NoErase, linestyle = 2, thick = 2.0
;plot  sou hemisphere SS sim03
cgPlot, year_series_ss, sou_ss_sim_03, color = 'brown', /overplot, /NoErase , thick = 2.0
;plot nor hemisphere NOAA sim04
cgPlot, year_series_noaa_sim, nor_noaa_sim_04, color = 'maroon', /overplot, /NoErase, linestyle = 2, thick = 2.0
;plot sou hemisphere NOAA sim04
cgPlot, year_series_noaa_sim, sou_noaa_sim_04, color = 'maroon', /overplot, /NoErase, thick = 2.0
;plot nor hemisphere OGI sim04
cgPlot, year_series_ogi, nor_ogi_sim_04, color = 'maroon', /overplot, /NoErase, linestyle = 2, thick = 1.5, psym = 2
;plot sou hemisphere OGI sim04
cgPlot, year_series_ogi, sou_ogi_sim_04, color = 'maroon', /overplot, /NoErase, psym = 2
;plot nor hemisphere SS sim04
cgPlot, year_series_ss, nor_ss_sim_04, color = 'maroon', /overplot, /NoErase, linestyle = 2, thick = 2.0
;plot  sou hemisphere SS sim04
cgPlot, year_series_ss, sou_ss_sim_04, color = 'maroon', /overplot, /NoErase , thick = 2.0
;plot nor hemisphere NOAA sim05
cgPlot, year_series_noaa_sim, nor_noaa_sim_05, color = 'gray', /overplot, /NoErase, linestyle = 2, thick = 2.0
;plot sou hemisphere NOAA sim05
cgPlot, year_series_noaa_sim, sou_noaa_sim_05, color = 'gray', /overplot, /NoErase, thick = 2.0
;plot nor hemisphere OGI sim05
cgPlot, year_series_ogi, nor_ogi_sim_05, color = 'gray', /overplot, /NoErase, linestyle = 2, thick = 1.5, psym = 2
;plot sou hemisphere OGI sim05
cgPlot, year_series_ogi, sou_ogi_sim_05, color = 'gray', /overplot, /NoErase, psym = 2
;plot nor hemisphere SS sim05
cgPlot, year_series_ss, nor_ss_sim_05, color = 'gray', /overplot, /NoErase, linestyle = 2, thick = 2.0
;plot  sou hemisphere SS sim05
cgPlot, year_series_ss, sou_ss_sim_05, color = 'gray', /overplot, /NoErase , thick = 2.0
;plot nor hemisphere Simpson data
cgPlot, year_series_ss, nor_ss, psym = 2, color = 'red', /overplot, err_yhigh = nor_err_ss, err_ylow = nor_err_ss, $
	/err_clip, /NoErase
;plot sou hemisphere Simpson data
cgPlot, year_series_ss, sou_ss, psym = 18, color = 'red', /overplot, err_yhigh = sou_err_ss, err_ylow = sou_err_ss, $
	/err_clip, /NoErase
;plot nor hemisphere NOAA
cgPlot, year_series_noaa, nor_noaa, psym = 2, color = 'blue', /overplot, err_yhigh = nor_err_noaa, err_ylow = nor_err_noaa, $
	/err_clip, /NoErase
;plot sou hemisphere NOAA
cgPlot, year_series_noaa, sou_noaa, psym = 18, color = 'blue', /overplot, err_yhigh = sou_err_noaa, err_ylow = sou_err_noaa, $
	/err_clip, /NoErase
;plot nor hemisphere OGI
cgPlot, year_series_ogi, nor_ogi, psym = 2, color = 'forest green', /overplot, err_yhigh = nor_err_ogi, err_ylow = nor_err_ogi, $
	/err_clip, /NoErase
;plot sou hemisphere OGI\
cgPlot, year_series_ogi, sou_ogi, psym = 18, color = 'forest green', /overplot, err_yhigh = sou_err_ogi, err_ylow = sou_err_ogi, $
	/err_clip, /NoErase


;calculating the northern and southern hemispheric ratio:
hem_ratio_noaa = nor_noaa[ * ] / sou_noaa[ * ]
hem_ratio_ss = nor_ss[ * ] / sou_ss[ * ]
hem_ratio_ogi = nor_ogi[ * ] / sou_ogi[ * ]
hem_ratio_noaa_sim = nor_noaa_sim / sou_noaa_sim
hem_ratio_ogi_sim = nor_ogi_sim / sou_ogi_sim 
hem_ratio_ss_sim = nor_ss_sim / sou_ss_sim
hem_ratio_noaa_sim01 = nor_noaa_sim_01 / sou_noaa_sim_01
hem_ratio_ogi_sim01 = nor_ogi_sim_01 / sou_ogi_sim_01
hem_ratio_ss_sim01 = nor_ss_sim_01 / sou_ss_sim_01
hem_ratio_ss_sim02 = nor_ss_sim_02 / sou_ss_sim_02
hem_ratio_noaa_sim02 = nor_noaa_sim_02 / sou_noaa_sim_02
hem_ratio_ogi_sim02 = nor_ogi_sim_02 / sou_ogi_sim_02

;calculating the northern and southern hemispheric differences
diff_noaa = nor_noaa - sou_noaa
diff_ss = nor_ss - sou_ss
diff_ogi = nor_ogi - sou_ogi
diff_noaa_sim = nor_noaa_sim - sou_noaa_sim
diff_ss_sim = nor_ss_sim - sou_ss_sim
diff_ogi_sim = nor_ogi_sim - sou_ogi_sim 
diff_noaa_sim01 = nor_noaa_sim_01 - sou_noaa_sim_01
diff_ss_sim01 = nor_ss_sim_01 - sou_ss_sim_01
diff_ogi_sim01 = nor_ogi_sim_01 - sou_ogi_sim_01
diff_noaa_sim02 = nor_noaa_sim_02 - sou_noaa_sim_02
diff_ss_sim02 = nor_ss_sim_02 - sou_ss_sim_02
diff_ogi_sim02 = nor_ogi_sim_02 - sou_ogi_sim_02
diff_noaa_sim03 = nor_noaa_sim_03 - sou_noaa_sim_03
diff_ss_sim03 = nor_ss_sim_03 - sou_ss_sim_03
diff_ogi_sim03 = nor_ogi_sim_03 - sou_ogi_sim_03
diff_noaa_sim04 = nor_noaa_sim_04 - sou_noaa_sim_04
diff_ss_sim04 = nor_ss_sim_04 - sou_ss_sim_04
diff_ogi_sim04 = nor_ogi_sim_04 - sou_ogi_sim_04
diff_noaa_sim05 = nor_noaa_sim_05 - sou_noaa_sim_05
diff_ss_sim05 = nor_ss_sim_05 - sou_ss_sim_05
diff_ogi_sim05 = nor_ogi_sim_05 - sou_ogi_sim_05

err_diff_noaa = sqrt( nor_err_noaa ^ 2 + sou_err_noaa ^ 2 )
err_diff_ss = sqrt( nor_err_ss ^ 2 + sou_err_ss ^ 2 )
err_diff_ogi = sqrt( nor_err_ogi ^ 2 + sou_err_ss ^ 2 ) 


;print, hem_ratio_noaa, hem_ratio_ss, hem_ratio_ogi, hem_ratio_noaa_sim, hem_ratio_ogi_sim, hem_ratio_ss_sim, hem_ratio_noaa_sim01, $
;	hem_ratio_ogi_sim01, hem_ratio_ss_sim01, hem_ratio_ss_sim02, hem_ratio_noaa_sim02, hem_ratio_ogi_sim02
;cgPlot, year_series_ogi, hem_ratio_ogi, psym = 14, ytitle = 'Northen / Southern Hemispheric Ratio', XTickformat = '(A1)', $
;	position = [0.15 , 0.52 , 0.9 , 0.9], yrange = [ 2 , 15 ], color = 'forest green', xrange = [ 1980 , 2016 ], /NoErase, $
;	err_ylow = ratio_err_ogi, err_yhigh = ratio_err_ogi, /err_clip
;cgPlot, year_series_noaa, hem_ratio_noaa, psym = 16, /NoErase, color = 'forest green', /overplot, err_ylow = ratio_err_noaa, $
;	err_yhigh = ratio_err_noaa, /err_clip
;cgPlot, year_series_ss, hem_ratio_ss, psym = 15, /overplot, /NoErase, color = 'forest green', err_ylow = ratio_err_ss, $
;	err_yhigh = ratio_err_ss, /err_clip
;cgPlot, year_series_noaa, hem_ratio_noaa_sim, psym = -3, /overplot, /NoErase, color = 'black', thick = 2.0
;cgPlot, year_series_ss, hem_ratio_ss_sim, psym = -3, /overplot, /NoErase, color = 'black', thick = 2.0
;cgPlot, year_series_ogi, hem_ratio_ogi_sim, psym = 2, /overplot, /NoErase, color = 'black', thick = 2.0
;cgPlot, year_series_noaa, hem_ratio_noaa_sim01, psym = -3, /overplot, /NoErase, color = 'violet', thick = 2.0
;cgPlot, year_series_ogi, hem_ratio_ogi_sim01, psym = 2, /overplot, /NoErase, color = 'violet', thick = 2.0
;cgPlot, year_series_ss, hem_ratio_ss_sim01, psym = -3, /overplot, /NoErase, color = 'violet', thick = 2.0
;cgPlot, year_series_noaa, hem_ratio_noaa_sim02, psym = -3, /overplot, /NoErase, color = 'maroon', thick = 2.0
;cgPlot, year_series_ogi, hem_ratio_ogi_sim02, psym = 2, /overplot, /NoErase, color = 'maroon', thick = 2.0
;cgPlot, year_series_ss, hem_ratio_ss_sim02, psym = -3, /overplot, /NoErase, color = 'maroon', thick = 2.0

cgPlot, year_series_ogi, diff_ogi, psym = 14, XTickformat = '(A1)', $
	position = [0.15 , 0.52 , 0.9 , 0.9], color = 'black', xrange = [ 1980 , 2016 ], /NoErase, yrange = [ 400, 1200], $
	err_ylow = err_diff_ogi, err_yhigh = err_diff_ogi, /err_clip, symsize = 1.5, charsize = 2, $
	ytitle = 'Ethane mixing ratio (pptv)' 

cgPlot, year_series_noaa_sim, diff_noaa_sim, color = 'violet', /overplot, /NoErase, thick = 2.0
cgPlot, year_series_ogi, diff_ogi_sim, color = 'violet', /overplot, /NoErase, thick = 2.0, psym = 1
cgPlot, year_series_ss, diff_ss_sim, color = 'violet', /overplot, /NoErase, thick = 2.0
cgPlot, year_series_noaa_sim, diff_noaa_sim01, color = 'forest green', /overplot, /NoErase, thick = 2.0
cgPlot, year_series_ogi, diff_ogi_sim01, color = 'forest green', /overplot, /NoErase, thick = 2.0, psym = 1
cgPlot, year_series_ss, diff_ss_sim01, color = 'forest green', /overplot, /NoErase, thick = 2.0
cgPlot, year_series_noaa_sim, diff_noaa_sim02, color = 'tomato', /overplot, /NoErase, thick = 2.0
cgPlot, year_series_ogi, diff_ogi_sim02, color = 'tomato', /overplot, /NoErase, thick = 2.0, psym = 1
cgPlot, year_series_ss, diff_ss_sim02, color = 'tomato', /overplot, /NoErase, thick = 2.0
cgPlot, year_series_noaa_sim, diff_noaa_sim03, color = 'brown', /overplot, /NoErase, thick = 2.0
cgPlot, year_series_ogi, diff_ogi_sim03, color = 'brown', /overplot, /NoErase, thick = 2.0, psym = 1
cgPlot, year_series_ss, diff_ss_sim03, color = 'brown', /overplot, /NoErase, thick = 2.0
cgPlot, year_series_noaa_sim, diff_noaa_sim04, color = 'maroon', /overplot, /NoErase, thick = 2.0
cgPlot, year_series_ogi, diff_ogi_sim04, color = 'maroon', /overplot, /NoErase, thick = 2.0, psym = 1
cgPlot, year_series_ss, diff_ss_sim04, color = 'maroon', /overplot, /NoErase, thick = 2.0
cgPlot, year_series_noaa_sim, diff_noaa_sim05, color = 'gray', /overplot, /NoErase, thick = 2.0
cgPlot, year_series_ogi, diff_ogi_sim05, color = 'gray', /overplot, /NoErase, thick = 2.0, psym = 1
cgPlot, year_series_ss, diff_ss_sim05, color = 'gray', /overplot, /NoErase, thick = 2.0
cgPlot, year_series_noaa, diff_noaa, psym = 16, err_yhigh = err_diff_noaa, err_ylow = err_diff_noaa, /err_clip, /overplot, /NoErase, $
	symsize = 1.5
cgPlot, year_series_ss, diff_ss, psym = 15, err_yhigh = err_diff_ss, err_ylow = err_diff_ss, /err_clip, /overplot, /NoErase, $
	symsize = 1.5

cgLEGEND, title = ['Base emission (Xiao et al.)', 'Rice et al. normalized', 'Rice et al.', 'Aydin et al. normalized', 'Aydin et al.', 'Simpson et al.'], $
	/Box, /Background, BG_Color = 'rose', location = [0.20, 0.90], $
	color = ['tomato', 'violet', 'gray', 'forest green', 'brown', 'maroon']
	
Plots, [0.15, 0.15], [0.50, 0.52], /Normal ; Fix left axis.
Plots, [0.90, 0.90], [0.50, 0.52], /Normal ; Fix right axis.


;cgLEGEND, title = ['NOAA Northern Hem. averages', 'NOAA Southern Hem. averages', 'UCI Northern Hem. averages', 'UCI Southern Hem. averages', $
;	'OGI Northern Hem. average', 'OGI Southern Hem. average'], psym = [ -16, -16, -15, -15, -14, -14], $
;	color = [ 'red', 'blue', 'red', 'blue', 'red', 'blue'], location = [ 0.20, 0.49], /Box, /Background, BG_Color = 'rose', Length = 0.0
;	
;cgLEGEND, title = ['NOAA Hemispheric Ratios', 'UCI Hemispheric Ratios', 'OGI Hemispheric Ratio'], psym = [-16, -15, -14], $
;	color = ['forest green', 'forest green', 'forest green'], location = [ 0.20, 0.88 ], /Box, /Background, BG_Color = 'rose', length = 0.0

;Make a JPEG for the graphs
READ, option, PROMPT= 'Do I need to save this graph as a jpeg so I can print it out later? Enter "1" for HELLYEAH or "0" for NOPE       '
name = strarr(1)
IF (option EQ 1) THEN BEGIN
	READ, name, PROMPT= 'Name of the JPEG file? (No need to include the jpeg extension, just name is fine)   '
	SCREEN2JPG, name
	print, 'JPEG has been saved to current directory'
ENDIF





END
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
