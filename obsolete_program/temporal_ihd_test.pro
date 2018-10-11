PRO TEMPORAL_IHD_TEST

;this is a test for the new algorithm to calculate the IHD of the sim and obs data.
;the main difference of this test and the normal version is that this script test out 
;the new way that the sim data is called out thus reduces run time.
;for the purpose of this test. To simplify the problem, only deal with the NOAA data for now.

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
;version 2.0: change the variable naming scheme to bin numbers instead of year numbers so if needed be, I can do a 2 - year or 5 - year mean
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
			bin04_de_idx = where( temp_de_year EQ 2012 OR temp_de_year EQ 2013 ) 
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

;-----------This new section using a new algorithm to call the sim data for NOAA-----------------

;specify file name 
filename = "/home/excluded-from-backup/data/C2H6/trac_avg.PSUSF_1981_2015.bpch"

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

;This new section is new for this new algorithm
read, lat_var, prompt = 'Enter lat: '
read, lon_var, prompt = 'Enter lon: '

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
	print, i
	help, data
	for k = 0, n_elements( data[ * , 0 ] ) - 1 do begin
		for j = 0, n_elements( data[ 0 , * ] ) - 1 do begin
			simArr[ i , k , j ] = data[ k , j ]
		endfor
	endfor
endfor

;this changes the lat and lon values to index values to call the data
;out from the sim array 
lat_idx = round( ( lat_var - (-90) ) / 2 )
lon_idx = round( ( lon_var - (-180) ) / 2.5 )
print, 'Lat and Lon indecies: ', lat_idx, lon_idx
print, simArr[ sim_bin01_idx_noaa, lon_idx, lat_idx ] / 2 * 1000


END

