PRO HEM_RATIO_INTEGRAL_METHOD

;this program calculate the hemispheric ratio by calculating weighted mean of the 
;integral of the poly fit

;it is built on average_all_data_lat_dis

COMPILE_OPT IDL2				;set compile options
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
;=====Finish reading in NOAA data files=====
;============================================================================================================
;============================================================================================================
;=====Start reading in OGI data files=====
;get the files
infile1= "/home/jakechung/IDL/myIDL/KhalilEtAl_data_set.dat"
infile2= "/home/jakechung/IDL/myIDL/KhalilEtAl_locations.dat"

;get number of lines in the data file
n1_ogi = file_lines(infile1)
n2_ogi = file_lines(infile2)

print, 'The input has ', n1_ogi, '   data lines with', n2_ogi, '  locations'

;create an empty array to store the data that was read in lon, lat, avg, date and locations
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

;Convert unix time to decimal year
time_ogi = (unix_date_ogi/3600 - 131496)/8760 + 1985
;=====Finish reading in OGI data files=====
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
;============================================================================================================
;=====Finish reading in all the data files for the three data sets
;============================================================================================================
;PART TWO===================================================================================================
;SUB PART TWO-ONE===========================================================================================
;Dealing with the NOAA data====================================================================================

;I'm not smart enough, so I'm going to do a cheat  In order to write the extracted measured data out to seperate arrays, it will write to a temporary ASCII file
;and then read that file in outside of the loop to extract it into an array. It is not very elegant but it should work just fine

;Set the link to the temp files
;These files will hold the extracted data
infile_jan = "/home/jakechung/IDL/myIDL/temp_folder/temp_jan_NOAA.dat"
infile_feb = "/home/jakechung/IDL/myIDL/temp_folder/temp_feb_NOAA.dat"
infile_mar = "/home/jakechung/IDL/myIDL/temp_folder/temp_mar_NOAA.dat"
infile_apr = "/home/jakechung/IDL/myIDL/temp_folder/temp_apr_NOAA.dat"
infile_may = "/home/jakechung/IDL/myIDL/temp_folder/temp_may_NOAA.dat"
infile_jun = "/home/jakechung/IDL/myIDL/temp_folder/temp_jun_NOAA.dat"
infile_jul = "/home/jakechung/IDL/myIDL/temp_folder/temp_jul_NOAA.dat"
infile_aug = "/home/jakechung/IDL/myIDL/temp_folder/temp_aug_NOAA.dat"
infile_sep = "/home/jakechung/IDL/myIDL/temp_folder/temp_sep_NOAA.dat"
infile_oct = "/home/jakechung/IDL/myIDL/temp_folder/temp_oct_NOAA.dat"
infile_nov = "/home/jakechung/IDL/myIDL/temp_folder/temp_nov_NOAA.dat"
infile_dec = "/home/jakechung/IDL/myIDL/temp_folder/temp_dec_NOAA.dat"

;open to write
openw , lun1 , infile_jan , /get_lun
openw , lun2 , infile_feb , /get_lun
openw , lun3 , infile_mar , /get_lun
openw , lun4 , infile_apr , /get_lun
openw , lun5 , infile_may , /get_lun
openw , lun6 , infile_jun , /get_lun
openw , lun7 , infile_jul , /get_lun
openw , lun8 , infile_aug , /get_lun
openw , lun9 , infile_sep , /get_lun
openw , lun10 , infile_oct , /get_lun
openw , lun11 , infile_nov , /get_lun
openw , lun12 , infile_dec , /get_lun

;loop through the latitude array to find matching locations indices and pull out values from the other arrays
for j = 0, n2_noaa - 1 do begin
	if ( j NE n2_noaa - 1 ) then b = j + 1 ELSE b = j
	if (strcmp(location_noaa[ j ], location_noaa[ b ], /FOLD_CASE) EQ 1) then begin
	;comparing the locations
		;varible "a" holds the indices ONE specific station in the loop.
		a = where(strmatch(location_noaa[ * ], location_noaa[ j ], /FOLD_CASE) EQ 1)
		if (n_elements(a) GT 20) then begin ;a check to filter out stations with not enough values to do calculations
		;this line is from the Simpson data files, but since NOAA data is richer than Simpson, so all stations should be able 
		;to clear this filter
			temp_avg = avg_noaa[ a ]
			;pulling out monthly data at a specific station
			jan = where( ( month_noaa[ a ] EQ 1 ), count )
				if ( count GT 0 ) then begin 
				;putting a 'if then' loop here to filter out false values, if the
				;where funtion cannot find match values, it will return -1 values which
				;will cause the temp_avg array to calculate incorrect values, could have used /NULL function here
				;but oh well!
					obs_jan_var = mean( temp_avg[ jan ] )
				endif else begin
							obs_jan_var = !Values.F_NAN
						endelse 
			feb = where( ( month_noaa[ a ] EQ 2 ), count )
				if ( count GT 0 ) then begin 
					obs_feb_var = mean( temp_avg[ feb ] )
				endif else begin
							obs_feb_var = !Values.F_NAN
						endelse 
			mar = where( ( month_noaa[ a ] EQ 3 ), count )
				if ( count GT 0 ) then begin
					obs_mar_var = mean( temp_avg[ mar ] )
				endif else begin
							obs_mar_var = !Values.F_NAN
						endelse 
			apr = where( ( month_noaa[ a ] EQ 4 ), count )
				if ( count GT 0 ) then begin
					obs_apr_var = mean( temp_avg[ apr ] )
				endif else begin
							obs_apr_var = !Values.F_NAN
						endelse 
			may = where( ( month_noaa[ a ] EQ 5 ), count )
				if ( count GT 0 ) then begin
					obs_may_var = mean( temp_avg[ may ] )
				endif else begin
							obs_may_var = !Values.F_NAN
						endelse 					
			jun = where( ( month_noaa[ a ] EQ 6 ), count )
				if ( count GT 0 ) then begin
					obs_jun_var = mean( temp_avg[ jun ] )
				endif else begin
							obs_jun_var = !Values.F_NAN
						endelse 
			jul = where( ( month_noaa[ a ] EQ 7 ), count )
				if ( count GT 0 ) then begin
					obs_jul_var = mean( temp_avg[ jul ] )
				endif else begin
							obs_jul_var = !Values.F_NAN
						endelse 
			aug = where( ( month_noaa[ a ] EQ 8 ), count )
				if ( count GT 0 ) then begin
					obs_aug_var = mean( temp_avg[ aug ] )
				endif else begin
							obs_aug_var = !Values.F_NAN
						endelse 
			sep = where( ( month_noaa[ a ] EQ 9 ), count )
				if ( count GT 0 ) then begin
					obs_sep_var = mean( temp_avg[ sep ] )
				endif else begin
							obs_sep_var = !Values.F_NAN
						endelse 
			oct = where( ( month_noaa[ a ] EQ 10 ), count )
				if ( count GT 0 ) then begin
					obs_oct_var = mean( temp_avg[ oct ] )
				endif else begin
							obs_oct_var = !Values.F_NAN
						endelse 
			nov = where( ( month_noaa[ a ] EQ 11 ), count )
				if ( count GT 0 ) then begin
					obs_nov_var = mean( temp_avg[ nov ] )
				endif else begin
							obs_nov_var = !Values.F_NAN
						endelse 
			dec = where( ( month_noaa[ a ] EQ 12 ), count )
				if ( count GT 0 ) then begin
					obs_dec_var = mean( temp_avg[ dec ] )
				endif else begin
							obs_dec_var = !Values.F_NAN
						endelse 
			;writing it to the temp files
			printf, lun1, obs_jan_var, lat1_noaa[ j ], lon1_noaa[ j ]
			printf, lun2, obs_feb_var, lat1_noaa[ j ], lon1_noaa[ j ]
			printf, lun3, obs_mar_var, lat1_noaa[ j ], lon1_noaa[ j ]
			printf, lun4, obs_apr_var, lat1_noaa[ j ], lon1_noaa[ j ]
			printf, lun5, obs_may_var, lat1_noaa[ j ], lon1_noaa[ j ]
			printf, lun6, obs_jun_var, lat1_noaa[ j ], lon1_noaa[ j ]
			printf, lun7, obs_jul_var, lat1_noaa[ j ], lon1_noaa[ j ]
			printf, lun8, obs_aug_var, lat1_noaa[ j ], lon1_noaa[ j ]
			printf, lun9, obs_sep_var, lat1_noaa[ j ], lon1_noaa[ j ]
			printf, lun10, obs_oct_var, lat1_noaa[ j ], lon1_noaa[ j ]
			printf, lun11, obs_nov_var, lat1_noaa[ j ], lon1_noaa[ j ]
			printf, lun12, obs_dec_var, lat1_noaa[ j ], lon1_noaa[ j ]
		endif
		j = a[n_elements(a) - 1]
		a = 0
	endif else if ( strcmp( location_noaa[ j ], location_noaa[ b ], /FOLD_CASE) EQ 0) then j = j + 1
endfor

;free then and open then again to refresh the lun values
free_lun, lun1, lun2, lun3, lun4, lun5, lun6, lun7, lun8, lun9, lun10, lun11, lun12

;reopen the files to get the values pulled out from the loop above
;these lun values now has the values of all jan or all feb, so on and so on classifiied station by station
openr, lun1, infile_jan, /get_lun
openr, lun2, infile_feb, /get_lun
openr, lun3, infile_mar, /get_lun
openr, lun4, infile_apr, /get_lun
openr, lun5, infile_may, /get_lun
openr, lun6, infile_jun, /get_lun
openr, lun7, infile_jul, /get_lun
openr, lun8, infile_aug, /get_lun
openr, lun9, infile_sep, /get_lun
openr, lun10, infile_oct, /get_lun
openr, lun11, infile_nov, /get_lun
openr, lun12, infile_dec, /get_lun

;get the sizes of the files
jan_size = file_lines( infile_jan )
feb_size = file_lines( infile_feb )
mar_size = file_lines( infile_mar )
apr_size = file_lines( infile_apr )
may_size = file_lines( infile_may )
jun_size = file_lines( infile_jun )
jul_size = file_lines( infile_jul )
aug_size = file_lines( infile_aug )
sep_size = file_lines( infile_sep )
oct_size = file_lines( infile_oct )
nov_size = file_lines( infile_nov )
dec_size = file_lines( infile_dec )

;creating monthly super arrays to store monthly data 
;again, all the stations are seperated by lon and lat values
;data structure as follow 
; [0 , *]: concentration
; [1 , *]: latitude
; [2, *]: longitude
obs_jan_arr_noaa = fltarr(3, jan_size)
obs_feb_arr_noaa = fltarr(3, feb_size)
obs_mar_arr_noaa = fltarr(3, mar_size)
obs_apr_arr_noaa = fltarr(3, apr_size)
obs_may_arr_noaa = fltarr(3, may_size)
obs_jun_arr_noaa = fltarr(3, jun_size)
obs_jul_arr_noaa = fltarr(3, jul_size)
obs_aug_arr_noaa = fltarr(3, aug_size)
obs_sep_arr_noaa = fltarr(3, sep_size)
obs_oct_arr_noaa = fltarr(3, oct_size)
obs_nov_arr_noaa = fltarr(3, nov_size)
obs_dec_arr_noaa = fltarr(3, dec_size)

;this section is a super loop to read in the temp files created above and seperate out all the months 
for t = 0, jan_size - 1 do begin
	readf, lun1, avg_var, lat_var, lon_var
	obs_jan_arr_noaa[ 0 , t ] = avg_var
	obs_jan_arr_noaa[ 1 , t ] = lat_var
	obs_jan_arr_noaa[ 2 , t ] = lon_var
endfor
for t = 0, feb_size - 1 do begin
	readf, lun2, avg_var, lat_var, lon_var
	obs_feb_arr_noaa[ 0 , t ] = avg_var
	obs_feb_arr_noaa[ 1 , t ] = lat_var
	obs_feb_arr_noaa[ 2 , t ] = lon_var
endfor
for t = 0, mar_size - 1 do begin
	readf, lun3, avg_var, lat_var, lon_var
	obs_mar_arr_noaa[ 0 , t ] = avg_var
	obs_mar_arr_noaa[ 1 , t ] = lat_var
	obs_mar_arr_noaa[ 2 , t ] = lon_var
endfor
for t = 0, apr_size - 1 do begin
	readf, lun4, avg_var, lat_var, lon_var
	obs_apr_arr_noaa[ 0 , t ] = avg_var
	obs_apr_arr_noaa[ 1 , t ] = lat_var
	obs_apr_arr_noaa[ 2 , t ] = lon_var
endfor
for t = 0, may_size - 1 do begin
	readf, lun5, avg_var, lat_var, lon_var
	obs_may_arr_noaa[ 0 , t ] = avg_var
	obs_may_arr_noaa[ 1 , t ] = lat_var
	obs_may_arr_noaa[ 2 , t ] = lon_var
endfor
for t = 0, jun_size - 1 do begin
	readf, lun6, avg_var, lat_var, lon_var
	obs_jun_arr_noaa[ 0 , t ] = avg_var
	obs_jun_arr_noaa[ 1 , t ] = lat_var
	obs_jun_arr_noaa[ 2 , t ] = lon_var
endfor
for t = 0, jul_size - 1 do begin
	readf, lun7, avg_var, lat_var, lon_var
	obs_jul_arr_noaa[ 0 , t ] = avg_var
	obs_jul_arr_noaa[ 1 , t ] = lat_var
	obs_jul_arr_noaa[ 2 , t ] = lon_var
endfor
for t = 0, aug_size - 1 do begin
	readf, lun8, avg_var, lat_var, lon_var
	obs_aug_arr_noaa[ 0 , t ] = avg_var
	obs_aug_arr_noaa[ 1 , t ] = lat_var
	obs_aug_arr_noaa[ 2 , t ] = lon_var
endfor
for t = 0, sep_size - 1 do begin
	readf, lun9, avg_var, lat_var, lon_var
	obs_sep_arr_noaa[ 0 , t ] = avg_var
	obs_sep_arr_noaa[ 1 , t ] = lat_var
	obs_sep_arr_noaa[ 2 , t ] = lon_var
endfor
for t = 0, oct_size - 1 do begin
	readf, lun10, avg_var, lat_var, lon_var
	obs_oct_arr_noaa[ 0 , t ] = avg_var
	obs_oct_arr_noaa[ 1 , t ] = lat_var
	obs_oct_arr_noaa[ 2 , t ] = lon_var
endfor
for t = 0, nov_size - 1 do begin
	readf, lun11, avg_var, lat_var, lon_var
	obs_nov_arr_noaa[ 0 , t ] = avg_var
	obs_nov_arr_noaa[ 1 , t ] = lat_var
	obs_nov_arr_noaa[ 2 , t ] = lon_var
endfor
for t = 0, dec_size - 1 do begin
	readf, lun12, avg_var, lat_var, lon_var
	obs_dec_arr_noaa[ 0 , t ] = avg_var
	obs_dec_arr_noaa[ 1 , t ] = lat_var
	obs_dec_arr_noaa[ 2 , t ] = lon_var
endfor

free_lun, lun1, lun2, lun3, lun4, lun5, lun6, lun7, lun8, lun9, lun10, lun11, lun12

;Finish SUB PART TWO-ONE
;Start SUB PART TWO-TWO===================================================================================================
;Dealing with the Simpson data
;The Simpson et al. data will go through the exact same process as the NOAA data with different variables and temp data files
;The Simpson et al. data only has 4 months data to deal with March, June, September and December

;Set the link to the temp files
;These files will hold the extracted data
infile_mar = "/home/jakechung/IDL/myIDL/temp_folder/temp_mar_Simpson.dat"
infile_jun = "/home/jakechung/IDL/myIDL/temp_folder/temp_jun_Simpson.dat"
infile_sep = "/home/jakechung/IDL/myIDL/temp_folder/temp_sep_Simpson.dat"
infile_dec = "/home/jakechung/IDL/myIDL/temp_folder/temp_dec_Simpson.dat"

;open to write
openw , lun3 , infile_mar , /get_lun
openw , lun6 , infile_jun , /get_lun
openw , lun9 , infile_sep , /get_lun
openw , lun12 , infile_dec , /get_lun

;loop through the latitude array to find matching locations indices and pull out values from the other arrays
for j = 0, n2_ss - 1 do begin
	if ( j NE n2_ss - 1 ) then b = j + 1 ELSE b = j
	if (strcmp(lat1_ss[ j ], lat1_ss[ b ], /FOLD_CASE) EQ 1) then begin
	;comparing the locations
		;varible "a" holds the indices ONE specific station in the loop.
		a = where(strmatch(lat1_ss[ * ], lat1_ss[ j ], /FOLD_CASE) EQ 1)
		if (n_elements(a) GT 30) then begin ;a check to filter out stations with not enough values to do calculations
			temp_avg = avg_ss[ a ]
			;pulling out monthly data at a specific station
			mar = where( ( tau_time_ss.month[ a ] EQ 3 ), count )
				;putting a 'if then' loop here to filter out false values, if the
				;where funtion cannot find match values, it will return -1 values which
				;will cause the temp_avg array to calculate incorrect values, could have used /NULL function here
				;but oh well!
				if ( count GT 0 ) then begin
					obs_mar_var = mean( temp_avg[ mar ] )
				endif else begin
							obs_mar_var = !Values.F_NAN
						endelse 				
			jun = where( ( tau_time_ss.month[ a ] EQ 6 ), count )
				if ( count GT 0 ) then begin
					obs_jun_var = mean( temp_avg[ jun ] )
				endif else begin
							obs_jun_var = !Values.F_NAN
						endelse 
			sep = where( ( tau_time_ss.month[ a ] EQ 9 ), count )
				if ( count GT 0 ) then begin
					obs_sep_var = mean( temp_avg[ sep ] )
				endif else begin
							obs_sep_var = !Values.F_NAN
						endelse  
			dec = where( ( tau_time_ss.month[ a ] EQ 12 ), count )
				if ( count GT 0 ) then begin
					obs_dec_var = mean( temp_avg[ dec ] )
				endif else begin
							obs_dec_var = !Values.F_NAN
						endelse 
			;writing it to the temp files
			printf, lun3, obs_mar_var, lat1_ss[ j ], lon1_ss[ j ]
			printf, lun6, obs_jun_var, lat1_ss[ j ], lon1_ss[ j ]
			printf, lun9, obs_sep_var, lat1_ss[ j ], lon1_ss[ j ]
			printf, lun12, obs_dec_var, lat1_ss[ j ], lon1_ss[ j ]
		endif
		j = a[n_elements(a) - 1]
		a = 0
	endif else if ( strcmp( lat1_ss[ j ], lat1_ss[ b ], /FOLD_CASE) EQ 0) then j = j + 1
endfor

;free the values and open then again to refresh the lun values
free_lun, lun3, lun6, lun9, lun12

;reopen the files to get the values pulled out from the loop above
;these lun values now has the values of all jan or all feb, so on and so on classifiied station by station
openr, lun3, infile_mar, /get_lun
openr, lun6, infile_jun, /get_lun
openr, lun9, infile_sep, /get_lun
openr, lun12, infile_dec, /get_lun

;get the sizes of the files
mar_size = file_lines( infile_mar )
jun_size = file_lines( infile_jun )
sep_size = file_lines( infile_sep )
dec_size = file_lines( infile_dec )

;creating monthly super arrays to store monthly data 
;again, all the stations are seperated by lon and lat values
;data structure as follow 
; [0 , *]: concentration
; [1 , *]: latitude
; [2, *]: longitude
obs_mar_arr_ss = fltarr(3, mar_size)
obs_jun_arr_ss = fltarr(3, jun_size)
obs_sep_arr_ss = fltarr(3, sep_size)
obs_dec_arr_ss = fltarr(3, dec_size)

;this section is a super loop to read in the temp files created above and seperate out all the months 
for t = 0, mar_size - 1 do begin
	readf, lun3, avg_var, lat_var, lon_var
	obs_mar_arr_ss[ 0 , t ] = avg_var
	obs_mar_arr_ss[ 1 , t ] = lat_var
	obs_mar_arr_ss[ 2 , t ] = lon_var
endfor
for t = 0, jun_size - 1 do begin
	readf, lun6, avg_var, lat_var, lon_var
	obs_jun_arr_ss[ 0 , t ] = avg_var
	obs_jun_arr_ss[ 1 , t ] = lat_var
	obs_jun_arr_ss[ 2 , t ] = lon_var
endfor
for t = 0, sep_size - 1 do begin
	readf, lun9, avg_var, lat_var, lon_var
	obs_sep_arr_ss[ 0 , t ] = avg_var
	obs_sep_arr_ss[ 1 , t ] = lat_var
	obs_sep_arr_ss[ 2 , t ] = lon_var
endfor
for t = 0, dec_size - 1 do begin
	readf, lun12, avg_var, lat_var, lon_var
	obs_dec_arr_ss[ 0 , t ] = avg_var
	obs_dec_arr_ss[ 1 , t ] = lat_var
	obs_dec_arr_ss[ 2 , t ] = lon_var
endfor

free_lun, lun3, lun6, lun9, lun12

;Finish SUB PART TWO-TWO
;Start SUB PART TWO-THREE===================================================================================================
;Dealing with  the OGI data
;The OGI data will go through the exact same process as the NOAA data with different variables and temp data files

;Set the link to the temp files
;These files will hold the extracted data
infile_jan = "/home/jakechung/IDL/myIDL/temp_folder/temp_jan_OGI.dat"
infile_feb = "/home/jakechung/IDL/myIDL/temp_folder/temp_feb_OGI.dat"
infile_mar = "/home/jakechung/IDL/myIDL/temp_folder/temp_mar_OGI.dat"
infile_apr = "/home/jakechung/IDL/myIDL/temp_folder/temp_apr_OGI.dat"
infile_may = "/home/jakechung/IDL/myIDL/temp_folder/temp_may_OGI.dat"
infile_jun = "/home/jakechung/IDL/myIDL/temp_folder/temp_jun_OGI.dat"
infile_jul = "/home/jakechung/IDL/myIDL/temp_folder/temp_jul_OGI.dat"
infile_aug = "/home/jakechung/IDL/myIDL/temp_folder/temp_aug_OGI.dat"
infile_sep = "/home/jakechung/IDL/myIDL/temp_folder/temp_sep_OGI.dat"
infile_oct = "/home/jakechung/IDL/myIDL/temp_folder/temp_oct_OGI.dat"
infile_nov = "/home/jakechung/IDL/myIDL/temp_folder/temp_nov_OGI.dat"
infile_dec = "/home/jakechung/IDL/myIDL/temp_folder/temp_dec_OGI.dat"

;open to write
openw , lun1 , infile_jan , /get_lun
openw , lun2 , infile_feb , /get_lun
openw , lun3 , infile_mar , /get_lun
openw , lun4 , infile_apr , /get_lun
openw , lun5 , infile_may , /get_lun
openw , lun6 , infile_jun , /get_lun
openw , lun7 , infile_jul , /get_lun
openw , lun8 , infile_aug , /get_lun
openw , lun9 , infile_sep , /get_lun
openw , lun10 , infile_oct , /get_lun
openw , lun11 , infile_nov , /get_lun
openw , lun12 , infile_dec , /get_lun

;loop through the latitude array to find matching locations indices and pull out values from the other arrays
for j = 0, n2_ogi - 1 do begin
	if ( j NE n2_ogi - 1 ) then b = j + 1 ELSE b = j
	if (strcmp(location_ogi[ j ], location_ogi[ b ], /FOLD_CASE) EQ 1) then begin
	;comparing the locations
		;varible "a" holds the indices ONE specific station in the loop.
		a = where(strmatch(location_ogi[ * ], location_ogi[ j ], /FOLD_CASE) EQ 1)
		if (n_elements(a) GT 20) then begin ;a check to filter out stations with not enough values to do calculations
		;this line is from the Simpson data files, but since OGI data is richer than Simpson, so all stations should be able 
		;to clear this filter
			temp_avg = avg_ogi[ a ]
			;pulling out monthly data at a specific station
			jan = where( ( tau_time_ogi.month[ a ] EQ 1 ), count )
				if ( count GT 0 ) then begin 
				;putting a 'if then' loop here to filter out false values, if the
				;where funtion cannot find match values, it will return -1 values which
				;will cause the temp_avg array to calculate incorrect values, could have used /NULL function here
				;but oh well!
					obs_jan_var = mean( temp_avg[ jan ] )
				endif else begin
							obs_jan_var = !Values.F_NAN
						endelse 
			feb = where( ( tau_time_ogi.month[ a ] EQ 2 ), count )
				if ( count GT 0 ) then begin 
					obs_feb_var = mean( temp_avg[ feb ] )
				endif else begin
							obs_feb_var = !Values.F_NAN
						endelse 
			mar = where( ( tau_time_ogi.month[ a ] EQ 3 ), count )
				if ( count GT 0 ) then begin
					obs_mar_var = mean( temp_avg[ mar ] )
				endif else begin
							obs_mar_var = !Values.F_NAN
						endelse 
			apr = where( ( tau_time_ogi.month[ a ] EQ 4 ), count )
				if ( count GT 0 ) then begin
					obs_apr_var = mean( temp_avg[ apr ] )
				endif else begin
							obs_apr_var = !Values.F_NAN
						endelse 
			may = where( ( tau_time_ogi.month[ a ] EQ 5 ), count )
				if ( count GT 0 ) then begin
					obs_may_var = mean( temp_avg[ may ] )
				endif else begin
							obs_may_var = !Values.F_NAN
						endelse 					
			jun = where( ( tau_time_ogi.month[ a ] EQ 6 ), count )
				if ( count GT 0 ) then begin
					obs_jun_var = mean( temp_avg[ jun ] )
				endif else begin
							obs_jun_var = !Values.F_NAN
						endelse 
			jul = where( ( tau_time_ogi.month[ a ] EQ 7 ), count )
				if ( count GT 0 ) then begin
					obs_jul_var = mean( temp_avg[ jul ] )
				endif else begin
							obs_jul_var = !Values.F_NAN
						endelse 
			aug = where( ( tau_time_ogi.month[ a ] EQ 8 ), count )
				if ( count GT 0 ) then begin
					obs_aug_var = mean( temp_avg[ aug ] )
				endif else begin
							obs_aug_var = !Values.F_NAN
						endelse 
			sep = where( ( tau_time_ogi.month[ a ] EQ 9 ), count )
				if ( count GT 0 ) then begin
					obs_sep_var = mean( temp_avg[ sep ] )
				endif else begin
							obs_sep_var = !Values.F_NAN
						endelse 
			oct = where( ( tau_time_ogi.month[ a ] EQ 10 ), count )
				if ( count GT 0 ) then begin
					obs_oct_var = mean( temp_avg[ oct ] )
				endif else begin
							obs_oct_var = !Values.F_NAN
						endelse 
			nov = where( ( tau_time_ogi.month[ a ] EQ 11 ), count )
				if ( count GT 0 ) then begin
					obs_nov_var = mean( temp_avg[ nov ] )
				endif else begin
							obs_nov_var = !Values.F_NAN
						endelse 
			dec = where( ( tau_time_ogi.month[ a ] EQ 12 ), count )
				if ( count GT 0 ) then begin
					obs_dec_var = mean( temp_avg[ dec ] )
				endif else begin
							obs_dec_var = !Values.F_NAN
						endelse 
			;writing it to the temp files
			printf, lun1, obs_jan_var, lat1_ogi[ j ], lon1_ogi[ j ]
			printf, lun2, obs_feb_var, lat1_ogi[ j ], lon1_ogi[ j ]
			printf, lun3, obs_mar_var, lat1_ogi[ j ], lon1_ogi[ j ]
			printf, lun4, obs_apr_var, lat1_ogi[ j ], lon1_ogi[ j ]
			printf, lun5, obs_may_var, lat1_ogi[ j ], lon1_ogi[ j ]
			printf, lun6, obs_jun_var, lat1_ogi[ j ], lon1_ogi[ j ]
			printf, lun7, obs_jul_var, lat1_ogi[ j ], lon1_ogi[ j ]
			printf, lun8, obs_aug_var, lat1_ogi[ j ], lon1_ogi[ j ]
			printf, lun9, obs_sep_var, lat1_ogi[ j ], lon1_ogi[ j ]
			printf, lun10, obs_oct_var, lat1_ogi[ j ], lon1_ogi[ j ]
			printf, lun11, obs_nov_var, lat1_ogi[ j ], lon1_ogi[ j ]
			printf, lun12, obs_dec_var, lat1_ogi[ j ], lon1_ogi[ j ]
		endif
		j = a[n_elements(a) - 1]
		a = 0
	endif else if ( strcmp( location_ogi[ j ], location_ogi[ b ], /FOLD_CASE) EQ 0) then j = j + 1
endfor

;free then and open then again to refresh the lun values
free_lun, lun1, lun2, lun3, lun4, lun5, lun6, lun7, lun8, lun9, lun10, lun11, lun12

;reopen the files to get the values pulled out from the loop above
;these lun values now has the values of all jan or all feb, so on and so on classifiied station by station
openr, lun1, infile_jan, /get_lun
openr, lun2, infile_feb, /get_lun
openr, lun3, infile_mar, /get_lun
openr, lun4, infile_apr, /get_lun
openr, lun5, infile_may, /get_lun
openr, lun6, infile_jun, /get_lun
openr, lun7, infile_jul, /get_lun
openr, lun8, infile_aug, /get_lun
openr, lun9, infile_sep, /get_lun
openr, lun10, infile_oct, /get_lun
openr, lun11, infile_nov, /get_lun
openr, lun12, infile_dec, /get_lun

;get the sizes of the files
jan_size = file_lines( infile_jan )
feb_size = file_lines( infile_feb )
mar_size = file_lines( infile_mar )
apr_size = file_lines( infile_apr )
may_size = file_lines( infile_may )
jun_size = file_lines( infile_jun )
jul_size = file_lines( infile_jul )
aug_size = file_lines( infile_aug )
sep_size = file_lines( infile_sep )
oct_size = file_lines( infile_oct )
nov_size = file_lines( infile_nov )
dec_size = file_lines( infile_dec )

;creating monthly super arrays to store monthly data 
;again, all the stations are seperated by lon and lat values
;data structure as follow 
; [0 , *]: concentration
; [1 , *]: latitude
; [2, *]: longitude
obs_jan_arr_ogi = fltarr(3, jan_size)
obs_feb_arr_ogi = fltarr(3, feb_size)
obs_mar_arr_ogi = fltarr(3, mar_size)
obs_apr_arr_ogi = fltarr(3, apr_size)
obs_may_arr_ogi = fltarr(3, may_size)
obs_jun_arr_ogi = fltarr(3, jun_size)
obs_jul_arr_ogi = fltarr(3, jul_size)
obs_aug_arr_ogi = fltarr(3, aug_size)
obs_sep_arr_ogi = fltarr(3, sep_size)
obs_oct_arr_ogi = fltarr(3, oct_size)
obs_nov_arr_ogi = fltarr(3, nov_size)
obs_dec_arr_ogi = fltarr(3, dec_size)

;this section is a super loop to read in the temp files created above and seperate out all the months 
for t = 0, jan_size - 1 do begin
	readf, lun1, avg_var, lat_var, lon_var
	obs_jan_arr_ogi[ 0 , t ] = avg_var
	obs_jan_arr_ogi[ 1 , t ] = lat_var
	obs_jan_arr_ogi[ 2 , t ] = lon_var
endfor
for t = 0, feb_size - 1 do begin
	readf, lun2, avg_var, lat_var, lon_var
	obs_feb_arr_ogi[ 0 , t ] = avg_var
	obs_feb_arr_ogi[ 1 , t ] = lat_var
	obs_feb_arr_ogi[ 2 , t ] = lon_var
endfor
for t = 0, mar_size - 1 do begin
	readf, lun3, avg_var, lat_var, lon_var
	obs_mar_arr_ogi[ 0 , t ] = avg_var
	obs_mar_arr_ogi[ 1 , t ] = lat_var
	obs_mar_arr_ogi[ 2 , t ] = lon_var
endfor
for t = 0, apr_size - 1 do begin
	readf, lun4, avg_var, lat_var, lon_var
	obs_apr_arr_ogi[ 0 , t ] = avg_var
	obs_apr_arr_ogi[ 1 , t ] = lat_var
	obs_apr_arr_ogi[ 2 , t ] = lon_var
endfor
for t = 0, may_size - 1 do begin
	readf, lun5, avg_var, lat_var, lon_var
	obs_may_arr_ogi[ 0 , t ] = avg_var
	obs_may_arr_ogi[ 1 , t ] = lat_var
	obs_may_arr_ogi[ 2 , t ] = lon_var
endfor
for t = 0, jun_size - 1 do begin
	readf, lun6, avg_var, lat_var, lon_var
	obs_jun_arr_ogi[ 0 , t ] = avg_var
	obs_jun_arr_ogi[ 1 , t ] = lat_var
	obs_jun_arr_ogi[ 2 , t ] = lon_var
endfor
for t = 0, jul_size - 1 do begin
	readf, lun7, avg_var, lat_var, lon_var
	obs_jul_arr_ogi[ 0 , t ] = avg_var
	obs_jul_arr_ogi[ 1 , t ] = lat_var
	obs_jul_arr_ogi[ 2 , t ] = lon_var
endfor
for t = 0, aug_size - 1 do begin
	readf, lun8, avg_var, lat_var, lon_var
	obs_aug_arr_ogi[ 0 , t ] = avg_var
	obs_aug_arr_ogi[ 1 , t ] = lat_var
	obs_aug_arr_ogi[ 2 , t ] = lon_var
endfor
for t = 0, sep_size - 1 do begin
	readf, lun9, avg_var, lat_var, lon_var
	obs_sep_arr_ogi[ 0 , t ] = avg_var
	obs_sep_arr_ogi[ 1 , t ] = lat_var
	obs_sep_arr_ogi[ 2 , t ] = lon_var
endfor
for t = 0, oct_size - 1 do begin
	readf, lun10, avg_var, lat_var, lon_var
	obs_oct_arr_ogi[ 0 , t ] = avg_var
	obs_oct_arr_ogi[ 1 , t ] = lat_var
	obs_oct_arr_ogi[ 2 , t ] = lon_var
endfor
for t = 0, nov_size - 1 do begin
	readf, lun11, avg_var, lat_var, lon_var
	obs_nov_arr_ogi[ 0 , t ] = avg_var
	obs_nov_arr_ogi[ 1 , t ] = lat_var
	obs_nov_arr_ogi[ 2 , t ] = lon_var
endfor
for t = 0, dec_size - 1 do begin
	readf, lun12, avg_var, lat_var, lon_var
	obs_dec_arr_ogi[ 0 , t ] = avg_var
	obs_dec_arr_ogi[ 1 , t ] = lat_var
	obs_dec_arr_ogi[ 2 , t ] = lon_var
endfor

free_lun, lun1, lun2, lun3, lun4, lun5, lun6, lun7, lun8, lun9, lun10, lun11, lun12

;This is a special section that is added in only for this script to calculate stations averages
obs_stt_avg_noaa = fltarr( 3 , n_elements( obs_mar_arr_noaa[ 0 , * ] ) )
obs_stt_avg_ogi = fltarr( 3 , n_elements( obs_mar_arr_ogi[ 0 , * ] ) )
obs_stt_avg_ss = fltarr( 3 , n_elements( obs_mar_arr_ss[ 0 , * ] ) ) 

for i = 0, n_elements( obs_mar_arr_noaa[ 0 , * ] ) - 1 do begin
	obs_stt_avg_noaa[ 0 , i ] = obs_mar_arr_noaa[ 1 , i ]
	obs_stt_avg_noaa[ 1 , i ] = obs_mar_arr_noaa[ 2 , i ]
	obs_stt_avg_noaa[ 2 , i ] = mean( [ obs_mar_arr_noaa[ 0 , i ], obs_jun_arr_noaa[ 0 , i ], obs_sep_arr_noaa[ 0 , i ], $
		obs_dec_arr_noaa[ 0 , i ] ], /NAN )
endfor

for i = 0, n_elements( obs_mar_arr_ogi[ 0 , * ] ) - 1 do begin
	obs_stt_avg_ogi[ 0 , i ] = obs_mar_arr_ogi[ 1 , i ]
	obs_stt_avg_ogi[ 1 , i ] = obs_mar_arr_ogi[ 2 , i ]
	obs_stt_avg_ogi[ 2 , i ] = mean( [ obs_mar_arr_ogi[ 0 , i ], obs_jun_arr_ogi[ 0 , i ], obs_sep_arr_ogi[ 0 , i ], $
		obs_dec_arr_ogi[ 0 , i ] ], /NAN )
endfor

for i = 0, n_elements( obs_mar_arr_ss[ 0 , * ] ) - 1 do begin
	obs_stt_avg_ss[ 0 , i ] = obs_mar_arr_ss[ 1 , i ]
	obs_stt_avg_ss[ 1 , i ] = obs_mar_arr_ss[ 2 , i ]
	obs_stt_avg_ss[ 2 , i ] = mean( [ obs_mar_arr_ss[ 0 , i ], obs_jun_arr_ss[ 0 , i ], obs_sep_arr_ss[ 0 , i ], $
		obs_dec_arr_ss[ 0 , i ] ], /NAN )
endfor
;END of the special section


;Finish SUB PART TWO-THREE==================================================================================================
;============================================================================================================================
;============================================================================================================================
;PART THREE================================================================================================================
;Read in the deseasonal data
;Start SUB PART THREE-ONE--------Dealing with NOAA deseasonal data================================================================

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
infile_noaa_de = "/home/jakechung/IDL/myIDL/temp_folder/temp_average_NOAA_deseason.dat"

;open to write in the temp files
openw , lun1 , infile_noaa_de , /get_lun

;loop to pull out the months, sort by latitude
for j = 0, deseasonal_lines - 1 do begin
	if (j NE deseasonal_lines - 1) then b = j + 1 else b = j
		if ( deseasonal_arr_noaa[ 1 , j ] EQ deseasonal_arr_noaa[ 1 , b ] ) then begin
			a = where( deseasonal_arr_noaa[ 1 , * ] EQ deseasonal_arr_noaa[ 1 , j ] )
			temp_de_var = deseasonal_arr_noaa[ 0 , a ]
			avg_de_var = stddev( temp_de_var ) / sqrt( n_elements( temp_de_var ) )
			printf, lun1, deseasonal_arr_noaa[ 1 , j ], avg_de_var
			j = a[n_elements(a) - 1]					
			a = 0
		endif else if ( deseasonal_arr_noaa[ 1 , j ] NE deseasonal_arr_noaa[ 1 , b ] ) then j = j + 1
endfor

free_lun, lun1

;refresh the temp files to write the data into arrays for plotting
openr , lun1 , infile_noaa_de , /get_lun

avg_de_size = file_lines( infile_noaa_de )

de_obs_avg_noaa = fltarr( 2 , avg_de_size)

for t = 0, avg_de_size - 1 do begin
	readf, lun1, lat_var_de, stderr
	de_obs_avg_noaa[ 0 , t ] = lat_var_de
	de_obs_avg_noaa[ 1 , t ] = stderr
endfor

free_lun, lun1

;Finish SUB PART THREE-ONE==================================================================================================
;Start SUB PART THREE-TWO-------Dealing with the Simpson data=====================================================================
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
infile_ogi_de = "/home/jakechung/IDL/myIDL/temp_folder/temp_avg_SS_deseason.dat"

;open to write in the temp files
openw , lun1 , infile_ogi_de , /get_lun


;loop to pull out the months, sort by latitude
for j = 0, deseasonal_lines - 1 do begin
	if (j NE deseasonal_lines - 1) then b = j + 1 else b = j
		if ( deseasonal_arr_ss[ 1 , j ] EQ deseasonal_arr_ss[ 1 , b ] ) then begin
			a = where( deseasonal_arr_ss[ 1 , * ] EQ deseasonal_arr_ss[ 1 , j ] )
			temp_de_var = deseasonal_arr_ss[ 0 , a ]
			avg_de_var = stddev( temp_de_var ) / sqrt( n_elements( temp_de_var ) )
			printf, lun1, deseasonal_arr_ss[ 1 , j ], avg_de_var
			j = a[n_elements(a) - 1]					
			a = 0
		endif else if ( deseasonal_arr_ss[ 1 , j ] NE deseasonal_arr_ss[ 1 , b ] ) then j = j + 1
endfor

free_lun, lun1
;refresh the temp files to write the data into arrays for plotting
openr , lun1 , infile_ogi_de , /get_lun

avg_de_size = file_lines( infile_ogi_de )

de_obs_avg_ss = fltarr( 2 , avg_de_size)

for t = 0, avg_de_size - 1 do begin
	readf, lun1, lat_var_de, stderr
	de_obs_avg_ss[ 0 , t ] = lat_var_de
	de_obs_avg_ss[ 1 , t ] = stderr
endfor

free_lun, lun1
;Finish SUB PART THREE-TWO==================================================================================================


;Start SUB PART THREE-ONE-------Dealing with the OGI data=====================================================================
;Create link to the data file
temp_deseason = "/home/jakechung/IDL/myIDL/temp_folder/temp_residual_OGI.dat"
openr, deseason, temp_deseason, /get_lun

deseasonal_lines = file_lines(temp_deseason)
deseasonal_arr_ogi = fltarr(5, deseasonal_lines)
for u = 0, deseasonal_lines - 1 do begin
	readf, deseason, deseasonal_arr0, lat0, lon0, year0, month0
	deseasonal_arr_ogi[ 0 , u ] = deseasonal_arr0
	deseasonal_arr_ogi[ 1 , u ] = lat0
	deseasonal_arr_ogi[ 2 , u ] = lon0
	deseasonal_arr_ogi[ 3 , u ] = year0
	deseasonal_arr_ogi[ 4 , u ] = month0
endfor

free_lun, deseason

;again, using the same method to seperate the months out from the deseasonal temp array from the deseason_all_noaa program
infile_avg_de = "/home/jakechung/IDL/myIDL/temp_folder/temp_avg_OGI_deseason.dat"

;open to write in the temp files
openw , lun1 , infile_avg_de , /get_lun

;loop to pull out the months, sort by latitude
for j = 0, deseasonal_lines - 1 do begin
	if (j NE deseasonal_lines - 1) then b = j + 1 else b = j
		if ( deseasonal_arr_ogi[ 1 , j ] EQ deseasonal_arr_ogi[ 1 , b ] ) then begin
			a = where( deseasonal_arr_ogi[ 1 , * ] EQ deseasonal_arr_ogi[ 1 , j ] )
			temp_de_var = deseasonal_arr_ogi[ 0 , a ]
			avg_de_var = stddev( temp_de_var ) / sqrt( n_elements( temp_de_var ) )
			printf, lun1, deseasonal_arr_ogi[ 1 , j ], avg_de_var
			j = a[n_elements(a) - 1]					
			a = 0
		endif else if ( deseasonal_arr_ogi[ 1 , j ] NE deseasonal_arr_ogi[ 1 , b ] ) then j = j + 1
endfor

free_lun, lun1

;refresh the temp files to write the data into arrays for plotting
openr , lun1 , infile_avg_de , /get_lun

avg_de_size = file_lines( infile_avg_de )

de_obs_avg_ogi = fltarr( 2 , avg_de_size)

for t = 0, avg_de_size - 1 do begin
	readf, lun1, lat_var_de, stderr
	de_obs_avg_ogi[ 0 , t ] = lat_var_de
	de_obs_avg_ogi[ 1 , t ] = stderr
endfor

free_lun, lun1
;Finish SUB PART THREE-TWO==================================================================================================
;Finish PART THREE===========================================================================================================
;============================================================================================================================
;============================================================================================================================
;PART FOUR==================================================================================================================

;First, need to compile the list of all the coordinates of the three data sets
;The three arrays that contain the coordinates of all the stations of the three data sets are:
;obs_[month]_arr_noaa
;obs_[month]_arr_ogi
;obs_[month]_arr_ss
;0 is avg, 1 is latitude and 2 is longitude

;array all_coordinates will store the coordinates of all stations of all three data sets.
all_coordinates = fltarr( 2, n_elements(obs_mar_arr_noaa[ 1 , * ]) + n_elements(obs_mar_arr_ogi[ 1 , * ]) + 	$
	n_elements(obs_mar_arr_ss[ 1 , * ]) )

for i = 0, n_elements(obs_mar_arr_noaa[ 1 , * ]) - 1 do begin
	all_coordinates[ 0 , i ] = obs_mar_arr_noaa[ 1 , i ]
	all_coordinates[ 1 , i ] = obs_mar_arr_noaa[ 2 , i ]
endfor

for i = n_elements(obs_mar_arr_noaa[ 1 , * ]), n_elements(obs_mar_arr_ogi[ 1 , * ]) +	$
		n_elements(obs_mar_arr_noaa[ 1 , * ]) - 1 do begin
	all_coordinates[ 0 , i ] = obs_mar_arr_ogi[ 1 , i - n_elements(obs_mar_arr_noaa[ 1 , * ]) - 1 ]
	all_coordinates[ 1 , i ] = obs_mar_arr_ogi[ 2 , i - n_elements(obs_mar_arr_noaa[ 1 , * ]) - 1 ]
endfor

for i = n_elements(obs_mar_arr_ogi[ 1 , * ]) + n_elements(obs_mar_arr_noaa[ 1 , * ]), 	$
		n_elements(all_coordinates[ 0 , * ]) - 1 do begin
	all_coordinates[ 0 , i ] = obs_mar_arr_ss[ 1 , i - 	$
		n_elements(obs_mar_arr_noaa[ 1 , * ]) - n_elements(obs_mar_arr_ogi[ 1 , * ]) - 1 ]
	all_coordinates[ 1 , i ] = obs_mar_arr_ss[ 2 , i - 	$
		n_elements(obs_mar_arr_noaa[ 1 , * ]) - n_elements(obs_mar_arr_ogi[ 1 , * ]) - 1 ]
endfor

;=====================================
;===Start extracting simulated data===
;=====================================
; Define the filename
;filename = "/home/excluded-from-backup/data/C2H6/trac_avg.spinup_1981_2015.198101010000"; default emissions
;filename = "/home/excluded-from-backup/data/C2H6/trac_avg.AydinSF_GFED4_MAVG_1981_2015.198101010000" ;Aydin et al. ems_sce
;filename = "/home/excluded-from-backup/data/C2H6/trac_avg.PSUSF_GFED4_MAVG_1981_2015.198101010000" ;PSU ems_sce
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

;locate indexes of all months of the simulated data 
;sim_jan = where( (sim_yymmdd.month EQ 1) AND (sim_yymmdd.year GT 1990) AND (sim_yymmdd.year LT 2017) )
;sim_feb = where( (sim_yymmdd.month EQ 2) AND (sim_yymmdd.year GT 1990) AND (sim_yymmdd.year LT 2017) )
;sim_mar = where( (sim_yymmdd.month EQ 3) AND (sim_yymmdd.year GT 1990) AND (sim_yymmdd.year LT 2017) )
;sim_apr = where( (sim_yymmdd.month EQ 4) AND (sim_yymmdd.year GT 1990) AND (sim_yymmdd.year LT 2017) )
;sim_may = where( (sim_yymmdd.month EQ 5) AND (sim_yymmdd.year GT 1990) AND (sim_yymmdd.year LT 2017) )
;sim_jun = where( (sim_yymmdd.month EQ 6) AND (sim_yymmdd.year GT 1990) AND (sim_yymmdd.year LT 2017) )
;sim_jul = where( (sim_yymmdd.month EQ 7) AND (sim_yymmdd.year GT 1990) AND (sim_yymmdd.year LT 2017) )
;sim_aug = where( (sim_yymmdd.month EQ 8) AND (sim_yymmdd.year GT 1990) AND (sim_yymmdd.year LT 2017) )
;sim_sep = where( (sim_yymmdd.month EQ 9) AND (sim_yymmdd.year GT 1990) AND (sim_yymmdd.year LT 2017) )
;sim_oct = where( (sim_yymmdd.month EQ 10) AND (sim_yymmdd.year GT 1990) AND (sim_yymmdd.year LT 2017) )
;sim_nov = where( (sim_yymmdd.month EQ 11) AND (sim_yymmdd.year GT 1990) AND (sim_yymmdd.year LT 2017) )
;sim_dec = where( (sim_yymmdd.month EQ 12) AND (sim_yymmdd.year GT 1990) AND (sim_yymmdd.year LT 2017) )

;extracting the simulated data for the corresponding months and years and locations
;sim_jan_var = fltarr(n_elements(all_coordinates[0 , *]))
;for i = 0, n_elements(all_coordinates[0 , *]) - 1 do begin
;	for t = 0, nt - 1 do begin
;		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
;		gridinfo= gridinfo, lat= [all_coordinates[0, i] - 0.5, all_coordinates[0, i] + 0.5], $
;		lon= [all_coordinates[1, i] - 0.5, all_coordinates[1, i] + 0.5], alrange = [1, 3], average=7)
;		globalavg[ t ]= data/2 * 1000
;	endfor
;	sim_jan_var[ i ] = mean(globalavg[ sim_jan ])
;endfor

;sim_feb_var = fltarr(n_elements(all_coordinates[0 , *]))
;for i = 0, n_elements(all_coordinates[0 , *]) - 1 do begin
;	for t = 0, nt - 1 do begin
;		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
;		gridinfo= gridinfo, lat= [all_coordinates[0, i] - 0.5, all_coordinates[0, i] + 0.5], $
;		lon= [all_coordinates[1, i] - 0.5, all_coordinates[1, i] + 0.5], alrange = [1, 3], average=7)
;		globalavg[ t ]= data/2 * 1000
;	endfor
;	sim_feb_var[ i ] = mean(globalavg[ sim_feb ])
;endfor

;sim_mar_var = fltarr(n_elements(all_coordinates[0 , *]))
;for i = 0, n_elements(all_coordinates[0 , *]) - 1 do begin
;	for t = 0, nt - 1 do begin
;		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
;		gridinfo= gridinfo, lat= [all_coordinates[0, i] - 0.5, all_coordinates[0, i] + 0.5], $
;		lon= [all_coordinates[1, i] - 0.5, all_coordinates[1, i] + 0.5], alrange = [1, 3], average=7)
;		globalavg[ t ]= data/2 * 1000
;	endfor
;	sim_mar_var[ i ] = mean(globalavg[ sim_mar ])
;endfor

;sim_apr_var = fltarr(n_elements(all_coordinates[0 , *]))
;for i = 0, n_elements(all_coordinates[0 , *]) - 1 do begin
;	for t = 0, nt - 1 do begin
;		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
;		gridinfo= gridinfo, lat= [all_coordinates[0, i] - 0.5, all_coordinates[0, i] + 0.5], $
;		lon= [all_coordinates[1, i] - 0.5, all_coordinates[1, i] + 0.5], alrange = [1, 3], average=7)
;		globalavg[ t ]= data/2 * 1000
;	endfor
;	sim_apr_var[ i ] = mean(globalavg[ sim_apr ])
;endfor

;sim_may_var = fltarr(n_elements(all_coordinates[0 , *]))
;for i = 0, n_elements(all_coordinates[0 , *]) - 1 do begin
;	for t = 0, nt - 1 do begin
;		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
;		gridinfo= gridinfo, lat= [all_coordinates[0, i] - 0.5, all_coordinates[0, i] + 0.5], $
;		lon= [all_coordinates[1, i] - 0.5, all_coordinates[1, i] + 0.5], alrange = [1, 3], average=7)
;		globalavg[ t ]= data/2 * 1000
;	endfor
;	sim_may_var[ i ] = mean(globalavg[ sim_may ])
;endfor

;sim_jun_var = fltarr(n_elements(all_coordinates[0 , *]))
;for i = 0, n_elements(all_coordinates[0 , *]) - 1 do begin
;	for t = 0, nt - 1 do begin
;		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
;		gridinfo= gridinfo, lat= [all_coordinates[0, i] - 0.5, all_coordinates[0, i] + 0.5], $
;		lon= [all_coordinates[1, i] - 0.5, all_coordinates[1, i] + 0.5], alrange = [1, 3], average=7)
;		globalavg[ t ]= data/2 * 1000
;	endfor
;	sim_jun_var[ i ] = mean(globalavg[ sim_jun ])
;endfor

;sim_jul_var = fltarr(n_elements(all_coordinates[0 , *]))
;for i = 0, n_elements(all_coordinates[0 , *]) - 1 do begin
;	for t = 0, nt - 1 do begin
;		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
;		gridinfo= gridinfo, lat= [all_coordinates[0, i] - 0.5, all_coordinates[0, i] + 0.5], $
;		lon= [all_coordinates[1, i] - 0.5, all_coordinates[1, i] + 0.5], alrange = [1, 3], average=7)
;		globalavg[ t ]= data/2 * 1000
;	endfor
;	sim_jul_var[ i ] = mean(globalavg[ sim_jul ])
;endfor

;sim_aug_var = fltarr(n_elements(all_coordinates[0 , *]))
;for i = 0, n_elements(all_coordinates[0 , *]) - 1 do begin
;	for t = 0, nt - 1 do begin
;		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
;		gridinfo= gridinfo, lat= [all_coordinates[0, i] - 0.5, all_coordinates[0, i] + 0.5], $
;		lon= [all_coordinates[1, i] - 0.5, all_coordinates[1, i] + 0.5], alrange = [1, 3], average=7)
;		globalavg[ t ]= data/2 * 1000
;	endfor
;	sim_aug_var[ i ] = mean(globalavg[ sim_aug ])
;endfor

;sim_sep_var = fltarr(n_elements(all_coordinates[0 , *]))
;for i = 0, n_elements(all_coordinates[0 , *]) - 1 do begin
;	for t = 0, nt - 1 do begin
;		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
;		gridinfo= gridinfo, lat= [all_coordinates[0, i] - 0.5, all_coordinates[0, i] + 0.5], $
;		lon= [all_coordinates[1, i] - 0.5, all_coordinates[1, i] + 0.5], alrange = [1, 3], average=7)
;		globalavg[ t ]= data/2 * 1000
;	endfor
;	sim_sep_var[ i ] = mean(globalavg[ sim_sep ])
;endfor

;sim_oct_var = fltarr(n_elements(all_coordinates[0 , *]))
;for i = 0, n_elements(all_coordinates[0 , *]) - 1 do begin
;	for t = 0, nt - 1 do begin
;		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
;		gridinfo= gridinfo, lat= [all_coordinates[0, i] - 0.5, all_coordinates[0, i] + 0.5], $
;		lon= [all_coordinates[1, i] - 0.5, all_coordinates[1, i] + 0.5], alrange = [1, 3], average=7)
;		globalavg[ t ]= data/2 * 1000
;	endfor
;	sim_oct_var[ i ] = mean(globalavg[ sim_oct ])
;endfor

;sim_nov_var = fltarr(n_elements(all_coordinates[0 , *]))
;for i = 0, n_elements(all_coordinates[0 , *]) - 1 do begin
;	for t = 0, nt - 1 do begin
;		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
;		gridinfo= gridinfo, lat= [all_coordinates[0, i] - 0.5, all_coordinates[0, i] + 0.5], $
;		lon= [all_coordinates[1, i] - 0.5, all_coordinates[1, i] + 0.5], alrange = [1, 3], average=7)
;		globalavg[ t ]= data/2 * 1000
;	endfor
;	sim_nov_var[ i ] = mean(globalavg[ sim_nov ])
;endfor

;sim_dec_var = fltarr(n_elements(all_coordinates[0 , *]))
;for i = 0, n_elements(all_coordinates[0 , *]) - 1 do begin
;	for t = 0, nt - 1 do begin
;		data= CTM_EXTRACT( *(datainfo[t].data), modelinfo=modelinfo, $
;		gridinfo= gridinfo, lat= [all_coordinates[0, i] - 0.5, all_coordinates[0, i] + 0.5], $
;		lon= [all_coordinates[1, i] - 0.5, all_coordinates[1, i] + 0.5], alrange = [1, 3], average=7)
;		globalavg[ t ]= data/2 * 1000
;	endfor
;	sim_dec_var[ i ] = mean(globalavg[ sim_dec ])
;endfor

;Finish PART FOUR=================================================================================================
;Start PART FIVE- All plotting routines==================================================================================
;Plot concentration vs latitude for the NOAA, Simpson, and OGI data
;Note: the Simpson et al. data is only available for March, June, September and December

;Setting up windows size, there will be 12 graphs total for 12 months
!P.Multi = 0
;cgDisplay, 1270, 1100
;!P.Multi = [0, 3, 4, 0, 0]

measure_errors = de_obs_avg_noaa[ 1 , * ]
p_fit_noaa = poly_fit( obs_stt_avg_noaa[ 0 , * ], obs_stt_avg_noaa[ 2 , * ], 6, measure_errors = measure_errors, sigma = sigma_noaa, $
	yfit = yfit_noaa, /double )

measure_errors = de_obs_avg_ss[ 1 , * ]
p_fit_ss = poly_fit( obs_stt_avg_ss[ 0 , * ], obs_stt_avg_ss[ 2 , * ], 6, measure_errors = measure_errors, sigma = sigma_ss, $
	yfit = yfit_ss, /double )	

measure_errors = de_obs_avg_ogi[ 1 , * ]
p_fit_ogi = poly_fit( obs_stt_avg_ogi[ 0 , * ], obs_stt_avg_ogi[ 2 , * ], 6, measure_errors = measure_errors, sigma = sigma_ogi, $
	yfit = yfit_ogi, /double )

cof_noaa0 = p_fit_noaa[ 0 ]
cof_noaa1 = p_fit_noaa[ 1 ]
cof_noaa2 = p_fit_noaa[ 2 ]
cof_noaa3 = p_fit_noaa[ 3 ]
cof_noaa4 = p_fit_noaa[ 4 ]
cof_noaa5 = p_fit_noaa[ 5 ]
cof_noaa6 = p_fit_noaa[ 6 ]

cof_ogi0 = p_fit_ogi[ 0 ]
cof_ogi1 = p_fit_ogi[ 1 ]
cof_ogi2 = p_fit_ogi[ 2 ]
cof_ogi3 = p_fit_ogi[ 3 ]
cof_ogi4 = p_fit_ogi[ 4 ]
cof_ogi5 = p_fit_ogi[ 5 ]
cof_ogi6 = p_fit_ogi[ 6 ]

cof_ss0 = p_fit_ss[ 0 ]
cof_ss1 = p_fit_ss[ 1 ]
cof_ss2 = p_fit_ss[ 2 ]
cof_ss3 = p_fit_ss[ 3 ]
cof_ss4 = p_fit_ss[ 4 ]
cof_ss5 = p_fit_ss[ 5 ]
cof_ss6 = p_fit_ss[ 6 ]
print, p_fit_ss
print, p_fit_noaa
cgplot, obs_stt_avg_ss[ 0 , sort( obs_stt_avg_ss[ 0 , * ] ) ], yfit_ss[ sort( obs_stt_avg_ss[ 0 , * ] ) ], xrange = [ -100 , 100 ], $
	color = 'red'
cgplot, obs_stt_avg_noaa[ 0 , sort( obs_stt_avg_noaa[ 0 , * ] ) ], yfit_noaa[ sort( obs_stt_avg_noaa[ 0 , * ] ) ], /overplot, $
	color = 'blue'
cgplot, obs_stt_avg_ogi[ 0 , sort( obs_stt_avg_ogi[ 0 , * ] ) ], yfit_ogi[ sort( obs_stt_avg_ogi[ 0 , * ] ) ], /overplot, $
	color = 'forest green'
;bin01 = -90 to -75
;bin02 = -75 to -60
;bin03 = -60 to -45
;bin04 = -45 to -30
;bin05 = -30 to -15
;bin06 = -15 to 0
;bin07 = 0 to 15
;bin08 = 15 to 30
;bin09 = 30 to 45
;bin10 = 45 to 60
;bin11 = 60 to 75
;bin12 = 75 tp 90

sou_var_noaa = fltarr( 6 )
;sou_var_noaa[ 0 ] = def_integrate_poly6( cof_noaa0, cof_noaa1, cof_noaa2, cof_noaa3, cof_noaa4, cof_noaa5, $
;	cof_noaa6, -75, -90 )
;sou_var_noaa[ 1 ] = def_integrate_poly6( cof_noaa0, cof_noaa1, cof_noaa2, cof_noaa3, cof_noaa4, cof_noaa5, $
;	cof_noaa6, -60, -75 )
;sou_var_noaa[ 2 ] = def_integrate_poly6( cof_noaa0, cof_noaa1, cof_noaa2, cof_noaa3, cof_noaa4, cof_noaa5, $
;	cof_noaa6, -45, -60 )
sou_var_noaa[ 3 ] = def_integrate_poly6( cof_noaa0, cof_noaa1, cof_noaa2, cof_noaa3, cof_noaa4, cof_noaa5, $
	cof_noaa6, -30, -45 )
sou_var_noaa[ 4 ] = def_integrate_poly6( cof_noaa0, cof_noaa1, cof_noaa2, cof_noaa3, cof_noaa4, cof_noaa5, $
	cof_noaa6, -15, -30 )
sou_var_noaa[ 5 ] = def_integrate_poly6( cof_noaa0, cof_noaa1, cof_noaa2, cof_noaa3, cof_noaa4, cof_noaa5, $
	cof_noaa6, 0, -15 )
	
nor_var_noaa = fltarr( 6 )
nor_var_noaa[ 0 ] = def_integrate_poly6( cof_noaa0, cof_noaa1, cof_noaa2, cof_noaa3, cof_noaa4, cof_noaa5, $
	cof_noaa6, 15, 0 )
nor_var_noaa[ 1 ] = def_integrate_poly6( cof_noaa0, cof_noaa1, cof_noaa2, cof_noaa3, cof_noaa4, cof_noaa5, $
	cof_noaa6, 30, 15 )
nor_var_noaa[ 2 ] = def_integrate_poly6( cof_noaa0, cof_noaa1, cof_noaa2, cof_noaa3, cof_noaa4, cof_noaa5, $
	cof_noaa6, 45, 30 )
nor_var_noaa[ 3 ] = def_integrate_poly6( cof_noaa0, cof_noaa1, cof_noaa2, cof_noaa3, cof_noaa4, cof_noaa5, $
	cof_noaa6, 60, 45 )
nor_var_noaa[ 4 ] = def_integrate_poly6( cof_noaa0, cof_noaa1, cof_noaa2, cof_noaa3, cof_noaa4, cof_noaa5, $
	cof_noaa6, 75, 60 )
;nor_var_noaa[ 5 ] = def_integrate_poly6( cof_noaa0, cof_noaa1, cof_noaa2, cof_noaa3, cof_noaa4, cof_noaa5, $
;	cof_noaa6, 90, 75 )
	
avg_sou_noaa = ( sou_var_noaa[ 0 ] * sin( 82.5 * !pi / 180 ) + sou_var_noaa[ 1 ] * sin( 67.5 * !pi / 180 ) + $
	sou_var_noaa[ 2 ] * sin( 52.5 * !pi / 180 ) + sou_var_noaa[ 3 ] * sin( 37.5 * !pi / 180 ) + sou_var_noaa[ 4 ] * $
	sin( 22.5 * !pi / 180 ) + sou_var_noaa[ 5 ] * sin( 7.5 * !pi / 180 ) ) / 3
	
avg_nor_noaa = ( nor_var_noaa[ 5 ] * sin( 82.5 * !pi / 180 ) + nor_var_noaa[ 4 ] * sin( 67.5 * !pi / 180 ) + $
	nor_var_noaa[ 3 ] * sin( 52.5 * !pi / 180 ) + nor_var_noaa[ 2 ] * sin( 37.5 * !pi / 180 ) + nor_var_noaa[ 1 ] * $
	sin( 22.5 * !pi / 180 ) + nor_var_noaa[ 0 ] * sin( 7.5 * !pi / 180 ) ) / 5
	
hem_noaa_ratio = avg_nor_noaa / avg_sou_noaa

sou_var_ss = fltarr( 6 )
;sou_var_ss[ 0 ] = def_integrate_poly6( cof_ss0, cof_ss1, cof_ss2, cof_ss3, cof_ss4, cof_ss5, $
;	cof_ss6, -75, -90 )
;sou_var_ss[ 1 ] = def_integrate_poly6( cof_ss0, cof_ss1, cof_ss2, cof_ss3, cof_ss4, cof_ss5, $
;	cof_ss6, -60, -75 )
;sou_var_ss[ 2 ] = def_integrate_poly6( cof_ss0, cof_ss1, cof_ss2, cof_ss3, cof_ss4, cof_ss5, $
;	cof_ss6, -45, -60 )
sou_var_ss[ 3 ] = def_integrate_poly6( cof_ss0, cof_ss1, cof_ss2, cof_ss3, cof_ss4, cof_ss5, $
	cof_ss6, -30, -45 )
sou_var_ss[ 4 ] = def_integrate_poly6( cof_ss0, cof_ss1, cof_ss2, cof_ss3, cof_ss4, cof_ss5, $
	cof_ss6, -15, -30 )
sou_var_ss[ 5 ] = def_integrate_poly6( cof_ss0, cof_ss1, cof_ss2, cof_ss3, cof_ss4, cof_ss5, $
	cof_ss6, 0, -15 )
	
nor_var_ss = fltarr( 6 )
nor_var_ss[ 0 ] = def_integrate_poly6( cof_ss0, cof_ss1, cof_ss2, cof_ss3, cof_ss4, cof_ss5, $
	cof_ss6, 15, 0 )
nor_var_ss[ 1 ] = def_integrate_poly6( cof_ss0, cof_ss1, cof_ss2, cof_ss3, cof_ss4, cof_ss5, $
	cof_ss6, 30, 15 )
nor_var_ss[ 2 ] = def_integrate_poly6( cof_ss0, cof_ss1, cof_ss2, cof_ss3, cof_ss4, cof_ss5, $
	cof_ss6, 45, 30 )
nor_var_ss[ 3 ] = def_integrate_poly6( cof_ss0, cof_ss1, cof_ss2, cof_ss3, cof_ss4, cof_ss5, $
	cof_ss6, 60, 45 )
nor_var_ss[ 4 ] = def_integrate_poly6( cof_ss0, cof_ss1, cof_ss2, cof_ss3, cof_ss4, cof_ss5, $
	cof_ss6, 75, 60 )
;nor_var_ss[ 5 ] = def_integrate_poly6( cof_ss0, cof_ss1, cof_ss2, cof_ss3, cof_ss4, cof_ss5, $
;	cof_ss6, 90, 75 )
	
avg_sou_ss = ( sou_var_ss[ 0 ] * sin( 82.5 * !pi / 180 ) + sou_var_ss[ 1 ] * sin( 67.5 * !pi / 180 ) + $
	sou_var_ss[ 2 ] * sin( 52.5 * !pi / 180 ) + sou_var_ss[ 3 ] * sin( 37.5 * !pi / 180 ) + sou_var_ss[ 4 ] * $
	sin( 22.5 * !pi / 180 ) + sou_var_ss[ 5 ] * sin( 7.5 * !pi / 180 ) ) / 3
	
avg_nor_ss = ( nor_var_ss[ 5 ] * sin( 82.5 * !pi / 180 ) + nor_var_ss[ 4 ] * sin( 67.5 * !pi / 180 ) + $
	nor_var_ss[ 3 ] * sin( 52.5 * !pi / 180 ) + nor_var_ss[ 2 ] * sin( 37.5 * !pi / 180 ) + nor_var_ss[ 1 ] * $
	sin( 22.5 * !pi / 180 ) + nor_var_ss[ 0 ] * sin( 7.5 * !pi / 180 ) ) / 5
	
hem_ss_ratio = avg_nor_ss / avg_sou_ss

sou_var_ogi = fltarr( 6 )
;sou_var_ogi[ 0 ] = def_integrate_poly6( cof_ogi0, cof_ogi1, cof_ogi2, cof_ogi3, cof_ogi4, cof_ogi5, $
;	cof_ogi6, -75, -90 )
;sou_var_ogi[ 1 ] = def_integrate_poly6( cof_ogi0, cof_ogi1, cof_ogi2, cof_ogi3, cof_ogi4, cof_ogi5, $
;	cof_ogi6, -60, -75 )
;sou_var_ogi[ 2 ] = def_integrate_poly6( cof_ogi0, cof_ogi1, cof_ogi2, cof_ogi3, cof_ogi4, cof_ogi5, $
;	cof_ogi6, -45, -60 )
sou_var_ogi[ 3 ] = def_integrate_poly6( cof_ogi0, cof_ogi1, cof_ogi2, cof_ogi3, cof_ogi4, cof_ogi5, $
	cof_ogi6, -30, -45 )
sou_var_ogi[ 4 ] = def_integrate_poly6( cof_ogi0, cof_ogi1, cof_ogi2, cof_ogi3, cof_ogi4, cof_ogi5, $
	cof_ogi6, -15, -30 )
sou_var_ogi[ 5 ] = def_integrate_poly6( cof_ogi0, cof_ogi1, cof_ogi2, cof_ogi3, cof_ogi4, cof_ogi5, $
	cof_ogi6, 0, -15 )
	
nor_var_ogi = fltarr( 6 )
nor_var_ogi[ 0 ] = def_integrate_poly6( cof_ogi0, cof_ogi1, cof_ogi2, cof_ogi3, cof_ogi4, cof_ogi5, $
	cof_ogi6, 15, 0 )
nor_var_ogi[ 1 ] = def_integrate_poly6( cof_ogi0, cof_ogi1, cof_ogi2, cof_ogi3, cof_ogi4, cof_ogi5, $
	cof_ogi6, 30, 15 )
nor_var_ogi[ 2 ] = def_integrate_poly6( cof_ogi0, cof_ogi1, cof_ogi2, cof_ogi3, cof_ogi4, cof_ogi5, $
	cof_ogi6, 45, 30 )
nor_var_ogi[ 3 ] = def_integrate_poly6( cof_ogi0, cof_ogi1, cof_ogi2, cof_ogi3, cof_ogi4, cof_ogi5, $
	cof_ogi6, 60, 45 )
nor_var_ogi[ 4 ] = def_integrate_poly6( cof_ogi0, cof_ogi1, cof_ogi2, cof_ogi3, cof_ogi4, cof_ogi5, $
	cof_ogi6, 75, 60 )
;nor_var_ogi[ 5 ] = def_integrate_poly6( cof_ogi0, cof_ogi1, cof_ogi2, cof_ogi3, cof_ogi4, cof_ogi5, $
;	cof_ogi6, 90, 75 )
	print, sou_var_ogi

avg_sou_ogi = ( sou_var_ogi[ 0 ] * sin( 82.5 * !pi / 180 ) + sou_var_ogi[ 1 ] * sin( 67.5 * !pi / 180 ) + $
	sou_var_ogi[ 2 ] * sin( 52.5 * !pi / 180 ) + sou_var_ogi[ 3 ] * sin( 37.5 * !pi / 180 ) + sou_var_ogi[ 4 ] * $
	sin( 22.5 * !pi / 180 ) + sou_var_ogi[ 5 ] * sin( 7.5 * !pi / 180 ) ) / 3
print, avg_sou_ogi
avg_nor_ogi = ( nor_var_ogi[ 5 ] * sin( 82.5 * !pi / 180 ) + nor_var_ogi[ 4 ] * sin( 67.5 * !pi / 180 ) + $
	nor_var_ogi[ 3 ] * sin( 52.5 * !pi / 180 ) + nor_var_ogi[ 2 ] * sin( 37.5 * !pi / 180 ) + nor_var_ogi[ 1 ] * $
	sin( 22.5 * !pi / 180 ) + nor_var_ogi[ 0 ] * sin( 7.5 * !pi / 180 ) ) / 5
	
hem_ogi_ratio = avg_nor_ogi / avg_sou_ogi

print, hem_noaa_ratio, hem_ss_ratio, hem_ogi_ratio

END




































