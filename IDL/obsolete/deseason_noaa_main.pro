PRO DESEASON_NOAA_MAIN

;Author: Jake Chung
;REU 2016/ AGU 2016
;Description: this program deseasonalize the NOAA data and then store the residual data at
;"/home/jakechung/IDL/myIDL/temp_folder/temp_residual_noaa.dat"

;Notes: plotting procedures can be added to plot the residual if needed be.
;The residual follows a Gaussian distribution (refer to figure 2, AGU poster) . 

COMPILE_OPT IDL2				;set compile options

;==========================
;<<<<<<Extract NOAA data>>>>>>
;==========================
;get the files
infile1= "/home/jakechung/IDL/myIDL/NOAA_data.dat"
infile2= "/home/jakechung/IDL/myIDL/NOAA_locations.dat"

;get number of lines in the data file
n1 = file_lines(infile1)
n2 = file_lines(infile2)

print, 'The NOAA data has ', n1, '   data lines with', n2, '  locations'

;create empty arrays to store the data that was read in lon(longitudes), lat(latitudes), 
;avg(concentration), date and locations

lat1= fltarr(n1)
lon1= fltarr(n1)
avg= fltarr(n1)
location= strarr(n2)
year= fltarr(n1)
month= fltarr(n1)
day= fltarr(n1)
hour= fltarr(n1)
alt= fltarr(n1)
time= fltarr(n1)

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
for i = 0, n1 - 1 do begin
	;read one line from the file
	readf, iunit1, sample_year, sample_month, sample_day, sample_hour, analysis_value,	$
				sample_latitude, sample_longitude, sample_altitude
	;store the values
	year[i] = sample_year
	month[i] = sample_month
	day[i] = sample_day
	hour[i] = sample_hour
	avg[i] = analysis_value
	lat1[i] = sample_latitude
	lon1[i] = sample_longitude
	alt[i] = sample_altitude
endfor

;read float data of the location file
for i = 0, n2-1 do begin
	readf, iunit2, sample_site_code
	location[i] = sample_site_code
endfor

;close input file
free_lun, iunit1, iunit2


;combine the year, month, day, hour into a decimal year values
;in this script this line is not necesary but keep it here anyway
for i = 0, n1 - 1 do begin	
	time[i] = year[i] + hour[i]/8544 + day[i]/365 + month[i]/12
endfor

;Remove the -999 values from the C2H6 arrays
;Note, the RemoveRows function is from the Coyote Library.
;The RemoveRows function removes rows but the arrays are in collums,
;so need to use the rotate function so that the RemoveRows function
;would work

neg999 = where(avg LT 0, count)
location = RemoveRows(rotate(location, 1), neg999)
time = RemoveRows(rotate(time, 1), neg999)
lat1 = RemoveRows(rotate(lat1, 1), neg999)
lon1 = RemoveRows(rotate(lon1, 1), neg999)
avg = RemoveRows(rotate(avg, 1), neg999)
year = RemoveRows(rotate(year, 1), neg999)
month = RemoveRows(rotate(month, 1), neg999)
day = RemoveRows(rotate(day, 1), neg999)

;check the new size of the two arrays, n1 and n2 should be the same.
n2 = n_elements(location)
n1 = n_elements(avg)
;print, n1, n2

;read in the stations ID file
;The station ID file contains the ID of all the stations in the NOAA data.
stations_ids = "/home/jakechung/IDL/myIDL/NOAA_stations_IDs.dat"
site_ids = ' '
n3 = file_lines(stations_ids)
station = strarr(n3)
openr, iunit3, stations_ids, /get_lun
for i = 0, n3 - 1 do begin
	readf, iunit3, site_ids
	station[ i ] = site_ids
endfor
free_lun, iunit3

;making the temp file to store the residual values
temp_residual = "/home/jakechung/IDL/myIDL/temp_folder/temp_residual_noaa.dat"
openw, iunit4, temp_residual, /get_lun

;------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
;start a giant loop right here to loop in, read all the sites IDs and deseasonalize the entire data set.
for r = 0, n3 - 1 do begin  ;<<<<<<<< loop starts here
z = where(strmatch(location[*], station[r], /FOLD_CASE) EQ 1) ;<<<calling out each station
print, 'Station ', station[r],' Index on the data set is ', z[0], ' to ', z[n_elements(z)-1]
avg_ice = avg[z]
time_ice = time[z]
lat1_ice = lat1[z]
lon1_ice = lon1[z]
month_ice = month[z]
day_ice = day[z]
year_ice = year[z]

;these arrays store the months of the measured data set. 
jan = where(month_ice[*] EQ 1)
feb = where(month_ice[*] EQ 2)
mar = where(month_ice[*] EQ 3)
apr = where(month_ice[*] EQ 4)
may = where(month_ice[*] EQ 5)
jun = where(month_ice[*] EQ 6)
jul = where(month_ice[*] EQ 7)
aug = where(month_ice[*] EQ 8)
sep = where(month_ice[*] EQ 9)
oct = where(month_ice[*] EQ 10)
nov = where(month_ice[*] EQ 11)
dec = where(month_ice[*] EQ 12)

jan_avg_var = mean( avg_ice [ jan ] ) 
feb_avg_var = mean( avg_ice [ feb ] ) 
mar_avg_var = mean( avg_ice [ mar ] ) 
apr_avg_var = mean( avg_ice [ apr ] ) 
may_avg_var = mean( avg_ice [ may ] ) 
jun_avg_var = mean( avg_ice [ jun ] ) 
jul_avg_var = mean( avg_ice [ jul ] ) 
aug_avg_var = mean( avg_ice [ aug ] ) 
sep_avg_var = mean( avg_ice [ sep ] ) 
oct_avg_var = mean( avg_ice [ oct ] ) 
nov_avg_var = mean( avg_ice [ nov ] ) 
dec_avg_var = mean( avg_ice [ dec ] ) 

residue_arr = fltarr(n_elements(avg_ice))
for j = 0, n_elements(avg_ice) - 1 do begin
	if (j NE n_elements(avg_ice) - 1) then b = j + 1 else b = j
	if ( year_ice [ j ] EQ year_ice [ b ] ) then begin
		a = where( year_ice EQ year_ice [ j ] )
		temp_avg_var = avg_ice [ a ]
		temp_month = month_ice [ a ] 
		for i = 0, n_elements(a) - 1 do begin
			if (temp_month[i] EQ 1) then $
				residue_arr[ i ] = jan_avg_var - temp_avg_var[ i ] else $
			if (temp_month[i] EQ 2) then $
				residue_arr[ i ] = feb_avg_var - temp_avg_var[ i ] else $
			if (temp_month[i] EQ 3) then $
				residue_arr[ i ] = mar_avg_var - temp_avg_var[ i ] else $
			if (temp_month[i] EQ 4) then $
				residue_arr[ i ] = apr_avg_var - temp_avg_var[ i ] else $
			if (temp_month[i] EQ 5) then $
				residue_arr[ i ] = may_avg_var - temp_avg_var[ i ] else $
			if (temp_month[i] EQ 6) then $
				residue_arr[ i ] = jun_avg_var - temp_avg_var[ i ] else $
			if (temp_month[i] EQ 7) then $
				residue_arr[ i ] = jul_avg_var - temp_avg_var[ i ] else $
			if (temp_month[i] EQ 8) then $
				residue_arr[ i ] = aug_avg_var - temp_avg_var[ i ] else $
			if (temp_month[i] EQ 9) then $
				residue_arr[ i ] = sep_avg_var - temp_avg_var[ i ] else $
			if (temp_month[i] EQ 10) then $
				residue_arr[ i ] = oct_avg_var - temp_avg_var[ i ] else $
			if (temp_month[i] EQ 11) then $
				residue_arr[ i ] = nov_avg_var - temp_avg_var[ i ] else $
			if (temp_month[i] EQ 12) then $
				residue_arr[ i ] = dec_avg_var - temp_avg_var[ i ] else $
				residue_arr[i] = 0
			if ( residue_arr[i] NE 0 ) then begin
				printf, iunit4, residue_arr[i], lat1_ice[5], lon1_ice[5], year_ice[ j ], temp_month[i]
			endif
		endfor
		j = a[n_elements(a) - 1]					
		a = 0
	endif 
endfor

endfor ;<<<<<<<< The big loop to loop in all the stations IDs ends here

free_lun, iunit4

END
