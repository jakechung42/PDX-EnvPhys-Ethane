PRO HEM_RATIO_ALLSTATION

;This program read in the all three data sets, calculate their hemispheric ratio and compare them with one another
;Further revision need to have an algorithm of some kind to remove the stations that fluctuate too much 
;New implemented method is to calculate the weighted means of the area under the curve of the poly_fit function
;from those averages, the hemispheric can be derived

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


;Start calculating the ratio of the three data sets.
;Start with the NOAA data set first
nor_noaa_idx = where( lat1_noaa GT 0 AND lat1_noaa NE 54.95 ) ; exclude the Alberto station because it bias the mean
nor_noaa_var = mean( avg_noaa[ nor_noaa_idx ] ) 

sou_noaa_idx = where( lat1_noaa LT 0 )
sou_noaa_var = mean( avg_noaa[ sou_noaa_idx ] )

hem_ratio_noaa = nor_noaa_var / sou_noaa_var
print, 'Hemispheric ratio of the NOAA data is', hem_ratio_noaa

;OGI data set
nor_ogi_idx = where( lat1_ogi GT 0 )
nor_ogi_var = mean( avg_ogi[ nor_ogi_idx ] )

sou_ogi_idx = where( lat1_ogi LT 0 )
sou_ogi_var = mean( avg_ogi[ sou_ogi_idx ] )

hem_ratio_ogi = nor_ogi_var / sou_ogi_var
print, 'Hemispheric ratio of the OGI data is', hem_ratio_ogi

;Now Simpson et al. data will be tricky because not all ss data is used in analysis because of the lack of data at some stations
;To keep it consistent with other programs, this program will pick out only the specific stations that are used in analysis
lat_var_GT20_ss = [ 71.3000, 64.5000, 64.4900, 60.7500, 57.8400, 57.8300, 57.8000, 46.9800, 45.9400, 44.6800, 42.8300, 37.9900, $
	30.4000, 29.9400, 23.4400, 22.8700, 22.1500, 20.2600, 15.2000, 13.3600, 9.52000, 7.42000, 7.05000, 6.92000, 5.32000, 1.42000, $
	-1.33000, -8.52000, -14.2300, -21.2300, -29.0200, -34.4200, -34.9000, -35.0300, -36.8200, -42.7200, -43.4200, -43.8800, -46.4400 ]

ss_temp_data = fltarr( 2 , n_elements( lat_var_GT20_ss ) )
for i = 0, n_elements( lat_var_GT20_ss ) - 1 do begin 
	a = where( lat1_ss EQ lat_var_GT20_ss[ i ] )
	ss_temp_data[ 0 , i ] = lat1_ss[ a[ 0 ] ] 
	ss_temp_data[ 1 , i ] = mean( avg_ss[ a ] )
endfor 
	
nor_ss_var = mean( ss_temp_data[ 1 , where( ss_temp_data[ 0 , * ] GT 0 ) ] )
sou_ss_var = mean( ss_temp_data[ 1 , where( ss_temp_data[ 0 , * ] LT 0 ) ] )

hem_ratio_ss = nor_ss_var / sou_ss_var 
print, 'Hemispheric ratio of the Simpson et al. data is', hem_ratio_ss







end

