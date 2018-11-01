PRO TEMPORAL_IHR_FULL

;	08/29/17 - Jake Chung
;	This program is an overhaul of the ihr_temporal_all_plot1
;	The new data structure is added along with data filtering
;and Aydin firn air data.
;	Would also try to use a better algorithm to deseasonalize
;the data

compile_opt idl2




;Read NOAA data
print, 'I-Read data'
print,  ' '
print, 'a-Read NOAA data'
infile1= "/home/excluded-from-backup/ethane/data/raw_data/NOAA/NOAA_data.dat"
infile2= "/home/excluded-from-backup/ethane/data/raw_data/NOAA/NOAA_locations.dat"

;get number of lines in the data file
n1_noaa = file_lines(infile1)
n2_noaa = file_lines(infile2)

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
;-
;-
;-
;-
;-
;Read UCI data
;read in the modded UCI data where the lon and lat are
;changed from uci_mod_lat_lon
print, 'b-Read UCI data (Notes: read in modded UCI data file from uci_mod_lat_lon.pro)'
infile1= "/home/excluded-from-backup/ethane/IDL/temp_file/modded_uci_data_v2.dat"
infile2= "/home/excluded-from-backup/ethane/data/raw_data/UCI/uci_stations_v2.dat"

n1_ss = file_lines(infile1)
n2_ss = file_lines(infile2)

;create an empty array to store the data that was read in from the data files
;since the time data in the UCI data set is in a different format,
;so the time data when read into
;UCI will be handled differently compared to NOAA data.
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


;Sort UCI data
print, 'i-Sort UCI data'
cen_lat = findgen(91) * 2 - 90
cen_lon = findgen(144) * 2.5 - 180
;change the first and last value of the cen_lat array to match  GMAO 2° x 2.5° lattitudes
cen_lat[ 0 ] = -89.5
cen_lat[ n_elements(cen_lat) - 1 ] = 89.5
new_time = fltarr( n_elements(avg_ss) )
new_time = !NULL
new_avg = fltarr( n_elements(avg_ss) )
new_avg = !NULL
new_lat = fltarr( n_elements(avg_ss) )
new_lon = fltarr( n_elements(avg_ss) )
;debug_x = fltarr( n_elements(avg_ss) )
;debug_x = !NULL

for i = 0, n_elements(cen_lat) - 1 do begin
	for j = 0, n_elements(cen_lon) - 1 do begin
		x = where(lat1_ss eq cen_lat[ i ] and lon1_ss eq cen_lon[ j ], count)
		if (count gt 0) then begin
			temp_time = unix_date_ss[ x ]
			temp_avg = avg_ss[ x ]
			a = sort(temp_time)
			temp_time = temp_time[ a ]
			temp_avg = temp_avg[ a ]
			new_lat[ n_elements(new_avg) : n_elements(new_avg) + $
				n_elements(temp_avg) - 1 ] = cen_lat[ i ]
			new_lon[ n_elements(new_avg) : n_elements(new_avg) + $
				n_elements(temp_avg) - 1 ] = cen_lon[ j ]
			new_time = [new_time, temp_time]
			new_avg = [new_avg, temp_avg]
;			debug_x = [debug_x , x]
		endif
	endfor
endfor

unix_date_ss = reverse( new_time )
avg_ss = reverse( new_avg )
lat1_ss = reverse( new_lat )
lon1_ss = reverse( new_lon )

;verifying the code - CODE IS GOOD
;infile = "/home/excluded-from-backup/ethane/IDL/temp_file/debug_temporal_ihr_full.dat"

;openw, lun, infile, /get_lun

;for i = 0, n_elements(avg_ss) - 1 do begin
;	printf, lun, unix_date_ss[ i ], lat1_ss[ i ], lon1_ss[ i ], avg_ss[ i ], format= "(4F15.2)"
;endfor

;free_lun, lun
;debug_x = debug_x[ sort(debug_x) ]
;for i = 0, n_elements(debug_x) - 1 do begin
;	if ( debug_x[i] ne i ) then print, i
;endfor


;Because time in the Simpson et al. data set denoted time as Unix time, these lines convert
;Unix time to Geos_chem time and stored in tau_time array
tau_time_ss = (unix_date_ss/3600 - 131496)
tau_time_ss = tau2yymmdd(tau_time_ss, /GEOS1)

;Convert unix time to decimal year
time_ss = (unix_date_ss/3600 - 131496)/8760 + 1985
;-
;-
;-
;-
;-
;Read OGI data
print, 'c-Read OGI data'
print,  ' '
infile1="/home/excluded-from-backup/ethane/data/raw_data/OGI/KhalilEtAl_data_set.dat"
infile2= "/home/excluded-from-backup/ethane/data/raw_data/OGI/KhalilEtAl_locations.dat"

;get number of lines in the data file
n1_ogi = file_lines(infile1)
n2_ogi = file_lines(infile2)

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

;Because time in the OGI data set denoted time as Unix time, these lines convert
;Unix time to Geos_chem time and stored in tau_time array
tau_time_ogi = (unix_date_ogi/3600 - 131496)
tau_time_ogi = tau2yymmdd(tau_time_ogi, /GEOS1)
;-
;-
;-
;-
;-
;Build Structure for the data
print, 'II-Build data structure from NOAA, UCI, OGI'
print,  ' '
noaa = {site : location_noaa, $
		ratio : avg_noaa, $
		lat : lat1_noaa, $
		lon : lon1_noaa, $
		year : year_noaa, $
		month : month_noaa, $
		residual : fltarr(n_elements(avg_noaa) ) }

uci = {site : location_ss, $
		ratio : avg_ss, $
		lat : lat1_ss, $
		lon: lon1_ss, $
		year : tau_time_ss.year, $
		month : tau_time_ss.month, $
		residual : fltarr(n_elements(avg_ss) ) }

ogi = {site : location_ogi, $
		ratio : avg_ogi, $
		lat : lat1_ogi, $
		lon : lon1_ogi, $
		year : tau_time_ogi.year, $
		month : tau_time_ogi.month, $
		residual : fltarr(n_elements(avg_ogi) ) }

print, 'a-NOAA Structure: '
help, noaa, /str
print, 'b-UCI Structure: '
help, uci, /str
print, 'c-OGI Structure: '
help, ogi, /str
;-
;-
;-
;-
;-
;Detrend all data
print, ' '
print, 'III-Detrend'
print, ' '

;Read NOAA station ID
print, 'a-Detrend NOAA'
print, '  i-Read NOAA station IDs'
;link to file
noaa_id_file = "/home/jakechung/IDL/myIDL/NOAA_stations_IDs.dat"
noaa_id_size = file_lines( noaa_id_file )
noaa_id = strarr( noaa_id_size )
id = ' '
;open to read
openr, id_lun, noaa_id_file, /get_lun

;read the data into the noaa_ids array
for i = 0, noaa_id_size - 1 do begin
	readf, id_lun, id
	noaa_id [ i ] = id
endfor

free_lun, id_lun

print, '  ii-Detrend process'
;loop through to each site to detrend
for i = 0, noaa_id_size - 1 do begin
	idx = where( strmatch ( noaa.site[ * ], noaa_id[ i ], /FOLD_CASE) EQ 1 )
	temp_ratio = noaa.ratio[ idx ]
	temp_month = noaa.month[ idx ]
	temp_year = noaa.year[ idx ]
	temp_residual = fltarr( n_elements(temp_ratio) )
	holder = fltarr( n_elements(temp_ratio) )
	b = fltarr( 13 )
	for j = 1, 12 do begin
		a = where( temp_month eq j, count )
		if count gt 0 then $
		b[ j ] = mean( temp_ratio[ a ], /NAN)
		holder[ a ] = b[ j ]
	endfor
	temp_residual = temp_ratio - holder
	noaa.residual[ idx ] = temp_residual
endfor

;detrend UCI
print, 'b-Detrend UCI'

test_arr = fltarr(5000)
test_arr = !NULL
for i = 0, n_elements(cen_lat) - 1 do begin
	for j = 0, n_elements(cen_lon) - 1 do begin
		x = where( uci.lat eq cen_lat[i] and uci.lon eq cen_lon[j], count)
		if (count gt 0) then begin
			;print, count
			;to be continued
		endif
	endfor
endfor

end
