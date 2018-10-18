FUNCTION read_uci

;compile_opt idl2
  ;this function read the UCI data set and put that data into a structure
  ;the function has no input


infile1= "/home/excluded-from-backup/ethane/IDL/temp_file/modded_uci_data_v2.dat"
infile2= "/home/excluded-from-backup/ethane/data/raw_data/UCI/uci_stations_v2.dat"

n1_ss = file_lines(infile1)
n2_ss = file_lines(infile2)

  ;create an empty array to store the data that was read in from the data files
  ;since the time data in the UCI data set is in a different format,
  ;the time data when read into
  ;UCI will be handled differently compared to NOAA data.
unix_date_ss = fltarr(n1_ss)
lat = fltarr(n1_ss)
lon = fltarr(n1_ss)
ratio = fltarr(n1_ss)
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
	lat[ i ] = lat0
	lon[ i ] = lon0
	ratio[ i ] = avg0
endfor

  ;read float data of the location file
for i = 0, n2_ss-1 do begin
	readf, iunit2, location0
	location_ss[ i ] = location0
endfor
  ;close input file
free_lun, iunit1, iunit2

;Because time in the Simpson et al. data set denoted time as Unix time, these lines convert
;Unix time to Geos_chem time and stored in tau_time array
tau_time_ss = (unix_date_ss/3600 - 131496)
tau_time_ss = tau2yymmdd(tau_time_ss, /GEOS1)

;convert to type float of the tau_time data array
year = fltarr(n_elements(tau_time_ss.year))
month = fltarr(n_elements(tau_time_ss.month))

year = float(tau_time_ss.year[*])
month = float(tau_time_ss.month[*])

;this is an old section on soting the UCI data by latitude and time
;but the new algorithm requires the data to be sorted in latitudes and
;longitudes, so this sorting algorithm is not used anymore, and also
;to be consistent with the read_uci8496 function
;   ;Sort UCI data
; cen_lat = findgen(91) * 2 - 90
; cen_lon = findgen(144) * 2.5 - 180
;   ;change the first and last value of the cen_lat array to match  GMAO 2° x 2.5° lattitudes
; cen_lat[ 0 ] = -89.5
; cen_lat[ n_elements(cen_lat) - 1 ] = 89.5
; new_time = fltarr( n_elements(avg_ss) )
; new_time = !NULL
; new_avg = fltarr( n_elements(avg_ss) )
; new_avg = !NULL
; new_lat = fltarr( n_elements(avg_ss) )
; new_lon = fltarr( n_elements(avg_ss) )
; ;debug_x = fltarr( n_elements(avg_ss) )
; ;debug_x = !NULL
;
; for i = 0, n_elements(cen_lat) - 1 do begin
; 	for j = 0, n_elements(cen_lon) - 1 do begin
; 		x = where(lat1_ss eq cen_lat[ i ] and lon1_ss eq cen_lon[ j ], count)
; 		if (count gt 0) then begin
; 			temp_time = unix_date_ss[ x ]
; 			temp_avg = avg_ss[ x ]
; 			a = sort(temp_time)
; 			temp_time = temp_time[ a ]
; 			temp_avg = temp_avg[ a ]
; 			new_lat[ n_elements(new_avg) : n_elements(new_avg) + $
; 				n_elements(temp_avg) - 1 ] = cen_lat[ i ]
; 			new_lon[ n_elements(new_avg) : n_elements(new_avg) + $
; 				n_elements(temp_avg) - 1 ] = cen_lon[ j ]
; 			new_time = [new_time, temp_time]
; 			new_avg = [new_avg, temp_avg]
; ;			debug_x = [debug_x , x]
; 		endif
; 	endfor
; endfor
;
; unix_date_ss = reverse( new_time )
; avg_ss = reverse( new_avg )
; lat1_ss = reverse( new_lat )
; lon1_ss = reverse( new_lon )

;this is the new algorithm of the for this function. Copied from read_uci8496
;algorithm to sort the data. Subsequence algorithms of the project requires the
;latitudes and longtitudes to be sorted in ascended order, so the arrays will be sorted and
;rebuilt before the data structure is built
;x contains the index of the latitudinally sorted array
x = sort(lat)
lat = lat[x]
ratio = ratio[x]
lon = lon[x]
year = year[x]
month = month[x]

;slat, slon, sratio, syear, smonth are the new arrays that will store the sorted
;values
slat = fltarr(n_elements(lat))
slat = lat
slon = fltarr(n_elements(lon))
slon = !NULL
sratio = fltarr(n_elements(ratio))
sratio = !NULL
syear = fltarr(n_elements(year))
syear = !NULL
smonth = fltarr(n_elements(month))
smonth = !NULL
;algorithm to sort longtides of each latitude value
for i = 0, n_elements(lat) - 1 do begin
  lat0 = lat[i]
  y = where(lat eq lat0, count)
  lon_temp = lon[y]
  ratio_temp = ratio[y]
  year_temp = year[y]
  month_temp = month[y]
  z = sort(lon_temp)
  lon_temp = lon_temp[z]
  ratio_temp = ratio_temp[z]
  year_temp = year_temp[z]
  month_temp = month_temp[z]
  ;rebuild the sorted data into the new arrays
  slon = [slon, lon_temp]
  sratio = [sratio, ratio_temp]
  syear = [syear, year_temp]
  smonth = [smonth, month_temp]
  i = i + count - 1
  if (i gt n_elements(lat)) then break
endfor


;build the data structure for UCI data
uci = { ratio : sratio, $
		lat : slat, $
		lon: slon, $
		year : syear, $
		month : smonth }

;print the structure to a text file to validate
;infile = '/home/excluded-from-backup/ethane/IDL/temp_file/uci_data_struct.dat'

;openw, lun, infile, /get_lun

;printf, lun, 'Month|Year|Latitude|Longitude|Ratio'
;for i = 0, n_elements(uci.ratio) - 1 do begin
;  printf, lun, uci.month[i], uci.year[i], uci.lat[i], uci.lon[i], uci.ratio[i]
;endfor

;free_lun, lun



return, uci

end
