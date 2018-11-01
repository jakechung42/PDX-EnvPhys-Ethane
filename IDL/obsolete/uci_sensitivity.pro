PRO UCI_SENSITIVITY

;09/21 - Jake Chung
;This program calculates the latitude bin average to determine the best way to average 
;the UCI data (whether to pick more than 30 or 20 data points)

;Read UCI data
;read in the modded UCI data where the lon and lat are 
;changed from uci_mod_lat_lon
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
		endif
	endfor
endfor

unix_date_ss = reverse( new_time )
avg_ss = reverse( new_avg )
lat1_ss = reverse( new_lat )
lon1_ss = reverse( new_lon )

;Because time in the Simpson et al. data set denoted time as Unix time, these lines convert
;Unix time to Geos_chem time and stored in tau_time array
tau_time_ss = (unix_date_ss/3600 - 131496)
tau_time_ss = tau2yymmdd(tau_time_ss, /GEOS1)

uci = {site : location_ss, $
		ratio : avg_ss, $
		lat : lat1_ss, $
		lon: lon1_ss, $
		year : tau_time_ss.year, $
		month : tau_time_ss.month }		

		
		
		
		
;Filter out sites with more than 10 data points
gt10_site_lat = fltarr(100)
gt10_site_lon = fltarr(100)
gt10_site_lat = !NULL
gt10_site_lon = !NULL
for i = 0, n_elements(cen_lat) - 1 do begin
	for j = 0, n_elements(cen_lon) - 1 do begin
		x = where( lat1_ss eq cen_lat[i] and lon1_ss eq cen_lon[j], count)
		if (count gt 10) then begin
			a = x[0]
			temp_lat = lat1_ss[a]
			temp_lon = lon1_ss[a]
			gt10_site_lat = [gt10_site_lat , temp_lat]
			gt10_site_lon = [gt10_site_lon , temp_lon]
		endif
	endfor
endfor
lat = [-90,-60,-30,0,30,60,90]
;isolate the data with more than 10 data points
gt10_avg = fltarr(5000)
gt10_lat = fltarr(5000)
gt10_lon = fltarr(5000)
gt10_avg = !NULL
gt10_lat = !NULL
gt10_lon = !NULL

for i = 0, n_elements(gt10_site_lat) - 1 do begin
	x = where(uci.lat eq gt10_site_lat[i] and uci.lon eq gt10_site_lon[i], count)
	lat_holder = fltarr(n_elements(x))
	lon_holder = fltarr(n_elements(x))
	lat_holder[*] = gt10_site_lat[i]
	lon_holder[*] = gt10_site_lon[i]
	gt10_avg = [gt10_avg , uci.ratio[x]]
	gt10_lat = [gt10_lat , lat_holder]
	gt10_lon = [gt10_lon , lon_holder]
endfor		

;Seperate into latitude bins of the data that has more than 10 data points
gt10_uci_ratio = fltarr(6)
print, 'GT 10: '
for i = 0, 5 do begin
	x = where(gt10_lat[*] gt lat[i] and gt10_lat[*] lt lat[i+1], count)
	if (count gt 0) then begin
		gt10_uci_ratio[i] = mean(gt10_avg[x])
	endif
endfor




;Filter out sites with more than 30 data points
gt30_site_lat = fltarr(100)
gt30_site_lon = fltarr(100)
gt30_site_lat = !NULL
gt30_site_lon = !NULL
for i = 0, n_elements(cen_lat) - 1 do begin
	for j = 0, n_elements(cen_lon) - 1 do begin
		x = where( lat1_ss eq cen_lat[i] and lon1_ss eq cen_lon[j], count)
		if (count gt 30) then begin
			a = x[0]
			temp_lat = lat1_ss[a]
			temp_lon = lon1_ss[a]
			gt30_site_lat = [gt30_site_lat , temp_lat]
			gt30_site_lon = [gt30_site_lon , temp_lon]
		endif
	endfor
endfor
lat = [-90,-60,-30,0,30,60,90]
;isolate the data with more than 10 data points
gt30_avg = fltarr(5000)
gt30_lat = fltarr(5000)
gt30_lon = fltarr(5000)
gt30_avg = !NULL
gt30_lat = !NULL
gt30_lon = !NULL

for i = 0, n_elements(gt30_site_lat) - 1 do begin
	x = where(uci.lat eq gt30_site_lat[i] and uci.lon eq gt30_site_lon[i], count)
	lat_holder = fltarr(n_elements(x))
	lon_holder = fltarr(n_elements(x))
	lat_holder[*] = gt30_site_lat[i]
	lon_holder[*] = gt30_site_lon[i]
	gt30_avg = [gt30_avg , uci.ratio[x]]
	gt30_lat = [gt30_lat , lat_holder]
	gt30_lon = [gt30_lon , lon_holder]
endfor		

;Seperate into latitude bins of the data that has more than 10 data points
gt30_uci_ratio = fltarr(6)
print, 'GT 10: '
for i = 0, 5 do begin
	x = where(gt30_lat[*] gt lat[i] and gt30_lat[*] lt lat[i+1], count)
	if (count gt 0) then begin
		gt30_uci_ratio[i] = mean(gt30_avg[x])
	endif
endfor
print, gt30_uci_ratio




;Filter out sites with more than 60 data points
gt60_site_lat = fltarr(100)
gt60_site_lon = fltarr(100)
gt60_site_lat = !NULL
gt60_site_lon = !NULL
for i = 0, n_elements(cen_lat) - 1 do begin
	for j = 0, n_elements(cen_lon) - 1 do begin
		x = where( lat1_ss eq cen_lat[i] and lon1_ss eq cen_lon[j], count)
		if (count gt 60) then begin
			a = x[0]
			temp_lat = lat1_ss[a]
			temp_lon = lon1_ss[a]
			gt60_site_lat = [gt60_site_lat , temp_lat]
			gt60_site_lon = [gt60_site_lon , temp_lon]
		endif
	endfor
endfor
lat = [-90,-60,-30,0,30,60,90]
;isolate the data with more than 10 data points
gt60_avg = fltarr(5000)
gt60_lat = fltarr(5000)
gt60_lon = fltarr(5000)
gt60_avg = !NULL
gt60_lat = !NULL
gt60_lon = !NULL

for i = 0, n_elements(gt60_site_lat) - 1 do begin
	x = where(uci.lat eq gt60_site_lat[i] and uci.lon eq gt60_site_lon[i], count)
	lat_holder = fltarr(n_elements(x))
	lon_holder = fltarr(n_elements(x))
	lat_holder[*] = gt60_site_lat[i]
	lon_holder[*] = gt60_site_lon[i]
	gt60_avg = [gt60_avg , uci.ratio[x]]
	gt60_lat = [gt60_lat , lat_holder]
	gt60_lon = [gt60_lon , lon_holder]
endfor		

;Seperate into latitude bins of the data that has more than 10 data points
gt60_uci_ratio = fltarr(6)
print, 'GT 10: '
for i = 0, 5 do begin
	x = where(gt60_lat[*] gt lat[i] and gt60_lat[*] lt lat[i+1], count)
	if (count gt 0) then begin
		gt60_uci_ratio[i] = mean(gt60_avg[x])
	endif
endfor
print, gt60_uci_ratio






;Seperate into latitude bins of all data
all_uci_ratio = fltarr(6)

print, 'All data: '
for i = 0, 5 do begin
	x = where(uci.lat gt lat[i] and uci.lat lt lat[i+1], count)
;	print, count
	if (count gt 0) then begin 
		all_uci_ratio[i] = mean(uci.ratio[x])
	endif
endfor

x_axis = [-75, -45, -15, 15, 45, 75]

cgplot, x_axis, all_uci_ratio, psym = 2, xrange = [-90, 90], xticks = 6
	cgplot, x_axis, gt10_uci_ratio, psym = 1, /overplot, color = 'red'

end

