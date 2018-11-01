PRO DESEASON_OGI_MAIN

;Author: Jake Chung
;REU 2016/ AGU 2016
;Description: this program deseasonalize the OGI data and then store the residual data at
;"/home/jakechung/IDL/myIDL/temp_folder/temp_residual_OGI.dat"

;Notes: plotting procedures can be added to plot the residual if needed be.
;The residual follows a Gaussian distribution (refer to figure 2, AGU poster) 

COMPILE_OPT IDL2				;set compile options

;=====Read in the files of OGI data set=====
;get the files
infile1= "/home/jakechung/IDL/myIDL/KhalilEtAl_data_set.dat"
infile2= "/home/jakechung/IDL/myIDL/KhalilEtAl_locations.dat"

;get number of lines in the data file
n1 = file_lines(infile1)
n2 = file_lines(infile2)

print, 'The input has ', n1, '   data lines with', n2, '  locations'

;create empty arrays to store the data that was read in lon(longitudes), lat(latitudes), 
;avg(concentration), date and locations
unix_date= fltarr(n1)
avg= fltarr(n1)
location= strarr(n2)
lat1= fltarr(n1)
lon1= fltarr(n1)

;FLOAT input variable
location0= ' '
date0= 0.0
avg0= 0.0
lat0= 0.0
lon0= 0.0

;open input files
openr, iunit1, infile1, /get_lun
openr, iunit2, infile2, /get_lun

;read float data of the data file
for i = 0, n1-1 do begin
	;read one line from the file
	readf, iunit1, date0, avg0, lat0, lon0
	;store the values
	unix_date[i] = date0
	avg[i] = avg0
	lat1[i] = lat0
	lon1[i] = lon0
endfor

;read float data of the location file
for i = 0, n2-1 do begin
	readf, iunit2, location0
	location[i] = location0
endfor

;close input file
free_lun, iunit1
free_lun, iunit2

;Because time in the OGI data set denoted time as Unix time, these lines convert Unix time to Geos_chem time and stored in tau_time array
tau_time = (unix_date/3600 - 131496)
tau_time = tau2yymmdd(tau_time, /GEOS1)

;Convert unix time to decimal year
time = (unix_date/3600 - 131496)/8760 + 1985
;=====Finish reading in the OGI data set=====


;=====Start deseasonal calculations=====
;making the temp file to store the residual values
temp_residual = "/home/jakechung/IDL/myIDL/temp_folder/temp_residual_OGI.dat"
openw, iunit4, temp_residual, /get_lun

station = [ 'Barrow', 'Cape Meares', 'Mauna Loa', 'Samoa', 'Cape Grim, Tasmania', 'South Pole' ]

;------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
;start a giant loop right here to loop in, read all the sites IDs and deseasonalize the entire data set.
for r = 0, 5 do begin  ;<<<<<<<< loop starts here
z = where(strmatch(location[*], station[r], /FOLD_CASE) EQ 1)
print, 'Station ', station[r],' Index on the data set is ', z[0], ' to ', z[n_elements(z)-1]
avg_ice = avg[z]
time_ice = time[z]
month_ice = tau_time.month[z]
day_ice = tau_time.day[z]
year_ice = tau_time.year[z]
lat_ice = lat1[z]
lon_ice = lon1[z]

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
					printf, iunit4, residue_arr[i], lat_ice[j], lon_ice[j], year_ice[ j ], temp_month[i]
				endif
			endfor
			j = a[n_elements(a) - 1]					
			a = 0
		endif else if ( year_ice [ j ] NE year_ice [ b ] ) then j = j + 1
endfor

endfor ;<<<<<<<< The big loop to loop in all the stations IDs ends here

free_lun, iunit4

end