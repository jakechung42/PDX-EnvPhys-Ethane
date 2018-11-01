PRO DESEASON_UCI_MAIN

;Author: Jake Chung
;REU 2016/ AGU 2016
;Description: this program deseasonalize the UCI data and then store the residual data at
;"/home/jakechung/IDL/myIDL/temp_folder/temp_residual_Simpson.dat"

;Notes: plotting procedures can be added to plot the residual if needed be.
;The residual follows a Gaussian distribution (refer to figure 2, AGU poster) 

COMPILE_OPT IDL2				;set compile options

;=====Read in the data files of the Simpson et al. data=====
;get the files
infile1= "/home/jakechung/IDL/myIDL/data_set.dat"
infile2= "/home/jakechung/IDL/myIDL/location_set.dat"

;get number of lines in the data file
n1 = file_lines(infile1)
n2 = file_lines(infile2)

print, 'The input has ', n1, '   data lines with', n2, '  locations'

;create empty arrays to store the data that was read in lon(longitudes), lat(latitudes), 
;avg(concentration), date and locations
unix_date= fltarr(n1)
lat1= fltarr(n1)
lon1= fltarr(n1)
avg= fltarr(n1)
location= strarr(n2)

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
for i = 0, n1-1 do begin
	;read one line from the file
	readf, iunit1, date0, lat0, lon0, avg0
	;store the values
	unix_date[i] = date0
	lat1[i] = lat0
	lon1[i] = lon0
	avg[i] = avg0
endfor

;read float data of the location file
for i = 0, n2-1 do begin
	readf, iunit2, location0
	location[i] = location0
	
endfor

;close input file
free_lun, iunit1, iunit2

;Because time in the Simpson et al. data set denoted time as Unix time, these lines convert Unix time to Geos_chem time and stored in tau_time array
tau_time = (unix_date/3600 - 131496)
tau_time = tau2yymmdd(tau_time, /GEOS1)

;Convert unix time to decimal year
time = (unix_date/3600 - 131496)/8760 + 1985
;=====Finish reading the Simpson et al. data file=====



;=====Start deseasonal calculations=====
;making the temp file to store the residual values
temp_residual = "/home/jakechung/IDL/myIDL/temp_folder/temp_residual_Simpson.dat"
openw, iunit4, temp_residual, /get_lun

;Cannot go by staion ID, Simpson data set has very inconsistent naming scheme, go by latitudes instead.
;These latitudes are the latitudes of the stations with more than 30 data points.
lat_var_GT30 = [       71.3000, $
      64.4900, $
      60.7500, $
      57.8300, $
      57.8000, $
      42.8300, $
      30.4000, $
      29.9400, $
      23.4400, $
      22.8700, $
      20.2600, $
      15.2000, $
      13.3600, $
      9.52000, $
      7.05000, $
      6.92000, $
      1.42000, $
     -8.52000, $
     -14.2300, $
     -21.2300, $
     -29.0200, $
     -34.9000, $
     -36.8200, $
     -42.7200, $
     -43.8800	]
;------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
;start a giant loop right here to loop in, read all latitudes and deseasonalize the entire data set.
for r = 0, n_elements(lat_var_GT30) - 1 do begin  ;<<<<<<<< loop starts here
z = where(lat1[*] EQ lat_var_GT30[r])

avg_ice = avg[z]
time_ice = time[z]
month_ice = tau_time.month[z]
day_ice = tau_time.day[z]
year_ice = tau_time.year[z]
lat_ice = lat1[z]
lon_ice = lon1[z]

;because the Simpson data dates are not in order, this will provide a fix for that by creating a new sorted array for the search to run
;This sort is only sorting the individual station that the for loop is calling out. It's not sorting the entire UCI array.
ordered_year = year_ice[ sort( year_ice ) ]
ordered_avg = avg_ice[ sort(year_ice) ]
ordered_month = month_ice[ sort(year_ice) ]

;these arrays store the months of the measured data set. 
;jan array has all the JAN of the measure data set
mar = where(ordered_month[*] EQ 3)
jun = where(ordered_month[*] EQ 6)
sep = where(ordered_month[*] EQ 9)
dec = where(ordered_month[*] EQ 12)

mar_avg_var = mean( ordered_avg [ mar ] ) 
jun_avg_var = mean( ordered_avg [ jun ] ) 
sep_avg_var = mean( ordered_avg [ sep ] ) 
dec_avg_var = mean( ordered_avg [ dec ] ) 

;Now find the residue of the measured and the simulated data by subtracting the simulated data from the measured data
;The monthly averages of simulated data subtracts the measured data for the corresponding month
residue_arr = fltarr(n_elements(avg_ice))
for j = 0, n_elements(avg_ice) - 1 do begin
	if (j NE n_elements(avg_ice) - 1) then b = j + 1 else b = j
		if ( ordered_year [ j ] EQ ordered_year [ b ] ) then begin
			a = where( ordered_year EQ ordered_year [ j ] )
			print, a
			temp_avg_var = ordered_avg [ a ]
			temp_month = ordered_month [ a ] 
			for i = 0, n_elements(a) - 1 do begin
				if (ordered_month[i] EQ 3) then $
					residue_arr[ i ] = mar_avg_var - temp_avg_var[ i ] else $
				if (ordered_month[i] EQ 6) then $
					residue_arr[ i ] = jun_avg_var - temp_avg_var[ i ] else $
				if (ordered_month[i] EQ 9) then $
					residue_arr[ i ] = sep_avg_var - temp_avg_var[ i ] else $
				if (ordered_month[i] EQ 12) then $
					residue_arr[ i ] = dec_avg_var - temp_avg_var[ i ] else $
					residue_arr[i] = 0
				if ( residue_arr[i] NE 0 ) then begin
					printf, iunit4, residue_arr[i], lat_ice[0], lon_ice[0], ordered_year[ j ], ordered_month[ i ]
				endif
			endfor
			j = a[n_elements(a) - 1]					
			a = 0
		endif else if ( ordered_year [ j ] NE ordered_year [ b ] ) then j = j + 1
endfor

endfor ;<<<<<<<< The big loop to loop in all the stations IDs ends here

free_lun, iunit4

end