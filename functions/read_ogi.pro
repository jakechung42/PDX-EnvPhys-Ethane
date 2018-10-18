function read_ogi

  ;this function read the ogi data and build the data structure

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

;convert to type float of the tau_time_ogi
year_ogi = fltarr(n_elements(tau_time_ogi.year))
month_ogi = fltarr(n_elements(tau_time_ogi.month))

year_ogi = float(tau_time_ogi.year[*])
month_ogi = float(tau_time_ogi.month[*])


ogi = {site : location_ogi, $
		ratio : avg_ogi, $
		lat : lat1_ogi, $
		lon : lon1_ogi, $
		year : year_ogi, $
		month : month_ogi }

return, ogi

end
