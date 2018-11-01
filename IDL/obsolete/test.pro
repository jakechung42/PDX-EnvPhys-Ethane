PRO TEST

; 07/05/18 - Jake Chung
; This is a test file to test the functionality of Atom editor
; Execute simple commands

;read in the DAT file for OGI and UCI data set.
infile = '/home/jakechung/IDL/myIDL/IHR_w2sites_OGI_UCI.dat'
;Note: because the String in the data file is problematic,
;to solve the problem, OGI is 42 and UCI is 27

;get number of lines in the data file
n1 = file_lines ( infile )

;create temporary arrays to store the values from the ASCII file
temp_time_ss_ogi = fltarr( n1 )
temp_ratio_ss_ogi = fltarr( n1 )
temp_lon_ss_ogi = fltarr( n1 )
temp_lat_ss_ogi = fltarr( n1 )
temp_network_ss_ogi = fltarr( n1 )

time = 0.0
ratio = 0.0
lon = 0.0
lat = 0.0
network = 0.0

openr, iunit, infile, /get_lun

for i = 0, n1 - 1 do begin
	readf, iunit, time, ratio, lon, lat, network
	temp_time_ss_ogi[ i ] = time
	temp_ratio_ss_ogi[ i ] = ratio
	temp_lon_ss_ogi[ i ] = lon
	temp_lat_ss_ogi[ i ] = lat
	temp_network_ss_ogi[ i ] = network
endfor

free_lun, iunit

print, ratio

end
