PRO DEBUG_SORT_MOD_UCI

;08/31/2017 - Jake Chung: This program compares the modified uci data from 
;uci_mod_lat_lon.pro with the debug output from temporal_ihr_full to make sure
;the sorting algorithm is valid.

;09/01/2017 - Jake Chung:
;	Also plot the timeseries for each site.

compile_opt idl2

;read in the original 
infile1= "/home/excluded-from-backup/ethane/IDL/temp_file/modded_uci_data_v2.dat"

;read in the sorted
infile2 = "/home/excluded-from-backup/ethane/IDL/temp_file/debug_temporal_ihr_full.dat"

n1 = file_lines(infile1)
n2 = file_lines(infile2)

openr, iunit1, infile1, /get_lun
openr, iunit2, infile2, /get_lun

uci_ini = fltarr( 4 , n1 )
uci_sort = fltarr( 4 , n2 )

for i = 0, n1 - 1 do begin
	readf, iunit1, date1, lat1, lon1, var1
	uci_ini[ 0 , i ] = date1
	uci_ini[ 1 , i ] = lat1
	uci_ini[ 2 , i ] = lon1
	uci_ini[ 3 , i ] = var1
endfor

for i = 0 , n2 - 1 do begin
	readf, iunit2, date2, lat2, lon2, var2
	uci_sort[ 0 , i ] = date2
	uci_sort[ 1 , i ] = lat2
	uci_sort[ 2 , i ] = lon2
	uci_sort[ 3 , i ] = var2
endfor

free_lun, iunit1,iunit2

cen_lat = findgen(91) * 2 - 90
;change the first and last value of the cen_lat array to match  GMAO 2° x 2.5° lattitudes
cen_lat[ 0 ] = -89.5
cen_lat[ n_elements(cen_lat) - 1 ] = 89.5
cen_lon = findgen(144) * 2.5 - 180

for i = 0, n_elements(uci_ini[ 1 , * ]) - 1 do begin
	if (uci_ini[ 1 , i ] ne uci_sort[ 1 , i ]) then print, 'problem'
	x = where(uci_sort[0 , *] eq uci_ini[0 , i] and $
		uci_sort[3 , *] eq uci_ini[3 , i] and $
		uci_sort[2 , *] eq uci_ini[2 , i] and $
		uci_sort[1 , *] eq uci_ini[1 , i], count)
	if (count eq 0) then print, 'problem'
endfor

;converting unix time to dates
tau_time_ss = (uci_sort[0 , *]/3600 - 131496)
tau_time_ss = tau2yymmdd(tau_time_ss, /GEOS1)

time_uci = (uci_sort[0 , *]/3600 - 131496)/8760 + 1985


open_device, /ps, /color, file='temp.eps', /landscape, margin=0.05, xsize = 10.0, ysize = 7.5
!x.thick=5
!y.thick=5
!p.font=0
!p.thick =3

!p.multi = [0, 3, 5, 0, 0]
;multiplot, /default    ;resets multiplot settings
;multiplot, [3,5], ygap=0.02, xgap=0.025;  sets up multiplot 
n = 0
d = 0s
for i = 0, n_elements(cen_lat) - 1 do begin
	for j = 0, n_elements(cen_lon) - 1 do begin
		x = where( uci_sort[1 , *] eq cen_lat[i] and uci_sort[2 , *] eq cen_lon[j], count)
		if (count gt 15) then begin
			n = n + 1
			d = d + 1s
			lat = string(cen_lat[i], format = '(F0.2)')
			lon = string(cen_lon[j], format = '(F0.2)')
			n_data = string(n_elements(x), format = '(I0)')
			id = string(d, format = '(I0)')
			cgplot, time_uci[x], uci_sort[3 , x], $
				xticks = 5, $
				char = 1.7, $
				title = id + '/ lat: ' + lat + ' lon: ' + lon + ' n: ' + n_data				
			if n eq 15 then begin
				close_device
				spawn, 'gv temp.eps'
				b = ' '
				read, b, prompt='Enter anything to continue: '
				n = 0
				open_device, /ps, /color, file='temp.eps', /landscape, margin=0.05, xsize = 10.0, ysize = 7.5
				!x.thick=5
				!y.thick=5
				!p.font=0
				!p.thick =3
				!p.multi = [0, 3, 5, 0, 0] 
			endif
		endif
	endfor
endfor

close_device
spawn, 'gv temp.eps'

end
	