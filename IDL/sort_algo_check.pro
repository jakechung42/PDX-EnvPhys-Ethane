PRO sort_algo_check

;this program takes the original file and the sorted file
;to validate the sorting algorithm
;mainly used for the 2 UCI data sets when the latitudes and longitudes
;are sorted.

compile_opt idl2

;read the original data
infile = '/home/excluded-from-backup/ethane/IDL/temp_file/modded_84_96_uci_data.dat'
openr, uci_lun, infile, /get_lun

size = file_lines(infile)

;create the arrays to store the data
utime = fltarr(size)
ratio = fltarr(size)
lon = fltarr(size)
lat = fltarr(size)

for i = 0, size - 1 do begin
  readf, uci_lun, c1, c2, c3, c4
  utime[i] = c1
  lat[i] = c2
  lon[i] = c3
  ratio[i] = c4
endfor

;release the lun
free_lun, uci_lun

;remove -999 data
neg999 = where(lon eq -999)
utime = RemoveRows(rotate(utime, 1), neg999)
lat = RemoveRows(rotate(lat, 1), neg999)
lon = RemoveRows(rotate(lon, 1), neg999)
ratio = RemoveRows(rotate(ratio, 1), neg999)

;convert the unix time to GEOS_CHEM time and to year and month using the
;tau2yymmdd routine
tau_time = (utime/3600 - 131496)
tau_time = tau2yymmdd(tau_time, /GEOS1)

;convert to type float of the tau_time data array
year = fltarr(n_elements(tau_time.year))
month = fltarr(n_elements(tau_time.month))

year = float(tau_time.year[*])
month = float(tau_time.month[*])



;read in the sorted data
infile = '/home/excluded-from-backup/ethane/IDL/temp_file/uci_8496_data_struct.dat'
openr, uci_lun, infile, /get_lun

size = file_lines(infile)

;create the arrays to store the data
month1 = fltarr(size)
year1 = fltarr(size)
ratio1 = fltarr(size)
lon1 = fltarr(size)
lat1 = fltarr(size)

for i = 0, size - 1 do begin
  readf, uci_lun, c1, c2, c3, c4, c5
  month1[i] = c1
  year1[i] = c2
  lat1[i] = c3
  lon1[i] = c4
  ratio1[i] = c5
endfor

;release the lun
free_lun, uci_lun



;check the data by comparing the dates, lat, lon, ratio of one
;each data point between the original and the sorted.

for i = 0, n_elements(ratio) - 1 do begin
  idx = where(year1 eq year[i] and $
    month1 eq month[i] and $
    lat1 eq lat[i] and $
    lon1 eq lon[i] and $
    ratio1 eq ratio[i], count)
  if (count lt 1) then begin
    print, 'Value does not match'
    print, i
  endif
endfor


end
