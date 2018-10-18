FUNCTION read_uci8496

;this function read in the  1984 - 1996 UCI data and build the data structure

;compile_opt idl2

;read 84_96 UCI data
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

;build structure
uci = { ratio : sratio, $
	lat : slat, $
	lon: slon, $
	year : syear, $
	month : smonth }

;print the structure to a text file to validate
;infile = '/home/excluded-from-backup/ethane/IDL/temp_file/uci_8496_data_struct.dat'

; openw, lun, infile, /get_lun
;
; printf, lun, 'Month|Year|Latitude|Longitude|Ratio'
; for i = 0, n_elements(sratio) - 1 do begin
;   printf, lun, smonth[i], syear[i], slat[i], slon[i], sratio[i]
; endfor
;
; free_lun, lun

;return the data structure
return, uci


end
