PRO WMAP_ALL_STATIONS

;plot all stations of the data sets on a map using the Gamap routines

COMPILE_OPT IDL2				;set compile options

;=======Start pulling in deseasonalized data values=======
temp_deseason = "/home/jakechung/IDL/myIDL/temp_folder/temp_residual_noaa.dat"
openr, deseason, temp_deseason, /get_lun

deseasonal_lines = file_lines(temp_deseason)
print, deseasonal_lines
deseasonal_arr = fltarr(5, deseasonal_lines)
for u = 0, deseasonal_lines - 1 do begin
	readf, deseason, deseasonal_arr0, lat0, lon0, year0, month0
	deseasonal_arr[ 0 , u ] = deseasonal_arr0
	deseasonal_arr[ 1 , u ] = lat0
	deseasonal_arr[ 2 , u ] = lon0
	deseasonal_arr[ 3 , u ] = year0
	deseasonal_arr[ 4 , u ] = month0
endfor

free_lun, deseason

;get the files
infile1= "/home/jakechung/IDL/myIDL/data_set.dat"
infile2= "/home/jakechung/IDL/myIDL/location_set.dat"

;get number of lines in the data file
n1 = file_lines(infile1)
n2 = file_lines(infile2)

print, 'The input has ', n1, '   data lines with', n2, '  locations'

;create an empty array to store the data that was read in lon, lat, avg, date and locations
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
	print, i
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

!P.Multi = 0
cgDisplay, 1380, 920
;Set up a world map for coordinates to be drawn on
MAP_SET, LIMIT=[ -90, -180, 90, 180 ], GRID=5, COLOR=!MYCT.BLUE, /CYL, /NOBORDER
MAP_CONTINENTS, /COUNTRIES, /COASTS, COLOR=!MYCT.BLACK

for j = 0, deseasonal_lines - 1 do begin
	if (j NE deseasonal_lines - 1) then b = j + 1 else b = j
		if ( deseasonal_arr[ 1 , j ] EQ deseasonal_arr[ 1 , b ] ) then begin
			a = where( deseasonal_arr[ 1 , * ] EQ deseasonal_arr[ 1 , j ] )
			temp_lat = deseasonal_arr[ 1 , a ]
			temp_lon = deseasonal_arr[ 2 , a ]
			print, temp_lon[ 0 ], temp_lat[ 0 ]
			DRAWDOTONMAP, temp_lon[ 5 ], temp_lat[ 5 ], 0, 0, ' ', !MYCT.red
			j = a[n_elements(a) - 1]
			a = 0
		endif else if ( deseasonal_arr[ 1 , j ] NE deseasonal_arr[ 1 , b ] ) then j = j + 1
endfor


;this loop pull out all the measured data of March, June, September and December
for j = 0, n2 - 1 do begin
	if ( j NE n2 - 1 ) then b = j + 1 ELSE b = j
	if ( lat1[ j ] EQ lat1[ b ] ) then begin
		;a is a index holder so the loop can skip the section where lat1[ j ] recurring
		a = where( lat1 EQ lat1[ j ] )
		if (n_elements(a) GT 30) then begin
		temp_lat = lat1[ a ]
		temp_lon = lon1[ a ]
		temp_location = location[ a ]
		DRAWDOTONMAP, temp_lon[ 0 ], temp_lat[ 0 ], 0, 0, ' ', !MYCT.green
		endif
		j = a[n_elements(a) - 1]
		a = 0
	endif else if ( lat1[ j ] NE lat1[ b ] ) then j = j + 1
endfor

DRAWDOTONMAP, -156.50, 71.16, 0, 0, ' ', !MYCT.black
DRAWDOTONMAP, -124, 45.50, 0, 0, ' ', !MYCT.black
DRAWDOTONMAP, -170.60, -14.10, 0, 0, ' ', !MYCT.black
DRAWDOTONMAP, -157.10, 21.08, 0, 0, ' ', !MYCT.black
DRAWDOTONMAP, 145, -42.00, 0, 0, ' ', !MYCT.black
DRAWDOTONMAP, 0, -90, 0, 0, ' ', !MYCT.black

cgLegend, Title= ['Simpson et al. data', 'NOAA data', 'OGI data'], PSym=[16, 16, 16], Color=['forest green', 'red', 'black'], location=[0.45,0.25],	$
	/Center_sym, /Box, /Background, BG_color= 'rose', Symsize= 2.0, Length=0.0, VSpace=2.0

;Make a JPEG for the graphs
READ, option, PROMPT= 'Do I need to save this graph as a jpeg so I can print it out later? Enter "1" for HELLYEAH or "0" for NOPE       '
name = strarr(1)
IF (option EQ 1) THEN BEGIN
	READ, name, PROMPT= 'Name of the JPEG file? (No need to include the jpeg extension, just name is fine)   '
	SCREEN2JPG, name
	print, 'JPEG has been saved to current directory'
ENDIF

end
