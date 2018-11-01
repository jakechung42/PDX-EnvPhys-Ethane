PRO UCI_84_96_MOD_LAT_LON

;09/19/2017 - Jake Chung
;	This program changes the sampling sites of the new UCI data from Blake 
;that covers the gap of the original UCI data and combine them
;into a single station represented by an coordinate. The coordinate grid follow the 
;GMAO 2° x 2.5° grid.
;	This program is also reading in the new UCI data with the lat and lon corrected.
;  GMAO 2° x 2.5° longitude edge to edge:
;-178.750     -176.250     -173.750     -171.250     -168.750     -166.250
;-163.750     -161.250     -158.750     -156.250     -153.750     -151.250
;-148.750     -146.250     -143.750     -141.250     -138.750     -136.250
;-133.750     -131.250     -128.750     -126.250     -123.750     -121.250
;-118.750     -116.250     -113.750     -111.250     -108.750     -106.250
;-103.750     -101.250     -98.7500     -96.2500     -93.7500     -91.2500
;-88.7500     -86.2500     -83.7500     -81.2500     -78.7500     -76.2500
;-73.7500     -71.2500     -68.7500     -66.2500     -63.7500     -61.2500
;-58.7500     -56.2500     -53.7500     -51.2500     -48.7500     -46.2500
;-43.7500     -41.2500     -38.7500     -36.2500     -33.7500     -31.2500
;-28.7500     -26.2500     -23.7500     -21.2500     -18.7500     -16.2500
;-13.7500     -11.2500     -8.75000     -6.25000     -3.75000     -1.25000
;1.25000      3.75000      6.25000      8.75000      11.2500      13.7500
;16.2500      18.7500      21.2500      23.7500      26.2500      28.7500
;31.2500      33.7500      36.2500      38.7500      41.2500      43.7500
;46.2500      48.7500      51.2500      53.7500      56.2500      58.7500
;61.2500      63.7500      66.2500      68.7500      71.2500      73.7500
;76.2500      78.7500      81.2500      83.7500      86.2500      88.7500
;91.2500      93.7500      96.2500      98.7500      101.250      103.750
;106.250      108.750      111.250      113.750      116.250      118.750
;121.250      123.750      126.250      128.750      131.250      133.750
;136.250      138.750      141.250      143.750      146.250      148.750
;151.250      153.750      156.250      158.750      161.250      163.750
;166.250      168.750      171.250      173.750      176.250      178.750

;  GMAO 2° x 2.5° latitude edge to edge:
;-91.0000     -89.0000     -87.0000     -85.0000     -83.0000     -81.0000
;-79.0000     -77.0000     -75.0000     -73.0000     -71.0000     -69.0000
;-67.0000     -65.0000     -63.0000     -61.0000     -59.0000     -57.0000
;-55.0000     -53.0000     -51.0000     -49.0000     -47.0000     -45.0000
;-43.0000     -41.0000     -39.0000     -37.0000     -35.0000     -33.0000
;-31.0000     -29.0000     -27.0000     -25.0000     -23.0000     -21.0000
;-19.0000     -17.0000     -15.0000     -13.0000     -11.0000     -9.00000
;-7.00000     -5.00000     -3.00000     -1.00000      1.00000      3.00000
;5.00000      7.00000      9.00000      11.0000      13.0000      15.0000
;17.0000      19.0000      21.0000      23.0000      25.0000      27.0000
;29.0000      31.0000      33.0000      35.0000      37.0000      39.0000
;41.0000      43.0000      45.0000      47.0000      49.0000      51.0000
;53.0000      55.0000      57.0000      59.0000      61.0000      63.0000
;65.0000      67.0000      69.0000      71.0000      73.0000      75.0000
;77.0000      79.0000      81.0000      83.0000      85.0000      87.0000
;89.0000      91.0000

compile_opt idl2

;Read UCI data

infile1= "/home/excluded-from-backup/ethane/data/raw_data/UCI/84_96_uci.dat"

n1_ss = file_lines(infile1)

unix_date_ss = fltarr(n1_ss)
lat1_ss = fltarr(n1_ss)
lon1_ss = fltarr(n1_ss)
avg_ss = fltarr(n1_ss)

;FLOAT input variable
date0= 0.0
lat0= 0.0
lon0= 0.0
avg0= 0.0

;open input files
openr, iunit1, infile1, /get_lun

;read float data of the data file
for i = 0, n1_ss - 1  do begin


	;read one line from the file
	readf, iunit1, date0, lat0, lon0, avg0
	;store the values
	unix_date_ss[ i ] = date0
	lat1_ss[ i ] = lat0
	lon1_ss[ i ] = lon0
	avg_ss[ i ] = avg0
endfor

;close input file
free_lun, iunit1

;create two arrays that contain the edge to edge lat and lon
edge_lat = findgen(92) * 2 - 91
edge_lon = findgen(146) * 2.5 - 181.25

;create two arrays that contain the center lat and lon which will replace the 
;lat and lon values in the UCI data set.
cen_lat = findgen(91) * 2 - 90
cen_lon = findgen(144) * 2.5 - 180

;change the first and last value of the cen_lat array to match  GMAO 2° x 2.5° lattitudes
cen_lat[ 0 ] = -89.5
cen_lat[ n_elements(cen_lat) - 1 ] = 89.5

for i = 0, n_elements(edge_lat) - 2 do begin
	for j = 0, n_elements(edge_lon) - 2 do begin
		x = where(lat1_ss GE edge_lat[ i ] AND $
				lat1_ss LE edge_lat[ i + 1 ] AND $
				lon1_ss GE edge_lon[ j ] AND $
				lon1_ss LE edge_lon[ j + 1], count)
		if (count gt 0) then begin
			lat1_ss[ x ] = ( edge_lat[ i ] + edge_lat[ i + 1 ] ) / 2
			lon1_ss[ x ] = ( edge_lon[ j ] + edge_lon[ j + 1 ] ) / 2
		endif
	endfor
endfor 

;The flaw in the algorithm is that the longitude 180 is not one of the center of the grid
;For this particular application, this flaw can be corrected by manually change the 180 
;longtitude to a corresponding center on the map, which is -180
x = where(lon1_ss eq 180)
lon1_ss[ x ] = -180

!P.Multi = 0
;cgDisplay, 1380, 920
;Set up a world map for coordinates to be drawn on
MAP_SET, LIMIT=[ -90, -180, 90, 180 ], /GRID, GLINESTYLE = 2, COLOR=!MYCT.BLUE, /CYL
MAP_CONTINENTS, /COASTS, COLOR=!MYCT.BLACK

;This section prints sites with more than 10 data points along with its lat and lon
for i = 0, n_elements(cen_lat) - 1 do begin
	for j = 0, n_elements(cen_lon) - 1 do begin
		x = where(lat1_ss eq cen_lat[ i ] and $
			lon1_ss eq cen_lon[ j ], count)
		if (n_elements(x) gt 10) then begin
			a = x[0]
			print, 'n: ', n_elements(x), ' lat: ', lat1_ss[ a ], ' lon: ', lon1_ss[ a ]
			DRAWDOTONMAP, lon1_ss[ a ], lat1_ss[ a ], 0, 0, ' ', !MYCT.red
		endif
	endfor
endfor 

;write to file
infile = "/home/excluded-from-backup/ethane/IDL/temp_file/modded_84_96_uci_data.dat"

openw, lun, infile, /get_lun

for i = 0, n_elements( lat1_ss ) - 1 do begin
	printf, lun, unix_date_ss[ i ], lat1_ss[ i ], lon1_ss[ i ], avg_ss[ i ], format = "(4F15.2)"
endfor

free_lun, lun

end
