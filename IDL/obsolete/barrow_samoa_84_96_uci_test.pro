PRO BARROW_SAMOA_84_96_UCI_TEST

;09/21 - Jake Chung
;	This program plots the 84_96 uci data set of the Barrow and Samoa sites 
;to determine the best way to average the data.

compile_opt idl2

infile = '/home/excluded-from-backup/ethane/data/raw_data/UCI/84_96_uci_BarrowSamoa.dat'

n = file_lines(infile)

;create an array to read in the 84-96 uci data
;uci_84_96[ 0 , * ] contains the UNIX time
;uci_84_96[ 1 , * ] contains the lattitudes
;uci_84_96[ 2 , * ] contains the longitudes
;uci_84_96[ 3 , * ] contains the mixing ratio

uci_84_96 = fltarr( 4 , n )

time = 0.0
ratio = 0.0
lat = 0.0
lon = 0.0

openr, iunit, infile, /get_lun

for i = 0, n - 1 do begin
	readf, iunit, time, lat, lon, ratio
	uci_84_96[ 0 , i ] = time
	uci_84_96[ 1 , i ] = lat
	uci_84_96[ 2 , i ] = lon
	uci_84_96[ 3 , i ] = ratio
endfor

free_lun, iunit

;separate the Barrow site and the Samoa site from the 84-96 uci data
x = where(uci_84_96[ 1 , * ] eq -14)
y = where(uci_84_96[ 1 , * ] eq 72)
samoa_new_uci = fltarr( 4 , n_elements(x) )
barrow_new_uci = fltarr( 4 , n_elements(y) )

samoa_new_uci[ 0 , * ] = uci_84_96[ 0 , x ]
samoa_new_uci[ 1 , * ] = uci_84_96[ 1 , x ]
samoa_new_uci[ 2 , * ] = uci_84_96[ 2 , x ]
samoa_new_uci[ 3 , * ] = uci_84_96[ 3 , x ]

barrow_new_uci[ 0 , * ] = uci_84_96[ 0 , y ]
barrow_new_uci[ 1 , * ] = uci_84_96[ 1 , y ]
barrow_new_uci[ 2 , * ] = uci_84_96[ 2 , y ]
barrow_new_uci[ 3 , * ] = uci_84_96[ 3 , y ]

x = sort(barrow_new_uci[ 0 , * ])
barrow_new_uci[ 0 , * ] = barrow_new_uci[ 0 , x ]
barrow_new_uci[ 1 , * ] = barrow_new_uci[ 1 , x ]
barrow_new_uci[ 2 , * ] = barrow_new_uci[ 2 , x ]
barrow_new_uci[ 3 , * ] = barrow_new_uci[ 3 , x ]

x = sort(samoa_new_uci[ 0 , * ])
samoa_new_uci[ 0 , * ] = samoa_new_uci[ 0 , x ]
samoa_new_uci[ 1 , * ] = samoa_new_uci[ 1 , x ]
samoa_new_uci[ 2 , * ] = samoa_new_uci[ 2 , x ]
samoa_new_uci[ 3 , * ] = samoa_new_uci[ 3 , x ]

;convert the Unix time to dates
barrow_new_uci[ 0 , * ] = ( barrow_new_uci[ 0 , * ]/3600 - 131496 )
barrow_new_uci_time = tau2yymmdd(barrow_new_uci[ 0 , * ], /GEOS1)

samoa_new_uci[ 0 , * ] = ( samoa_new_uci[ 0 , * ]/3600 - 131496 )
samoa_new_uci_time = tau2yymmdd(samoa_new_uci[ 0 , * ], /GEOS1)

;convert tau date to decimal dates
;1 year = 8765.76 hours
barrow_new_uci[ 0 , * ] = barrow_new_uci[ 0 , * ] / 8765.76 + 1985
samoa_new_uci[ 0 , * ] = samoa_new_uci[ 0 , * ] / 8765.76 + 1985

;deseasonalize the new Barrow data<<<09/21 revision
;find and then average the data to get the mean of 12 months
detrend_new_barrow_uci = fltarr( n_elements( barrow_new_uci[ 0 , * ] ) )
b = fltarr( 13 )
holder = fltarr( n_elements( barrow_new_uci[ 0 , * ] ) )

for i = 1, 12 do begin
	a = where( barrow_new_uci_time.month eq i, count )
	if count gt 0 then $
	b[ i ] = mean( barrow_new_uci[ 3 , a ], /NAN)
	holder[ a ] = b[ i ]
endfor 

detrend_new_barrow_uci = barrow_new_uci[ 3 , * ] - holder

;deseasonalize the new Samoa data<<<09/21 revision
;find and then average the data to get the mean of 12 months
detrend_new_samoa_uci = fltarr( n_elements( samoa_new_uci[ 0 , * ] ) )
b = fltarr( 13 )
holder = fltarr( n_elements( samoa_new_uci[ 0 , * ] ) )

for i = 1, 12 do begin
	a = where( samoa_new_uci_time.month eq i, count )
	if count gt 0 then $
	b[ i ] = mean( samoa_new_uci[ 3 , a ], /NAN)
	holder[ a ] = b[ i ]
endfor 

detrend_new_samoa_uci = samoa_new_uci[ 3 , * ] - holder

open_device, /ps, /color, file='temp.eps', margin=0.05, xsize = 10.0, ysize = 7.5
!x.thick=5
!y.thick=5
!p.font=0
!p.thick =3

!P.MULTI = [ 0 , 1 , 2 , 0 , 0 ]

cgplot, barrow_new_uci[ 0 , * ], barrow_new_uci[ 3 , * ], xticks = 12, xrange = [1984,1996]

cgplot, samoa_new_uci[ 0 , * ], samoa_new_uci[ 3 , * ], xticks = 11, xrange = [1985,1996]

close_device
spawn, 'gv temp.eps'
end
