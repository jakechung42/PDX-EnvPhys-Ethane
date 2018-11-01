PRO GAUSSIAN_FILTER_3SIGMA

;Author: Jake Chung
;08/25/2017
;	This program is built from PLOT_NOAA_BARROW_SUMMIT. It filter out data
;that is more or less than 3 standard deviation from the residual data.
;	The standard deviation is calculated by using the chHistoplot and GaussFit
;routines on the residual.

COMPILE_OPT IDL2	

;>>>GET NOAA DATA<<<

infile1 = "/home/excluded-from-backup/ethane/data/raw_data/NOAA/NOAA_data.dat"
infile2 = "/home/excluded-from-backup/ethane/data/raw_data/NOAA/NOAA_locations.dat"

n1 = file_lines(infile1)
n2 = file_lines(infile2)

lat1= fltarr(n1)
lon1= fltarr(n1)
avg= fltarr(n1)
location= strarr(n2)
year= fltarr(n1)
month= fltarr(n1)
day= fltarr(n1)
hour= fltarr(n1)
alt= fltarr(n1)
time= fltarr(n1)

sample_site_code= ' '
sample_latitude= 0.0
sample_longitude= 0.0
analysis_value= 0.0
sample_year= 0.0
sample_month= 0.0
sample_day= 0.0
sample_hour= 0.0
sample_altitude= 0.0

openr, iunit1, infile1, /get_lun
openr, iunit2, infile2, /get_lun

for i = 0, n1 - 1 do begin
	;read one line from the file
	readf, iunit1, sample_year, sample_month, sample_day, sample_hour, analysis_value,	$
				sample_latitude, sample_longitude, sample_altitude
	;store the values
	year[i] = sample_year
	month[i] = sample_month
	day[i] = sample_day
	hour[i] = sample_hour
	avg[i] = analysis_value
	lat1[i] = sample_latitude
	lon1[i] = sample_longitude
	alt[i] = sample_altitude
endfor

;read float data of the location file
for i = 0, n2-1 do begin
	readf, iunit2, sample_site_code
	location[i] = sample_site_code
endfor

;close input file
free_lun, iunit1, iunit2

;combine the year, month, day, hour into a decimal year values
for i = 0, n1 - 1 do begin	
	time[i] = year[i] + hour[i]/8544 + day[i]/365 + month[i]/12
endfor

;Remove the -999 values from the C2H6 arrays
;Note, the RemoveRows function is from the Coyote Library.
;The RemoveRows function removes rows but the arrays are in collums,
;so need to use the rotate function so that the RemoveRows function
;would work

neg999 = where(avg LT 0, count)
location = RemoveRows(rotate(location, 1), neg999)
time = RemoveRows(rotate(time, 1), neg999)
lat1 = RemoveRows(rotate(lat1, 1), neg999)
lon1 = RemoveRows(rotate(lon1, 1), neg999)
avg = RemoveRows(rotate(avg, 1), neg999)
year = RemoveRows(rotate(year, 1), neg999)
month = RemoveRows(rotate(month, 1), neg999)
day = RemoveRows(rotate(day, 1), neg999)

;combine the NOAA data into one single array called noaa_data
noaa_data = fltarr( 5 , n_elements(location) )
location_noaa = fltarr( n_elements(location) )
location_noaa = location
noaa_data[ 0 , * ] = time 
noaa_data[ 1 , * ] = lat1
noaa_data[ 2 , * ] = lon1
noaa_data[ 3 , * ] = avg
noaa_data[ 4 , * ] = month

;call the Summit data
summit_idx = where( strmatch(location_noaa, 'SUM', /FOLD_CASE) EQ 1 )

;call the Barrow data
barrow_idx = where( strmatch(location_noaa, 'BRW', /FOLD_CASE) EQ 1 )

;call the South Pole data
southpole_idx = where( strmatch(location_noaa, 'SPO', /FOLD_CASE) EQ 1 )

;call the Samoa data
samoa_idx = where( strmatch(location_noaa, 'SMO', /FOLD_CASE) EQ 1 )

;call the Cape Grim data
capegrim_idx =  where( strmatch(location_noaa, 'CGO', /FOLD_CASE) EQ 1 )

;create new arrays for each site
capegrim = fltarr( 5 , n_elements(location_noaa[capegrim_idx]))
summit = fltarr( 5 , n_elements(location_noaa[summit_idx]))
barrow = fltarr( 5 , n_elements(location_noaa[barrow_idx]))
southpole = fltarr( 5 , n_elements(location_noaa[southpole_idx]))
samoa = fltarr( 5 , n_elements(location_noaa[samoa_idx]))

;combine the data into one single array
capegrim[ 0 , * ] = noaa_data[ 0 , capegrim_idx ]
capegrim[ 1 , * ] = noaa_data[ 1 , capegrim_idx ]
capegrim[ 2 , * ] = noaa_data[ 2 , capegrim_idx ]
capegrim[ 3 , * ] = noaa_data[ 3 , capegrim_idx ]
capegrim[ 4 , * ] = noaa_data[ 4 , capegrim_idx ]

summit[ 0 , * ] = noaa_data[ 0 , summit_idx ]
summit[ 1 , * ] = noaa_data[ 1 , summit_idx ]
summit[ 2 , * ] = noaa_data[ 2 , summit_idx ]
summit[ 3 , * ] = noaa_data[ 3 , summit_idx ]
summit[ 4 , * ] = noaa_data[ 4 , summit_idx ]

barrow[ 0 , * ] = noaa_data[ 0 , barrow_idx ]
barrow[ 1 , * ] = noaa_data[ 1 , barrow_idx ]
barrow[ 2 , * ] = noaa_data[ 2 , barrow_idx ]
barrow[ 3 , * ] = noaa_data[ 3 , barrow_idx ]
barrow[ 4 , * ] = noaa_data[ 4 , barrow_idx ]

southpole[ 0 , * ] = noaa_data[ 0 , southpole_idx ]
southpole[ 1 , * ] = noaa_data[ 1 , southpole_idx ]
southpole[ 2 , * ] = noaa_data[ 2 , southpole_idx ]
southpole[ 3 , * ] = noaa_data[ 3 , southpole_idx ]
southpole[ 4 , * ] = noaa_data[ 4 , southpole_idx ]

samoa[ 0 , * ] = noaa_data[ 0 , samoa_idx ]
samoa[ 1 , * ] = noaa_data[ 1 , samoa_idx ]
samoa[ 2 , * ] = noaa_data[ 2 , samoa_idx ]
samoa[ 3 , * ] = noaa_data[ 3 , samoa_idx ]
samoa[ 4 , * ] = noaa_data[ 4 , samoa_idx ]


;calculate the residual of the data 
detrend_capegrim = fltarr( n_elements( capegrim[ 0 , * ] ) )
b = fltarr( 13 )
holder = fltarr( n_elements( capegrim[ 0 , * ] ) )

for i = 1, 12 do begin
	a = where( capegrim[ 4 , * ] eq i, count )
	if count gt 0 then $
	b[ i ] = mean( capegrim[ 3 , a ], /NAN)
	holder[ a ] = b[ i ]
endfor 

detrend_capegrim = capegrim[ 3 , * ] - holder

detrend_summit = fltarr( n_elements( summit[ 0 , * ] ) )
b = fltarr( 13 )
holder = fltarr( n_elements( summit[ 0 , * ] ) )

for i = 1, 12 do begin
	a = where( summit[ 4 , * ] eq i, count )
	if count gt 0 then $
	b[ i ] = mean( summit[ 3 , a ], /NAN)
	holder[ a ] = b[ i ]
endfor 

detrend_summit = summit[ 3 , * ] - holder

detrend_barrow = fltarr( n_elements( barrow[ 0 , * ] ) )
b = fltarr( 13 )
holder = fltarr( n_elements( barrow[ 0 , * ] ) )

for i = 1, 12 do begin
	a = where( barrow[ 4 , * ] eq i, count )
	if count gt 0 then $
	b[ i ] = mean( barrow[ 3 , a ], /NAN)
	holder[ a ] = b[ i ]
endfor 

detrend_barrow = barrow[ 3 , * ] - holder

detrend_southpole = fltarr( n_elements( southpole[ 0 , * ] ) )
b = fltarr( 13 )
holder = fltarr( n_elements( southpole[ 0 , * ] ) )

for i = 1, 12 do begin
	a = where( southpole[ 4 , * ] eq i, count )
	if count gt 0 then $
	b[ i ] = mean( southpole[ 3 , a ], /NAN)
	holder[ a ] = b[ i ]
endfor 

detrend_southpole = southpole[ 3 , * ] - holder

detrend_samoa = fltarr( n_elements( samoa[ 0 , * ] ) )
b = fltarr( 13 )
holder = fltarr( n_elements( samoa[ 0 , * ] ) )

for i = 1, 12 do begin
	a = where( samoa[ 4 , * ] eq i, count )
	if count gt 0 then $
	b[ i ] = mean( samoa[ 3 , a ], /NAN)
	holder[ a ] = b[ i ]
endfor 

detrend_samoa = samoa[ 3 , * ] - holder



for i = 0 , 4 do begin
	delvar, binsize, coeff
	if i eq 0 then detrend_data = detrend_southpole else if $
		i eq 1 then detrend_data = detrend_samoa else if $
		i eq 2 then detrend_data = detrend_barrow else if $
		i eq 3 then detrend_data = detrend_capegrim else if $
		i eq 4 then detrend_data = detrend_summit
	;use the cgHistplot to calculate the Gaussian fit
	open_device, /ps, /color, file='temp.eps', margin=0.05, xsize = 10.0, ysize = 7.5
	!x.thick=5
	!y.thick=5
	!p.font=0
	!p.thick =3

	;Use the cgHistoplot to generate the data needed for GaussFit routine
	cgHistoplot, detrend_data, /fillpoly, LOCATIONS=loc, BINSIZE=binsize, HISTDATA=h
	binCenters = loc + (binsize / 2.0)
	yfit = GaussFit(binCenters, h, coeff, NTERMS=3)

	close_device
	
	if i eq 0 then coeff_southpole = coeff else if $
		i eq 1 then coeff_samoa = coeff else if $
		i eq 2 then coeff_barrow = coeff else if $
		i eq 3 then coeff_capegrim = coeff else if $
		i eq 4 then coeff_summit = coeff
		
endfor
print, 'Before filtering'
print, 'South Pole: ', coeff_southpole
print, 'Samoa: ', coeff_samoa
print, 'Barrow: ', coeff_barrow
print, 'Cape Grim: ', coeff_capegrim
print, 'Summit: ', coeff_summit
;removing the outliers

	
outlier_idx_up = where(detrend_southpole GT coeff_southpole[2] * 3)
outlier_idx_low = where(detrend_southpole LT - coeff_southpole[2] * 3)
;concatnate the two arrays to perform RemoveRows function
outlier_idx = [outlier_idx_up , outlier_idx_low]
southpole = RemoveRows(southpole, outlier_idx)
detrend_southpole = RemoveRows(detrend_southpole, outlier_idx)
	
outlier_idx_up = where(detrend_samoa GT coeff_samoa[2] * 3)
outlier_idx_low = where(detrend_samoa LT - coeff_samoa[2] * 3)
;concatnate the two arrays to perform RemoveRows function
outlier_idx = [outlier_idx_up , outlier_idx_low]
samoa = RemoveRows(samoa, outlier_idx)
detrend_samoa = RemoveRows(detrend_samoa, outlier_idx)

outlier_idx_up = where(detrend_barrow GT coeff_barrow[2] * 3)
outlier_idx_low = where(detrend_barrow LT - coeff_barrow[2] * 3)
;concatnate the two arrays to perform RemoveRows function
outlier_idx = [outlier_idx_up , outlier_idx_low]
barrow = RemoveRows(barrow, outlier_idx)
detrend_barrow = RemoveRows(detrend_barrow, outlier_idx)

outlier_idx_up = where(detrend_capegrim GT coeff_capegrim[2] * 3)
outlier_idx_low = where(detrend_capegrim LT - coeff_capegrim[2] * 3)
;concatnate the two arrays to perform RemoveRows function
outlier_idx = [outlier_idx_up , outlier_idx_low]
capegrim = RemoveRows(capegrim, outlier_idx)
detrend_capegrim = RemoveRows(detrend_capegrim, outlier_idx)

outlier_idx_up = where(detrend_summit GT coeff_summit[2] * 3)
outlier_idx_low = where(detrend_summit LT - coeff_summit[2] * 3)
;concatnate the two arrays to perform RemoveRows function
outlier_idx = [outlier_idx_up , outlier_idx_low]
summit = RemoveRows(summit, outlier_idx)
detrend_summit = RemoveRows(detrend_summit, outlier_idx)

;calculate the uncertainty of the five sites
for i = 0 , 4 do begin
	delvar, binsize, coeff
	if i eq 0 then detrend_data = detrend_southpole else if $
		i eq 1 then detrend_data = detrend_samoa else if $
		i eq 2 then detrend_data = detrend_barrow else if $
		i eq 3 then detrend_data = detrend_capegrim else if $
		i eq 4 then detrend_data = detrend_summit
	;use the cgHistplot to calculate the Gaussian fit
	open_device, /ps, /color, file='temp.eps', margin=0.05, xsize = 10.0, ysize = 7.5
	!x.thick=5
	!y.thick=5
	!p.font=0
	!p.thick =3

	;Use the cgHistoplot to generate the data needed for GaussFit routine
	cgHistoplot, detrend_data, /fillpoly, LOCATIONS=loc, BINSIZE=binsize, HISTDATA=h
	binCenters = loc + (binsize / 2.0)
	yfit = GaussFit(binCenters, h, coeff, NTERMS=3)

	close_device
	
	if i eq 0 then coeff_southpole = coeff else if $
		i eq 1 then coeff_samoa = coeff else if $
		i eq 2 then coeff_barrow = coeff else if $
		i eq 3 then coeff_capegrim = coeff else if $
		i eq 4 then coeff_summit = coeff
		
endfor
print, 'After filtering'
print, 'South Pole: ', coeff_southpole
print, 'Samoa: ', coeff_samoa
print, 'Barrow: ', coeff_barrow
print, 'Cape Grim: ', coeff_capegrim
print, 'Summit: ', coeff_summit
;print the uncertainty calculated from the Gaussian fit
print, 'From Gaussian fit'
print, 'South Pole: ', coeff_southpole[2] / sqrt(n_elements(detrend_southpole))
print, n_elements(detrend_southpole)
print, 'Summit: ', coeff_summit[2] / sqrt(n_elements(detrend_summit))
print, n_elements(detrend_summit)
print, 'Barrow: ', coeff_barrow[2] / sqrt(n_elements(detrend_barrow))
print, n_elements(detrend_barrow)
print, 'Cape Grim: ', coeff_capegrim[2] / sqrt(n_elements(detrend_capegrim))
print, n_elements(detrend_capegrim)
print, 'Samoa: ', coeff_samoa[2] / sqrt(n_elements(detrend_samoa))
print, n_elements(detrend_samoa)
print, 'From the filtered residual'
print, 'South Pole: ', stddev(detrend_southpole) / sqrt(n_elements(detrend_southpole))
print, 'Summit: ', stddev(detrend_summit) / sqrt(n_elements(detrend_summit))
print, 'Barrow: ', stddev(detrend_barrow) / sqrt(n_elements(detrend_barrow))
print, 'Cape Grim: ', stddev(detrend_capegrim) / sqrt(n_elements(detrend_capegrim))
print, 'Samoa: ', stddev(detrend_samoa) / sqrt(n_elements(detrend_samoa))


open_device, /ps, /color, file='temp.eps', margin=0.05, xsize = 10.0, ysize = 7.5
!x.thick=5
!y.thick=5
!p.font=0
!p.thick =3

multiplot, /default    ; resets multiplot settings
multiplot, [1,2], ygap=0.002, xgap=0;  sets up multiplot 

;Plot the Northern Hemisphere
cgplot, noaa_data[ 0 , summit_idx], noaa_data[ 3 , summit_idx], /nodata, ytitle = 'NH (pptv)', $
	xrange = [ 2005 , 2016 ], yrange = [500 , 3500]

cgplot, summit[ 0 , * ], summit[ 3 , * ], /overplot, color = 'red'
cgplot, barrow[ 0 , * ], barrow[ 3 , * ], /overplot, color = 'forest green'

multiplot, /doyaxis, /doxaxis

;Plot the Southern Hemisphere
cgplot, noaa_data[ 0 , summit_idx], noaa_data[ 3 , summit_idx], /nodata, ytitle = 'SH (pptv)', xtitle = 'Time', $
	xrange = [ 2005 , 2016 ], yrange = [ 0 , 700]
cgplot, southpole[ 0 , * ], southpole[ 3 , * ], /overplot, color = 'brown'
cgplot, samoa[ 0 , * ], samoa[ 3 , * ], /overplot, color = 'blue'
cgplot, capegrim[ 0 , * ], capegrim[ 3 , * ], /overplot, color = 'black'

AL_Legend, ['Barrow', 'Summit', 'South Pole', 'Samoa', 'Cape Grim'], $	
	color = ['forest green', 'red', 'brown', 'blue', 'black'], $
	linestyle = [ 0, 0, 0, 0, 0], box=0, $
	position = [ 2009, 600 ], charsize = 0.9
	
close_device
;spawn, 'gv temp.eps'<<<< can be uncommented to produce plot.



end

