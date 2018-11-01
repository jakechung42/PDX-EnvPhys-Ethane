PRO PLOT_NOAA_BARROW_SUMMIT

;Author: Jake Chung
;08/18/2017
;Descrition: This program plot the Barrow, Summit, South Pole and Samoa sites from the
;NOAA data set to validate with the longitudinal distribution

;08/24/17 - Jake Chung
;Compare the original data with the data with all outliers removed

;08/25/17 - Jake Chung
;Adding histogram plot option.

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

;create new arrays that specific to the sites
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


;calculate the sigma values of each site 
sigma_capegrim = stddev(detrend_capegrim) 
sigma_summit = stddev(detrend_summit) 
sigma_barrow = stddev(detrend_barrow)
sigma_southpole = stddev(detrend_southpole)
sigma_samoa = stddev(detrend_samoa)


;removing the outliers
help, capegrim
outlier_idx_up = where(detrend_capegrim GT sigma_capegrim * 3)
outlier_idx_low = where(detrend_capegrim LT - sigma_capegrim * 3)
;concatnate the two arrays to perform RemoveRows function
outlier_idx = [outlier_idx_up , outlier_idx_low]
capegrim = RemoveRows(capegrim, outlier_idx)
detrend_capegrim2 = RemoveRows(detrend_capegrim, outlier_idx)
help, capegrim


help, summit
outlier_idx_up = where(detrend_summit GT sigma_summit * 3)
outlier_idx_low = where(detrend_summit LT - sigma_summit * 3)
;concatnate the two arrays to perform RemoveRows function
outlier_idx = [outlier_idx_up , outlier_idx_low]
summit = RemoveRows(summit, outlier_idx)
detrend_summit2 = RemoveRows(detrend_summit, outlier_idx)
help, summit


help, barrow
outlier_idx_up = where(detrend_barrow GT sigma_barrow * 3)
outlier_idx_low = where(detrend_barrow LT - sigma_barrow * 3)
;concatnate the two arrays to perform RemoveRows function
outlier_idx = [outlier_idx_up , outlier_idx_low]
barrow = RemoveRows(barrow, outlier_idx)
detrend_barrow2 = RemoveRows(detrend_barrow, outlier_idx)
help, barrow


help, southpole
outlier_idx_up = where(detrend_southpole GT sigma_southpole * 3)
outlier_idx_low = where(detrend_southpole LT - sigma_southpole * 3)
;concatnate the two arrays to perform RemoveRows function
outlier_idx = [outlier_idx_up , outlier_idx_low]
southpole = RemoveRows(southpole, outlier_idx)
detrend_southpole2 = RemoveRows(detrend_southpole, outlier_idx)
help, southpole


help, samoa
outlier_idx_up = where(detrend_samoa GT sigma_samoa * 3)
outlier_idx_low = where(detrend_samoa LT - sigma_samoa * 3)
;concatnate the two arrays to perform RemoveRows function
outlier_idx = [outlier_idx_up , outlier_idx_low]
samoa = RemoveRows(samoa, outlier_idx)
detrend_samoa2 = RemoveRows(detrend_samoa, outlier_idx)
help, samoa

;recalculate the sigma values of each site 
sigma_capegrim2 = stddev(detrend_capegrim2) 
sigma_summit2 = stddev(detrend_summit2) 
sigma_barrow2 = stddev(detrend_barrow2)
sigma_southpole2 = stddev(detrend_southpole2)
sigma_samoa2 = stddev(detrend_samoa2)



open_device, /ps, /color, file='temp.eps', margin=0.05, xsize = 10.0, ysize = 7.5
!x.thick=5
!y.thick=5
!p.font=0
!p.thick =3

!p.multi = [ 0 , 1 , 2 ]
;multiplot, /default    ; resets multiplot settings
;multiplot, [1,2], ygap=0.002, xgap=0;  sets up multiplot 

;Tinker with histogam plot

cgHistoplot, detrend_southpole, /fillpoly, LOCATIONS=loc, BINSIZE=binsize1, HISTDATA=h

binCenters = loc + (binsize1 / 2.0)
yfit = GaussFit(binCenters, h, coeff, NTERMS=3)
cgPlot, binCenters, yfit, COLOR='dodger blue', THICK=3, /OVERPLOT
print, coeff
maxfit = String(coeff[0], FORMAT='(F0.2)')
centerfit = String(coeff[1], FORMAT='(F0.2)')
fwhm = String(2 * SQRT(2 * ALOG(2)) * coeff[2], FORMAT='(F0.2)')
gauss_sigma = String(coeff[2], FORMAT='(F0.2)')
sigma_southpole = String(sigma_southpole, FORMAT='(F0.2)')
cgText, 0.5, 0.85, /NORMAL, 'Maximum: ' + maxfit, COLOR='navy'
cgText, 0.5, 0.80, /NORMAL, 'Center: ' + centerfit, COLOR='navy'
cgText, 0.5, 0.75, /NORMAL, 'FWHM: ' + fwhm, COLOR='navy'
cgText, 0.5, 0.70, /NORMAL, 'Fit Sigma: ' + gauss_sigma, COLOR='navy'
cgText, 0.5, 0.65, /NORMAL, 'Data Sigma: ' + sigma_southpole, COLOR='navy'

print, h
;cgText, 0.5, 0.75, /NORMAL, 'FWHM: ' sigma_southpole, COLOR='navy'
;Plot the Northern Hemisphere
;cgplot, noaa_data[ 0 , summit_idx], noaa_data[ 3 , summit_idx], /nodata, ytitle = 'NH (pptv)', $
;	xrange = [ 2005 , 2016 ], yrange = [500 , 3500]
;cgplot, noaa_data[ 0 , summit_idx], noaa_data[ 3 , summit_idx], /overplot, color = 'red'
;cgplot, noaa_data[ 0 , barrow_idx], noaa_data[ 3 , barrow_idx], /overplot, color = 'forest green'

;cgplot, summit[ 0 , * ], summit[ 3 , * ], /overplot, color = 'red'
;cgplot, barrow[ 0 , * ], barrow[ 3 , * ], /overplot, color = 'forest green'

;multiplot, /doyaxis, /doxaxis

cgHistoplot, detrend_southpole2, /fillpoly, LOCATIONS=loc, BINSIZE=binsize2, HISTDATA=h

binCenters = loc + (binsize2 / 2.0)
yfit = GaussFit(binCenters, h, coeff, NTERMS=3)
cgPlot, binCenters, yfit, COLOR='dodger blue', THICK=3, /OVERPLOT
print, coeff
maxfit = String(coeff[0], FORMAT='(F0.2)')
centerfit = String(coeff[1], FORMAT='(F0.2)')
fwhm = String(2 * SQRT(2 * ALOG(2)) * coeff[2], FORMAT='(F0.2)')
gauss_sigma = String(coeff[2], FORMAT='(F0.2)')
sigma_southpole2 = String(sigma_southpole2, FORMAT='(F0.2)')
cgText, 0.5, 0.35, /NORMAL, 'Maximum: ' + maxfit, COLOR='navy'
cgText, 0.5, 0.30, /NORMAL, 'Center: ' + centerfit, COLOR='navy'
cgText, 0.5, 0.25, /NORMAL, 'FWHM: ' + fwhm, COLOR='navy'
cgText, 0.5, 0.20, /NORMAL, 'Fit Sigma: ' + gauss_sigma, COLOR='navy'
cgText, 0.5, 0.15, /NORMAL, 'Data Sigma: ' + sigma_southpole2, COLOR='navy'
;cgText, 0.5, 0.25, /NORMAL, 'FWHM: ' sigma_southpole, COLOR='navy'


;Plot the Southern Hemisphere
;cgplot, noaa_data[ 0 , summit_idx], noaa_data[ 3 , summit_idx], /nodata, ytitle = 'SH (pptv)', xtitle = 'Time', $
;	xrange = [ 2005 , 2016 ], yrange = [ 0 , 700]
;cgplot, noaa_data[ 0 , southpole_idx], noaa_data[ 3 , southpole_idx], /overplot, color = 'brown'
;cgplot, noaa_data[ 0 , samoa_idx], noaa_data[ 3 , samoa_idx], /overplot, color = 'blue'
;cgplot, noaa_data[ 0 , capegrim_idx], noaa_data[ 3 , capegrim_idx], /overplot, color = 'black'

;cgplot, southpole[ 0 , * ], southpole[ 3 , * ], /overplot, color = 'brown'
;cgplot, samoa[ 0 , * ], samoa[ 3 , * ], /overplot, color = 'blue'
;cgplot, capegrim[ 0 , * ], capegrim[ 3 , * ], /overplot, color = 'black'

;AL_Legend, ['Barrow', 'Summit', 'South Pole', 'Samoa', 'Cape Grim'], $	
;	color = ['forest green', 'red', 'brown', 'blue', 'black'], $
;	linestyle = [ 0, 0, 0, 0, 0], box=0, $
;	position = [ 2009, 1000 ], charsize = 0.8
close_device
spawn, 'gv temp.eps'

end

