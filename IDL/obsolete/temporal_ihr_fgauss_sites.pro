PRO TEMPORAL_IHR_FGAUSS_SITES


;08/29/2017 - Jake Chung:
;	This program is built from temporal_ihr_rm_outliers_FiveSites and
;gaussian_filter_3sigma
;	It plots the temporal ihr of the Sites that represent NH and SH
;but using the Gaussian Fit on the detrend data to remove the outliers instead
;of the standard deviation of the detrend data.
;	It also has the numerical firn air data provided by Aydin et al.
;	First attempt on constructing a data structure



;read in the DAT file for OGI and UCI data set.
infile = '/home/excluded-from-backup/ethane/data/raw_data/Barrow_Samoa/IHR_w2sites_OGI_UCI.dat'
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

;08/17 revision: since the UCI and OGI data are together, cannot modify OGI data without modifying the 
;UCI data. So just keep this here but not use the arrays contain ogi_Samoa stations in the code.
;---------------------------------------------------------------------------------------------------
;read in the DAT file of NOAA data set for the Samoa and Barrow sites
infile = '/home/excluded-from-backup/ethane/data/raw_data/Barrow_Samoa/NOAA_Barrow_Samoa.dat'

;get the number of lines
n2 = file_lines( infile )

;create temporary arrays to store the values from the ASCII file
temp_year_noaa = fltarr( n2 )
temp_month_noaa = fltarr( n2 )
temp_day_noaa = fltarr( n2 )
temp_ratio_noaa = fltarr( n2 )
temp_lon_noaa = fltarr( n2 )
temp_lat_noaa = fltarr( n2 )

year = 0.0
month = 0.0
day = 0.0
ratio = 0.0
lon = 0.0
lat = 0.0

openr, iunit, infile, /get_lun

for i = 0, n2 - 1 do begin
	readf, iunit, year, month, day, ratio, lon, lat
	temp_year_noaa[ i ] = year
	temp_month_noaa[ i ] = month
	temp_day_noaa[ i ] = day
	temp_ratio_noaa[ i ] = ratio
	temp_lon_noaa[ i ] = lon
	temp_lat_noaa[ i ] = lat
endfor

free_lun, iunit

;remove the -999 values in the NOAA data 
a = where( temp_ratio_noaa lt 0 )
temp_year_noaa = RemoveRows( rotate(temp_year_noaa, 1), a )
temp_month_noaa = RemoveRows( rotate(temp_month_noaa, 1), a )
temp_day_noaa = RemoveRows( rotate(temp_day_noaa, 1), a )
temp_ratio_noaa = RemoveRows( rotate(temp_ratio_noaa, 1), a )
temp_lon_noaa = RemoveRows( rotate(temp_lon_noaa, 1), a )
temp_lat_noaa = RemoveRows( rotate(temp_lat_noaa, 1), a )


;08/17 revision: the next section read in the Cape Grim OGI data.<<<<<<<<<<<<<
infile = '/home/excluded-from-backup/ethane/data/raw_data/OGI/CapeGrim_OGI.dat'

n3 = file_lines(infile)
;create an array to read in the OGI data
;capegrim_ogi[ * , 0 ] contains the UNIX time
;capegrim_ogi[ 0 , * ] contains the mixing ratio
capegrim_ogi = fltarr( 2 , n3 )

time = 0.0
ratio = 0.0

openr, iunit, infile, /get_lun

for i = 0, n3 - 1 do begin
	readf, iunit, time, ratio
	capegrim_ogi[ 0 , i ] = time
	capegrim_ogi[ 1 , i ] = ratio
endfor

free_lun, iunit


;separate the data from the main arrays
;42 is OGI and 27 is UCI
;isolate the Samoa data from the OGI network
a = where( temp_network_ss_ogi eq 42 and temp_lat_ss_ogi eq -14.1 )
;isolate the Samoa data 
samoa_ogi = fltarr( 4, n_elements(a) )
samoa_ogi[ 0 , * ] = temp_ratio_ss_ogi[ a ]
samoa_ogi[ 1 , * ] = temp_lat_ss_ogi[ a ]
samoa_ogi[ 2 , * ] = temp_lon_ss_ogi[ a ]
samoa_ogi[ 3 , * ] = temp_time_ss_ogi[ a ]

;isolate the Barrow data from the OGI network
a = where( temp_network_ss_ogi eq 42 and temp_lat_ss_ogi eq 71.16 )
;isolate the Barrow data from the OGI network
barrow_ogi = fltarr( 4 , n_elements( a ) )
barrow_ogi[ 0 , * ] = temp_ratio_ss_ogi[ a ]
barrow_ogi[ 1 , * ] = temp_lat_ss_ogi[ a ] 
barrow_ogi[ 2 , * ] = temp_lon_ss_ogi[ a ]
barrow_ogi[ 3 , * ] = temp_time_ss_ogi[ a ]

;isolate the Samoa data from the UCI network
a = where( temp_network_ss_ogi eq 27 and temp_lat_ss_ogi eq -14.1 )
samoa_uci = fltarr( 4 , n_elements(a) )
samoa_uci[ 0 , * ] = temp_ratio_ss_ogi[ a ]
samoa_uci[ 1 , * ] = temp_lat_ss_ogi[ a ]
samoa_uci[ 2 , * ] = temp_lon_ss_ogi[ a ]
samoa_uci[ 3 , * ] = temp_time_ss_ogi[ a ]

;isolate the Barrow data from the UCI network
a = where( temp_network_ss_ogi eq 27 and temp_lat_ss_ogi eq 71.16 )
barrow_uci = fltarr( 4 , n_elements(a) )
barrow_uci[ 0 , * ] = temp_ratio_ss_ogi[ a ]
barrow_uci[ 1 , * ] = temp_lat_ss_ogi[ a ]
barrow_uci[ 2 , * ] = temp_lon_ss_ogi[ a ]
barrow_uci[ 3 , * ] = temp_time_ss_ogi[ a ]

;isolate the Samoa data from the NOAA data set
a = where( temp_lat_noaa eq -14.1 )
samoa_noaa = fltarr( 6 , n_elements(a) )
samoa_noaa[ 0 , * ] = temp_ratio_noaa[ a ]
samoa_noaa[ 1 , * ] = temp_lat_noaa[ a ]
samoa_noaa[ 2 , * ] = temp_lon_noaa[ a ]
samoa_noaa[ 3 , * ] = temp_year_noaa[ a ]
samoa_noaa[ 4 , * ] = temp_month_noaa[ a ]
samoa_noaa[ 5 , * ] = temp_day_noaa[ a ]

;isolate the Barrow data from the NOAA data set
a = where( temp_lat_noaa eq 71.16 )
barrow_noaa = fltarr( 6 , n_elements(a) )
barrow_noaa[ 0 , * ] = temp_ratio_noaa[ a ]
barrow_noaa[ 1 , * ] = temp_lat_noaa[ a ]
barrow_noaa[ 2 , * ] = temp_lon_noaa[ a ]
barrow_noaa[ 3 , * ] = temp_year_noaa[ a ]
barrow_noaa[ 4 , * ] = temp_month_noaa[ a ]
barrow_noaa[ 5 , * ] = temp_day_noaa[ a ]
;---------------------------------------------------------------------------------------------------
;convert Unix time in the UCI and OGI data to day and time arrays
samoa_ogi[ 3 , * ] = ( samoa_ogi[ 3 , * ]/3600 - 131496 )
samoa_ogi_time = tau2yymmdd(samoa_ogi[ 3 , * ], /GEOS1)

barrow_ogi[ 3 , * ] = ( barrow_ogi[ 3 , * ]/3600 - 131496 )
barrow_ogi_time = tau2yymmdd(barrow_ogi[ 3 , * ], /GEOS1)

samoa_uci[ 3 , * ] = ( samoa_uci[ 3 , * ]/3600 - 131496 )
samoa_uci_time = tau2yymmdd(samoa_uci[ 3 , * ], /GEOS1)

barrow_uci[ 3 , * ] = ( barrow_uci[ 3 , * ]/3600 - 131496 )
barrow_uci_time = tau2yymmdd(barrow_uci[ 3 , * ], /GEOS1)

capegrim_ogi[ 0 , * ] = ( capegrim_ogi[ 0 , * ]/3600 - 131496 );<<<<08/17 revision
capegrim_ogi_time =  tau2yymmdd(capegrim_ogi[ 0 , * ], /GEOS1)



;detrend the data
;find and average the data to get the mean of 12 months
detrend_capegrim_ogi = fltarr( n_elements( capegrim_ogi[ 0 , * ] ) )
b = fltarr( 13 )
holder = fltarr( n_elements( capegrim_ogi[ 0 , * ] ) )

for i = 1, 12 do begin
	a = where( capegrim_ogi_time.month eq i, count )
	if count gt 0 then $
	b[ i ] = mean( capegrim_ogi[ 1 , a ], /NAN)
	holder[ a ] = b[ i ]
endfor 

detrend_capegrim_ogi = capegrim_ogi[ 1 , * ] - holder

;deseasonalize the samoa OGI data
;find and average the data to get the mean of 12 months
detrend_samoa_ogi = fltarr( n_elements( samoa_ogi[ 0 , * ] ) )
b = fltarr( 13 )
holder = fltarr( n_elements( samoa_ogi[ 0 , * ] ) )

for i = 1, 12 do begin
	a = where( samoa_ogi_time.month eq i, count )
	if count gt 0 then $
	b[ i ] = mean( samoa_ogi[ 0 , a ], /NAN)
	holder[ a ] = b[ i ]
endfor 

detrend_samoa_ogi = samoa_ogi[ 0 , * ] - holder

;deseasonalize the Barrow OGI data
;find and then average the data to get the mean of 12 months
detrend_barrow_ogi = fltarr( n_elements( barrow_ogi[ 0 , * ] ) )
b = fltarr( 13 )
holder = fltarr( n_elements( barrow_ogi[ 0 , * ] ) )

for i = 1, 12 do begin
	a = where( barrow_ogi_time.month eq i, count )
	if count gt 0 then $
	b[ i ] = mean( barrow_ogi[ 0 , a ], /NAN)
	holder[ a ] = b[ i ]
endfor 

detrend_barrow_ogi = barrow_ogi[ 0 , * ] - holder

;deseasonalize the samoa UCI data
;find and then average the data to get the mean of 12 months
detrend_samoa_uci = fltarr( n_elements( samoa_uci[ 0 , * ] ) )
b = fltarr( 13 )
holder = fltarr( n_elements( samoa_uci[ 0 , * ] ) )

for i = 1, 12 do begin
	a = where( samoa_uci_time.month eq i, count )
	if count gt 0 then $
	b[ i ] = mean( samoa_uci[ 0 , a ], /NAN)
	holder[ a ] = b[ i ]

endfor 

detrend_samoa_uci = samoa_uci[ 0 , * ] - holder


;deseasonalize the Barrow UCI data
;find and then average the data to get the mean of 12 months
detrend_barrow_uci = fltarr( n_elements( barrow_uci[ 0 , * ] ) )
b = fltarr( 13 )
holder = fltarr( n_elements( barrow_uci[ 0 , * ] ) )

for i = 1, 12 do begin
	a = where( barrow_uci_time.month eq i, count )
	if count gt 0 then $
	b[ i ] = mean( barrow_uci[ 0 , a ], /NAN)
	holder[ a ] = b[ i ]
endfor 

detrend_barrow_uci = barrow_uci[ 0 , * ] - holder

;deseasonalize the Samoa NOAA data
;find and then average the data to get the mean of 12 months
detrend_samoa_noaa = fltarr( n_elements( samoa_noaa[ 0 , * ] ) )
b = fltarr( 13 )
holder = fltarr( n_elements( samoa_noaa[ 0 , * ] ) )

for i = 1, 12 do begin
	a = where( samoa_noaa[ 4 , * ] eq i, count )
	if count gt 0 then $
	b[ i ] = mean( samoa_noaa[ 0 , a ], /NAN)
	holder[ a ] = b[ i ]
endfor 

detrend_samoa_noaa = samoa_noaa[ 0 , * ] - holder

;deseasonalize the Barrow NOAA data
;find and then average the data to get the mean of 12 months
detrend_barrow_noaa = fltarr( n_elements( barrow_noaa[ 0 , * ] ) )
b = fltarr( 13 )
holder = fltarr( n_elements( barrow_noaa[ 0 , * ] ) )

for i = 1, 12 do begin
	a = where( barrow_noaa[ 4 , * ] eq i, count )
	if count gt 0 then $
	b[ i ] = mean( barrow_noaa[ 0 , a ], /NAN)
	holder[ a ] = b[ i ]
endfor 

detrend_barrow_noaa = barrow_noaa[ 0 , * ] - holder

;Gaussian fit to get the standard deviation

for i = 0 , 6 do begin
	delvar, binsize, coeff
	if i eq 0 then detrend_data = detrend_capegrim_ogi else if $
		i eq 1 then detrend_data = detrend_samoa_ogi else if $
		i eq 2 then detrend_data = detrend_barrow_ogi else if $
		i eq 3 then detrend_data = detrend_samoa_uci else if $
		i eq 4 then detrend_data = detrend_barrow_uci else if $
		i eq 5 then detrend_data = detrend_samoa_noaa else if $
		i eq 6 then detrend_data = detrend_barrow_noaa
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
	
	if i eq 0 then coeff_capegrim_ogi = coeff else if $
		i eq 1 then coeff_samoa_ogi = coeff else if $
		i eq 2 then coeff_barrow_ogi = coeff else if $
		i eq 3 then coeff_samoa_uci = coeff else if $
		i eq 4 then coeff_barrow_uci = coeff else if $
		i eq 5 then coeff_samoa_noaa = coeff else if $
		i eq 6 then coeff_barrow_noaa = coeff
		
endfor


;remove outliers from Gaussian Fit
;Cape Grim - OGI
outlier_idx_up = where(detrend_capegrim_ogi GT coeff_capegrim_ogi[2] * 3)
outlier_idx_low = where(detrend_capegrim_ogi LT - coeff_capegrim_ogi[2] * 3)
;concatnate the two arrays to perform RemoveRows function
outlier_idx = [outlier_idx_up , outlier_idx_low]
capegrim_ogi = RemoveRows(capegrim_ogi, outlier_idx)
detrend_capegrim_ogi = RemoveRows(detrend_capegrim_ogi, outlier_idx)


;Samoa - OGI
outlier_idx_up = where(detrend_samoa_ogi GT coeff_samoa_ogi[2] * 3)
outlier_idx_low = where(detrend_samoa_ogi LT - coeff_samoa_ogi[2] * 3)
;concatnate the two arrays to perform RemoveRows function
outlier_idx = [outlier_idx_up , outlier_idx_low]
samoa_ogi = RemoveRows(samoa_ogi, outlier_idx)
detrend_samoa_ogi = RemoveRows(detrend_samoa_ogi, outlier_idx)


;Barrow - OGI
outlier_idx_up = where(detrend_barrow_ogi GT coeff_barrow_ogi[2] * 3)
outlier_idx_low = where(detrend_barrow_ogi LT - coeff_barrow_ogi[2] * 3)
;concatnate the two arrays to perform RemoveRows function
outlier_idx = [outlier_idx_up , outlier_idx_low]
barrow_ogi = RemoveRows(barrow_ogi, outlier_idx)
detrend_barrow_ogi = RemoveRows(detrend_barrow_ogi, outlier_idx)


;Samoa - UCI
outlier_idx_up = where(detrend_samoa_uci GT coeff_samoa_uci[2] * 3)
outlier_idx_low = where(detrend_samoa_uci LT - coeff_samoa_uci[2] * 3)
;concatnate the two arrays to perform RemoveRows function
outlier_idx = [outlier_idx_up , outlier_idx_low]
samoa_uci = RemoveRows(samoa_uci, outlier_idx)
detrend_samoa_uci = RemoveRows(detrend_samoa_uci, outlier_idx)


;Barrow - UCI
outlier_idx_up = where(detrend_barrow_uci GT coeff_barrow_uci[2] * 3)
outlier_idx_low = where(detrend_barrow_uci LT - coeff_barrow_uci[2] * 3)
;concatnate the two arrays to perform RemoveRows function
outlier_idx = [outlier_idx_up , outlier_idx_low]
barrow_uci = RemoveRows(barrow_uci, outlier_idx)
detrend_barrow_uci = RemoveRows(detrend_barrow_uci, outlier_idx)


;Samoa - NOAA
outlier_idx_up = where(detrend_samoa_noaa GT coeff_samoa_noaa[2] * 3)
outlier_idx_low = where(detrend_samoa_noaa LT - coeff_samoa_noaa[2] * 3)
;concatnate the two arrays to perform RemoveRows function
outlier_idx = [outlier_idx_up , outlier_idx_low]
samoa_noaa = RemoveRows(samoa_noaa, outlier_idx)
detrend_samoa_noaa = RemoveRows(detrend_samoa_noaa, outlier_idx)


;Barrow - NOAA
outlier_idx_up = where(detrend_barrow_noaa GT coeff_barrow_noaa[2] * 3)
outlier_idx_low = where(detrend_barrow_noaa LT - coeff_barrow_noaa[2] * 3)
;concatnate the two arrays to perform RemoveRows function
outlier_idx = [outlier_idx_up , outlier_idx_low]
barrow_noaa = RemoveRows(barrow_noaa, outlier_idx)
detrend_barrow_noaa = RemoveRows(detrend_barrow_noaa, outlier_idx)


;recalculate the time array
;Reconvert Unix time in the UCI and OGI data to day and time arrays

samoa_ogi_time = tau2yymmdd(samoa_ogi[ 3 , * ], /GEOS1)

barrow_ogi_time = tau2yymmdd(barrow_ogi[ 3 , * ], /GEOS1)

samoa_uci_time = tau2yymmdd(samoa_uci[ 3 , * ], /GEOS1)

barrow_uci_time = tau2yymmdd(barrow_uci[ 3 , * ], /GEOS1)

capegrim_ogi_time =  tau2yymmdd(capegrim_ogi[ 0 , * ], /GEOS1)



;Re - Detrend the data
;find and average the data to get the mean of 12 months
detrend_capegrim_ogi = fltarr( n_elements( capegrim_ogi[ 0 , * ] ) )
b = fltarr( 13 )
holder = fltarr( n_elements( capegrim_ogi[ 0 , * ] ) )

for i = 1, 12 do begin
	a = where( capegrim_ogi_time.month eq i, count )
	if count gt 0 then $
	b[ i ] = mean( capegrim_ogi[ 1 , a ], /NAN)
	holder[ a ] = b[ i ]
endfor 

detrend_capegrim_ogi = capegrim_ogi[ 1 , * ] - holder

;deseasonalize the Barrow OGI data
;find and then average the data to get the mean of 12 months
detrend_barrow_ogi = fltarr( n_elements( barrow_ogi[ 0 , * ] ) )
b = fltarr( 13 )
holder = fltarr( n_elements( barrow_ogi[ 0 , * ] ) )

for i = 1, 12 do begin
	a = where( barrow_ogi_time.month eq i, count )
	if count gt 0 then $
	b[ i ] = mean( barrow_ogi[ 0 , a ], /NAN)
	holder[ a ] = b[ i ]
endfor 

detrend_barrow_ogi = barrow_ogi[ 0 , * ] - holder

;deseasonalize the samoa UCI data
;find and then average the data to get the mean of 12 months
detrend_samoa_uci = fltarr( n_elements( samoa_uci[ 0 , * ] ) )
b = fltarr( 13 )
holder = fltarr( n_elements( samoa_uci[ 0 , * ] ) )

for i = 1, 12 do begin
	a = where( samoa_uci_time.month eq i, count )
	if count gt 0 then $
	b[ i ] = mean( samoa_uci[ 0 , a ], /NAN)
	holder[ a ] = b[ i ]
endfor 

detrend_samoa_uci = samoa_uci[ 0 , * ] - holder

;deseasonalize the Barrow UCI data
;find and then average the data to get the mean of 12 months
detrend_barrow_uci = fltarr( n_elements( barrow_uci[ 0 , * ] ) )
b = fltarr( 13 )
holder = fltarr( n_elements( barrow_uci[ 0 , * ] ) )

for i = 1, 12 do begin
	a = where( barrow_uci_time.month eq i, count )
	if count gt 0 then $
	b[ i ] = mean( barrow_uci[ 0 , a ], /NAN)
	holder[ a ] = b[ i ]
endfor 

detrend_barrow_uci = barrow_uci[ 0 , * ] - holder

;deseasonalize the Samoa NOAA data
;find and then average the data to get the mean of 12 months
detrend_samoa_noaa = fltarr( n_elements( samoa_noaa[ 0 , * ] ) )
b = fltarr( 13 )
holder = fltarr( n_elements( samoa_noaa[ 0 , * ] ) )

for i = 1, 12 do begin
	a = where( samoa_noaa[ 4 , * ] eq i, count )
	if count gt 0 then $
	b[ i ] = mean( samoa_noaa[ 0 , a ], /NAN)
	holder[ a ] = b[ i ]
endfor 

detrend_samoa_noaa = samoa_noaa[ 0 , * ] - holder

;deseasonalize the Barrow NOAA data
;find and then average the data to get the mean of 12 months
detrend_barrow_noaa = fltarr( n_elements( barrow_noaa[ 0 , * ] ) )
b = fltarr( 13 )
holder = fltarr( n_elements( barrow_noaa[ 0 , * ] ) )

for i = 1, 12 do begin
	a = where( barrow_noaa[ 4 , * ] eq i, count )
	if count gt 0 then $
	b[ i ] = mean( barrow_noaa[ 0 , a ], /NAN)
	holder[ a ] = b[ i ]
endfor 

detrend_barrow_noaa = barrow_noaa[ 0 , * ] - holder

;After the outliers have been removed, putting all data sets and all sites
;in a data structure

;construct the NOAA network
barrow = {ratio : barrow_noaa[ 0 , * ], $
		lat : barrow_noaa[ 1 , * ], $
		lon : barrow_noaa[ 2 , * ], $
		year : barrow_noaa[ 3 , * ], $
		month : barrow_noaa[ 4 , * ], $
		residual : detrend_barrow_noaa }
		
samoa = {ratio : samoa_noaa[ 0 , * ], $
		lat : samoa_noaa[ 1 , * ], $
		lon : samoa_noaa[ 2 , * ], $
		year : samoa_noaa[ 3 , * ], $
		month : samoa_noaa[ 4 , * ], $
		residual : detrend_samoa_noaa }

noaa = {barrow : barrow, $
		samoa : samoa}
		
;construct the UCI network
barrow = {ratio : barrow_uci[ 0 , * ], $
		lat : barrow_uci[ 1 , * ], $
		lon : barrow_uci[ 2 , * ], $
		year : barrow_uci_time.year, $
		month : barrow_uci_time.month, $
		residual : detrend_barrow_uci }
		
samoa = {ratio : samoa_uci[ 0 , * ], $
		lat : samoa_uci[ 1 , * ], $
		lon : samoa_uci[ 2 , * ], $
		year : samoa_uci_time.year, $
		month : samoa_uci_time.month, $
		residual : detrend_samoa_uci }

uci = {barrow : barrow, $
		samoa : samoa}

;construc the OGI network
barrow = {ratio : barrow_ogi[ 0 , * ], $
		lat : barrow_ogi[ 1 , * ], $
		lon : barrow_ogi[ 2 , * ], $
		year : barrow_ogi_time.year, $
		month : barrow_ogi_time.month, $
		residual : detrend_barrow_ogi }
		
capegrim = {ratio : capegrim_ogi[ 1 , * ], $
		lat : -42, $
		lon : -145, $
		year : capegrim_ogi_time.year, $
		month : capegrim_ogi_time.month, $
		residual : detrend_capegrim_ogi }

ogi = {barrow : barrow, $
		capegrim : capegrim}
		
		
c2h6data = {noaa : noaa, $
			uci : uci, $
			ogi : ogi }

;calculations for IHR
;for OGI 
ihr_ogi = mean( c2h6data.ogi.barrow.ratio ) / mean( c2h6data.ogi.capegrim.ratio )

;for UCI, samoa station
;seperate the data into bins
a = 0 
a = where( samoa_uci_time.year eq 1996 or $
	samoa_uci_time.year eq 1997 or $
	samoa_uci_time.year eq 1998 or $
	samoa_uci_time.year eq 1999 or $
	samoa_uci_time.year eq 2000 )
samoa_bin01_uci = fltarr( 2 , n_elements( a ) )
samoa_bin01_uci[ 0 , * ] = samoa_uci[ 0 , a ]
samoa_bin01_uci[ 1 , * ] = samoa_uci[ 1 , a ]

a = 0
a = where( samoa_uci_time.year eq 2001 or $
	samoa_uci_time.year eq 2002 or $
	samoa_uci_time.year eq 2003 or $
	samoa_uci_time.year eq 2004 or $
	samoa_uci_time.year eq 2005 )
samoa_bin02_uci = fltarr( 2 , n_elements( a ) )
samoa_bin02_uci[ 0 , * ] = samoa_uci[ 0 , a ]
samoa_bin02_uci[ 1 , * ] = samoa_uci[ 1 , a ]

a = 0
a = where( samoa_uci_time.year eq 2006 or $
	samoa_uci_time.year eq 2007 or $
	samoa_uci_time.year eq 2008 or $
	samoa_uci_time.year eq 2009 )
samoa_bin03_uci = fltarr( 2 , n_elements( a ) )
samoa_bin03_uci[ 0 , * ] = samoa_uci[ 0 , a ]
samoa_bin03_uci[ 1 , * ] = samoa_uci[ 1 , a ]

;for UCI, Barrow station
;separate the data into bins
a = 0
a = where( barrow_uci_time.year eq 1996 or $
	barrow_uci_time.year eq 1997 or $
	barrow_uci_time.year eq 1998 or $
	barrow_uci_time.year eq 1999 or $
	barrow_uci_time.year eq 2000 )
barrow_bin01_uci = fltarr( 2 , n_elements( a ) )
barrow_bin01_uci[ 0 , * ] = barrow_uci[ 0 , a ]
barrow_bin01_uci[ 1 , * ] = barrow_uci[ 1 , a ]

a = 0
a = where( barrow_uci_time.year eq 2001 or $
	barrow_uci_time.year eq 2002 or $
	barrow_uci_time.year eq 2003 or $
	barrow_uci_time.year eq 2004 or $
	barrow_uci_time.year eq 2005 )
barrow_bin02_uci = fltarr( 2 , n_elements( a ) )
barrow_bin02_uci[ 0 , * ] = barrow_uci[ 0 , a ]
barrow_bin02_uci[ 1 , * ] = barrow_uci[ 1 , a ]

a = 0
a = where( barrow_uci_time.year eq 2006 or $
	barrow_uci_time.year eq 2007 or $
	barrow_uci_time.year eq 2008 or $
	barrow_uci_time.year eq 2009 )
barrow_bin03_uci = fltarr( 2 , n_elements( a ) )
barrow_bin03_uci[ 0 , * ] = barrow_uci[ 0 , a ]
barrow_bin03_uci[ 1 , * ] = barrow_uci[ 1 , a ]

;for NOAA, Samoa station
;separate the data into bins
a = 0
a = where( samoa_noaa[ 3 , * ] eq 2006 or $
	samoa_noaa[ 3 , * ] eq 2007 )
samoa_bin01_noaa = fltarr( 2 , n_elements( a ) )
samoa_bin01_noaa[ 0 , * ] = samoa_noaa[ 0 , a ]
samoa_bin01_noaa[ 1 , * ] = samoa_noaa[ 1 , a ]

a = 0
a = where( samoa_noaa[ 3 , * ] eq 2008 or $
	samoa_noaa[ 3 , * ] eq 2009 )
samoa_bin02_noaa = fltarr( 2 , n_elements( a ) )
samoa_bin02_noaa[ 0 , * ] = samoa_noaa[ 0 , a ]
samoa_bin02_noaa[ 1 , * ] = samoa_noaa[ 1 , a ]

a = 0
a = where( samoa_noaa[ 3 , * ] eq 2010 or $
	samoa_noaa[ 3 , * ] eq 2011 )
samoa_bin03_noaa = fltarr( 2 , n_elements( a ) )
samoa_bin03_noaa[ 0 , * ] = samoa_noaa[ 0 , a ]
samoa_bin03_noaa[ 1 , * ] = samoa_noaa[ 1 , a ]

a = 0
a = where( samoa_noaa[ 3 , * ] eq 2012 or $
	samoa_noaa[ 3 , * ] eq 2013 )
samoa_bin04_noaa = fltarr( 2 , n_elements( a ) )
samoa_bin04_noaa[ 0 , * ] = samoa_noaa[ 0 , a ]
samoa_bin04_noaa[ 1 , * ] = samoa_noaa[ 1 , a ]

a = 0
a = where( samoa_noaa[ 3 , * ] eq 2014 )
samoa_bin05_noaa = fltarr( 2 , n_elements( a ) )
samoa_bin05_noaa[ 0 , * ] = samoa_noaa[ 0 , a ]
samoa_bin05_noaa[ 1 , * ] = samoa_noaa[ 1 , a ]

;for NOAA, Barrow station
a = 0
a = where( barrow_noaa[ 3 , * ] eq 2006 or $
	barrow_noaa[ 3 , * ] eq 2007 )
barrow_bin01_noaa = fltarr( 2 , n_elements( a ) )
barrow_bin01_noaa[ 0 , * ] = barrow_noaa[ 0 , a ]
barrow_bin01_noaa[ 1 , * ] = barrow_noaa[ 1 , a ]

a = 0
a = where( barrow_noaa[ 3 , * ] eq 2008 or $
	barrow_noaa[ 3 , * ] eq 2009 )
barrow_bin02_noaa = fltarr( 2 , n_elements( a ) )
barrow_bin02_noaa[ 0 , * ] = barrow_noaa[ 0 , a ]
barrow_bin02_noaa[ 1 , * ] = barrow_noaa[ 1 , a ]

a = 0
a = where( barrow_noaa[ 3 , * ] eq 2010 or $
	barrow_noaa[ 3 , * ] eq 2011 )
barrow_bin03_noaa = fltarr( 2 , n_elements( a ) )
barrow_bin03_noaa[ 0 , * ] = barrow_noaa[ 0 , a ]
barrow_bin03_noaa[ 1 , * ] = barrow_noaa[ 1 , a ]

a = 0
a = where( barrow_noaa[ 3 , * ] eq 2012 or $
	barrow_noaa[ 3 , * ] eq 2013 )
barrow_bin04_noaa = fltarr( 2 , n_elements( a ) )
barrow_bin04_noaa[ 0 , * ] = barrow_noaa[ 0 , a ]
barrow_bin04_noaa[ 1 , * ] = barrow_noaa[ 1 , a ]

a = 0
a = where( barrow_noaa[ 3 , * ] eq 2014 )
barrow_bin05_noaa = fltarr( 2 , n_elements( a ) )
barrow_bin05_noaa[ 0 , * ] = barrow_noaa[ 0 , a ]
barrow_bin05_noaa[ 1 , * ] = barrow_noaa[ 1 , a ]
;calculate the ihr for the UCI network
ihr_uci = fltarr( 3 ) 
ihr_uci[ 0 ] = mean( barrow_bin01_uci[ 0 , * ] ) / mean( samoa_bin01_uci[ 0 , * ] ) 
ihr_uci[ 1 ] = mean( barrow_bin02_uci[ 0 , * ] ) / mean( samoa_bin02_uci[ 0 , * ] ) 
ihr_uci[ 2 ] = mean( barrow_bin03_uci[ 0 , * ] ) / mean( samoa_bin03_uci[ 0 , * ] ) 

;calculate the ihr for the NOAA network
ihr_noaa = fltarr( 5 )
ihr_noaa[ 0 ] = mean( barrow_bin01_noaa[ 0 , * ] ) / mean( samoa_bin01_noaa[ 0 , * ] )
ihr_noaa[ 1 ] = mean( barrow_bin02_noaa[ 0 , * ] ) / mean( samoa_bin02_noaa[ 0 , * ] )
ihr_noaa[ 2 ] = mean( barrow_bin03_noaa[ 0 , * ] ) / mean( samoa_bin03_noaa[ 0 , * ] )
ihr_noaa[ 3 ] = mean( barrow_bin04_noaa[ 0 , * ] ) / mean( samoa_bin04_noaa[ 0 , * ] )
ihr_noaa[ 4 ] = mean( barrow_bin05_noaa[ 0 , * ] ) / mean( samoa_bin05_noaa[ 0 , * ] )
;---------------------------------------------------------------------------------------------------
;separate the detrend data of the UCI for Samoa station
;for UCI, samoa station
a = 0 
a = where( samoa_uci_time.year eq 1996 or $
	samoa_uci_time.year eq 1997 or $
	samoa_uci_time.year eq 1998 or $
	samoa_uci_time.year eq 1999 or $
	samoa_uci_time.year eq 2000 )
samoa_detre_bin01_uci = fltarr( 1 , n_elements( a ) )
samoa_detre_bin01_uci[ 0 , * ] = detrend_samoa_uci[ 0 , a ]

a = 0
a = where( samoa_uci_time.year eq 2001 or $
	samoa_uci_time.year eq 2002 or $
	samoa_uci_time.year eq 2003 or $
	samoa_uci_time.year eq 2004 or $
	samoa_uci_time.year eq 2005 )
samoa_detre_bin02_uci = fltarr( 1 , n_elements( a ) )
samoa_detre_bin02_uci[ 0 , * ] = detrend_samoa_uci[ 0 , a ]

a = 0
a = where( samoa_uci_time.year eq 2006 or $
	samoa_uci_time.year eq 2007 or $
	samoa_uci_time.year eq 2008 or $
	samoa_uci_time.year eq 2009 )
samoa_detre_bin03_uci = fltarr( 1 , n_elements( a ) )
samoa_detre_bin03_uci[ 0 , * ] = detrend_samoa_uci[ 0 , a ]

;detrend data of the Barrow station from UCI network

a = 0
a = where( barrow_uci_time.year eq 1996 or $
	barrow_uci_time.year eq 1997 or $
	barrow_uci_time.year eq 1998 or $
	barrow_uci_time.year eq 1999 or $
	barrow_uci_time.year eq 2000 )
barrow_detre_bin01_uci = fltarr( 1 , n_elements( a ) )
barrow_detre_bin01_uci[ 0 , * ] = detrend_barrow_uci[ 0 , a ]

a = 0
a = where( barrow_uci_time.year eq 2001 or $
	barrow_uci_time.year eq 2002 or $
	barrow_uci_time.year eq 2003 or $
	barrow_uci_time.year eq 2004 or $
	barrow_uci_time.year eq 2005 )
barrow_detre_bin02_uci = fltarr( 1 , n_elements( a ) )
barrow_detre_bin02_uci[ 0 , * ] = detrend_barrow_uci[ 0 , a ]

a = 0
a = where( barrow_uci_time.year eq 2006 or $
	barrow_uci_time.year eq 2007 or $
	barrow_uci_time.year eq 2008 or $
	barrow_uci_time.year eq 2009 )
barrow_detre_bin03_uci = fltarr( 1 , n_elements( a ) )
barrow_detre_bin03_uci[ 0 , * ] = detrend_barrow_uci[ 0 , a ]

;separate the detrend data of the Samoa station from the NOAA network
a = 0
a = where( samoa_noaa[ 3 , * ] eq 2006 or $
	samoa_noaa[ 3 , * ] eq 2007 )
samoa_detre_bin01_noaa = fltarr( 1 , n_elements( a ) )
samoa_detre_bin01_noaa[ 0 , * ] = detrend_samoa_noaa[ 0 , a ]

a = 0
a = where( samoa_noaa[ 3 , * ] eq 2008 or $
	samoa_noaa[ 3 , * ] eq 2009 )
samoa_detre_bin02_noaa = fltarr( 1 , n_elements( a ) )
samoa_detre_bin02_noaa[ 0 , * ] = detrend_samoa_noaa[ 0 , a ]

a = 0
a = where( samoa_noaa[ 3 , * ] eq 2010 or $
	samoa_noaa[ 3 , * ] eq 2011 )
samoa_detre_bin03_noaa = fltarr( 1 , n_elements( a ) )
samoa_detre_bin03_noaa[ 0 , * ] = detrend_samoa_noaa[ 0 , a ]

a = 0
a = where( samoa_noaa[ 3 , * ] eq 2012 or $
	samoa_noaa[ 3 , * ] eq 2013 )
samoa_detre_bin04_noaa = fltarr( 1 , n_elements( a ) )
samoa_detre_bin04_noaa[ 0 , * ] = detrend_samoa_noaa[ 0 , a ]

a = 0
a = where( samoa_noaa[ 3 , * ] eq 2014 )
samoa_detre_bin05_noaa = fltarr( 1 , n_elements( a ) )
samoa_detre_bin05_noaa[ 0 , * ] = detrend_samoa_noaa[ 0 , a ]

;barrow station
a = 0
a = where( barrow_noaa[ 3 , * ] eq 2006 or $
	barrow_noaa[ 3 , * ] eq 2007 )
barrow_detre_bin01_noaa = fltarr( 1 , n_elements( a ) )
barrow_detre_bin01_noaa[ 0 , * ] = detrend_barrow_noaa[ 0 , a ]

a = 0
a = where( barrow_noaa[ 3 , * ] eq 2008 or $
	barrow_noaa[ 3 , * ] eq 2009 )
barrow_detre_bin02_noaa = fltarr( 1 , n_elements( a ) )
barrow_detre_bin02_noaa[ 0 , * ] = detrend_barrow_noaa[ 0 , a ]

a = 0
a = where( barrow_noaa[ 3 , * ] eq 2010 or $
	barrow_noaa[ 3 , * ] eq 2011 )
barrow_detre_bin03_noaa = fltarr( 1 , n_elements( a ) )
barrow_detre_bin03_noaa[ 0 , * ] = detrend_barrow_noaa[ 0 , a ]

a = 0
a = where( barrow_noaa[ 3 , * ] eq 2012 or $
	barrow_noaa[ 3 , * ] eq 2013 )
barrow_detre_bin04_noaa = fltarr( 1 , n_elements( a ) )
barrow_detre_bin04_noaa[ 0 , * ] = detrend_barrow_noaa[ 0 , a ]

a = 0
a = where( barrow_noaa[ 3 , * ] eq 2014 )
barrow_detre_bin05_noaa = fltarr( 1 , n_elements( a ) )
barrow_detre_bin05_noaa[ 0 , * ] = detrend_barrow_noaa[ 0 , a ]

;calculate the uncertainties for the OGI hemispheric means and ihr data
samoa_ogi_err = std_dev_fgauss( detrend_samoa_ogi ) / sqrt( n_elements( detrend_samoa_ogi ) )
barrow_ogi_err = std_dev_fgauss( detrend_barrow_ogi ) / sqrt( n_elements( detrend_barrow_ogi ) )
capegrim_ogi_err = std_dev_fgauss( detrend_capegrim_ogi ) / sqrt( n_elements( detrend_capegrim_ogi ) )
nor_ogi = mean( barrow_ogi[ 0 , * ] )
sou_ogi = mean( samoa_ogi[ 0 , * ] )
sou_ogi_w_capegrim = mean( capegrim_ogi[ 1 , * ] )
ihr_ogi_err = sqrt( (barrow_ogi_err/sou_ogi_w_capegrim)^2 + (nor_ogi*capegrim_ogi_err/(sou_ogi_w_capegrim^2))^2 )

;calculate the uncertainties for the UCI hemispheric means and ihr data
samoa_uci_err = fltarr( 3 )
barrow_uci_err = fltarr( 3 )
samoa_uci_err[ 0 ] = std_dev_fgauss( samoa_detre_bin01_uci ) / sqrt( n_elements( samoa_detre_bin01_uci ) )
samoa_uci_err[ 1 ] = std_dev_fgauss( samoa_detre_bin02_uci ) / sqrt( n_elements( samoa_detre_bin02_uci ) )
samoa_uci_err[ 2 ] = std_dev_fgauss( samoa_detre_bin03_uci ) / sqrt( n_elements( samoa_detre_bin03_uci ) )
barrow_uci_err[ 0 ] = std_dev_fgauss( barrow_detre_bin01_uci ) / sqrt( n_elements( barrow_detre_bin01_uci ) )
barrow_uci_err[ 1 ] = std_dev_fgauss( barrow_detre_bin02_uci ) / sqrt( n_elements( barrow_detre_bin02_uci ) )
barrow_uci_err[ 2 ] = std_dev_fgauss( barrow_detre_bin03_uci ) / sqrt( n_elements( barrow_detre_bin03_uci ) )

nor_uci = fltarr( 3 )
sou_uci = fltarr( 3 )
nor_uci[ 0 ] = mean( barrow_bin01_uci[ 0 , * ] )
nor_uci[ 1 ] = mean( barrow_bin02_uci[ 0 , * ] )
nor_uci[ 2 ] = mean( barrow_bin03_uci[ 0 , * ] )

sou_uci[ 0 ] = mean( samoa_bin01_uci[ 0 , * ] )
sou_uci[ 1 ] = mean( samoa_bin02_uci[ 0 , * ] )
sou_uci[ 2 ] = mean( samoa_bin03_uci[ 0 , * ] )

ihr_uci_err = sqrt( (barrow_uci_err/sou_uci)^2 + (nor_uci*samoa_uci_err/(sou_uci^2))^2 )

;calculate the uncertainties for the NOAA hemispheric means and ihr data
samoa_noaa_err = fltarr( 5 )
barrow_noaa_err = fltarr( 5 )
samoa_noaa_err[ 0 ] = std_dev_fgauss( samoa_detre_bin01_noaa ) / sqrt( n_elements( samoa_detre_bin01_noaa ) )
samoa_noaa_err[ 1 ] = std_dev_fgauss( samoa_detre_bin02_noaa ) / sqrt( n_elements( samoa_detre_bin02_noaa ) )
samoa_noaa_err[ 2 ] = std_dev_fgauss( samoa_detre_bin03_noaa ) / sqrt( n_elements( samoa_detre_bin03_noaa ) )
samoa_noaa_err[ 3 ] = std_dev_fgauss( samoa_detre_bin04_noaa ) / sqrt( n_elements( samoa_detre_bin04_noaa ) )
samoa_noaa_err[ 4 ] = std_dev_fgauss( samoa_detre_bin05_noaa ) / sqrt( n_elements( samoa_detre_bin05_noaa ) )

barrow_noaa_err[ 0 ] = std_dev_fgauss( barrow_detre_bin01_noaa ) / sqrt( n_elements( barrow_detre_bin01_noaa ) )
barrow_noaa_err[ 1 ] = std_dev_fgauss( barrow_detre_bin02_noaa ) / sqrt( n_elements( barrow_detre_bin02_noaa ) )
barrow_noaa_err[ 2 ] = std_dev_fgauss( barrow_detre_bin03_noaa ) / sqrt( n_elements( barrow_detre_bin03_noaa ) )
barrow_noaa_err[ 3 ] = std_dev_fgauss( barrow_detre_bin04_noaa ) / sqrt( n_elements( barrow_detre_bin04_noaa ) )
barrow_noaa_err[ 4 ] = std_dev_fgauss( barrow_detre_bin05_noaa ) / sqrt( n_elements( barrow_detre_bin05_noaa ) )

nor_noaa = fltarr( 5 )
sou_noaa = fltarr( 5 )
nor_noaa[ 0 ] = mean( barrow_bin01_noaa[ 0 , * ] )
nor_noaa[ 1 ] = mean( barrow_bin02_noaa[ 0 , * ] )
nor_noaa[ 2 ] = mean( barrow_bin03_noaa[ 0 , * ] )
nor_noaa[ 3 ] = mean( barrow_bin04_noaa[ 0 , * ] )
nor_noaa[ 4 ] = mean( barrow_bin05_noaa[ 0 , * ] )

sou_noaa[ 0 ] = mean( samoa_bin01_noaa[ 0 , * ] )
sou_noaa[ 1 ] = mean( samoa_bin02_noaa[ 0 , * ] )
sou_noaa[ 2 ] = mean( samoa_bin03_noaa[ 0 , * ] )
sou_noaa[ 3 ] = mean( samoa_bin04_noaa[ 0 , * ] )
sou_noaa[ 4 ] = mean( samoa_bin05_noaa[ 0 , * ] )

ihr_noaa_err = sqrt( (barrow_noaa_err/sou_noaa)^2 + (nor_noaa*samoa_noaa_err/(sou_noaa^2))^2 )
;------------------------------------------------------------------------------------------------------
;New section for the 08/16 revision to add in the Aydin et al firn air data

infile_firn_nor = "/home/excluded-from-backup/ethane/data/raw_data/Aydin_et_al_firn_air/Summit_Firn.dat"
infile_firn_sou = "/home/excluded-from-backup/ethane/data/raw_data/Aydin_et_al_firn_air/SouthPole_Firn.dat"

openr, lun_nor, infile_firn_nor, /get_lun
openr, lun_sou, infile_firn_sou, /get_lun

firn_nor = fltarr( 2 , file_lines(infile_firn_nor) )
firn_sou = fltarr( 2 , file_lines(infile_firn_sou) )

;Read in Northern Firn Air
for i = 0, file_lines(infile_firn_nor) - 1 do begin
	readf, lun_nor, time_year, mix_ratio
	firn_nor[ 0 , i ] = time_year
	firn_nor[ 1 , i ] = mix_ratio
endfor

;Read in Southern Firn Air
for i = 0, file_lines(infile_firn_sou) - 1 do begin
	readf, lun_sou, time_year, mix_ratio
	firn_sou[ 0 , i ] = time_year
	firn_sou[ 1 , i ] = mix_ratio
endfor

free_lun, lun_nor, lun_sou

;ihr_firn contains the interhemispheric ratio of the Aydin et al firn air
ihr_firn = firn_nor[ 1 , *] / firn_sou[ 1 , * ]



;------------------------------------------------------------------------------------------------------
;dealing with the simulated data

;sim_title = 'PSU emissions scaled to Xiao et al over 1996-2003'
filename1 = "/home/excluded-from-backup/data/C2H6/trac_avg.PSUSF_1981_2015.bpch"

filename2 = "/home/excluded-from-backup/data/C2H6/trac_avg.AydinSF_1981_2015.bpch"

;sim_title = 'Constant default base emissions' 
filename3 = "/home/excluded-from-backup/data/C2H6/trac_avg.spinup_GFED4_MAVG_1981_2015.198101010000"

;sim_title = 'Aydin et al. no normalization to Xiao et al.'
filename4 = "/home/excluded-from-backup/data/C2H6/trac_avg.AydinAbsSF_1981_2015.bpch"

;sim_title = 'Simpson et al. no normalization to Xiao et al.'
filename5 = "/home/excluded-from-backup/data/C2H6/trac_avg.SimpsonSF_1981_2015.bpch"

;sim_title = 'Unscaled PSU emission, MER3BB, MER18FF'
filename6 = "/home/excluded-from-backup/data/C2H6/trac_avg.PSU_MER3BB_MER18FF_1981_2015.bpch"



for x = 0, 5 do begin ;<<<<<< Beginning of a big loop to loop through all 6 emission scenarios.

print, 'Processing Emission scenario', x
if x eq 0 then tempfilename = filename1 else if $
	x eq 1 then tempfilename = filename2 else if $
	x eq 2 then tempfilename = filename3 else if $
	x eq 3 then tempfilename = filename4 else if $
	x eq 4 then tempfilename = filename5 else if $
	x eq 5 then tempfilename = filename6
	

;read in the file
ctm_get_data, datainfo, 'IJ-AVG-$', tracer=1, filename=tempfilename, tau0=tau0

;Get MODELINFO and GRIDINFO structures, xmid, ymid hold lon/lat centers
getmodelandgridinfo, datainfo[0], modelinfo, gridinfo

;nt is the number of data points in the simulation.
nt = n_elements(datainfo)
globalavg= fltarr(nt)
tempSimArr = fltarr( nt, 144, 91 )
;sim_yymmdd stores the entire time element of the input simulation
sim_yymmdd= tau2yymmdd(tau0, /GEOS1)

;this section call out the data blocks of the simulated data and store it in simArr
;simArr is a 3-D array with the following attributes:
;simArr[ *, 0, 0 ]: the time dimension
;simArr[ 0, *, 0 ]: longitudes
;simArr[ 0, 0, * ]: latitudes

for i = 0, nt - 1 do begin
	data= CTM_EXTRACT( *(datainfo[ i ].data), modelinfo= modelinfo, $
	gridinfo= gridinfo, lat= [-90, 90], lon= [-180, 180], alrange = [1, 2], $
	average= 4 )
	for k = 0, n_elements( data[ * , 0 ] ) - 1 do begin
		for j = 0, n_elements( data[ 0 , * ] ) - 1 do begin
			tempSimArr[ i , k , j ] = data[ k , j ]
		endfor
	endfor
endfor
	
;------------------------------------------------------------------------------------------
;this following is a special loop to average the simArr annually and then
;separate the data into 6 latitudinal bands
;creating an array to store the annual values
temp_global_sim = fltarr( 34 , 6 )
temp_annualSimArr = fltarr( 34, 144 , 91 )
temp_global_ihr = fltarr(34)
temp_global_nor = fltarr(34)
temp_global_sou = fltarr(34)
for i = 0, 33 do begin	;there are 34 years in the simulation, so from 0 to 33
	j = ( i ) * 12 
	j0 = findgen( 12 ) + j
	temp_annualSimArr[ i , * , * ] = mean( tempSimArr[ j0 , * , * ], 1 ) / 2 * 1000
	
	temp_global_nor[i] = temp_annualSimArr[ i , 9 , 81 ] ;9 and 81 are the indecies of the Barrow station
	temp_global_sou[i] = temp_annualSimArr[ i , 4 , 38 ] ;130 and 25 are index for Cape Grim.
endfor
temp_global_ihr = temp_global_nor/temp_global_sou



;remove the first year of the simulation, spin-up year
temp_global_nor[ 0 ] = !VALUES.F_NAN
temp_global_sou[ 0 ] = !VALUES.F_NAN
temp_global_ihr[ 0 ] = !VALUES.F_NAN



if x eq 0 then global_nor1 = temp_global_nor else if $
	x eq 1 then global_nor2 = temp_global_nor else if $
	x eq 2 then global_nor3 = temp_global_nor else if $
	x eq 3 then global_nor4 = temp_global_nor else if $
	x eq 4 then global_nor5 = temp_global_nor else if $
	x eq 5 then global_nor6 = temp_global_nor

if x eq 0 then global_sou1 = temp_global_sou else if $
	x eq 1 then global_sou2 = temp_global_sou else if $
	x eq 2 then global_sou3 = temp_global_sou else if $
	x eq 3 then global_sou4 = temp_global_sou else if $
	x eq 4 then global_sou5 = temp_global_sou else if $
	x eq 5 then global_sou6 = temp_global_sou
	
if x eq 0 then global_ihr1 = temp_global_ihr else if $
	x eq 1 then global_ihr2 = temp_global_ihr else if $
	x eq 2 then global_ihr3 = temp_global_ihr else if $
	x eq 3 then global_ihr4 = temp_global_ihr else if $
	x eq 4 then global_ihr5 = temp_global_ihr else if $
	x eq 5 then global_ihr6 = temp_global_ihr
	
	
endfor

;------------------------------------------------------------------------------------------------------
;chi-square test to determine goodness-of-fit
ihr_obs = fltarr( 9 )
ihr_sim = fltarr( 9 )
sigma_obs_ihr = fltarr( 9 )
chi_sq_var = fltarr( 6 )

for x = 0, 5 do begin

if x eq 0 then temp_global_ihr = global_ihr1 else if $
	x eq 1 then temp_global_ihr = global_ihr2 else if $
	x eq 2 then temp_global_ihr = global_ihr3 else if $
	x eq 3 then temp_global_ihr = global_ihr4 else if $
	x eq 4 then temp_global_ihr = global_ihr5 else if $
	x eq 5 then temp_global_ihr = global_ihr6
	
;OGI period, simulation
ihr_sim[ 0 ] = mean( temp_global_ihr[ 1 : 5] )

;UCI period, simulation
ihr_sim[ 1 ] = mean( temp_global_ihr[ 15 : 19 ] )
ihr_sim[ 2 ] = mean( temp_global_ihr[ 20 : 24 ] )
ihr_sim[ 3 ] = mean( temp_global_ihr[ 25 : 28 ] )

;NOAA period, simulation
ihr_sim[ 4 ] = mean( temp_global_ihr[ 25 : 26 ] )
ihr_sim[ 5 ] = mean( temp_global_ihr[ 27 : 28 ] )
ihr_sim[ 6 ] = mean( temp_global_ihr[ 29 : 30 ] )
ihr_sim[ 7 ] = mean( temp_global_ihr[ 31 : 32 ] )
ihr_sim[ 8 ] = mean( temp_global_ihr[ 33 ] )

if x eq 0 then ihr_sim1 = ihr_sim else if $
	x eq 1 then ihr_sim2 = ihr_sim else if $
	x eq 2 then ihr_sim3 = ihr_sim else if $
	x eq 3 then ihr_sim4 = ihr_sim else if $
	x eq 4 then ihr_sim5 = ihr_sim else if $
	x eq 5 then ihr_sim6 = ihr_sim
	
endfor 
;OGI obs
ihr_obs[ 0 ] = ihr_ogi

;UCI obs
ihr_obs[ 1 ] = ihr_uci[ 0 ]
ihr_obs[ 2 ] = ihr_uci[ 1 ]
ihr_obs[ 3 ] = ihr_uci[ 2 ]

;NOAA obs
ihr_obs[ 4 ] = ihr_noaa[ 0 ]
ihr_obs[ 5 ] = ihr_noaa[ 1 ]
ihr_obs[ 6 ] = ihr_noaa[ 2 ]
ihr_obs[ 7 ] = ihr_noaa[ 3 ]
ihr_obs[ 8 ] = ihr_noaa[ 4 ]

;Uncertainty of the obs data
;for OGI data
sigma_obs_ihr[ 0 ] = ihr_ogi_err

;for UCI data
sigma_obs_ihr[ 1 ] = ihr_uci_err[ 0 ]
sigma_obs_ihr[ 2 ] = ihr_uci_err[ 1 ]
sigma_obs_ihr[ 3 ] = ihr_uci_err[ 2 ]

;for NOAA data
sigma_obs_ihr[ 4 ] = ihr_noaa_err[ 0 ]
sigma_obs_ihr[ 5 ] = ihr_noaa_err[ 1 ]
sigma_obs_ihr[ 6 ] = ihr_noaa_err[ 2 ]
sigma_obs_ihr[ 7 ] = ihr_noaa_err[ 3 ]
sigma_obs_ihr[ 8 ] = ihr_noaa_err[ 4 ]

chi_sq_var[ 0 ] = total( ( ihr_obs - ihr_sim1 )^2/( sigma_obs_ihr ^ 2) )
chi_sq_var[ 1 ] = total( ( ihr_obs - ihr_sim2 )^2/( sigma_obs_ihr ^ 2) )
chi_sq_var[ 2 ] = total( ( ihr_obs - ihr_sim3 )^2/( sigma_obs_ihr ^ 2) )
chi_sq_var[ 3 ] = total( ( ihr_obs - ihr_sim4 )^2/( sigma_obs_ihr ^ 2) )
chi_sq_var[ 4 ] = total( ( ihr_obs - ihr_sim5 )^2/( sigma_obs_ihr ^ 2) )
chi_sq_var[ 5 ] = total( ( ihr_obs - ihr_sim6 )^2/( sigma_obs_ihr ^ 2) )





;------------------------------------------------------------------------------------------------------
;plotting procedure
ogi_time = [ 1985 ]
ss_time = [ 1998 , 2003 , 2008 ]
noaa_time = [2007 , 2009 , 2011 , 2013 , 2014]

all_time = findgen( 34 ) + 1981

;set up plots and all plotting procedures of this program



open_device, /ps, /color, file='temp.eps', margin=0.05, xsize = 10.0, ysize = 7.5
!x.thick=5
!y.thick=5
!p.font=0
!p.thick =3

multiplot, /default    ; resets multiplot settings
multiplot, [1,3], ygap=0.002, xgap=0;  sets up multiplot 

;Plot for the Northern Hemispheric Means
cgplot, ogi_time, nor_ogi, /nodata, ytitle = 'NHM (pptv)', xrange = [ 1982, 2015 ], charsize = 1.3, $
	symsize = 0.7, yrange = [ 1000, 2200 ]
cgplot, ogi_time, nor_ogi, /overplot, err_ylow = barrow_ogi_err, err_yhigh = barrow_ogi_err, psym = 4, color = 'black'

cgplot, ss_time, nor_uci, /overplot, err_ylow = barrow_uci_err, err_yhigh = barrow_uci_err, psym = 5, color = 'black'

cgplot, noaa_time, nor_noaa, /overplot, err_ylow = barrow_noaa_err, err_yhigh = barrow_noaa_err, psym = 9, color = 'black'

cgplot, all_time, global_nor1, /overplot, color = 'red', linestyle = 2

cgplot, all_time, global_nor2, /overplot, color = 'blue'

cgplot, all_time, global_nor3, /overplot, color = 'violet'

cgplot, all_time, global_nor4, /overplot, color = 'blue', linestyle = 2

cgplot, all_time, global_nor5, /overplot, color = 'forest green'

cgplot, all_time, global_nor6, /overplot, color = 'red'

cgplot, firn_nor[ 0 , * ], firn_nor[ 1 , * ], /overplot, color = 'black'

multiplot, /doyaxis, /doxaxis

;Plot for the Southern Hemispheric Means
cgplot, ogi_time, nor_ogi, /nodata, ytitle = 'SHM (pptv)', xrange = [ 1982, 2015 ], charsize = 1.3, $
	symsize = 0.7, yrange = [ 200, 400 ], /NoErase, XTickformat='(A1)'
	
cgplot, ogi_time, sou_ogi_w_capegrim, /overplot, err_ylow = samoa_ogi_err, err_yhigh = samoa_ogi_err, psym = 14, color = 'black'

cgplot, ss_time, sou_uci, /overplot, err_ylow = samoa_uci_err, err_yhigh = samoa_uci_err, psym = 17, color = 'black'

cgplot, noaa_time, sou_noaa, /overplot, err_ylow = samoa_noaa_err, err_yhigh = samoa_noaa_err, psym = 16, color = 'black'

cgplot, all_time, global_sou1, /overplot, color = 'red', linestyle = 2

cgplot, all_time, global_sou2, /overplot, color = 'blue'

cgplot, all_time, global_sou3, /overplot, color = 'violet'

cgplot, all_time, global_sou4, /overplot, color = 'blue', linestyle = 2

cgplot, all_time, global_sou5, /overplot, color = 'forest green'

cgplot, all_time, global_sou6, /overplot, color = 'red'

cgplot, firn_sou[ 0 , * ], firn_sou[ 1 , * ], /overplot, color = 'black'

multiplot, /doyaxis, /doxaxis

;AL_Legend, ['Scenario A', 'Scenario B', 'Scenario C', 'Scenario D', 'Scenario E', 'Scenario F'], $	
;	color = ['violet', 'forest green', 'blue', 'blue', 'red', 'red'], $
;	psym = [ -3, -3, -3, -3, -3, -3], linestyle = [ 0, 0, 2, 0, 2, 0 ], box = 0, $
;	position = [ 1983, 1200 ], charsize = 1.0
;AL_Legend, ['Barrow, NOAA', 'Barrow, UCI', 'Barrow, OGI', 'Samoa, NOAA', 'Samoa, UCI', 'Samoa, OGI'], $
;	psym = [ 9, 5, 4, 16, 17, 14 ], box=0, position = [ 1990, 1200 ], charsize = 1.0
	
	


;Interhemispheric Ratios
cgplot, ogi_time, ihr_ogi, /nodata, xrange = [ 1982, 2015 ], ytitle = 'Barrow/Samoa ratios', $
	yrange = [ 3.5 , 7.3 ], charsize = 1.2
cgplot, ogi_time, ihr_ogi, /overplot, err_ylow = ihr_ogi_err, err_yhigh = ihr_ogi_err, psym = 4, color = 'black'
cgplot, ss_time, ihr_uci, /overplot, err_ylow = ihr_uci_err, err_yhigh = ihr_uci_err, psym = 5, color = 'black'
cgplot, noaa_time, ihr_noaa, /overplot, err_ylow = ihr_noaa_err, err_yhigh = ihr_noaa_err, psym = 9, color = 'black'
cgplot, all_time, global_ihr1, /overplot, color = 'red', linestyle = 2
cgplot, all_time, global_ihr2, /overplot, color = 'blue'
cgplot, all_time, global_ihr3, /overplot, color = 'violet'
cgplot, all_time, global_ihr4, /overplot, color = 'blue', linestyle = 2
cgplot, all_time, global_ihr5, /overplot, color = 'forest green'
cgplot, all_time, global_ihr6, /overplot, color = 'red'
cgplot, firn_nor[ 0 , * ], ihr_firn, /overplot, color = 'black'

AL_Legend, ['Scenario A', 'Scenario B', 'Scenario C', $
	'Scenario D', 'Scenario E', 'Scenario F'], $	
	color = ['violet', 'forest green', 'blue', 'blue', 'red', 'red'], $
	psym = [ -3, -3, -3, -3, -3, -3], linestyle = [ 0, 0, 2, 0, 2, 0 ], box=0, $
	position = [ 1983, 5.45 ], charsize = 1.0
;AL_Legend, ['NOAA network', 'UCI network', 'OGI network'], $
;	psym = [ 9, 5, 4 ], box=0, position = [ 1988.3, 5.4 ], charsize = 1.0

close_device
spawn, 'gv temp.eps'

;------------------------------------------------------------------------------------------------------
;------------------------------------------------------------------------------------------------------
;------------------------------------------------------------------------------------------------------
;------------------------------------------------------------------------------------------------------
;------------------------------------------------------------------------------------------------------
;------------------------------------------------------------------------------------------------------
;------------------------------------------------------------------------------------------------------
;------------------------------------------------------------------------------------------------------
;------------------------------------------------------------------------------------------------------
;------------------------------------------------------------------------------------------------------
;------------------------------------------------------------------------------------------------------
;------------------------------------------------------------------------------------------------------
;------------------------------------------------------------------------------------------------------



end
