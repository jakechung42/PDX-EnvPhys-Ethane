PRO TEMPORAL_CAPEGRIM_LATITUDE

;this program plots the Cape Grim (cg) latitude +- 2 degrees from all data

compile_opt idl2

;read in the 1984 1996 UCI data. uci_new is the data structure of the UCI data
uci_new = read_uci8496()
print, 'UCI 1984 - 1996 data: '
help, uci_new, /str

;read in NOAA data
noaa = read_noaa()
print, 'NOAA data: '
help, noaa, /str

;read in UCI data
uci = read_uci()
print, 'UCI data: '
help, uci, /str

;read in OGI data
ogi = read_ogi()
print, 'OGI data: '
help, ogi, /str


;prepare the data
;get Cape Grim (cg) data from the OGI and NOAA data
;cg data for NOAA is CGO
;for OGI, the coordinate of cg is longitude 145, latitude -42
;NOAA:
cg_noaa_idx = where(strmatch(noaa.site[*], 'CGO', /FOLD_CASE) EQ 1)
cg_noaa = fltarr(2, n_elements(cg_noaa_idx))
cg_noaa[0,*] = noaa.ratio[cg_noaa_idx]
cg_noaa[1,*] = noaa.year[cg_noaa_idx] + noaa.month[cg_noaa_idx]/12

;OGI:
cg_ogi_idx = where(ogi.lat eq -42 and ogi.lon eq 145)
cg_ogi = fltarr(2, n_elements(cg_ogi_idx))
cg_ogi[0,*] = ogi.ratio[cg_ogi_idx]
cg_ogi[1,*] = ogi.year[cg_ogi_idx] + ogi.month[cg_ogi_idx]/12

;process the UCI data, isolate the sites on latitude band -40 to -44
;this algorithm will treat all the data between latitude -40 and -44 as 1 site
;ignoring longitude
;update: now will also pull the longitudes data to differentiate between the sites

;change the up_bound and low_bound values to change the latitude band boundaries.
up_bound = -40
low_bound = -44
cg_band_uci_idx = where(uci.lat ge low_bound and uci.lat le up_bound)
cg_band_uci = fltarr(4, n_elements(cg_band_uci_idx))
cg_band_uci[0,*] = uci.ratio[cg_band_uci_idx]
cg_band_uci[1,*] = uci.year[cg_band_uci_idx] + uci.month[cg_band_uci_idx]/12
cg_band_uci[2,*] = uci.lat[cg_band_uci_idx]
cg_band_uci[3,*] = uci.lon[cg_band_uci_idx]

;process UCI new data (1984-1996)
cg_band_new_uci_idx = where(uci_new.lat ge low_bound and uci_new.lat le up_bound)
cg_band_new_uci = fltarr(4, n_elements(cg_band_new_uci_idx))
cg_band_new_uci[0,*] = uci_new.ratio[cg_band_new_uci_idx]
cg_band_new_uci[1,*] = uci_new.year[cg_band_new_uci_idx] + uci_new.month[cg_band_new_uci_idx]/12
cg_band_new_uci[2,*] = uci_new.lat[cg_band_new_uci_idx]
cg_band_new_uci[3,*] = uci_new.lon[cg_band_new_uci_idx]


;algorithm to categorize each site of the UCI data set
;change the variables' names to something quick and easy to type that is only used
;within this algorithm
ratio = cg_band_uci[0,*]
year = cg_band_uci[1,*]
lat = cg_band_uci[2,*]
lon = cg_band_uci[3,*]
size = n_elements(ratio)
uq_lat_uci = !NULL ;this variable stores the latitude value of the unique site of UCI data
uq_lon_uci = !NULL ;this varable stores the longitude value of the unique site of UCI data
uq_n_uci = !NULL ;this variable stores the number of data points of each unique site'

;this section finds the unique coordinates and how many time each coordinate
;repeats
for i = 0, size - 1 do begin
  ;this algorithm will not be able to detect the last site if the last site
  ;only has 1 data point so a check for uniqueness has to be done for the last site
  if (i eq size-1) then begin
    uq_lat_uci = [uq_lat_uci, lat[i]]
    uq_lon_uci = [uq_lon_uci, lon[i]]
    uq_n_uci = [uq_n_uci, 1]
    break
  endif
  ;does the next lat same as current lat?
  if (lat[i] eq lat[i+1]) then begin
    ;case where next lon is the same as current lon
    if (lon[i] eq lon[i+1]) then begin
      ;if both next lat-lon eq to current lat-lon, we have a unique site.
      x = where(lat eq lat[i] and lon eq lon[i], count)
      uq_lat_uci = [uq_lat_uci, lat[i]]
      uq_lon_uci = [uq_lon_uci, lon[i]]
      uq_n_uci = [uq_n_uci, count]
      i = i + count - 1
      if (i ge size-1 ) then break
    endif else begin
      ;case where next lon is not the same as current lon
      uq_lat_uci = [uq_lat_uci, lat[i]]
      uq_lon_uci = [uq_lon_uci, lon[i]]
      uq_n_uci = [uq_n_uci, 1]
    endelse
  endif else begin
    ;case where next lat is not the same as current lat
    uq_lat_uci = [uq_lat_uci, lat[i]]
    uq_lon_uci = [uq_lon_uci, lon[i]]
    uq_n_uci = [uq_n_uci, 1]
  endelse
endfor

print, 'Latitude band: ', low_bound, ' to ', up_bound
print, 'UCI data set:'
print, 'Latitude|Longitude|Number of data points of that coordinate'
for i = 0, n_elements(uq_n_uci) - 1 do begin
  print, uq_lat_uci[i], uq_lon_uci[i], uq_n_uci[i]
endfor
;need double precision type on uq_n_uci
uq_n_uci = float(uq_n_uci)


;same algorithm to check for the unique sites for the UCI 84 96 data set.
;change the variables' names to something quick and easy to type that is only used
;within this algorithm
ratio = cg_band_new_uci[0,*]
year = cg_band_new_uci[1,*]
lat = cg_band_new_uci[2,*]
lon = cg_band_new_uci[3,*]
size = n_elements(ratio)
uq_lat_new_uci = !NULL ;this variable stores the latitude value of the unique site of the 1984-1996 UCI data
uq_lon_new_uci = !NULL ;this varable stores the longitude value of the unique site of the 1984-1996 UCI data
uq_n_new_uci = !NULL ;this variable stores the number of data points of each unique site
for i = 0, size - 1 do begin
  ;this algorithm will not be able to detect the last site if the last site
  ;only has 1 data point so a check for uniqueness has to be done for the last site
  if (i eq size-1) then begin
    uq_lat_new_uci = [uq_lat_new_uci, lat[i]]
    uq_lon_new_uci = [uq_lon_new_uci, lon[i]]
    uq_n_new_uci = [uq_n_new_uci, 1]
    break
  endif
  ;does the next lat same as current lat?
  if (lat[i] eq lat[i+1]) then begin
    ;case where next lon is the same as current lon
    if (lon[i] eq lon[i+1]) then begin
      ;if both next lat-lon eq to current lat-lon, we have a unique site.
      x = where(lat eq lat[i] and lon eq lon[i], count)
      uq_lat_new_uci = [uq_lat_new_uci, lat[i]]
      uq_lon_new_uci = [uq_lon_new_uci, lon[i]]
      uq_n_new_uci = [uq_n_new_uci, count]
      i = i + count - 1
      if (i ge size-1 ) then break
    endif else begin
      ;case where next lon is not the same as current lon
      uq_lat_new_uci = [uq_lat_new_uci, lat[i]]
      uq_lon_new_uci = [uq_lon_new_uci, lon[i]]
      uq_n_new_uci = [uq_n_new_uci, 1]
    endelse
  endif else begin
    ;case where next lat is not the same as current lat
    uq_lat_new_uci = [uq_lat_new_uci, lat[i]]
    uq_lon_new_uci = [uq_lon_new_uci, lon[i]]
    uq_n_new_uci = [uq_n_new_uci, 1]
  endelse
endfor

print, 'Latitude band: ', low_bound, ' to ', up_bound
print, '1984 to 1996 UCI data set:'
print, 'Latitude|Longitude|Number of data points of that coordinate'
for i = 0, n_elements(uq_n_new_uci) - 1 do begin
  print, uq_lat_new_uci[i], uq_lon_new_uci[i], uq_n_new_uci[i]
endfor
;need double precision type on uq_n_uci
uq_n_new_uci = float(uq_n_new_uci)

;constructing a color table
color_tb1 = strarr(16)
color_tb1[0] = 'Pale Green'
color_tb1[1] = 'Aquamarine'
color_tb1[2] = 'Spring Green'
color_tb1[3] = 'Cyan'
color_tb1[4] = 'Turquoise'
color_tb1[5] = 'Sea Green'
color_tb1[6] = 'Forest Green'
color_tb1[7] = 'Green Yellow'
color_tb1[8] = 'Chartreuse'
color_tb1[9] = 'Lawn Green'
color_tb1[10] = 'Green'
color_tb1[11] = 'Lime Green'
color_tb1[12] = 'Olive Drab'
color_tb1[13] = 'Olive'
color_tb1[14] = 'Dark Green'
color_tb1[15] = 'Pale Goldenrod'

;constructing another color table
color_tb2 = strarr(23)
color_tb2[0] = 'Peru'
color_tb2[1] = 'Indian Red'
color_tb2[2] = 'Chocolate'
color_tb2[3] = 'Sienna'
color_tb2[4] = 'Dark Salmon'
color_tb2[5] = 'Salmon'
color_tb2[6] = 'Light Salmon'
color_tb2[7] = 'Orange'
color_tb2[8] = 'Coral'
color_tb2[9] = 'Light Coral'
color_tb2[10] = 'Firebrick'
color_tb2[11] = 'Brown'
color_tb2[12] = 'Hot Pink'
color_tb2[13] = 'Deep Pink'
color_tb2[14] = 'Magenta'
color_tb2[15] = 'Tomato'
color_tb2[16] = 'Orange Red'
color_tb2[17] = 'Red'
color_tb2[18] = 'Violet Red'
color_tb2[19] = 'Maroon'
color_tb2[20] = 'Thistle'
color_tb2[21] = 'Plum'
color_tb2[22] = 'Violet'


;plotting procedure
open_device, /ps, /color, /landscape, file='temp.eps', margin=0.01, xsize = 10.0, ysize = 7.5
!x.thick=1
!y.thick=1
!p.font=0
!p.thick =1

;set up plot
cgPlot, cg_noaa[1,*], cg_noaa[0,*], xrange = [1980, 2016], charsize = 1.2, $
  /noerase, /nodata, ytitle = 'Mixing Ratio (pptv)', xtitle = 'Years', $
  title = 'Time series Cape Grim NOAA & OGI, lat band -44 to -40 UCI', xgridstyle = 1.0, $
  ygridstyle = 1.0, xticklen = 1.0, yticklen = 1.0, xticks = 12, yticks = 10, yrange = [0,1000]

cgPlot, cg_noaa[1,*], cg_noaa[0,*], psym = 11, color = 'steelblue', $
  symsize = 0.9, /overplot

cgplot, cg_ogi[1,*], cg_ogi[0,*], psym = 9, color = 'steelblue', /overplot, symsize = 0.9

;uncomment this plotting procedure if plotting all the sites within the lattitude band
;cgPlot, cg_band_uci[1,*], cg_band_uci[0,*], psym = 5, color = 'blue', /overplot, symsize = 1.1

;algorithm to plot the sites within a latitude band in different colors
;plot the data in different symbol sizes
symsize = fltarr(n_elements(uq_n_uci))

;this is a simple function to find the symbol size for the number of data points.
symsize = (1.7-0.7)/(70-1)*(uq_n_uci-1)+0.7
for i = 0, n_elements(uq_n_uci) - 1 do begin
  x = where(cg_band_uci[2,*] eq uq_lat_uci[i] and cg_band_uci[3,*] eq uq_lon_uci[i], count)
  cgPlot, cg_band_uci[1,x], cg_band_uci[0,x], psym = 5, /overplot, symsize = symsize[i], $
    color = color_tb1[i+8] ;the constant that is added to the color_tb1 index is just to
    ;scale the colors to find a set that can be represented well on the graph. Can just simply
    ;take it out if need to.
endfor

;plotting procedure for plotting multiple different sites of the 84 96 data set
symsize = fltarr(n_elements(uq_n_new_uci))
symsize = (1.7-0.7)/(33-1)*(uq_n_new_uci-1)+0.7
for i = 0, n_elements(uq_n_new_uci) - 1 do begin
  x = where(cg_band_new_uci[2,*] eq uq_lat_new_uci[i] and cg_band_new_uci[3,*] eq uq_lon_new_uci[i], count)
  cgPlot, cg_band_new_uci[1,x], cg_band_new_uci[0,x], psym = 4, /overplot, symsize = symsize[i], $
    color = color_tb2[i+7]
endfor

;uncomment this if I want to plot the 84 96 data as one single site.
;cgPlot, cg_band_new_uci[1,*], cg_band_new_uci[0,*], psym = 4, color = 'red', /overplot, $
;  symsize = 1.1, xgridstyle = 1, ygridstyle = 1

cgLegend, title = ['OGI', 'UCI 84-96', 'UCI', 'NOAA'], psyms = [9, 4, 5, 11], /box, /background, $
  bg_color = 'rose', SymColors = ['steelblue', 'red', 'forest green', 'steelblue'], /center_sym, length = 0, $
  location = [0.3, 0.7], charsize = 0.9, symsize = 1.1

close_device
spawn, 'gv temp.eps'
end
