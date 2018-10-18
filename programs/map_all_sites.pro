PRO MAP_ALL_SITES

  ;this programThis program will plot the locations of all the sites from the 3 data networks and the new UCI on to a world map.
  ;UCI data and UCI_new (1984-1996) data uses the modified coordinates.
  ;The modified coordinates are generated from the programs uci_mod_lat_lon.pro and uci_84_96_mod_lat_lon.pro.
  ;The coordinates of OGI data is input manually.

COMPILE_OPT IDL2


;get NOAA data
noaa = read_noaa()
print, 'NOAA data: '
help, noaa, /str
lat_noaa = noaa.lat
lon_noaa = noaa.lon


;get UCI data
uci = read_uci()
print, 'UCI data: '
help, uci, /str
lat_uci = uci.lat
lon_uci = uci.lon


;get UCI 1984 - 1996 data
uci_new = read_uci8496()
print, 'UCI 1984 - 1996 data: '
help, uci_new, /str
lat_uci8496 = uci_new.lat
lon_uci8496 = uci_new.lon


;get OGI data
ogi = read_ogi()
print, 'OGI data: '
help, ogi, /str


;algorithm to filter the UCI and UCI 84-96 data
;display how many data points for each coordinate and write to a text file
;generate the text file
infile = '/home/excluded-from-backup/ethane/IDL/temp_file/avail_coordn_all_network.dat'

openw, lun, infile, /get_lun

;output data for the UCI new 1984 - 1996
printf, lun, '  Available 1984-1996 UCI data: '
printf, lun, 'Latitudes | Longitudes | Number of data points'
cen_lat = findgen(91) * 2 - 90
cen_lon = findgen(144) * 2.5 - 180
cen_lon[0] = 180
for i = 0, n_elements(cen_lat) - 1 do begin
  for j = 0, n_elements(cen_lon) - 1 do begin
    a = where(lat_uci8496 eq cen_lat[i] AND $
      lon_uci8496 eq cen_lon[j], count)
    if count ne 0 then begin
      ;output the coordinates and number of data points to the text file
      printf, lun, cen_lat[i], cen_lon[j], count
    endif
  endfor
endfor

;output data for the UCI network
printf, lun, '  Available UCI data: '
printf, lun, 'Latitudes | Longitudes | Number of data points'
cen_lat = findgen(91) * 2 - 90
cen_lon = findgen(144) * 2.5 - 180
cen_lon[0] = 180
for i = 0, n_elements(cen_lat) - 1 do begin
  for j = 0, n_elements(cen_lon) - 1 do begin
    a = where(lat_uci eq cen_lat[i] AND $
      lon_uci eq cen_lon[j], count)
    if count ne 0 then begin
      ;output the coordinates and number of data points to the text file
      printf, lun, cen_lat[i], cen_lon[j], count
    endif
  endfor
endfor

;output data for the NOAA network
;the algorithm for outputting the NOAA data is different from the UCI data

;read in a file that contains the NOAA sites' IDs
infile2= "/home/excluded-from-backup/ethane/data/raw_data/NOAA/NOAA_stations_IDs.dat"
openr, id_lun, infile2, /get_lun

;noaa_id contains the content of NOAA_stations_IDs.dat
noaa_id = strarr(file_lines(infile2))
sample_site_code = ' '
for i = 0, file_lines(infile2)-1 do begin
	readf, id_lun, sample_site_code
	noaa_id[i] = sample_site_code
endfor
free_lun, id_lun

;print out the coordinates of all the sites to the text file
printf, lun, '  Available NOAA data: '
printf, lun, 'Latitudes | Longitudes | Number of data points'
size_noaa = n_elements(noaa.ratio)
for i = 0, n_elements(noaa_id) - 1 do begin
  z = where(strmatch(noaa.site[*], noaa_id[i], /FOLD_CASE) EQ 1) ;<<<calling out each station
  x = z[0] ;<<<< get the first index to call out values
  printf, lun, noaa.lat[x], noaa.lon[x], n_elements(z)
endfor


;output the OGI coordinates and number of data points to the text file
lat_ogi = [71.16, 45.5, -14.1, 21.08, -42, -90]
lon_ogi = [-156.5, -124, -170.6, -157.1, 145, 0]

printf, lun, '  Available OGI data: '
printf, lun, 'Latitudes | Longitudes | Number of data points'
for i = 0, 5 do begin
  a = []
  a = where(ogi.lat eq lat_ogi[i] and ogi.lon eq lon_ogi[i], count)
  printf, lun, lat_ogi[i], lon_ogi[i], count
endfor


free_lun, lun
;specify coordinates of the OGI sites
lat_ogi = [71.16, 45.5, -14.1, 21.08, -42, -90]
lon_ogi = [-156.5, -124, -170.6, -157.1, 145, 0]

;plotting procedure:
;there are 2 plotting functions:
;change plot_optn to 1 to plot all networks on one map
;change plot_optn to 2 to plot NOAA/OGI on one map and UCI/UCI8496 on another map
plot_optn = 3
print, 'change plot_optn to 1 to plot all networks on one map'
print, 'change plot_optn to 2 to plot NOAA/OGI on one map and UCI/UCI8496 on another map'
print, 'plot option: ', plot_optn

if (plot_optn eq 1) then begin
  ;setup the Coyote Graphics display
  cgDisplay, 1584, 864, Title='Map Plot'
  ;set up the map projection
  xrange = [-180,180] ;longitudes
  yrange = [-90,90] ;latitudes
  center_lon = 0
  map = Obj_New('cgMap', 'Equirectangular', Ellipsoid=19, $
     XRange=xrange, YRange=yrange, /LatLon_Ranges, CENTER_LON=center_lon, $
     Position=[0.1, 0.1, 0.9, 0.8], Limit=[-90, -180, 90, 180])
  map -> Draw

  ;plot the OGI sites
  cgPlot, lon_ogi, lat_ogi, color = 'red', psym = 2, symsize = 1.5, Position=[0.1, 0.1, 0.9, 0.8], $
    xstyle = 4, ystyle = 4, xrange = [-180, 180], yrange = [-90, 90]
  cgPlot, lon_noaa, lat_noaa, color = 'forest green', psym = 7, symsize = 1.5, Position=[0.1, 0.1, 0.9, 0.8], $
    xstyle = 4, ystyle = 4, xrange = [-180, 180], yrange = [-90, 90], /overplot
  cgPlot, lon_uci, lat_uci, color = 'blue', psym = 5, symsize = 1.5, Position=[0.1, 0.1, 0.9, 0.8], $
    xstyle = 4, ystyle = 4, xrange = [-180, 180], yrange = [-90, 90], /overplot
  cgPlot, lon_uci8496, lat_uci8496, color = 'violet', psym = 9, symsize = 1.5, Position=[0.1, 0.1, 0.9, 0.8], $
    xstyle = 4, ystyle = 4, xrange = [-180, 180], yrange = [-90, 90]   , /overplot
  ;add map annotations
  cgMap_Grid, Map=map, /Box, charsize = 2
  cgMap_Continents, Map=map, Color='black'
  ;add map legend
  cgLegend, SymColors = ['red', 'forestgreen', 'blue', 'violet'], PSyms = [2, 7, 5, 9], $
    Symsize = 2, Location = [0.7, 0.7], titles = ['OGI', 'NOAA', 'UCI', 'UCI 1984-1996'], $
    /Box, /Background, BG_Color = 'rose', /center_sym, length = 0
endif
if (plot_optn eq 2) then begin
;----plotting the OGI and NOAA data----
;setup the Coyote Graphics display
cgDisplay, 990, 540, Title='Map Plot 1'
;set up the map projection
xrange = [-180,180] ;longitudes
yrange = [-90,90] ;latitudes
center_lon = 0
map = Obj_New('cgMap', 'Equirectangular', Ellipsoid=19, $
   XRange=xrange, YRange=yrange, /LatLon_Ranges, CENTER_LON=center_lon, $
   Position=[0.1, 0.1, 0.9, 0.8], Limit=[-90, -180, 90, 180])
map -> Draw
;plot the OGI coordinates
cgPlot, lon_ogi, lat_ogi, color = 'red', psym = 2, symsize = 1.2, Position=[0.1, 0.1, 0.9, 0.8], $
  xstyle = 4, ystyle = 4, xrange = [-180, 180], yrange = [-90, 90]
;plot the NOAA coordinates
cgPlot, lon_noaa, lat_noaa, color = 'forestgreen', psym = 7, symsize = 1.2, Position=[0.1, 0.1, 0.9, 0.8], $
  /overplot
;add map annotations
cgMap_Grid, Map=map, /Box, charsize = 2
cgMap_Continents, Map=map, Color='black'
;add map legend
cgLegend, SymColors = ['red', 'forestgreen'], PSyms = [2, 7], $
  Symsize = 1.5, Location = [0.7, 0.7], titles = ['OGI', 'NOAA'], $
  /Box, /Background, BG_Color = 'rose', /center_sym, length = 0


;----plotting the UCI and UCI 84-96 data----
;create second display for the UCI data
cgDisplay, 990, 540, Title='Map Plot 2', /free
;set up the map projection
xrange = [-180,180] ;longitudes
yrange = [-90,90] ;latitudes
center_lon = 0
map = Obj_New('cgMap', 'Equirectangular', Ellipsoid=19, $
   XRange=xrange, YRange=yrange, /LatLon_Ranges, CENTER_LON=center_lon, $
   Position=[0.1, 0.1, 0.9, 0.8], Limit=[-90, -180, 90, 180])
map -> Draw
;plot the UCI coordinates (old and new)
cgPlot, lon_uci, lat_uci, color = 'red', psym = 5, symsize = 1.2, Position=[0.1, 0.1, 0.9, 0.8], $
  xstyle = 4, ystyle = 4, xrange = [-180, 180], yrange = [-90, 90]
cgPlot, lon_uci8496, lat_uci8496, color = 'forestgreen', psym = 9, symsize = 1.2, Position=[0.1, 0.1, 0.9, 0.8], $
  /overplot
;add map annotations
cgMap_Grid, Map=map, /Box, charsize = 2
cgMap_Continents, Map=map, Color='black'
;add map legend
cgLegend, SymColors = ['red', 'forestgreen'], PSyms = [5, 9], $
  Symsize = 1.5, Location = [0.7, 0.7], titles = ['UCI', 'UCI 84-96'], $
  /Box, /Background, BG_Color = 'rose', /center_sym, length = 0

endif

if (plot_optn eq 3) then begin

  ;setup the Coyote Graphics display
  cgDisplay, 1584, 864, Title='Map Plot'
  ;set up the map projection
  xrange = [-180,180] ;longitudes
  yrange = [-90,90] ;latitudes
  center_lon = 0

  ;OGI sites
  map = Obj_New('cgMap', 'Equirectangular', Ellipsoid=19, $
     XRange=xrange, YRange=yrange, /LatLon_Ranges, CENTER_LON=center_lon, $
     Position=[0.05, 0.05, 0.475, 0.475], Limit=[-90, -180, 90, 180])
  map -> Draw
  ;add map annotations
  cgMap_Grid, Map=map, /Box, charsize = 1.1
  cgMap_Continents, Map=map, Color='black'
  ;plot the OGI sites
  cgPlot, lon_ogi, lat_ogi, color = 'red', psym = 2, symsize = 1.5, Position=[0.05, 0.05, 0.475, 0.475], $
    xstyle = 4, ystyle = 4, xrange = [-180, 180], yrange = [-90, 90], label = 'OGI sites', charsize = 2

  ;NOAA sites
  map = Obj_New('cgMap', 'Equirectangular', Ellipsoid=19, $
     XRange=xrange, YRange=yrange, /LatLon_Ranges, CENTER_LON=center_lon, $
     Position=[0.525, 0.525, 0.95, 0.95], Limit=[-90, -180, 90, 180])
  map -> Draw
  ;add map annotations
  cgMap_Grid, Map=map, /Box, charsize = 1.1
  cgMap_Continents, Map=map, Color='black'
  ;plot the OGI sites
  cgPlot, lon_noaa, lat_noaa, color = 'blue', psym = 2, symsize = 1.5, Position=[0.525, 0.525, 0.95, 0.95], $
    xstyle = 4, ystyle = 4, xrange = [-180, 180], yrange = [-90, 90], label = 'NOAA sites', charsize = 2

endif

end
