PRO time_vs_lat_plot

  ;this program plots the time vs latitude of the new UCI data

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

;convert the time to decimal year for all data
noaa_dtime = fltarr(n_elements(noaa.year))
for i = 0, n_elements(noaa.year) - 1 do begin
  noaa_dtime[i] = noaa.year[i] + noaa.month[i]/12
endfor

uci_dtime = fltarr(n_elements(uci.year))
for i = 0, n_elements(uci.year) - 1 do begin
  uci_dtime[i] = uci.year[i] + uci.month[i]/12
endfor

uci_new_dtime = fltarr(n_elements(uci_new.year))
for i = 0, n_elements(uci_new.year) - 1 do begin
  uci_new_dtime[i] = uci_new.year[i] + uci_new.month[i]/12
endfor

ogi_dtime = fltarr(n_elements(ogi.year))
for i = 0, n_elements(ogi.year) - 1 do begin
  ogi_dtime[i] = ogi.year[i] + ogi.month[i]/12
endfor

;plotting procedure
open_device, /ps, /color, file='temp.eps', margin=0.05, xsize = 10.0, ysize = 7.5
!x.thick=2
!y.thick=2
!p.font=0
!p.thick =1

multiplot, /default    ; resets multiplot settings
;multiplot, [2,2], gap = 0.06;  sets up multiplot

;set up plot
cgPlot, ogi_dtime, ogi.lat, /nodata, yrange = [-100, 100], xrange = [1977, 2016], $
  ytitle = 'Latitudes', xtitle = 'Years', xticks = 10, xticklen = 1, $
  xgridstyle = 1, ygridstyle = 1

cgPlot, uci_new_dtime, uci_new.lat, psym = 4, color = 'red', /overplot

cgPlot, uci_dtime, uci.lat, psym = 5, color = 'blue', /overplot

cgPlot, noaa_dtime, noaa.lat, psym = 7, color = 'forest green', /overplot, $
  symsize = 0.7

cgPlot, ogi_dtime, ogi.lat, psym = 9, color = 'black', /overplot

cgLegend, SymColors = ['red', 'blue', 'forest green', 'black'], psyms = [4,5,7,9], $
  symsize = 0.8, /box, /background, bg_color = 'rose', /center_sym, length = 0, $
  titles = ['84 96 UCI', 'UCI', 'NOAA', 'OGI'], location = [0.45,0.3], charsize = 0.8

close_device
spawn, 'gv temp.eps'
end
