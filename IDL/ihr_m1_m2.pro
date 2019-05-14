;this program plot the IHR of mehod 1 and method 2 
PRO ihr_m1_m2

compile_opt idl2

;read method 1 IHR data
infile1 = '/home/excluded-from-backup/ethane/IDL/temp_file/method1_IHR.dat'
infile2 = '/home/excluded-from-backup/ethane/IDL/temp_file/method2_IHR.dat'


openr, lun1, infile1, /get_lun
openr, lun2, infile2, /get_lun

n1_line = file_lines(infile1)
n2_line = file_lines(infile2)

;create varible that stores method 1 and method 2 ihr
m1 = fltarr(4, n1_line)
m2 = fltarr(4, n2_line)
;read file
readf, lun1, m1
readf, lun2, m2

free_lun, lun1, lun2

;calculate linear fit for the data
;remove NAN years
not_nan = where(finite(m1[2, *]))
m1_nonan = m1[1 : 3, not_nan]
not_nan = where(finite(m2[2, *]))
m2_nonan = m2[1 : 3, not_nan]

print, m2_nonan
;seperate IHR after year 2000
m1_ge2k_idx = where(m1_nonan[0, *] ge 2000)
m1_ge2k = m1_nonan[*, m1_ge2k_idx]
m2_ge2k_idx = where(m2_nonan[0, *] ge 2000)
m2_ge2k = m2_nonan[*, m2_ge2k_idx]
;seperate IHR before year 2000
m1_lt2k_idx = where(m1_nonan[0, *] lt 2000)
m1_lt2k = m1_nonan[*, m1_lt2k_idx]
m2_lt2k_idx = where(m2_nonan[0, *] lt 2000)
m2_lt2k = m2_nonan[*, m2_lt2k_idx]
;linear fit the < 2000 and > 2000 data
m1_ge2k_fit = linfit(m1_ge2k[0, *], m1_ge2k[1, *], MEASURE_ERRORS=m1_ge2k[2, *], YFIT = m1_ge2k_y)
m1_lt2k_fit = linfit(m1_lt2k[0, *], m1_lt2k[1, *], MEASURE_ERRORS=m1_lt2k[2, *], YFIT = m1_lt2k_y)
m2_ge2k_fit = linfit(m2_ge2k[0, *], m2_ge2k[1, *], MEASURE_ERRORS=m2_ge2k[2, *], YFIT = m2_ge2k_y)
m2_lt2k_fit = linfit(m2_lt2k[0, *], m2_lt2k[1, *], MEASURE_ERRORS=m2_lt2k[2, *], YFIT = m2_lt2k_y)


;seperate the network from the main array
m1_noaa_idx = where(m1[0, *] eq 1)
m1_noaa = m1[1 : 3, m1_noaa_idx]
m1_uci_idx = where(m1[0, *] eq 2)
m1_uci = m1[1 : 3, m1_uci_idx]
m1_ogi_idx = where(m1[0, *] eq 3)
m1_ogi = m1[1 : 3, m1_ogi_idx]

m2_noaa_idx = where(m2[0, *] eq 1)
m2_noaa = m2[1 : 3, m2_noaa_idx]
m2_uci_idx = where(m2[0, *] eq 2)
m2_uci = m2[1 : 3, m2_uci_idx]
m2_ogi_idx = where(m2[0, *] eq 3)
m2_ogi = m2[1 : 3, m2_ogi_idx]

;plot the IHR of method 1 and method 2
;set up plot
open_device, /ps, /color, file='temp.eps', margin=0.05, xsize = 10.0, ysize = 7.5
!x.thick=1
!y.thick=1
!p.font=1
!p.thick = 2

cgPlot, m1_noaa[0, *], m1_noaa[1, *], xrange = [1982, 2016], xticklen = 1, xgridstyle = 1, $
	xticks = 17, yrange = [-3,3], ytitle = 'Residual', err_yhigh = m1_noaa[2, *], $
	err_ylow = m1_noaa[2, *], psym = 5
cgPlot, m1_uci[0, *], m1_uci[1, *], err_yhigh = m1_uci[2, *], psym = 2, err_ylow = m1_uci[2, *], /overplot
cgPlot, m1_ogi[0, *], m1_ogi[1, *], err_yhigh = m1_ogi[2, *], psym = 4, err_ylow = m1_ogi[2, *], /overplot


cgPlot, m2_noaa[0, *], m2_noaa[1, *], err_yhigh = m2_noaa[2, *], psym = 5, err_ylow = m2_noaa[2, *], /overplot, $
	color = 'forest green'
cgPlot, m2_uci[0, *], m2_uci[1, *], err_yhigh = m2_uci[2, *], psym = 2, err_ylow = m2_uci[2, *], /overplot, $
	color = 'forest green'
cgPlot, m2_ogi[0, *], m2_ogi[1, *], err_yhigh = m2_ogi[2, *], psym = 4, err_ylow = m2_ogi[2, *], /overplot, $
	color = 'forest green'

;plot the linear fit function
cgPlot, m1_ge2k[0, *], m1_ge2k_y, /overplot
cgPlot, m1_lt2k[0, *], m1_lt2k_y, /overplot
cgPlot, m2_ge2k[0, *], m2_ge2k_y, /overplot, color = 'forest green'
cgPlot, m2_lt2k[0, *], m2_lt2k_y, /overplot, color = 'forest green'
close_device

spawn, 'gv temp.eps'

end
