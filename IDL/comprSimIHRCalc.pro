FUNCTION get_cat, sce13

;function that separates the different categories out from the input files
;and outputs a structure that contains all the different categories

m1_idx = where(sce13[0, *] eq 1, count)
if count eq 0 then begin
	print, 'probably miss-matching veriable name, check input for get_cat function'
	stop
endif
m3_idx = where(sce13[0, *] eq 3, count)
if count eq 0 then begin
	print, 'probably miss-matching veriable name, check input for get_cat function'
	stop
endif

out = {m1: sce13[1:4, m1_idx], $
		m3: sce13[1:4, m3_idx] }
return, out

end

;======================================================================
FUNCTION sep_network, sce2

;this function separates the 3 networks out from the sce method 2 input files
;the output is a structure for each network

ogi_idx = where(sce2[0, *] eq 1, count)
if count eq 0 then begin
	print, 'probably miss-matching veriable name, check input for sep_network function'
	stop
endif
uci_idx = where(sce2[0, *] eq 2, count)
if count eq 0 then begin
	print, 'probably miss-matching veriable name, check input for sep_network function'
	stop
endif
noaa_idx = where(sce2[0, *] eq 3, count)
if count eq 0 then begin
	print, 'probably miss-matching veriable name, check input for sep_network function'
	stop
endif

out = {noaa: sce2[1:4, noaa_idx], $
		uci: sce2[1:4, uci_idx], $
		ogi: sce2[1:4, ogi_idx] }
return, out

end

;======================================================================
PRO comprSimIHRCalc, $
	method1 = method1, $
	method2 = method2, $
	method3 = method3, $
	normalize_plot = normalize_plot

;different plots will be plotted depending on the input method argument
;this program imports the IHR data from all 6 scenarios calculated using different methods
;the 3 main methods:
;	Using all global sim data
;	Using the coordinates only of the obs data to retrieve the sim data
;	Using the coordinates and time period of each location from the obs data to retrive the sim data

;import all data files
infileA13 = "/home/excluded-from-backup/ethane/IDL/temp_file/SceA_IHR_Method_1_3.dat"
infileB13 = "/home/excluded-from-backup/ethane/IDL/temp_file/SceB_IHR_Method_1_3.dat"
infileC13 = "/home/excluded-from-backup/ethane/IDL/temp_file/SceC_IHR_Method_1_3.dat"
infileD13 = "/home/excluded-from-backup/ethane/IDL/temp_file/SceD_IHR_Method_1_3.dat"
infileE13 = "/home/excluded-from-backup/ethane/IDL/temp_file/SceE_IHR_Method_1_3.dat"
infileF13 = "/home/excluded-from-backup/ethane/IDL/temp_file/SceF_IHR_Method_1_3.dat"

infileA2 = "/home/excluded-from-backup/ethane/IDL/temp_file/SceA_SenStudy_allNetworks.dat"
infileB2 = "/home/excluded-from-backup/ethane/IDL/temp_file/SceB_SenStudy_allNetworks.dat"
infileC2 = "/home/excluded-from-backup/ethane/IDL/temp_file/SceC_SenStudy_allNetworks.dat"
infileD2 = "/home/excluded-from-backup/ethane/IDL/temp_file/SceD_SenStudy_allNetworks.dat"
infileE2 = "/home/excluded-from-backup/ethane/IDL/temp_file/SceE_SenStudy_allNetworks.dat"
infileF2 = "/home/excluded-from-backup/ethane/IDL/temp_file/SceF_SenStudy_allNetworks.dat"

openr, lunA13, infileA13, /get_lun
openr, lunB13, infileB13, /get_lun
openr, lunC13, infileC13, /get_lun
openr, lunD13, infileD13, /get_lun
openr, lunE13, infileE13, /get_lun
openr, lunF13, infileF13, /get_lun
openr, lunA2, infileA2, /get_lun
openr, lunB2, infileB2, /get_lun
openr, lunC2, infileC2, /get_lun
openr, lunD2, infileD2, /get_lun
openr, lunE2, infileE2, /get_lun
openr, lunF2, infileF2, /get_lun

sceA13 = fltarr(5, file_lines(infileA13))
sceB13 = fltarr(5, file_lines(infileB13))
sceC13 = fltarr(5, file_lines(infileC13))
sceD13 = fltarr(5, file_lines(infileD13))
sceE13 = fltarr(5, file_lines(infileE13))
sceF13 = fltarr(5, file_lines(infileF13))

sceA2 = fltarr(5, file_lines(infileA2))
sceB2 = fltarr(5, file_lines(infileB2))
sceC2 = fltarr(5, file_lines(infileC2))
sceD2 = fltarr(5, file_lines(infileD2))
sceE2 = fltarr(5, file_lines(infileE2))
sceF2 = fltarr(5, file_lines(infileF2))


readf, lunA13, sceA13 
readf, lunB13, sceB13 
readf, lunC13, sceC13 
readf, lunD13, sceD13 
readf, lunE13, sceE13 
readf, lunF13, sceF13 
readf, lunA2, sceA2 
readf, lunB2, sceB2 
readf, lunC2, sceC2 
readf, lunD2, sceD2 
readf, lunE2, sceE2 
readf, lunF2, sceF2 

free_lun, lunA13, lunB13, lunC13, lunD13, lunE13, lunF13, lunA2, lunB2, lunC2, lunD2, lunE2, lunF2

sceA = get_cat(sceA13)
sceB = get_cat(sceB13)
sceC = get_cat(sceC13)
sceD = get_cat(sceD13)
sceE = get_cat(sceE13)
sceF = get_cat(sceF13)

sceA2 = sep_network(sceA2)
sceB2 = sep_network(sceB2)
sceC2 = sep_network(sceC2)
sceD2 = sep_network(sceD2)
sceE2 = sep_network(sceE2)
sceF2 = sep_network(sceF2)
if KEYWORD_SET(normalize_data) then begin
case 1 of
	KEYWORD_SET(method1): begin
		;set up plot
		open_device, /ps, /color, file='temp.eps', margin=0.05, xsize = 10.0, ysize = 7.5
		!x.thick=1
		!y.thick=1
		!p.font=1
		!p.thick =0.5

		multiplot, /default    ; resets multiplot settings
		multiplot, [1,3], ygap=0.002, xgap=0;  sets up multiplot 
		
		;plot norther hemisphere of method1
		cgPlot, sceA.m1[0, *], sceA.m1[1, *], /nodata, xrange = [1981, 2014], xticklen = 1, xgridstyle = 1, $
			xticks = 16, ytitle = 'Mixing ratio(pptv)', yrange = [500, 1500], $
			title = 'Sensitivity study of the simulated IHR (no obs data)'
		cgPlot, sceA.m1[0, *], sceA.m1[1, *], /overplot, color = 'blue violet'
		cgPlot, sceB.m1[0, *], sceB.m1[1, *], /overplot, color = 'blue'
		cgPlot, sceC.m1[0, *], sceC.m1[1, *], /overplot, color = 'forest green'
		cgPlot, sceD.m1[0, *], sceD.m1[1, *], /overplot, color = 'lime green'
		cgPlot, sceE.m1[0, *], sceE.m1[1, *], /overplot, color = 'red'
		cgPlot, sceF.m1[0, *], sceF.m1[1, *], /overplot, color = 'magenta'
		multiplot, /doyaxis, /doxaxis
		;plot southern hemisphere method1
		cgPlot, sceA.m1[0, *], sceA.m1[2, *], /nodata, xrange = [1981, 2014], xticklen = 1, xgridstyle = 1, $
			xticks = 16, ytitle = 'Mixing ratio(pptv)', XTickformat='(A1)', yrange = [250, 500]
		cgPlot, sceA.m1[0, *], sceA.m1[2, *], /overplot, color = 'blue violet'
		cgPlot, sceB.m1[0, *], sceB.m1[2, *], /overplot, color = 'blue'
		cgPlot, sceC.m1[0, *], sceC.m1[2, *], /overplot, color = 'forest green'
		cgPlot, sceD.m1[0, *], sceD.m1[2, *], /overplot, color = 'lime green'
		cgPlot, sceE.m1[0, *], sceE.m1[2, *], /overplot, color = 'red'
		cgPlot, sceF.m1[0, *], sceF.m1[2, *], /overplot, color = 'magenta'		
		multiplot, /doyaxis, /doxaxis
		cgPlot, sceA.m1[0, *], sceA.m1[2, *], /nodata, xrange = [1981, 2014], xticklen = 1, xgridstyle = 1, $
			xticks = 16, ytitle = 'IHR', yrange = [1, 4]
		cgPlot, sceA.m1[0, *], sceA.m1[3, *], /overplot, color = 'blue violet'
		cgPlot, sceB.m1[0, *], sceB.m1[3, *], /overplot, color = 'blue'
		cgPlot, sceC.m1[0, *], sceC.m1[3, *], /overplot, color = 'forest green'
		cgPlot, sceD.m1[0, *], sceD.m1[3, *], /overplot, color = 'lime green'
		cgPlot, sceE.m1[0, *], sceE.m1[3, *], /overplot, color = 'red'
		cgPlot, sceF.m1[0, *], sceF.m1[3, *], /overplot, color = 'magenta'	
		
		close_device

		spawn, 'gv temp.eps'

		multiplot, /default
	end
	KEYWORD_SET(method2): begin
		;set up plot
		open_device, /ps, /color, file='temp.eps', margin=0.05, xsize = 10.0, ysize = 7.5
		!x.thick=1
		!y.thick=1
		!p.font=1
		!p.thick =0.5

		multiplot, /default    ; resets multiplot settings
		multiplot, [1,3], ygap=0.002, xgap=0;  sets up multiplot 
		
		;plot norther hemisphere of method1
		cgPlot, sceA2.noaa[0, *], sceA2.noaa[1, *], /nodata, xrange = [1981, 2014], xticklen = 1, xgridstyle = 1, $
			xticks = 16, ytitle = 'Mixing ratio(pptv)', yrange = [500, 1500], $
			title = 'Sensitivity study of the simulated IHR (complete obs data)'
		cgPlot, sceA2.noaa[0, *], sceA2.noaa[1, *], /overplot, color = 'blue violet', psym = 5
		cgPlot, sceA2.uci[0, *], sceA2.uci[1, *], /overplot, color = 'blue violet', psym = 2
		cgPlot, sceA2.ogi[0, *], sceA2.ogi[1, *], /overplot, color = 'blue violet', psym = 4	
		cgPlot, sceB2.noaa[0, *], sceB2.noaa[1, *], /overplot, color = 'blue', psym = 5
		cgPlot, sceB2.uci[0, *], sceB2.uci[1, *], /overplot, color = 'blue', psym = 2
		cgPlot, sceB2.ogi[0, *], sceB2.ogi[1, *], /overplot, color = 'blue', psym = 4	
		cgPlot, sceC2.noaa[0, *], sceC2.noaa[1, *], /overplot, color = 'forest green', psym = 5
		cgPlot, sceC2.uci[0, *], sceC2.uci[1, *], /overplot, color = 'forest green', psym = 2
		cgPlot, sceC2.ogi[0, *], sceC2.ogi[1, *], /overplot, color = 'forest green', psym = 4	
		cgPlot, sceD2.noaa[0, *], sceD2.noaa[1, *], /overplot, color = 'lime green', psym = 5
		cgPlot, sceD2.uci[0, *], sceD2.uci[1, *], /overplot, color = 'lime green', psym = 2
		cgPlot, sceD2.ogi[0, *], sceD2.ogi[1, *], /overplot, color = 'lime green', psym = 4	
		cgPlot, sceE2.noaa[0, *], sceE2.noaa[1, *], /overplot, color = 'red', psym = 5
		cgPlot, sceE2.uci[0, *], sceE2.uci[1, *], /overplot, color = 'red', psym = 2
		cgPlot, sceE2.ogi[0, *], sceE2.ogi[1, *], /overplot, color = 'red', psym = 4	
		cgPlot, sceF2.noaa[0, *], sceF2.noaa[1, *], /overplot, color = 'magenta', psym = 5
		cgPlot, sceF2.uci[0, *], sceF2.uci[1, *], /overplot, color = 'magenta', psym = 2
		cgPlot, sceF2.ogi[0, *], sceF2.ogi[1, *], /overplot, color = 'magenta', psym = 4	
		multiplot, /doyaxis, /doxaxis
		;plot southern hemisphere method1
		cgPlot, sceA2.noaa[0, *], sceA2.noaa[2, *], /nodata, xrange = [1981, 2014], xticklen = 1, xgridstyle = 1, $
			xticks = 16, ytitle = 'Mixing ratio(pptv)', XTickformat='(A1)', yrange = [250, 500]
		cgPlot, sceA2.noaa[0, *], sceA2.noaa[2, *], /overplot, color = 'blue violet', psym = 5
		cgPlot, sceA2.uci[0, *], sceA2.uci[2, *], /overplot, color = 'blue violet', psym = 2
		cgPlot, sceA2.ogi[0, *], sceA2.ogi[2, *], /overplot, color = 'blue violet', psym = 4	
		cgPlot, sceB2.noaa[0, *], sceB2.noaa[2, *], /overplot, color = 'blue', psym = 5
		cgPlot, sceB2.uci[0, *], sceB2.uci[2, *], /overplot, color = 'blue', psym = 2
		cgPlot, sceB2.ogi[0, *], sceB2.ogi[2, *], /overplot, color = 'blue', psym = 4	
		cgPlot, sceC2.noaa[0, *], sceC2.noaa[2, *], /overplot, color = 'forest green', psym = 5
		cgPlot, sceC2.uci[0, *], sceC2.uci[2, *], /overplot, color = 'forest green', psym = 2
		cgPlot, sceC2.ogi[0, *], sceC2.ogi[2, *], /overplot, color = 'forest green', psym = 4	
		cgPlot, sceD2.noaa[0, *], sceD2.noaa[2, *], /overplot, color = 'lime green', psym = 5
		cgPlot, sceD2.uci[0, *], sceD2.uci[2, *], /overplot, color = 'lime green', psym = 2
		cgPlot, sceD2.ogi[0, *], sceD2.ogi[2, *], /overplot, color = 'lime green', psym = 4	
		cgPlot, sceE2.noaa[0, *], sceE2.noaa[2, *], /overplot, color = 'red', psym = 5
		cgPlot, sceE2.uci[0, *], sceE2.uci[2, *], /overplot, color = 'red', psym = 2
		cgPlot, sceE2.ogi[0, *], sceE2.ogi[2, *], /overplot, color = 'red', psym = 4	
		cgPlot, sceF2.noaa[0, *], sceF2.noaa[2, *], /overplot, color = 'magenta', psym = 5
		cgPlot, sceF2.uci[0, *], sceF2.uci[2, *], /overplot, color = 'magenta', psym = 2
		cgPlot, sceF2.ogi[0, *], sceF2.ogi[2, *], /overplot, color = 'magenta', psym = 4			
		multiplot, /doyaxis, /doxaxis
		cgPlot, sceA2.noaa[0, *], sceA2.noaa[3, *], /nodata, xrange = [1981, 2014], xticklen = 1, xgridstyle = 1, $
			xticks = 16, ytitle = 'IHR', yrange = [1, 4]
		cgPlot, sceA2.noaa[0, *], sceA2.noaa[3, *], /overplot, color = 'blue violet', psym = 5
		cgPlot, sceA2.uci[0, *], sceA2.uci[3, *], /overplot, color = 'blue violet', psym = 2
		cgPlot, sceA2.ogi[0, *], sceA2.ogi[3, *], /overplot, color = 'blue violet', psym = 4	
		cgPlot, sceB2.noaa[0, *], sceB2.noaa[3, *], /overplot, color = 'blue', psym = 5
		cgPlot, sceB2.uci[0, *], sceB2.uci[3, *], /overplot, color = 'blue', psym = 2
		cgPlot, sceB2.ogi[0, *], sceB2.ogi[3, *], /overplot, color = 'blue', psym = 4	
		cgPlot, sceC2.noaa[0, *], sceC2.noaa[3, *], /overplot, color = 'forest green', psym = 5
		cgPlot, sceC2.uci[0, *], sceC2.uci[3, *], /overplot, color = 'forest green', psym = 2
		cgPlot, sceC2.ogi[0, *], sceC2.ogi[3, *], /overplot, color = 'forest green', psym = 4	
		cgPlot, sceD2.noaa[0, *], sceD2.noaa[3, *], /overplot, color = 'lime green', psym = 5
		cgPlot, sceD2.uci[0, *], sceD2.uci[3, *], /overplot, color = 'lime green', psym = 2
		cgPlot, sceD2.ogi[0, *], sceD2.ogi[3, *], /overplot, color = 'lime green', psym = 4	
		cgPlot, sceE2.noaa[0, *], sceE2.noaa[3, *], /overplot, color = 'red', psym = 5
		cgPlot, sceE2.uci[0, *], sceE2.uci[3, *], /overplot, color = 'red', psym = 2
		cgPlot, sceE2.ogi[0, *], sceE2.ogi[3, *], /overplot, color = 'red', psym = 4	
		cgPlot, sceF2.noaa[0, *], sceF2.noaa[3, *], /overplot, color = 'magenta', psym = 5
		cgPlot, sceF2.uci[0, *], sceF2.uci[3, *], /overplot, color = 'magenta', psym = 2
		cgPlot, sceF2.ogi[0, *], sceF2.ogi[3, *], /overplot, color = 'magenta', psym = 4		
		
		close_device

		spawn, 'gv temp.eps'

		multiplot, /default
	end
	KEYWORD_SET(method3): begin
		;set up plot
		open_device, /ps, /color, file='temp.eps', margin=0.05, xsize = 10.0, ysize = 7.5
		!x.thick=1
		!y.thick=1
		!p.font=1
		!p.thick =0.5

		multiplot, /default    ; resets multiplot settings
		multiplot, [1,3], ygap=0.002, xgap=0;  sets up multiplot 
		
		;plot norther hemisphere of method1
		cgPlot, sceA.m3[0, *], sceA.m3[1, *], /nodata, xrange = [1981, 2014], xticklen = 1, xgridstyle = 1, $
			xticks = 16, ytitle = 'Mixing ratio(pptv)', yrange = [500, 1500], $
			title = 'Sensitivity study of the simulated IHR (only obs spatial data)'
		cgPlot, sceA.m3[0, *], sceA.m3[1, *], /overplot, color = 'blue violet'
		cgPlot, sceB.m3[0, *], sceB.m3[1, *], /overplot, color = 'blue'
		cgPlot, sceC.m3[0, *], sceC.m3[1, *], /overplot, color = 'forest green'
		cgPlot, sceD.m3[0, *], sceD.m3[1, *], /overplot, color = 'lime green'
		cgPlot, sceE.m3[0, *], sceE.m3[1, *], /overplot, color = 'red'
		cgPlot, sceF.m3[0, *], sceF.m3[1, *], /overplot, color = 'magenta'
		multiplot, /doyaxis, /doxaxis
		;plot southern hemisphere method1
		cgPlot, sceA.m3[0, *], sceA.m3[2, *], /nodata, xrange = [1981, 2014], xticklen = 1, xgridstyle = 1, $
			xticks = 16, ytitle = 'Mixing ratio(pptv)', XTickformat='(A1)', yrange = [250, 500]
		cgPlot, sceA.m3[0, *], sceA.m3[2, *], /overplot, color = 'blue violet'
		cgPlot, sceB.m3[0, *], sceB.m3[2, *], /overplot, color = 'blue'
		cgPlot, sceC.m3[0, *], sceC.m3[2, *], /overplot, color = 'forest green'
		cgPlot, sceD.m3[0, *], sceD.m3[2, *], /overplot, color = 'lime green'
		cgPlot, sceE.m3[0, *], sceE.m3[2, *], /overplot, color = 'red'
		cgPlot, sceF.m3[0, *], sceF.m3[2, *], /overplot, color = 'magenta'		
		multiplot, /doyaxis, /doxaxis
		cgPlot, sceA.m3[0, *], sceA.m3[2, *], /nodata, xrange = [1981, 2014], xticklen = 1, xgridstyle = 1, $
			xticks = 16, ytitle = 'IHR', yrange = [1, 4]
		cgPlot, sceA.m3[0, *], sceA.m3[3, *], /overplot, color = 'blue violet'
		cgPlot, sceB.m3[0, *], sceB.m3[3, *], /overplot, color = 'blue'
		cgPlot, sceC.m3[0, *], sceC.m3[3, *], /overplot, color = 'forest green'
		cgPlot, sceD.m3[0, *], sceD.m3[3, *], /overplot, color = 'lime green'
		cgPlot, sceE.m3[0, *], sceE.m3[3, *], /overplot, color = 'red'
		cgPlot, sceF.m3[0, *], sceF.m3[3, *], /overplot, color = 'magenta'	
		
		close_device

		spawn, 'gv temp.eps'

		multiplot, /default
	end
endcase

end
