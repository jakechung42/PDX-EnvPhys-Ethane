FUNCTION sind, angle

;this function takes sine value in degrees
angle_rad = angle*!PI/180
out = sin(angle_rad)
return, out

end

;======================================================================
FUNCTION MakeLatBin, bin_bound, lat, year, ratio

;a modified version of MakeLatBin subfunction from the program simGlobalDataAnalysis to run 
;for the simulated data set. The difference is that the simulated data only has annual data 
;so don't need the month vector.

;this function takes in the above inputs, then run through an algorithm to bin the latitudes
;as specicfied by the bin_bound input. Latitude bands that is not specified in the bin_bound 
;array will be discarded.
;band id will be set as the mid latitude of that lat band.
;the longitudes data is irrelevant in this calculation. So we can ignore the longitudes.
;we are going to make a new array and convert the latitudes to their appropriate bands

bin_lat = []
bin_year = []
bin_ratio = []
;algorithm to move the latitude into bins
;n_elements(bin_bound)-2 because we don't need it to run to the end. It just how the algorithm works
for i = 0, n_elements(bin_bound)-2 do begin
	x = where(lat ge bin_bound[i] and lat lt bin_bound[i + 1], count)
	if (count gt 0) then begin
		temp0 = (bin_bound[i] + bin_bound[i + 1])/2 ;id the different bins
		temp1 = fltarr(n_elements(x))
		temp1[*] = temp0
		bin_lat = [bin_lat, temp1]
		bin_year = [bin_year,year[x]]
		bin_ratio = [bin_ratio, ratio[x]]
	endif
endfor

bin_out = fltarr(3, n_elements(bin_ratio))
bin_out[0, *] = bin_lat
bin_out[1, *] = bin_year
bin_out[2, *] = bin_ratio

return, bin_out
end

;======================================================================
;this program studies the sensitivity of the simulated data when sampled with different methods
;refer to issue #10 of the Github repo

PRO simDataSensitivity

compile_opt idl2

;specify the latitude bins
;restriction for the bin_bound array:
;	Southern Hemisphere must be negative
;	Latitude 0 must be in the array
bin_bound = [-50, -30, 0, 30, 50, 75]

;read in sim data first. 
;use default emission scenario
ems_dir = "/home/excluded-from-backup/data/C2H6/trac_avg.PSUSF_1981_2015.bpch"
print, ems_dir

;extract the 3D array that contains the global mixing ratio
;read in the file
ctm_get_data, datainfo, 'IJ-AVG-$', tracer=1, filename=ems_dir, tau0=tau0
;Get MODELINFO and GRIDINFO structures, xmid, ymid hold lon/lat centers
getmodelandgridinfo, datainfo[0], modelinfo, gridinfo
;nt is the number of data points in the simulation.
nt = n_elements(datainfo)
simArr = fltarr(nt, 144, 91)
;sim_ymd stores the entire time information of the input simulation
sim_ymd = tau2yymmdd(tau0, /GEOS1)
;this section call out the data blocks of the simulated data and store it in simArr
;simArr is a 3-D array with the following attributes:
;simArr[ *, 0, 0 ]: the time dimension
;simArr[ 0, *, 0 ]: longitudes
;simArr[ 0, 0, * ]: latitudes
;the mixing ratio is the dependent variable and the time, lon, lat are the independent variables
for i = 0, nt - 1 do begin
	data= CTM_EXTRACT(*(datainfo[i].data), modelinfo= modelinfo, $
	gridinfo= gridinfo, lat= [-90, 90], lon= [-180, 180], alrange = [1, 2], $
	average= 4)
	for k = 0, n_elements(data[* ,0 ]) - 1 do begin
		for j = 0, n_elements(data[0 ,* ]) - 1 do begin
			simArr[i, k, j] = data[k, j]/2 * 1000
		endfor
	endfor
endfor
;clean up memory
CTM_CLEANUP 

;simYear contains the years that the simulation has
simYear = sim_ymd.year[sort(sim_ymd.year)]
simYear = simYear[uniq(simYear)]

;need to calculate the annual average of the sim data
;there are 34 years in the simulated data 
;create another 3D array that will contain the annual mean of the simulated data
annualSimGlobal = fltarr(34, 144, 91)
for i = 0, nt - 1 do begin
	if (i eq 407) then break
	annualSimGlobal[i/12, *, *] = mean(simArr[i : i + 11, *, *], 1)
	i = i + 10 ;jump 10 indices, not 11 because the for loop will add 1
endfor
;the annualSimGlobal variable contains the simulated annual global data 
;this variable can be used in all 3 methods of calculating the simulated data


;======================================================================
;the first method of calculating the IHR from the sim data is to use the complete data 
;and not consider the observed data

;average across longitude to make latitudinal profiles for each year thus reduces to 2 independent variables
annualLatProfile = fltarr(34, 91)
for i = 0, 33 do begin
	annualLatProfile[i, *] = mean(annualSimGlobal[i, *, *], 2)
endfor 

;convert the bin_bound to index of the simulated data output
bin_bound_idx = fltarr(n_elements(bin_bound))
for i = 0, n_elements(bin_bound) - 1 do begin
	lon = 1 ;give a random longitude value since we don't care about the longitude in this case
	CTM_INDEX, CTM_TYPE('GEOS1', RESOLUTION= 2), IDX_lon, IDX_lat, CENTER= [bin_bound[i], lon], /non_interactive
	bin_bound_idx[i] = IDX_lat
endfor

;CTM_INDEX doesn't account for IDL starting from 0 so need to shift the bin_bound_idx by 1
bin_bound_idx[*] = bin_bound_idx[*] - 1
;group the data into latitudinal bins as specified above in bin_bound variable
;the number of bin is determined by n_elements(bin_bound)-1
;all the different sites in the same latitudinal bin will be averaged to get an unified value
annualBinArr = fltarr(n_elements(bin_bound_idx)-1, n_elements(simYear))
for i = 0, 33 do begin
	for j = 0, n_elements(bin_bound_idx) -2 do begin
		annualBinArr[j, i] = mean(annualLatProfile[i, bin_bound_idx[j] : bin_bound_idx[j +1]])
	endfor
endfor
;annualBinArr contains the time series averaged of each bin
	;bin ---->
;time
;	|
;	|
;	|
;	V
      ; 271.783      356.959      729.085      1298.60      1495.73
      ; 266.137      346.695      733.309      1333.07      1492.77
      ; 270.524      361.553      734.369      1293.23      1474.48
      ; 271.667      361.351      735.280      1290.20      1492.03
      ; 278.804      359.126      732.173      1311.05      1491.47
      ; 279.984      355.169      734.358      1309.11      1485.26
      ; 281.931      354.071      730.357      1312.66      1481.43
      ; 278.257      359.854      730.640      1295.70      1499.27
      ; 276.357      357.693      729.350      1289.04      1497.99
      ; 279.750      360.397      725.935      1296.37      1517.48
      ; 271.323      353.321      731.940      1318.87      1494.20
      ; 270.098      361.074      735.077      1300.25      1475.76
      ; 276.441      361.725      734.450      1285.89      1488.58
      ; 279.308      358.691      732.405      1313.80      1498.73
      ; 279.538      353.889      734.379      1308.96      1481.17
      ; 281.881      354.801      732.025      1309.03      1484.52
      ; 278.871      359.305      729.530      1294.61      1494.63
      ; 275.850      358.387      726.589      1292.75      1508.74
      ; 279.591      361.673      726.157      1299.34      1510.67
      ; 268.921      349.955      733.593      1334.83      1498.87
      ; 269.321      352.831      731.556      1309.64      1487.99
      ; 269.458      361.639      735.450      1295.71      1493.60
      ; 280.286      359.892      733.219      1312.93      1501.03
      ; 279.887      352.775      732.961      1311.83      1473.14
      ; 280.699      356.700      733.079      1305.64      1484.97
      ; 277.105      358.810      727.935      1297.50      1499.95
      ; 276.460      357.477      727.362      1293.74      1514.27
      ; 278.876      361.770      727.004      1297.99      1510.17
      ; 268.145      347.650      736.318      1340.56      1492.70
      ; 270.683      356.255      732.098      1298.82      1477.31
      ; 268.575      361.786      738.374      1296.31      1490.98
      ; 276.529      359.320      735.026      1295.51      1487.05
      ; 281.373      361.384      733.593      1302.38      1499.19
      ; 280.199      351.951      733.789      1311.64      1470.20

;determine bins in the southern and bins in northern hemispheres 
sou_count = where(bin_bound lt 0, count)
sou_count = count
nor_count = where(bin_bound gt 0, count)
nor_count = count

;separate the north and south
souArr = annualBinArr[0 : (sou_count - 1), *]
norArr = annualBinArr[sou_count : sou_count + nor_count - 1, *] 

;calculate the weights of each lat bin:
;prep arrays to calculate weights
neg_idx = where(bin_bound lt 0)
neg = bin_bound[neg_idx]
pos_idx = where(bin_bound gt 0)
pos = bin_bound[pos_idx]
zero = 0
neg = [neg, zero]
pos = [zero, pos]

sou_weight = fltarr(n_elements(neg_idx))
nor_weight = fltarr(n_elements(pos_idx))

for i = 0, n_elements(neg) - 2 do begin
	sou_weight[i] = abs((sind(neg[i+1]) - sind(neg[i]))/sind(neg[0]))
endfor

for i = 0, n_elements(pos) - 2 do begin
	nor_weight[i] = abs((sind(pos[i+1]) - sind(pos[i]))/sind(pos[n_elements(pos)-1]))
endfor

annualSouMean = fltarr(n_elements(souArr[0, *]))
;apply weights to the hemispheric means
for i = 0, n_elements(souArr[0, *])-1 do begin
	annualSouMean[i] = total(souArr[*, i]*sou_weight)
endfor

annualNorMean = fltarr(n_elements(norArr[0, *]))
for i = 0, n_elements(norArr[0, *])-1 do begin
	annualNorMean[i] = total(norArr[*, i]*nor_weight)
endfor

;calculate the annual IHR 
annual_IHR = annualNorMean/annualSouMean

;======================================================================
;method 2 of calculating IHR of simulated data 
;import .dat file from the simgGlobalDataAnalysis program since it already had one of the IHR calculated
infile = "/home/excluded-from-backup/ethane/IDL/temp_file/SceE_SenStudy_allNetworks.dat"

;open file to read
openr, lun, infile, /get_lun
;initiate varibale to store the data
masterVar = fltarr(5, file_lines(infile))
;read file
readf, lun, masterVar
;release memory
free_lun, lun

;separate the different networks
ogi_idx = where(masterVar[0, *] eq 1)
ogi2 = masterVar[*, ogi_idx]

uci_idx = where(masterVar[0, *] eq 2)
uci2 = masterVar[*, uci_idx]

noaa_idx = where(masterVar[0, *] eq 3)
noaa2 = masterVar[*, noaa_idx]

;======================================================================
;method 3 of calculating the IHR of the simulated data
;calculate the IHR using the coordinates of all the available sites
;this method differs from the second one in the point that the simulated data 
;and the observed data only have spatial relationship but
;not temporal, since the time period of a site is ignored.

;import the coordinates of all sites from all networks
infile = "/home/excluded-from-backup/ethane/IDL/temp_file/network_coordinates.dat"

;read file
openr, lun, infile, /get_lun
allCoor = fltarr(2, file_lines(infile))
readf, lun, allCoor
;free up memory
free_lun, lun
;allCoor contains coordinate lat, lon for all sites from all 3 networks

idx = sort(allCoor[0, *])
allCoor = allCoor[*, idx]

;need to remove duplicate coordinates
;get unique latitudes
lat_idx = uniq(allCoor[0, *])
u_lat = allCoor[0, lat_idx]
all_lat = []
all_lon = []
;re-build the allCoor matrix to remove duplicate values
for i = 0, n_elements(u_lat)-1 do begin
	;retrieve latitude individually
	a = where(allCoor[0, *] eq u_lat[i])
	temp = allCoor[*, a] 
	;to get unique longitudes, need to sort it first
	lon_idx = sort(temp[1, *])
	temp = temp[*, lon_idx]
	lon_idx = uniq(temp[1, *])
	u_lon = temp[1, lon_idx]
	for j = 0, n_elements(u_lon)-1 do begin
		;stacking up the values
		all_lat = [all_lat, u_lat[i]]
		all_lon = [all_lon, u_lon[j]]
	endfor
endfor

allCoor = [rotate(all_lat, 1), rotate(all_lon, 1)]

;convert the coordinates to index values
allCoor_idx = fltarr(2, n_elements(allCoor[0, *]))
for i = 0, n_elements(allCoor[0, *])-1 do begin
	CTM_INDEX, CTM_TYPE('GEOS1', RESOLUTION= 2), a, b, CENTER= [allCoor[0, i], allCoor[1, i]], /non_interactive
	allCoor_idx[1, i] = a ;longitude
	allCoor_idx[0, i] = b ;latitude
endfor

;shift allCoor_idx by 1 index because CTM_INDEX starts from 1
allCoor_idx[*, *] = allCoor_idx[*, *] - 1
;define some variables to build the ratio table
t = [] ;time
la = [] ;lat
lo = [] ;lon
r = [] ;ratio
for i = 0, n_elements(allCoor_idx[0, *])-1 do begin
	ratio = annualSimGlobal[*, allCoor_idx[1, i], allCoor_idx[0, i]] ;get the ratio from coordinate indecies 
	for j = 0, n_elements(ratio)-1 do begin ;loop through ratio to fill in the data for lat, lon and time
		t = [t, j + 1981]
		la = [la, allCoor[0, i]]
		lo = [lo, allCoor[1, i]]
		r = [r, ratio[j]]
	endfor
endfor

;annualRawSim3 contains a table that has the following headers
;Time, Lat, Lon, Ratio which is from the simulated data
;sampled at the coordinates from the obs data.
annualRawSim3 = [rotate(t, 1), rotate(la, 1), rotate(lo, 1), rotate(r, 1)]

idx = sort(annualRawSim3[0, *])
annualRawSim3 = annualRawSim3[*, idx]

;separate the simulated data into bins based on latitudes 
annualBinSim3 = MakeLatBin(bin_bound, annualRawSim3[1, *], annualRawSim3[0, *], annualRawSim3[3, *])

;obtain the bin id, which is the mid latitude of bin_bound
temp0 = annualBinSim3[0, sort(annualBinSim3[0, *])]
mid_bin_bound = annualBinSim3[0, uniq(temp0)]

;prep variables
lat_id = []
year1 = []
ratio1 = []

;calculate the annual bin average
for i = 0, n_elements(mid_bin_bound)-1 do begin ;loop through each bin 
	y = where(annualBinSim3[0, *] eq mid_bin_bound[i], count);pull the the segment of the data with the appropriate bin
	if count eq 0 then begin ;checkpoit
		print, 'There must be at elast one annualBinSim3 eq to the bin id or else there is problem'
		stop
	endif
	temp0 = annualBinSim3[*, y] 
	for j = 0, n_elements(simYear)-1 do begin;loop through each year to calculate the mean of each year
		x = where(temp0[1, *] eq simYear[j])
		;rebuild the array
		lat_id = [lat_id, mid_bin_bound[i]]
		year1 = [year1, simYear[j]]
		ratio1 = [ratio1, mean(temp0[2, x])]
	endfor
endfor

annualAllBin = [rotate(lat_id, 1), rotate(year1, 1), rotate(ratio1, 1)]

;applying weights to the ratio
for i = 0, n_elements(annualAllBin[0, *])-1 do begin
	loc = where(annualAllBin[0, i] eq mid_bin_bound)
	if (annualAllBin[0, i] gt 0) then begin
		weight = abs((sind(bin_bound[loc+1]) - sind(bin_bound[loc]))/sind(bin_bound[n_elements(bin_bound)-1]))
		annualAllBin[2, i] = annualAllBin[2, i]*weight
	endif else begin
		weight = abs((sind(bin_bound[loc+1]) - sind(bin_bound[loc]))/sind(bin_bound[0]))
		annualAllBin[2, i] = annualAllBin[2, i]*weight
	endelse
endfor

;calculating hemispheric means
sou3_idx = where(annualAllBin[0, *] lt 0)
sou3 = annualAllBin[*, sou3_idx]
nor3_idx = where(annualAllBin[0, *] gt 0)
nor3 = annualAllBin[*, nor3_idx]

annualNor3 = fltarr(2, n_elements(simYear))
for i = 0, n_elements(simYear)-1 do begin
	x = where(nor3[1, *] eq simYear[i])
	annualNor3[0, i] = simYear[i]
	annualNor3[1, i] = total(nor3[2, x])
endfor

annualSou3 = fltarr(2, n_elements(simYear))
for i = 0, n_elements(simYear)-1 do begin
	x = where(sou3[1, *] eq simYear[i])
	annualSou3[0, i] = simYear[i]
	annualSou3[1, i] = total(sou3[2, x])
endfor


annual_IHR_3 = annualNor3[1, *]/ annualSou3[1, *]


;plotting procedure
;set up plot
open_device, /ps, /color, file='temp.eps', margin=0.05, xsize = 10.0, ysize = 7.5
!x.thick=1
!y.thick=1
!p.font=1
!p.thick =0.5

multiplot, /default    ; resets multiplot settings
multiplot, [1,3], ygap=0.002, xgap=0;  sets up multiplot 

;plot northern hemisphere of all 3 methods
cgPlot, simYear, annualNorMean, xrange = [1980, 2015], xticklen = 1, xgridstyle = 1, $
	xticks = 17, /nodata, ytitle = 'Mixing ratio (pptv)', $
	title = 'Sensitivity study of the simulated IHR', yrange = [600, 1200]
cgPlot, simYear, annualNorMean, /overplot, color = 'red'
cgPlot, simYear, annualNor3[1, *], /overplot, color = 'blue'
cgPlot, ogi2[1, *], ogi2[2, *], psym = 4, color = 'black', /overplot
cgPlot, uci2[1, *], uci2[2, *], psym = 2, color = 'black', /overplot
cgPlot, noaa2[1, *], noaa2[2, *], psym = 5, color = 'black', /overplot

multiplot, /doyaxis, /doxaxis

;plot southern hemisphere of all 3 methods
cgPlot, simYear, annualSouMean, xrange = [1980, 2015], xticklen = 1, xgridstyle = 1, $
	xticks = 17, /nodata, ytitle = 'Mixing ratio (pptv)', XTickformat='(A1)', $
	yrange = [200, 400]
cgPlot, simYear, annualSouMean, /overplot, color = 'red'
cgPlot, simYear, annualSou3[1, *], /overplot, color = 'blue'
cgPlot, ogi2[1, *], ogi2[3, *], psym = 4, color = 'black', /overplot
cgPlot, uci2[1, *], uci2[3, *], psym = 2, color = 'black', /overplot
cgPlot, noaa2[1, *], noaa2[3, *], psym = 5, color = 'black', /overplot

multiplot, /doyaxis, /doxaxis

;plot IHR of all 3 methods
cgPlot, simYear, annual_IHR, xrange = [1980, 2015], xticklen = 1, xgridstyle = 1, $
	xticks = 17, /nodata, ytitle = 'IHR', yrange = [2, 4]
cgPlot, simYear, annual_IHR, /overplot, color = 'red'
cgPlot, simYear, annual_IHR_3, /overplot, color = 'blue'
cgPlot, ogi2[1, *], ogi2[4, *], psym = 4, color = 'black', /overplot
cgPlot, uci2[1, *], uci2[4, *], psym = 2, color = 'black', /overplot
cgPlot, noaa2[1, *], noaa2[4, *], psym = 5, color = 'black', /overplot

cgLegend, title = ['All Global Data', 'All coordinates, ignore the temporal data'], $
	psym = [3, 3], linestyle = [0, 0], $
	color = ['red', 'blue'], $
	location = [1984, 2.5], charsize = 1, /data, vspace = 1
	
close_device

spawn, 'gv temp.eps'

multiplot, /default


end
