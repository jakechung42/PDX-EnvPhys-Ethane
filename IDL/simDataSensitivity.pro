FUNCTION sind, angle

;this function takes sine value in degrees

angle_rad = angle*!PI/180

out = sin(angle_rad)

return, out

end

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
ems_dir = "/home/excluded-from-backup/data/C2H6/trac_avg.spinup_GFED4_MAVG_1981_2015.198101010000"
print, ems_dir

;extract the 3D array that contains the global mixing ratio
;read in the file
ctm_get_data, datainfo, 'IJ-AVG-$', tracer=1, filename=ems_dir, tau0=tau0
;Get MODELINFO and GRIDINFO structures, xmid, ymid hold lon/lat centers
getmodelandgridinfo, datainfo[0], modelinfo, gridinfo
;nt is the number of data points in the simulation.
nt = n_elements(datainfo)
simArr = fltarr(nt, 144, 91)
;sim_yymmdd stores the entire time information of the input simulation
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
;method 2 of calculating IHR of simulated data (copy from simGlobalDataAnalysis.pro


end 