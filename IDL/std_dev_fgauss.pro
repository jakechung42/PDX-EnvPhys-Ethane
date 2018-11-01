function std_dev_fgauss, data

data = double( data )

;use the cgHistplot to calculate the Gaussian fit
open_device, /ps, /color, file='temp.eps', margin=0.05, xsize = 10.0, ysize = 7.5
!x.thick=1
!y.thick=1
!p.font=0
!p.thick =1

;Use the cgHistoplot to generate the data needed for GaussFit routine
cgHistoplot, data, /fillpoly, LOCATIONS=loc, BINSIZE=binsize, HISTDATA=h
binCenters = loc + (binsize / 2.0)
yfit = GaussFit(binCenters, h, coeff, NTERMS=3)

close_device

return, coeff[ 2 ]

end
