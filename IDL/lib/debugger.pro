FUNCTION debugger, input

;this function stops the program and print the input to a ASCII file 
;to resume program use command .continue in the IDL command prompt
;there is no output to this function


openw, lun, 'temp.dat', /get_lun
printf, lun, input

free_lun, lun

spawn, 'gedit temp.dat'

stop

end
