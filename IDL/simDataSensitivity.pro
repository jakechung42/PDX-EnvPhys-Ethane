;this program studies the sensitivity of the simulated data when sampled with different methods
;refer to issue #10 of the Github repo

PRO simDataSensitivity

compile_opt idl2


;specify the latitude bins
bin_bound = [-50, -30, 0, 30, 50, 75]