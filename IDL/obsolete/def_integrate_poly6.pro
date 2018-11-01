function def_integrate_poly6, cof0, $
						cof1, $
						cof2, $
						cof3, $
						cof4, $
						cof5, $
						cof6, $
						up_lim, $
						down_lim

cof0 = double( cof0 )
cof1 = double( cof1 )
cof2 = double( cof2 )
cof3 = double( cof3 )
cof4 = double( cof4 )
cof5 = double( cof5 )
cof6 = double( cof6 )
up_lim = double( up_lim )
down_lim = double( down_lim )

result_down = cof0 * (down_lim) + ( cof1 / 2 ) * (down_lim) ^ 2 + ( cof2 / 3 ) * (down_lim) ^ 3 + ( cof3 / 4 ) * (down_lim) ^ 4 + $
	( cof4 / 5 ) * (down_lim) ^ 5 + ( cof5 / 6 ) * (down_lim) ^ 6 + ( cof6 / 7 ) * (down_lim) ^ 7

result_up = cof0 * (up_lim) + ( cof1 / 2 ) * (up_lim) ^ 2 + ( cof2 / 3 ) * (up_lim) ^ 3 + ( cof3 / 4 ) * (up_lim) ^ 4 + $
	( cof4 / 5 ) * (up_lim) ^ 5 + ( cof5 / 6 ) * (up_lim) ^ 6 + ( cof6 / 7 ) * (up_lim) ^ 7

f_result = (result_up) - (result_down) 

return, f_result

end