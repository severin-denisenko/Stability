set pm3d map
set pm3d map interpolate 0,0
set xlabel 'a, au'
set ylabel 'e'
splot "result.dat" using ($1/100 * 8700 / 2 ):(0.5 - $2/200):3 matrix w pm3d
pause -1
