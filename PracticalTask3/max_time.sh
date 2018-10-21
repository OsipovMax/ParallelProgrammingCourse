gnuplot << EOP
set yrange[0:10]
set terminal jpeg size 640,480
set output "maxTime.jpg"
set title "max Time"
plot 'result' u 1:3 w l
EOP
