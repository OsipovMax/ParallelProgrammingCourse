gnuplot << EOP
set yrange[0:100]
set terminal jpeg size 640,480
set output "totalTime.jpg"
set title "total time"
plot 'result' u 1:2 w l
EOP
