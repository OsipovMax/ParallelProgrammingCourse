#!/bin/bash

gnuplot << EOP
set yrange[0:50]
set terminal jpeg size 640,480
set output "MaxTime.jpg"
set title "MaxTime"
set multiplot
set key font ",10" spacing 0.5
plot "512" u 1:3 w l title "512x512" lt rgb 'purple'
set key font ",10" spacing 4.5
plot "1024" u 1:3 w l title "1024x1024" lt rgb 'red'
set key font ",10" spacing 8.5
plot "2048" u 1:3 w l title "2048x2048" lt rgb 'blue'
set key font ",10" spacing 12.5
plot "4096" u 1:3 w l title "4096x4096" lt rgb 'green'
EOP

gnuplot << EOP
set yrange[400:512]
set terminal jpeg size 640,480
set output "Speedup.jpg"
set title "Speedup"
set multiplot
set key font ",10" spacing 0.5
plot "512" u 1:4 w l title "512x512" lt rgb 'purple'
set key font ",10" spacing 4.5
plot "1024" u 1:4 w l title "1024x1024" lt rgb 'red'
set key font ",10" spacing 8.5
plot "2048" u 1:4 w l title "2048x2048" lt rgb 'blue'
set key font ",10" spacing 12.5
plot "4096" u 1:4 w l title "4096x4096" lt rgb 'green'
EOP

gnuplot << EOP
set yrange[0.95:1.05]
set terminal jpeg size 640,480
set output "Efficiency.jpg"
set title "Efficiency"
set multiplot
set key font ",10" spacing 0.5
plot "512" u 1:5 w l title "512x512" lt rgb 'purple'
set key font ",10" spacing 4.5
plot "1024" u 1:5 w l title "1024x1024" lt rgb 'red'
set key font ",10" spacing 8.5
plot "2048" u 1:5 w l title "2048x2048" lt rgb 'blue'
set key font ",10" spacing 12.5
plot "4096" u 1:5 w l title "4096x4096" lt rgb 'green'
EOP

