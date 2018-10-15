#!/bin/bash

gnuplot << EOP
set yrange[0:3e+08]
set terminal jpeg size 640,480
set output "L1.jpg"
set title "L1 DCM"
set multiplot
set key font ",10" spacing 0.5
plot "L1_DCM0" w l lt rgb 'purple'
set key font ",10" spacing 1.75
plot "L1_DCM1" w l lt rgb 'red'
set key font ",10" spacing 3
plot "L1_DCM2" w l lt rgb 'blue'
EOP

gnuplot << EOP
set yrange[0:3e+07]
set terminal jpeg size 640,480
set output "L2.jpg"
set title "L2 DCM"
set multiplot
set key font ",10" spacing 0.5
plot "L2_DCM0" w l lt rgb 'purple'
set key font ",10" spacing 1.75
plot "L2_DCM1" w l lt rgb 'red'
set key font ",10" spacing 3
plot "L2_DCM2" w l lt rgb 'blue'
EOP

gnuplot << EOP
set yrange[350:480]
set terminal jpeg size 640,480
set output "Mflops.jpg"
set title "Mflops"
set multiplot
set key font ",10" spacing 0.5
plot "Mflops0" w l lt rgb 'purple'
set key font ",10" spacing 1.75
plot "Mflops1" w l lt rgb 'red'
set key font ",10" spacing 3
plot "Mflops2" w l lt rgb 'blue'
EOP

gnuplot << EOP
set yrange[3e+9:2e+11]
set terminal jpeg size 640,480
set output "RealCycl.jpg"
set title "RealCycl"
set multiplot
set key font ",10" spacing 0.5
plot "Cycles0" w l lt rgb 'purple'
set key font ",10" spacing 1.75
plot "Cycles1" w l lt rgb 'red'
set key font ",10" spacing 3
plot "Cycles2" w l lt rgb 'blue'
EOP

gnuplot << EOP
set yrange[0:45]
set terminal jpeg size 640,480
set output "RealTime.jpeg"
set title "RealTime"
set multiplot
set key font ",10" spacing 0.5
plot "RealTime0" w l lt rgb 'purple'
set key font ",10" spacing 1.75
plot "RealTime1" w l lt rgb 'red'
set key font ",10" spacing 3
plot "RealTime2" w l lt rgb 'blue'
EOP
