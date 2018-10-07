#!/bin/bash

gnuplot << EOP
set yrange[0.5:1.5]
set terminal jpeg size 640,480
set output "averaging.jpg"
plot "rep.txt" with lines
EOP
