#!/bin/bash

exec 1>>report.txt


for mode in 0 1 2 3 4 5
do
\echo -n "$mode "
( ./main _rA _rB _rC $mode; ) 2>> report.txt
\echo -n " "
done

gnuplot << EOP
#set yrange[0:7]
set terminal jpeg size 640,480
set output "temp.jpg"
plot "report.txt" with lines
EOP
