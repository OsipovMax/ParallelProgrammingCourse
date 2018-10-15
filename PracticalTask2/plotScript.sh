#!/bin/bash

for fg in 0 1
do
for mode in 0 1 2 
do
for size in  700 1000 1500 2000 2500
do
(./genMatrix f $size $size _rA ; ./genMatrix f $size $size _rB; ./PPT_prac _rA _rB _rC $mode $fg)
done 
done 
done

