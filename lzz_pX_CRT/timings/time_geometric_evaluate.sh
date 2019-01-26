#!/bin/sh

filename=../thresholds_geometric_evaluate_FFT.h
rm -f $filename
echo "#ifndef __THRESHOLDS_GEOMETRIC_EVALUATE_FFT_H" >> $filename
echo "#define __THRESHOLDS_GEOMETRIC_EVALUATE_FFT_H" >> $filename
awk -f find_thresholds.awk -v p=23068673 -v type=EVALUATE_SMALL time_geometric_evaluate.res >> $filename
awk -f find_thresholds.awk -v p=288230376151711813 -v type=EVALUATE_LARGE time_geometric_evaluate.res >> $filename
awk -f find_thresholds.awk -v p=0 -v type=EVALUATE_FFT time_geometric_evaluate.res >> $filename
echo "#endif" >> $filename
