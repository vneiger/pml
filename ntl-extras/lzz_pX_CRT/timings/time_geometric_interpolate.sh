#!/bin/sh

filename=../thresholds_geometric_interpolate_FFT.h
rm -f $filename
echo "#ifndef __THRESHOLDS_GEOMETRIC_INTERPOLATE_FFT_H" >> $filename
echo "#define __THRESHOLDS_GEOMETRIC_INTERPOLATE_FFT_H" >> $filename
awk -f find_thresholds.awk -v p=23068673 -v type=INTERPOLATE_SMALL time_geometric_interpolate.res >> $filename
awk -f find_thresholds.awk -v p=288230376151711813 -v type=INTERPOLATE_LARGE time_geometric_interpolate.res >> $filename
awk -f find_thresholds.awk -v p=0 -v type=INTERPOLATE_FFT time_geometric_interpolate.res >> $filename
echo "#endif" >> $filename
