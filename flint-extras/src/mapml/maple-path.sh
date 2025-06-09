#!/bin/tcsh

\mv mapml.mpl.bak mapml.mpl 

sed -e "s|path1|$1|" -i.bak mapml.mpl  

sed -e "s|path2|$2|" -i.bak2  mapml.mpl  

\rm mapml.mpl.bak2 