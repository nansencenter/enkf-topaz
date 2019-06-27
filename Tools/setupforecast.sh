#!/bin/bash

#KAL -- This script sets up the forecast for the EnKF. Input 
#KAL -- is the "base" of the analysis file name

[ $# -ne 1 ] && { echo "No base file name supplied" ; exit 1 ;}

# HYCOM files
if ! ls ${1}.[ab] > /dev/null || ! ls ${1}_mem???.[ab] > /dev/null  ; then
   echo "Could not find files with base $1"
   exit 1
fi

for i in ${1}.[ab] ${1}_mem???.[ab] ; do
   tailpart1=$(echo $i | tail -c9 )
   tailpart2=$(echo $i | tail -c9 | cut -c1-3)
   post=$(echo $i | sed "s/.*\.//")

   if [ ! "$tailpart2" == "mem" ] ; then
      tailpart1="mem001.${post}"
   fi

   finalname=$( echo $tailpart1 | sed "s/mem/forecast/")
   echo "$finalname -- > $i"
   ln -s $i $finalname  
done

# Ice file
if [ -f ${1}ICE.uf ] ; then
   finalname=forecastICE.uf
   echo "$finalname -- > ${1}ICE.uf"
   ln -s ${1}ICE.uf  forecastICE.uf ; 
fi
