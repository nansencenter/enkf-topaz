#!/bin/bash

#KAL -- This script sets up the forecast for the EnKF. Input 
#KAL -- is the "base" of the analysis file name

[ $# -ne 2 ] && { echo "No base file name and target supplied" ; exit 1 ;}

# HYCOM files
if ! ls ${1}[0-9][0-9][0-9].[ab] > /dev/null  ; then
   echo "Could not find files with base $1"
   exit 1
fi
#ls ${1}[0-9][0-9][0-9].[ab] 



for i in ${1}[0-9][0-9][0-9].[ab] ; do

   numpart=$(echo $i | sed "s/.*$1//" | sed "s/\..*//")
   abpart=$(echo $i | sed "s/.*\.//")
   #echo $i $numpart $abpart

   if [ $numpart -gt 1 ] ; then
      #echo "yes"
      tailpart="_mem$numpart.$abpart"
   else
      tailpart=".$abpart"
   fi


   newfile=${2}$tailpart
   echo "$newfile -> $i"
   ln -sf $i $newfile  
done

# Ice file
if [ -f ${1}ICE.uf ] ; then
   finalname=${2}ICE.uf
   echo "$finalname -- >  ${1}ICE.uf "
   ln -sf ${1}ICE.uf  $finalname ; 
fi
