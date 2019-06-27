#!/bin/bash
usage="####################################################################\n\
This routine collects the analysed fields from the EnKF, and assembles \n\
them into files of name analysisXXX.[ab] analysisICE.uf - these files can be \n\
used as restart files by HYCOM. The unassembled fields have names of type\n\
analysisXXX_procXX.[ab] \n\n\
Usage:   $(basename $0)  restart_template ice_template ensemble_member nproc \n\
Where: \n\
\t\"restart_template\" is an already working restart file  \n\
\t\"ice_template\"     is an already working ice restart file  \n\
\t\"ensemble member\"  is the ensemble member number of this restart file \n\
\t\"nproc\"            is the number of MPI threads used when running the EnKF  \n\
\n\
Example:\n\
\t EnKF_assemble.sh Forecast/ENSrestart2007_289_00.a Forecast/ENSrestart2007_289_00ICE.uf  3 4
\n\n\
NB:\n\
Note that the templates are needed to \"fill in\" what hasn't been updated by the \n\
EnKF in the final analysis file. \n\
####################################################################\n"


prog="$(dirname $0)/EnKF_assemble"


[ $# -ne 4 ] && { echo -e $usage ; exit 1 ; }

# Run EnKF postprocess -- This processes the analysis files,
# puts the files a final analyzed file, having the correct
# order of the restart fields
$prog $@
[ $? -ne 0 ]  && { echo EnKF_assemble failed ; exit 1 ; }

restartbase=$(echo $1 | sed "s/\.[ab]$//")

# The above only reorders the analysis fields into correct order, but
# the ".b" file will lack a header. This copies the header from the
# template file

cmem=$(echo 00$3 | tail -c4)
ppfile=analysis$cmem

myrandom="$RANDOM_$RANDOM"
head -n2 ${restartbase}.b > tmp$myrandom.b
cat ${ppfile}.b >> tmp$myrandom.b
mv tmp$myrandom.b ${ppfile}.b


# We should now have a analysis file named "analysisXXX.[ab]" and analysisICE.uf

echo $?
