#! /bin/bash

#============== Script 4DEnVar for the LOTOS-EUROS  =======================
# Andres Yarce Botero-Santiago Lopez-Restrepo
# Based in the Paper: "A review of operational methods of variational 
# and ensemble-variational data assimilation" R.N. Bannister
#
#==========================================================================
#  The true run was created with    Emis_Fact=0.5
#                                   Depo_Vd=1
#  stored in no2_column_True.mat  (in the code with recent the true is one ensamble member not perturbed)
#==========================================================================
#  First background created with    Emis_Fact= 1.5
#                                   Depo_Vd=1
#==========================================================================
#  To choice:
#
#  Parameter perturbation             pertu=1   (with Vd=1 only emissions perturbed)
#  Initial condition  perturbation    pertu=0
# =========================================================================

#===Path where the program is running===
mydir='/home/kalman/GITHUB/Ensemble_Data_Assimilation/FORTRAN/LOTOS-EUROS'
#===Path LOTOS-EUROS MODEL===
LE='/home/kalman/GITHUB/dataAssimilation/src/LOTOS-EUROS_CODE/arjo/lotos-euros/LOTOS-EUROSV2.2/lotos-euros/4DEnVAR'

#===Path where the Input Emission are===
emioriginal='/disk/kalman/inputdata_Colombia/emissions/EDGAR/v4.3/v4.3.2/data/NOx_Original/2012'
emi='/disk/kalman/inputdata_Colombia/emissions/EDGAR/v4.3/v4.3.2/data/NOx/2012'
#===Ojo igual que emi pero sin el primer / para evitar probleas en fortran===

emi2='disk/kalman/inputdata_Colombia/emissions/EDGAR/v4.3/v4.3.2/data/NOx/2012'
#===Path where netcdf-fortran and netcdf is installed

NETCDF_FORTRAN_HOME='/home/kalman/opt/netcdf/4.4.1.1/netcdf-fortran/4.4.4/'

NETCDF_HOME='/home/kalman/opt/netcdf/4.4.1.1'


#===Remove all previous run====
rm -rf /disk/kalman/scratch/projects/4DEnVar

#===Remove all temporal files====
rm ${mydir}/temp/*
#===Run Real State Simulation===
echo 'Running Real State'

#===Remove all the content on the Nox Folder===
rm -rf ${emi}
#== next step copy the content of the backup original NOx folder
cp -a ${emioriginal}  ${emi}

#==Name of all the files in the path and write it in /temp/files.dat===
cd ${emi}
echo ${emi} >> ${mydir}/temp/files.dat
for i in $(ls -C1)
do
echo $i >> ${mydir}/temp/files.dat
done 

cd ${mydir}
#===Compile and Execute FORTRAN run_lotos_euros_real to modifie emissions for Real State===

gfortran -o execute_run_lotos_euros_real  RUN_LOTOS-EUROS_REAL.F95 -lblas -llapack  -I${NETCDF_FORTRAN_HOME}/include -I${mydir}/MODULES -L${NETCDF_FORTRAN_HOME}/lib -lnetcdff -Wl,-rpath -Wl,${NETCDF_FORTRAN_HOME}/lib -L${NETCDF_HOME}/lib -lnetcdf -Wl,-rpath -Wl,${NETCDF_HOME}/lib

./execute_run_lotos_euros_real
rm execute_run_lotos_euros_real

#====Run LOTOS-EUROS MODEL===

#=Write Ensemble-id 0==
echo  'ensemble.id		:  15'> ${LE}/ensemble-settings.rc
cd ${LE}
${LE}/setup_le ${LE}/lotos-euros-2.rc -s -b > ${LE}/ensemble-15.out 2>&1 

cd ${mydir}
echo 'LE Compiled and Runnig'

