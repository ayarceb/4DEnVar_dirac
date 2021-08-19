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


#==================================================================================
# 					Modified by user
#==================================================================================

#===Path where the program is running===
mydir='/home/slopezr2/GITHUB/Ensemble_Data_Assimilation/FORTRAN/LOTOS-EUROS'

#===Path LOTOS-EUROS MODEL===
LE='/home/slopezr2/lekf_4DEnVAR/lekf/v3.0.003-beta'

#===Path where netcdf-fortran and netcdf is installed===
NETCDF_FORTRAN_HOME='/home/slopezr2/opt/netcdf/4.7.3/netcdf-fortran/4.5.2'
NETCDF_HOME='/home/slopezr2/opt/netcdf/4.7.3'

#===Run ID====
runid='Proof_Ensemble_Output_V2'

#===Date of simulations====
rm ${LE}/proj/eafit/000/rc/timerange.rc
start_date=20160331
echo 'timerange.start     :  2019-02-00 00:00:00'>>${LE}/proj/eafit/000/rc/timerange.rc
echo 'timerange.end       :  2019-04-00 00:00:00'>>${LE}/proj/eafit/000/rc/timerange.rc

#===Number of Ensembles===
Nens=25

#==================================================================================
# 			      End Modified by user
#==================================================================================

#===Remove all temporal files====
rm ${mydir}/temp/*

#===Path LOTOS-EUROS Ensemble Outputs===
LE_Outputs='/home/slopezr2/scratch/projects/4DEnVAR/'${runid}'/output'

#===Write runid, timerange and Ensemble files====
echo 'run.id             : '${runid}>${LE}/proj/eafit/000/rc/runid.rc



echo 'kf.nmodes             : '${Nens}> ${LE}/proj/eafit/000/rc/N_Ensembles.rc
echo ${LE_Outputs}>>${mydir}/temp/Ensembles.in
echo ${Nens}>>${mydir}/temp/Ensembles.in

#===Run Real State Simulation===
echo 'Running Model Real and Ensemble'

#====Run LOTOS-EUROS MODEL====
cd ${LE}

./launcher

#====Read LE Ensemble outputs====

#==Merging LE outputs for each ensemble member==
#echo 'Merging LE Ensembles'
#cd ${LE_Outputs}
#let "j=0"
#for i in $(ls LE_${runid}_conc-sfc_${start_date}_xi**a.nc)
#	do
#	let "j=j+1"
#	if [ $j -lt 10 ]
#	then
#		ncks -O -h --mk_rec_dmn time LE_${runid}_conc-sfc_${start_date}_xi0${j}a.nc  Merge_x0${j}.nc
#		mv LE_${runid}_conc-sfc_${start_date}_xi0${j}a.nc ..
#		ncrcat -h Merge_x0${j}.nc LE_${runid}_conc-sfc_2*_xi0${j}a.nc Ens_x0${j}.nc

#	else
#			ncks -O -h --mk_rec_dmn time LE_${runid}_conc-sfc_${start_date}_xi${j}a.nc  Merge_x${j}.nc
#		mv LE_${runid}_conc-sfc_${start_date}_xi${j}a.nc ..
#		ncrcat -h Merge_x${j}.nc LE_${runid}_conc-sfc_2*_xi${j}a.nc Ens_x${j}.nc
#	fi
#done 

#==Merging Observations==
#echo 'Merging Observations'

#ncks -O -h --mk_rec_dmn time LE_${runid}_conc-sfc_${start_date}_xb.nc  Merge_xb.nc
#mv LE_${runid}_conc-sfc_${start_date}_xb.nc ..
#ncrcat -h Merge_xb.nc LE_${runid}_conc-sfc_2*_xb.nc y.nc

cd ${mydir}

#==FORTRAN READ NC FILES==
#echo 'Reading LE outputs'
#gfortran -o read_le_ensemble_output  	READ_LE_ENSEMBLE_OUTPUTS.F95 -lblas -llapack  -I${NETCDF_FORTRAN_HOME}/include -I${mydir}/MODULES -L${NETCDF_FORTRAN_HOME}/lib -lnetcdff -Wl,-rpath -Wl,${NETCDF_FORTRAN_HOME}/lib -L${NETCDF_HOME}/lib -lnetcdf -Wl,-rpath -Wl,${NETCDF_HOME}/lib
#./read_le_ensemble_output
#rm read_le_ensemble_output

