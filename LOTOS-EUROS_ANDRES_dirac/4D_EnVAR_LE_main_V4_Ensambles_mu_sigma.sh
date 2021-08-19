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
mydir='/run/media/dirac/Datos/Reciente_Dropbox/users/arjo/lotos-euros/Repositorio_Personal_Slopez/Personal/FORTRAN/LOTOS-EUROS_ANDRES'

#===Path LOTOS-EUROS MODEL (OJO carpeta de LEKF)===
LE='/run/media/dirac/Datos/Reciente_Dropbox/users/arjo/lotos-euros/Version_WRF_04_2020/lekf_4DEnVAR/lekf/v3.0.003-beta'

#===Path where netcdf-fortran and netcdf is installed===
OPT=${HOME}'/opt'
NETCDF_FORTRAN_HOME='/usr/lib64'
#NETCDF_HOME=${OPT}'/netcdf/4.4.0'
NETCDF_HOME=${OPT}'/home/dirac/miniconda3/pkgs/libnetcdf-4.7.3-nompi_h9f9fd6a_101'
#===Run ID====
runid='Proof_Ensemble_Expermients_mu_1_disturbf_1'

#===Date of simulations====
if [ -f ${LE}/proj/eafit/000/rc/timerange.rc ]
then
	rm ${LE}/proj/eafit/000/rc/timerange.rc
	
fi	
start_date=20190202
#=====Days to be simulated======

days_simulation=3


echo 'timerange.start     :  2019-02-01 00:00:00'>>${LE}/proj/eafit/000/rc/timerange.rc
echo 'timerange.end       :  2019-02-04 00:00:00'>>${LE}/proj/eafit/000/rc/timerange.rc

echo ${start_date}>${mydir}/DATA_4DEnVAR/startdate.in
echo ${days_simulation}>>${mydir}/DATA_4DEnVAR/startdate.in

echo ${runid}>${mydir}/DATA_4DEnVAR/runid.in
#===Number of Ensembles===
Nens=10


#===Remove all temporal files====
if [ -d ${mydir}/temp ]
then
	rm ${mydir}/temp/*
fi
#===Path LOTOS-EUROS Ensemble Outputs===
LE_Outputs='/run/media/dirac/Datos/scratch/projects/4DEnVAR/'${runid}'/output'

#==================================================================================
# 			      End Modified by user
#==================================================================================


#===Write runid, timerange and Ensemble files====
echo 'run.id             : '${runid}>${LE}/proj/eafit/000/rc/runid.rc



echo 'kf.nmodes             : '${Nens}> ${LE}/proj/eafit/000/rc/N_Ensembles.rc
echo ${LE_Outputs}>>${mydir}/temp/Ensembles.in
echo ${Nens}>>${mydir}/temp/Ensembles.in


if [ -f ${LE}/proj/eafit/000/rc/Restart.rc ]
then
	rm ${LE}/proj/eafit/000/rc/Restart.rc
fi
	
echo 'kf.restart.path               : /run/media/dirac/Datos/scratch/projects/4DEnVAR/'${runid}'/restart'>>${LE}/proj/eafit/000/rc/Restart.rc
echo 'kf.restart.key                :  model=LEKF;expid='${runid}>>${LE}/proj/eafit/000/rc/Restart.rc



#===Run Real State Simulation===
echo 'Running Model Real and Ensemble'

#====Run LOTOS-EUROS MODEL====
cd ${LE}

#./launcher

#====Read LE Ensemble outputs====

#==Merging LE DC for each ensemble member ==
#echo 'Merging LE Ensembles DC'
#cd ${LE_Outputs}
#let "j=0"
#for i in $(ls LE_${runid}_dc_${start_date}_xi**a.nc)
#	do
#	let "j=j+1"
#	if [ $j -lt 10 ]
#	then
#		ncks -O -h --mk_rec_dmn time LE_${runid}_dc_${start_date}_xi0${j}a.nc  Merge_x0${j}.nc
#		mv LE_${runid}_dc_${start_date}_xi0${j}a.nc ..
#		ncrcat -O -h Merge_x0${j}.nc LE_${runid}_dc_2*_xi0${j}a.nc Ens_x0${j}.nc
#	else
#			ncks -O -h --mk_rec_dmn time LE_${runid}_dc_${start_date}_xi${j}a.nc  Merge_x${j}.nc
#		mv LE_${runid}_dc_${start_date}_xi${j}a.nc ..
#		ncrcat -O -h Merge_x${j}.nc LE_${runid}_dc_2*_xi${j}a.nc Ens_x${j}.nc
#	fi
#done 

##==Merging LE outputs for each ensemble member, con CDO y revisar formato salidas satelite simulado==
#echo 'Merging LE Ensembles Outputs'

#let "j=0"
#for i in $(ls LE_${runid}_column_${start_date}_xi**a.nc)
#	do
#	let "j=j+1"
#	if [ $j -lt 10 ]
#	then
#		ncks -O -h --mk_rec_dmn time LE_${runid}_column_${start_date}_xi0${j}a.nc  Y_Merge_x0${j}.nc
#		mv LE_${runid}_column_${start_date}_xi0${j}a.nc ..
#		ncrcat -O -h Y_Merge_x0${j}.nc LE_${runid}_column_2*_xi0${j}a.nc Y_Ens_x0${j}.nc
#
#	else
#			ncks -O -h --mk_rec_dmn time LE_${runid}_column_${start_date}_xi${j}a.nc  Y_Merge_x${j}.nc
#		mv LE_${runid}_column_${start_date}_xi${j}a.nc ..
#		ncrcat -O -h Y_Merge_x${j}.nc LE_${runid}_column_2*_xi${j}a.nc Y_Ens_x${j}.nc
#	fi
#done 




###==Merging Real State, No se necesita==
#echo 'Merging Real State'

#ncks -O -h --mk_rec_dmn time LE_${runid}_dc_${start_date}_xb.nc  Merge_xb.nc
#mv LE_${runid}_dc_${start_date}_xb.nc ..
#ncrcat -O -h Merge_xb.nc LE_${runid}_dc_2*_xb.nc X_real.nc

##==Merging Observations, Cambiar para formato nuevas salidas de satelite hacer con CDO, ya no es Xb si no el archivo de las observaciones reales==
#echo 'Merging Observations'

#ncks -O -h --mk_rec_dmn time LE_${runid}_column_${start_date}_xb.nc  Y_Merge_xb.nc
#mv LE_${runid}_column_${start_date}_xb.nc ..
#ncrcat -O -h Y_Merge_xb.nc LE_${runid}_column_2*_xb.nc Y.nc



#==FORTRAN READ NC FILES==

cd ${mydir}
echo 'Reading LE outputs'
gfortran -o read_le_ensemble_output  READ_LE_ENSEMBLE_OUTPUTS_V2.F95 -I${mydir}/MODULES -lblas -llapack  -I${NETCDF_FORTRAN_HOME}/include -I${mydir}/MODULES -L${NETCDF_FORTRAN_HOME}/lib -lnetcdff -Wl,-rpath -#Wl,${NETCDF_FORTRAN_HOME}/lib -L${NETCDF_HOME}/lib -lnetcdf -Wl,-rpath -Wl,${NETCDF_HOME}/lib  -I/usr/lib64/gfortran/modules
./read_le_ensemble_output
rm read_le_ensemble_output
#==FORTRAN 4DEnVAR Method



gfortran -c 4DEnVAR_Method_V2.F95 modulo_distribucion_normal.F95 module_matrix.F95 module_EnKF.F95  -lblas -llapack  -I${NETCDF_FORTRAN_HOME}/include -L${NETCDF_FORTRAN_HOME}/lib -lnetcdff -Wl,-rpath -Wl,${NETCDF_FORTRAN_HOME}/lib -L${NETCDF_HOME}/lib -lnetcdf -Wl,-rpath -Wl,${NETCDF_HOME}/lib -I/usr/lib64/gfortran/modules

gfortran -o 4DEnVAR_method 4DEnVAR_Method_V2.o modulo_distribucion_normal.o module_matrix.o module_EnKF.o  -lblas -llapack  -I${NETCDF_FORTRAN_HOME}/include -L${NETCDF_FORTRAN_HOME}/lib -lnetcdff -Wl,-rpath -Wl,${NETCDF_FORTRAN_HOME}/lib -L${NETCDF_HOME}/lib -lnetcdf -Wl,-rpath -Wl,${NETCDF_HOME}/lib -I/usr/lib64/gfortran/modules


./4DEnVAR_method
rm 4DEnVAR_method



