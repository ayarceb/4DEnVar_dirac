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


# remove folders old runs in scratch folder

rm -rf /run/media/dirac/Datos/scratch/projects/4DEnVAR_inner
mv /run/media/dirac/Datos/scratch/projects/4DEnVAR/40_ensembles_15_enero/LE_40_ensembles*.nc /run/media/dirac/Datos/scratch/projects/4DEnVAR/40_ensembles_15_enero/output
rm -rf /run/media/dirac/Datos/scratch/projects/4DEnVAR/40_ensembles_15_enero/output/Ens_x*.nc
rm -rf /run/media/dirac/Datos/scratch/projects/4DEnVAR/40_ensembles_15_enero/output/Merge_x*.nc
rm /home/dirac/4DEnVar_dirac/LOTOS-EUROS_ANDRES_dirac/data/*

#==================================================================================
# 					Modified by user
#==================================================================================


#===Path where the program is running===
mydir='/home/dirac/4DEnVar_dirac/LOTOS-EUROS_ANDRES_dirac'

#===Path LOTOS-EUROS MODEL (OJO carpeta de LEKF)===
LE='/run/media/dirac/Datos/Reciente_Dropbox/users/arjo/lotos-euros/Version_WRF_04_2020/lekf_4DEnVAR/lekf/v3.0.003-beta'

#===Path where netcdf-fortran and netcdf is installed===
OPT=${HOME}'/opt'
NETCDF_FORTRAN_HOME='/usr/lib64'
NETCDF_HOME='/lib64/gfortran/modules'
#NETCDF_HOME=${OPT}'/home/dirac/miniconda3/pkgs/libnetcdf-4.7.3-nompi_h9f9fd6a_101'

#===Run ID====
runid='40_ensembles_15_enero'

#===Date of simulations====
if [ -f ${LE}/proj/eafit/000/rc/timerange.rc ]
then
	rm ${LE}/proj/eafit/000/rc/timerange.rc
	
fi	

if [ -f ${LE}/proj/eafit/000/rc/timerange_inner.rc ]
then
	rm ${LE}/proj/eafit/000/rc/timerange_inner.rc
	
fi

if [ -f ${mydir}/DATA_4DEnVAR/parameters.in ]
then
	rm ${mydir}/DATA_4DEnVAR/parameters.in
	
fi


start_date=20190115 # Recordar modificar
#=====Days to be simulated======

days_simulation=5   # Recordar modificar


#===== Data assimilation parameters=====
dawindows=1    #Assimilation Window
st_ass=1   #Assimilation window start
inner=15    #Number of inner step
eps=0.0002  #Tolerance Inner Loop
sigma_obs=0.01  #Sigma Observations

echo ${dawindows}>>${mydir}/DATA_4DEnVAR/parameters.in
echo ${st_ass}>>${mydir}/DATA_4DEnVAR/parameters.in
echo ${inner}>>${mydir}/DATA_4DEnVAR/parameters.in
echo ${eps}>>${mydir}/DATA_4DEnVAR/parameters.in
echo ${sigma_obs}>>${mydir}/DATA_4DEnVAR/parameters.in


echo 'timerange.start     :  2019-01-19 00:00:00'>>${LE}/proj/eafit/000/rc/timerange.rc
echo 'timerange.end       :  2019-01-20 00:00:00'>>${LE}/proj/eafit/000/rc/timerange.rc

echo 'timerange.start     :  2019-01-16 00:00:00'>>${LE}/proj/eafit/000/rc/timerange_inner.rc
echo 'timerange.end       :  2019-01-17 00:00:00'>>${LE}/proj/eafit/000/rc/timerange_inner.rc



echo ${start_date}>${mydir}/DATA_4DEnVAR/startdate.in
echo ${days_simulation}>>${mydir}/DATA_4DEnVAR/startdate.in

echo ${runid}>${mydir}/DATA_4DEnVAR/runid.in
#===Number of Ensembles===
Nens=40


#===Parameter rho  esquema de optimizacion====
rho=0.1


echo ${rho}>${mydir}/DATA_4DEnVAR/rho.in


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

#cd ${LE}
#./launcher

#====Read LE Ensemble outputs====

#==Merging LE DC for each ensemble member ==
echo 'Merging LE Ensembles DC'
cd ${LE_Outputs}
let "j=0"
for i in $(ls LE_${runid}_dc_${start_date}_xi**a.nc)
	do
	let "j=j+1"
	if [ $j -lt 10 ]
	then
		ncks -O -h --mk_rec_dmn time LE_${runid}_dc_${start_date}_xi0${j}a.nc  Merge_x0${j}.nc
		mv LE_${runid}_dc_${start_date}_xi0${j}a.nc ..
		ncrcat -O -h Merge_x0${j}.nc LE_${runid}_dc_2*_xi0${j}a.nc Ens_x0${j}.nc

	else
			ncks -O -h --mk_rec_dmn time LE_${runid}_dc_${start_date}_xi${j}a.nc  Merge_x${j}.nc
		mv LE_${runid}_dc_${start_date}_xi${j}a.nc ..
		ncrcat -O -h Merge_x${j}.nc LE_${runid}_dc_2*_xi${j}a.nc Ens_x${j}.nc
	fi
	
	echo $j
	echo $i
	

done 




#==FORTRAN READ NC FILES==

cd ${mydir}
echo 'Reading LE outputs'
gfortran -o read_le_ensemble_output  READ_LE_ENSEMBLE_OUTPUTS_V2.F95 -I${mydir}/MODULES -lblas -llapack  -I${NETCDF_FORTRAN_HOME}/include -I${mydir}/MODULES -L${NETCDF_FORTRAN_HOME}/lib -lnetcdff -Wl,-rpath -Wl,${NETCDF_FORTRAN_HOME}/lib -L${NETCDF_HOME}/lib -lnetcdf -Wl,-rpath -Wl,${NETCDF_HOME}/lib  -I/usr/lib64/gfortran/modules -ffree-line-length-none -ffixed-line-length-none -fimplicit-none
./read_le_ensemble_output
rm read_le_ensemble_output
#==FORTRAN 4DEnVAR Method


gfortran -c module_EnKF.F95  -lblas -llapack  -I${NETCDF_HOME} -lnetcdff -Wl,-rpath, -lnetcdff -Wl,-rpath -ffree-line-length-none -ffixed-line-length-none -fimplicit-none -g -fcheck=all -Wall -fbacktrace



gfortran -c 4DEnVAR_Method_V3.F95  -lblas -llapack  -I${NETCDF_FORTRAN_HOME}/include  -L${NETCDF_FORTRAN_HOME}/lib -lnetcdff -Wl,-rpath -Wl,${NETCDF_FORTRAN_HOME}/lib -L${NETCDF_HOME}/lib -lnetcdf -Wl,-rpath -Wl,${NETCDF_HOME}/lib   -ffree-line-length-none -ffixed-line-length-none -fimplicit-none -g -fcheck=all -Wall -fbacktrace

gfortran -o 4DEnVAR_method 4DEnVAR_Method_V3.o modulo_distribucion_normal.o module_matrix.o module_EnKF.o  -I${mydir}/MODULES -lblas -llapack  -I${NETCDF_FORTRAN_HOME}/include -L${NETCDF_FORTRAN_HOME}/lib -lnetcdff -Wl,-rpath -Wl,${NETCDF_FORTRAN_HOME}/lib -L${NETCDF_HOME}/lib -lnetcdf -Wl,-rpath -Wl,${NETCDF_HOME}/lib  -ffree-line-length-none -ffixed-line-length-none -fimplicit-none -g -fcheck=all -Wall -fbacktrace


./4DEnVAR_method
rm 4DEnVAR_method



