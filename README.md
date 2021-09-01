# 4DEnVar_dirac



Instruction for the operation of the code 4DENVAR_LOTOS-EUROS


The code assimilates TROPOMI NO2 Column column as observations to modify the dc correction factors estimating it through the assimilation procedure. The dc correction factors are the parameter values that multiply the current emission inventories


The GITHUB path of the working code is https://github.com/ayarceb/4DEnVar_dirac


	Main file:  4D_EnVAR_LE_main_V4.sh (Here put number of ensembles, Simulation date, )

	In folder DATA_4DEnVar the file parameteres.in   (Assimilation Window, Spread, Sigma Errors)
                                       

	Launch file for the ensemble  lekf.rc   (Generate the ensembles, mod # ensembles) (lotos)


    	launck lekf_inner



4DEnVar_method.f95 (read parameters from the folder DATA_4DEnVar)
module_enkf.f95 (Create new dc files in the folder data)
 


----------------------------------------------------------------------------------------------------


Para modificar el dc noise se deben modificar con atención dos archivos en dirac

Uno es el READ_LE_ENSEMBLE_OUTPUTS_V2.F95   (/home/dirac/4DEnVar_dirac/LOTOS-EUROS_ANDRES_dirac)

y otro el lekf_noise_solo_ciudades.90      (/run/media/dirac/Datos/Reciente_Dropbox/users/arjo/lotos-euros/Version_WRF_04_2020/lekf_4DEnVAR/lekf/v3.0.003-beta/base/002/src)
