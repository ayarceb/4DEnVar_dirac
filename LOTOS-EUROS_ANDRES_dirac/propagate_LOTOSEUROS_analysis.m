% Function LOTOSEUROS propagate ensemble 
% Change of Emission and Deposition parameters to generate a parallel LOTOS-EUROS ensemble running
% Andres Yarce Botero-Arjo Segers 24-01-2019

% NOTE : the initial conditions are taken from the folder MEDIAS previously
% generated from the mean of the ensemble member initial conditions


function [le,ld,c5]=propagate_LOTOSEUROS_analysis(nens,le,ld)
nens=1;
%function [ens,le,ld,c5]=propagate_LOTOSEUROS_ensemble(base_emi,nens,le,ld)
%% Ensemble LOTOS-EUROS generation for emission and deposition parameter perturbation
 

 fprintf( '\n                                                        ' );  % nens is the number of ensembles to generate
 fprintf( 'Run the LOTOS-EUROS model from the mean of the ensembles updating\n' );

 emis_Fact=normrnd(6,sqrt(le),[1,nens]);                                  % Gaussian distributed parameter for emissions
 emis_Fact(emis_Fact<0.4)= 0.8;
 %emis_Fact=[1 1 1];
 %depo_Vd=normrnd(1,sqrt(ld),[1,nens]);                                    % Gaussian distributed parameter for emissions
 %depo_Vd(depo_Vd<0.4)=0.8;
 depo_Vd=[1 1 1];
 csvwrite('emis_Fact.csv',emis_Fact);csvwrite('depo_Fact.csv',depo_Vd);  % write csv files with emis and depo factors generated
 
%%  Cycle to generate the ensemble running with the modified parameters




% Part to uncomment for initial condition perturbation
 d='/run/media/dirac/Datos/Dropbox/users/arjo/lotos-euros/4DEnVar_LE/LE_2018_ensembles';
     cd(d); pwd
     % next step remove all the content on the NOx folder
     !rm -rf Boundary_perturbed/* 

 for j=1:nens  % In the Dirac PC, the process manager of the processors generate the wait list
 
 %% Part to change the Emission input folder multiplying for a random factor gaussianly generated 
 
 % ----------- Uncomment this first part: perturb parameter ---------------
     
%      d='/run/media/dirac/Datos/inputdata_Colombia/emissions/EDGAR/v42';
%      cd(d); pwd;
%      % next step remove all the content on the NOx folder
%      !rm -rf NOx/*     
%      % next step copy the content of the backup original NOx folder
%      !cp -a NOx-Original_copia/* NOx
%      cd NOx;pwd;
%      t='/run/media/dirac/Datos/inputdata_Colombia/emissions/EDGAR/v42/NOx';
%      tipo='*.nc';
%      Nombres=get_list_files(t,tipo);pwd;
%      j
%      for i=1:length(Nombres)
%      f=ncread(Nombres{i},'emi_nox')*emis_Fact(j);
%      ncwrite(Nombres{i},'emi_nox',f);
%      end
% 
%      pause(2)
%    
% -------------------------------------------------------------------------

% Uncomment this part to perturb the initial condition

    cd(d); 
        
     % next step copy the content of the backup original NOx folder
     cd Boundary_perturbed
     bc='ens-';
     bc1=num2str(j);
     bcstr=strcat(bc,bc1);
     mkdir(bcstr)
     cd(bcstr)
     
     !cp -a /run/media/dirac/Datos/Dropbox/users/arjo/lotos-euros/4DEnVar_LE/LE_2018_ensembles/Boundary_Conditions_original/* .
     
%      cd Boundary_perturbed;pwd
%      t='/run/media/dirac/Datos/Dropbox/users/arjo/lotos-euros/4DEnVar_LE/LE_2018_ensembles/Boundary_perturbed';
     tipo='*.nc';
     
     nn=50;nm=60;
%      figure
%      no2=ncread('LE_20180626_20180708_conc-3d_20180707.nc','no');
%       plot(squeeze(no2(nn,nm,1,:)));
%       hold on
     
     Nombres=get_list_files(pwd,tipo);pwd;
     
     for i=1:length(Nombres)
     f=ncread(Nombres{i},'no')*emis_Fact(j);
     ncwrite(Nombres{i},'no',f);
     f=ncread(Nombres{i},'no2')*emis_Fact(j);
     ncwrite(Nombres{i},'no2',f);
     f=ncread(Nombres{i},'no3a_f')*emis_Fact(j);
     ncwrite(Nombres{i},'no3a_f',f);
     f=ncread(Nombres{i},'no3a_c')*emis_Fact(j);
     ncwrite(Nombres{i},'no3a_c',f);
     f=ncread(Nombres{i},'n2o5')*emis_Fact(j);
     ncwrite(Nombres{i},'n2o5',f);
     f=ncread(Nombres{i},'hno3')*emis_Fact(j);
     ncwrite(Nombres{i},'hno3',f);
     f=ncread(Nombres{i},'nh3')*emis_Fact(j);
     ncwrite(Nombres{i},'nh3',f);
     end
     
     no2=ncread('LE_20180626_20180708_conc-3d_20180707.nc','no');
%       plot(squeeze(no2(nn,nm,1,:)));
% %      
          
     pause(2)
      d='/run/media/dirac/Datos/Dropbox/users/arjo/lotos-euros/4DEnVar_LE/LE_2018_ensembles';
     cd(d); 

% -------------------------------------------------------------------------

%    %% Run LOTOS-EUROS model
% 
%cd /run/media/dirac/Datos/users/arjo/lotos-euros/Modelo_test_4D  % for when the model to run is in DATOS no in dropbox
cd /run/media/dirac/Datos/Dropbox/users/arjo/lotos-euros/4DEnVar_LE/LE_2018_ensembles
% 
% 
 fileID = fopen('ensemble-settings.rc','w');
 fprintf(fileID, '! current ensemble number:\n');
 fprintf(fileID, 'ensemble.id  : %i\n',j);
 fprintf(fileID, '\n');
 fprintf(fileID, '!Factor Applied to depositions\n');
 fprintf(fileID,'deposition.adhoc.factor.nox  :  %f \n',depo_Vd(j));
 fclose(fileID);

% 
% % check with -h for an overview of the flags in the compiling and
% % submission.  -s argument is for change to run folder and submit inmediatelly
% 
% % print info ...

 fprintf( 'launching ensemble member %i \n', j );
 fprintf( 'emision factor %d \n', emis_Fact(j) );
 fprintf( 'deposition factor %d \n', depo_Vd(j) );
 % compile and start LOTOS-EUROS, run it in the background;
% % write messages to "ensemble-NN.out" file with NN the ensemble member,
% % and also write error messages to the same file
 iret = system( sprintf( './setup_le lotos-euros_analysis.rc -s -n -b > ensemble_analysis-%i.out 2>&1 \n', j ) );
 if iret ~= 0,
     disp('error from calling setup LE');
     break;
 end
 end 
 fprintf( '\n %i Ensemble members submited \n',j );
 %% check if all ensemble member finished

c1='/run/media/dirac/Datos/scratch/projects/LOTOS-EUROS/4DEnVar_ANALYSIS/';
c2='ens-';
c3=num2str(nens);
c4='/run';
c5=strcat(c1,c2,c3,c4);
cd(c5)
a=isfile('le.ok');
while a==0
pause(5)
a=isfile('le.ok');
end
clearvars a


 fprintf( '\n %i Ensemble members finished \n',j );
end


