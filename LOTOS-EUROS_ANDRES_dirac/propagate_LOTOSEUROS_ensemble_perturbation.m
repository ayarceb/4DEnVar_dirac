% Function LOTOSEUROS propagate ensemble 
% Change of Emission and Deposition parameters to generate a parallel LOTOS-EUROS ensemble running
% Andres Yarce Botero-Arjo Segers 24-01-2019
%=====================================================================================================
% nens   number of ensembles
% le     deviation perturbation noise to the emission parameter
% ld     deviation perturbation noise to the deposition velocity parameter
% pertu  flag to choice between parameter or initial state perturbation:
%        Parameter perturbation             pertu=1   (with Vd=1 only emissions perturbed)
%        Initial condition  perturbation    pertu=0


function [le,ld,c5]=propagate_LOTOSEUROS_ensemble(nens,le,ld,pertu,meane,infl)

%% Ensemble LOTOS-EUROS generation for emission and deposition parameter perturbation
 
 fprintf( '                                                        ' );   % nens is the number of ensembles to generate
 fprintf( 'Code to run %i ensembles of the LOTOS-EUROS model \n',nens );
 fprintf( 'le is =%i',le );
 fprintf( '                                                        ' ); 
  stdv=le*0.1
 emis_Fact=normrnd(meane,sqrt(stdv),[1,nens]);                                  % Gaussian distributed parameter for emissions and 
 %emis_Fact(end)=0.1;
 %emis_Fact(emis_Fact<0.4)= 0.8;
 %emis_Fact=[3 3 3];
 depo_Vd=normrnd(1,sqrt(ld),[1,nens]);                                    % Gaussian distributed parameter for emissions
 %depo_Vd(depo_Vd<0.4)=0.8;
 %depo_Vd=[1,1,1,1,1];
 %stdv=le*0.1
% xinit_Fact=normrnd(meane,sqrt(stdv),[1,nens]);                              % Gaussian distributed parameter for initial conditions 
 xinit_Fact=normrnd(0,sqrt(1),[1,nens]);                              % Gaussian distributed parameter for initial conditions 
 %xinit_Fact(end)=1;
 %xinit_Fact(end)=0.1;
 csvwrite('emis_Fact.csv',emis_Fact);csvwrite('depo_Fact.csv',depo_Vd);   % write csv files with emis and depo factors generated
 csvwrite('xinit_Fact.csv',xinit_Fact);
%  nuevoArchivo = "xinit_Fact" + le + ".csv";
%  movefile('xinit_Fact.csv',nuevoArchivo)

%%  Cycle to generate the ensemble running with the modified parameters
 for j=1:nens  % In the Dirac PC, the process manager of the processors generate the wait list
 
 %% Part to change the Emission input folder multiplying for a random factor gaussianly generated 
 
 % ----------- Uncomment this first part: perturb parameter ---------------
 
 if pertu==1
     
     d='/run/media/dirac/Datos/inputdata_Colombia/emissions/EDGAR/v42';
     cd(d); pwd;
     % next step remove all the content on the NOx folder
     !rm -rf NOx/*     
     % next step copy the content of the backup original NOx folder
     !cp -a NOx-Original_copia/* NOx
     cd NOx;pwd;
     t='/run/media/dirac/Datos/inputdata_Colombia/emissions/EDGAR/v42/NOx';
     tipo='*.nc';
     Nombres=get_list_files(t,tipo);pwd;
     j
     for i=1:length(Nombres)
     f=ncread(Nombres{i},'emi_nox')*emis_Fact(j);
     ncwrite(Nombres{i},'emi_nox',f);
     end
     sprintf('Parameters perturbed')
     pause(2)
   
% -------------------------------------------------------------------------

% Uncomment this part to perturb the initial condition
elseif pertu==0

    % Part to uncomment for initial condition perturbation
     d='/run/media/dirac/Datos/Dropbox/users/arjo/lotos-euros/4DEnVar_LE/LE_2018_ensembles';
     cd(d);
     % next step remove all the content on the NOx folder
 %    !rm -rf Boundary_perturbed/* 
    
     cd(d); 
        
     % next step copy the content of the backup original NOx folder
     cd Boundary_perturbed
     bc='ens-';
     bc1=num2str(j);
     bcstr=strcat(bc,bc1);
     mkdir(bcstr)
     cd(bcstr)
     
     !cp -a /run/media/dirac/Datos/Dropbox/users/arjo/lotos-euros/4DEnVar_LE/LE_2018_ensembles/Boundary_Conditions_original/* .
     
%    cd Boundary_perturbed;pwd
%    t='/run/media/dirac/Datos/Dropbox/users/arjo/lotos-euros/4DEnVar_LE/LE_2018_ensembles/Boundary_perturbed';
     tipo='*.nc';
     
     nn=50;nm=60;
  
     no2_orig=ncread('LE_20180626_20180708_conc-3d_20180707.nc','no2');
    
     
     % In this verison the perturbation is performed
     % x^pert=x(0)+0.05*normrand*.x(0)
     
     Nombres=get_list_files(pwd,tipo);pwd;
     j
     for i=1:length(Nombres)
     f=ncread(Nombres{i},'no')+infl*xinit_Fact(j).*ncread(Nombres{i},'no');
     ncwrite(Nombres{i},'no',f);
     f=ncread(Nombres{i},'no2')+infl*xinit_Fact(j).*ncread(Nombres{i},'no2');
     ncwrite(Nombres{i},'no2',f);
     f=ncread(Nombres{i},'no3a_f')+infl*xinit_Fact(j).*ncread(Nombres{i},'no3a_f');
     ncwrite(Nombres{i},'no3a_f',f);
     f=ncread(Nombres{i},'no3a_c')+infl*xinit_Fact(j).*ncread(Nombres{i},'no3a_c');
     ncwrite(Nombres{i},'no3a_c',f);
     f=ncread(Nombres{i},'n2o5')+infl*xinit_Fact(j).*ncread(Nombres{i},'n2o5');
     ncwrite(Nombres{i},'n2o5',f);
     f=ncread(Nombres{i},'hno3')+infl*xinit_Fact(j).*ncread(Nombres{i},'hno3');
     ncwrite(Nombres{i},'hno3',f);
     f=ncread(Nombres{i},'nh3')+infl*xinit_Fact(j).*ncread(Nombres{i},'nh3');
     ncwrite(Nombres{i},'nh3',f);
     end
     
     
%%  ==== part to plot to check if boundary initial condition is already perturbed ======     
     no2=ncread('LE_20180626_20180708_conc-3d_20180707.nc','no2');
    figure
    a1=plot(squeeze(no2_orig(nn,nm,1,:)));
    hold on
    a2=plot(squeeze(no2(nn,nm,1,:)));
    title('Plot to test that boundary condition were perturbed')  
    legend([a1 a2],'original initial conditions','perturbed initial conditions')
     ylabel('Concentration [mol no_2 mol-1 air]','FontSize',10,'Interpreter','latex')
     xlabel('time (hours)','FontSize',10,'Interpreter','latex')
    %%     
     sprintf('Initial conditions perturbed')     
     pause(2)
     d='/run/media/dirac/Datos/Dropbox/users/arjo/lotos-euros/4DEnVar_LE/LE_2018_ensembles';
     j
     cd(d); 
     
end
% -------------------------------------------------------------------------

%% ======================= Run LOTOS-EUROS model ==========================
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
pwd
 fprintf( 'launching ensemble member %i \n', j );
 fprintf( 'emision factor %d \n', emis_Fact(j) );
 fprintf( 'deposition factor %d \n', depo_Vd(j) );
 fprintf( 'initial condition factor %d \n', xinit_Fact(j) );
 % compile and start LOTOS-EUROS, run it in the background;
% % write messages to "ensemble-NN.out" file with NN the ensemble member,
% % and also write error messages to the same file
 iret = system( sprintf( './setup_le lotos-euros.rc -s -n -b > ensemble-%i.out 2>&1 \n', j ) );
 if iret ~= 0,
     disp('error from calling setup LE');
     break;
 end
 end 
 fprintf( '\n %i Ensemble members submited \n',j );
 %% check if all ensemble member finished

c1='/run/media/dirac/Datos/scratch/projects/LOTOS-EUROS/4DEnVar/';
c2='ens-';
c3=num2str(nens);
c4='/run';
c5=strcat(c1,c2,c3,c4);
cd(c5)
pwd
a=isfile('le.ok');
while a==0
pause(5)
a=isfile('le.ok');
end
clearvars a


 fprintf( '\n %i Ensemble members finished \n',j );
end


