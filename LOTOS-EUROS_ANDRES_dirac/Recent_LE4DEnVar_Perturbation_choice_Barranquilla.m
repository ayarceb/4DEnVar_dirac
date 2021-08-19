%============== Script 4DEnVar for the LOTOS-EUROS  =======================
% Andres Yarce Botero
% Based in the Paper: "A review of operational methods of variational 
% and ensemble-variational data assimilation" R.N. Bannister
%
%==========================================================================
%  The true run was created with    Emis_Fact=0.5
%                                   Depo_Vd=1
%  stored in no2_column_True.mat  (in the code with recent the true is one ensamble member not perturbed)
%==========================================================================
%  First background created with    Emis_Fact= 1.5
%                                   Depo_Vd=1
%==========================================================================
%  To choice:
%
%  Parameter perturbation             pertu=1   (with Vd=1 only emissions perturbed)
%  Initial condition  perturbation    pertu=0
% =========================================================================
clc;clear all;close all
 system('rm -rf /run/media/dirac/Datos/scratch/projects/LOTOS-EUROS/4DEnVar')

%% Results from previous runs

%load no2_column_True.mat                 % Load Truth
%load run_Perturbation_Initial_state.mat  %load the 25 ensemble init state
%load run_Perturb_Parameter.mat           %load the 25 ensemble parameter p

%% ===============  Some general info for the plots  ======================
x0=3000;y0=600;             % size and position of the graphic to generates
width=1800;height=700;      
Fx=15;                     % Fontsize Xlabel Ylabel
Ft=20;                     % Fontsize title
k=11;                      % Ensemble color grey for ensemble members
%% =========================== Parameters =================================
% model parameters


kkk=1;                  % this index is used when the experiment of inflation is performed 
gridx=11;               % Lateral dimension from the output LE grid
gridy=11;                                               
gridz=5;                % (CONSIDERAR EJE Z?) si por trabajar con los estados no salidas
pertu=0;                % Modify (0-initial conditions) or (1-parameters)
infl=0.8;               % Inflation factor
nest=gridx*gridy;                 % Number of outputs column 11 x 11 grid 
tsim=96;                          % Time simulation steps for 4 days

% [no2_column_T]=mat2vect(no2_column_True,tsim,nest);     % Convert output matrix true to vector (evaluar si necesito esto o con el true de recent omitir) 
lat=linspace(5.83,6.82,12);
lon=linspace(-76.15,-75.16,12);

%% =================== Draw grid study domain==============================
% set(gcf,'color','w');
% set(gcf,'position',[x0,y0,width,height]);
% set(gcf,'defaultTextInterpreter','latex');
% [X,Y]=meshgrid(lon,lat);
% plot(X,Y,X',Y','linewidth', 2,'color', [0 0 0])
% set(gca,'xtick',lon)
% set(gca,'ytick',lat)
% xlabel('Longitude ','FontSize',20,'Interpreter','latex')
% ylabel('Latitude ','FontSize',20,'Interpreter','latex')
% title('Study domain','FontSize',20,'Interpreter','latex')

% figure
% 
% Nx=11;
% Ny=11;
% Nz=6;
% clf
% hold on
% [i,j]=meshgrid(1:Nx,1:Ny);
% k=zeros(Ny,Nx)+Nz;
% surf(i,j,k)
% [i,k]=meshgrid(1:Nx,1:Nz);
% j=zeros(Nz,Nx)+Ny;
% surf(i,j,k)
% [j,k]=meshgrid(1:Ny,1:Nz);
% i=zeros(Nz,Ny)+Nx;
% surf(i,j,k)
% [i,j]=meshgrid(1:Nx,1:Ny);
% k=zeros(Ny,Nx)+1;
% surf(i,j,k)
% [i,k]=meshgrid(1:Nx,1:Nz);
% j=zeros(Nz,Nx)+1;
% surf(i,j,k)
% [j,k]=meshgrid(1:Ny,1:Nz);
% i=zeros(Nz,Ny)+1;
% surf(i,j,k)
% view(30,30)
%% =================== Observation parameters =============================
% load HDay.mat         % Observation matrix spots of regular pattern   (11x11x7)
% load HBase.mat        % Observation matrix spots of random pattern     (11x11)
load H.mat              % Observation matrix operator H  (m X nest) (22,121)
load observedstate.mat    % Vector that save the number of the state observed 
m=22;                     % Number of observations
sigma=0.01;               % Observation error
M=6;           % Assimilation window time  (1 Day) form start_ass 24 to stop_ass 28
R=sigma^2*eye(m,m);       % Covariance error measurements
inR=inv(R);
frequency=1;              % Frequency of observations  TROPOMI frequency must be 24 (Once per day)
%% =================== DA Parameters ======================================
inner=5;                     % number of inner steps at the inner loop
tolerance_inner=2e-12;         % Criteria to get out the inner cycle
le=1;                         % perturbed sigma deviation for emission factor
ld=0;                         % perturbed sigma deviation for deposition factor
nens=40;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     ;
nens=nens+1;                  % put here number of ensembles (the extra is for the true)
meane=5;                      % mean of the ensemble
St_Ass=24;                              % (2,St_Ass = Start Assimilation)
muestreo=St_Ass:frequency:St_Ass+M;     % Frequency to have observations into the window of assimilation
%% ================= Ensemble generation and propagation===================                                                                                                                  
%====Simulation Initialization===  
xb=zeros(nest,tsim,nens-1); % Background State. Dimensions=(states,time,ensemble member)
xa=zeros(nest,tsim,nens-1); % Analysis State.   Dimensions=(states,time,ensemble member)
Xb=zeros(nest,tsim,nens-1); % Mean deviation matrix
Yx=zeros(m,tsim,nens-1);
meanyb=zeros(m,tsim);     % multiplication H*mean_column;

EnseMem_column=zeros(nest,tsim,nens);      % used to import the column output of NO2 for each ensemble member  
EnseMem_state=zeros(nest*5,tsim,nens);     % used to import the concentration of NO2 for each ensemble member
tic 
sprintf('Running forward ensemble')
[le,ld,c5]=propagate_LOTOSEUROS_ensemble_perturbation(nens,kkk,ld,pertu,meane,infl);      % propagate nens the model (Indicate inside the function is the perturbation is applied to initial conditions or parameters)
% %load 16members_IC.mat
toc
%% ====================== Ensemble arrays =============================== 
%========  This part is important  to stay comented if the function propagate_LOTO...

%load UnderWorkspace.mat    % load all the workspace generated from the initial condition perturbation with synthetic observation 12 
% load OverWorkspace.mat      % load all the workspace generated from the initial condition perturbation with synthetic observation 1 
%load OverWorkspace1.mat
% M=48;                       % this is because the load workspace save by error M=24 and just take one day
%##########################  part for pertu==1 ########################    
% Parameters
if pertu==1
    for i=1:nens
        cd /run/media/dirac/Datos/scratch/projects/LOTOS-EUROS/4DEnVar                  % folder of the principal output output
        %cd /run/media/dirac/Datos/scratch/projects/LOTOS-EUROS/4DEnVar_ensembles        % folder of outputs organize info in each ensemble
        s1='ens-';                  
        s2=num2str(i);
        s3='/output';
        s4=strcat(s1,s2,s3);
        cd(s4) 
        
        system('cdo select,name=no2_column LE_4DEnVar_column_2018*.nc no2_column_20180701_20180708.nc');           % Concatenate for columns
        system('cdo select,name=no2 LE_4DEnVar_conc-3d_2018*.nc no2_20180701_20180708.nc');                        % Concatenate for states
        no2_column=ncread('no2_column_20180701_20180708.nc','no2_column') ;     % columns   OUTPUTS   (nx,ny,tsim)
        no2=ncread('no2_20180701_20180708.nc','no2') ;                          % concentrations STATES   (nx,ny,nz,tsim)
        cd /run/media/dirac/Datos/Dropbox/users/arjo/lotos-euros/4DEnVar_LE/LE_2018_ensembles           % folder of the main file
        [no2_column_b] = mat2vect(no2_column,tsim,nest);    % mat2vect function to convert to new dimensions (nest,tsim)
        [no2_state_vector] = mat3D2vect(no2,tsim,nest);     % mat3Dvect function to convert to new dimensionsconvert to (nest,nz,tsim)
        EnseMem_column(:,:,i)=no2_column_b;                 % Organice in an array the column vectors                   
        EnseMem_state(:,:,i)=no2_state_vector;              % Organice in an array the state vectors
    end


    mean_column=mean(EnseMem_column,3);    % Ensemble mean column
    mean_state=mean(EnseMem_state,3);    % Ensemble mean state

    for d=St_Ass:M
        meanyb(:,d)=H*mean_column(:,d);         % Model Observation based on the ensemble mean 
    end

    for i=1:nens-1  %Pregunta Santiago, Por que aca es nens y en linea 208 es nens-1? ojo que eso puede danar las dimensiones
        Xb(:,:,i)=EnseMem_column(:,:,i)-mean_column;   % Deviation matrix from the columns
        Yx(:,:,i)=H*EnseMem_column(:,:,i)-meanyb;      %Aproximation of H*M*Xb according with equation (38)   
        XXb(:,:,i)=EnseMem_state(:,:,i)-mean_state;    % Deviation matrix from the states
    end
    
    
    
    
    
    Xb=Xb*sqrt(1/(nens-2));
    Yx=Yx*sqrt(1/(nens-2));
    XXb=XXb*sqrt(1/(nens-2));
    i=1;

    
%##########################  part for pertu==0 ########################    
% States     
elseif pertu==0
    for i=1:nens
        cd /run/media/dirac/Datos/scratch/projects/LOTOS-EUROS/4DEnVar            % folder of outputs organize info in each ensemble
        %cd /run/media/dirac/Datos/scratch/projects/LOTOS-EUROS/4DEnVar_ensembles        % folder of %outputs organize info in each ensemble
        
        s1='ens-';                  
        s2=num2str(i);
        s3='/output';
        s4=strcat(s1,s2,s3);
        cd(s4) 
        
        system('cdo select,name=no2_column LE_4DEnVar_column_2018*.nc no2_column_20180701_20180708.nc');           % Concatenate for columns
        system('cdo select,name=no2 LE_4DEnVar_conc-3d_2018*.nc no2_20180701_20180708.nc');                        % Concatenate for states
        no2_column=ncread('no2_column_20180701_20180708.nc','no2_column') ;     % columns   OUTPUTS   (nx,ny,tsim)
        no2=ncread('no2_20180701_20180708.nc','no2') ;                          % concentrations STATES   (nx,ny,nz,tsim)
%         no2_column_save(:,:,:,i)=no2_column;               %to save the value in matrix for eac ensamble       
%         no2_save(:,:,:,i)=no2;                             %to save the value in matrix for eac ensamble       
        cd /run/media/dirac/Datos/Dropbox/users/arjo/lotos-euros/4DEnVar_LE/LE_2018_ensembles           % folder of the main file
        [no2_column_b] = mat2vect(no2_column,tsim,nest);                                                % convert to (nest,tsim)
        [no2_state_vector] = mat3D2vect(no2,tsim,nest);                                                 % convert to (nest,nz,tsim)
        EnseMem_column(:,:,i)=no2_column_b;                 % Organice in an array the column vectors                   
        EnseMem_state(:,:,i)=no2_state_vector;              % Organice in an array the state vectors
   
        EnseMem_column1(:,:,i)=no2_column_b;                 % Organice in an array the column vectors                   
        EnseMem_state1(:,:,i)=no2_state_vector;  
    end

    
    mean_column=mean(EnseMem_column(:,:,1:end-1),3);    % Ensemble mean column
    mean_state=mean(EnseMem_state(:,:,1:end-1),3);    % Ensemble mean state

    for d=St_Ass:M                              % This for is in the DA window
        meanyb(:,d)=H*mean_column(:,d);         % Model Observation based on the ensemble mean 
    end

    for i=1:nens-1
        Xb(:,:,i)=EnseMem_column(:,:,i)-mean_column;   % Deviation matrix from the columns
        Yx(:,:,i)=H*EnseMem_column(:,:,i)-meanyb;      %Aproximation of H*M*Xb according with equation (38)   
        XXb(:,:,i)=EnseMem_state(:,:,i)-mean_state;    % Deviation matrix from the states
    end

    True=EnseMem_column(:,:,end);       % True generated is the last member of the ensemble
    
    Xb=Xb*sqrt(1/(nens-2));
    Yx=Yx*sqrt(1/(nens-2));
    XXb=XXb*sqrt(1/(nens-2));2,
    i=1;
    
    
%%  Small section to ploth the deviation matrix plot 
% this small section works to plot the deviation member rest with the mean
% figure
% 
% hold on
% set(gcf,'color','w');
% set(gcf,'position',[x0,y0,width,height]);
% set(gcf,'defaultTextInterpreter','latex');
% 
% for i=1:nens-1
% 
% subplot(1,3,1)
% imagesc(EnseMem_column(:,:,i));
% subplot(1,3,2)
% imagesc(mean_column);
% subplot(1,3,3)
% c2=imagesc(abs(EnseMem_column(:,:,i)-mean_column)); 
% c2 = colormap(gca,'jet');
% suptitle(sprintf('Ensemble %i',i))    
% pause(0.5)
% 
% end
end

save mean_column5 mean_column
save mean_state5 mean_state
save EnseMem_column5 EnseMem_column
save EnseMem_state5 EnseMem_state
% 
% 
% 
% %% ========================== Analysis step ===============================
% % %Inner Loop
%      disp(['Analysis Step. Porcentage Analysis: ',num2str(1/inner*100),'%'])
%      dx=zeros(nens-1,inner);   % incremental direction in the ensemble space
%      ysinthetic=zeros(nest,tsim);
%      
%      xmin=-1;
%      xmax=1;
%      n=nens-1;
%      x=xmin+rand(1,n)*(xmax-xmin);
%      dx(:,1)=x';
% %    dx(:,1)=rand(nens,1);
%      xk=zeros(nest*5,tsim-St_Ass,inner);                                                %time of assimilation + time free propagation, the 5 is because the level take in acount states  
%      xk(:,St_Ass,1)=max(mean_state(:,St_Ass)+squeeze(XXb(:,St_Ass,:))*dx(:,1),0);
%      new_no2=reshape(xk(:,St_Ass,1),11,11,5,1);   
%      
%      mean_state_matrix=max(reshape(mean_state(:,St_Ass),11,11,5),0);
%      first_deviation= reshape(squeeze(XXb(:,St_Ass,:))*dx(:,1),11,11,5),0;
% %% ---------   plot of the initial step to perform the assimilation step
%      
%      %      imagesc(Mean_state_matrix(:,:,1))
% % suptitle('Mean surface concentration start assimilation window (24 hours)')     
% 
% %-----------------------------------------------------------------------------------------------------------------------
% % hold on
% % set(gcf,'color','w');
% % set(gcf,'position',[x0,y0,width,height]);
% % set(gcf,'defaultTextInterpreter','latex');
% % for i=1:5
% % subplot1=subplot(1,5,i)
% % 
% % imagesc(mean_state_matrix(:,:,i))
% % h=colorbar;
% % caxis([2e-10 18e-10])
% % h.Location='southoutside' 
% % set(subplot1,'CLim',[2e-10 1.8e-09],'Layer','top','XTickLabel',...
% %     {'-76.06','-75.88','-75.7','-75.52','-75.34'},'YTickLabel',...
% %     {'6.82','6.73','6.64','6.55','6.46','6.37','6.28','6.19','6.1','6.01','5.92'});
% % ylabel(h, '[mol no_2 mol-1 air]')
% % xlabel('Longitude ','FontSize',20,'Interpreter','latex')
% % ylabel('Latitude ','FontSize',20,'Interpreter','latex')
% % end
% %-------------------------------------------------------------------------------------------------------------------  
% figure
% 
% sb1=subplot(1,3,1);
% imagesc(new_no2(:,:,1))
% set(sb1,'CLim',[2e-10 1.8e-09],'Layer','top','XTickLabel',...
%     {'-76.06','-75.88','-75.7','-75.52','-75.34'},'YTickLabel',...
%     {'6.82','6.73','6.64','6.55','6.46','6.37','6.28','6.19','6.1','6.01','5.92'});
% h=colorbar;
% title('new NO_2')
% ylabel(h, '[mol no_2 mol-1 air]')
% xlabel('Longitude ','FontSize',20,'Interpreter','latex')
% ylabel('Latitude ','FontSize',20,'Interpreter','latex')
% 
% sb2=subplot(1,3,2);
% imagesc(mean_state_matrix(:,:,1))
% set(sb2,'CLim',[2e-10 1.8e-09],'Layer','top','XTickLabel',...
%     {'-76.06','-75.88','-75.7','-75.52','-75.34'},'YTickLabel',...
%     {'6.82','6.73','6.64','6.55','6.46','6.37','6.28','6.19','6.1','6.01','5.92'});
% h=colorbar;
% ylabel(h, '[mol no_2 mol-1 air]')
% xlabel('Longitude ','FontSize',20,'Interpreter','latex')
% ylabel('Latitude ','FontSize',20,'Interpreter','latex')
% 
% sb3=subplot(1,3,3);
% imagesc(first_deviation(:,:,1))
% set(sb3,'Layer','top','XTickLabel',...
%     {'-76.06','-75.88','-75.7','-75.52','-75.34'},'YTickLabel',...
%     {'6.82','6.73','6.64','6.55','6.46','6.37','6.28','6.19','6.1','6.01','5.92'});
% h=colorbar;
% ylabel(h, '[mol no_2 mol-1 air]')
% xlabel('Longitude ','FontSize',20,'Interpreter','latex')
% ylabel('Latitude ','FontSize',20,'Interpreter','latex')     
% %% =================== Analysis step:  transform the mean no2 to the matrix form (nx,ny,nz,tsim) ========================
% !mkdir /run/media/dirac/Datos/scratch/projects/LOTOS-EUROS/4DEnVar_ensembles/ens-1 /run/media/dirac/Datos/scratch/projects/LOTOS-EUROS/4DEnVarEns1
% !cp -a /run/media/dirac/Datos/scratch/projects/LOTOS-EUROS/4DEnVar_ensembles/ens-1 /run/media/dirac/Datos/scratch/projects/LOTOS-EUROS/4DEnVarEns1/ens-1
% 
% cd /run/media/dirac/Datos/scratch/projects/LOTOS-EUROS/4DEnVarEns1/ens-1/restart;      % go to folder where the restart file is going to consider
%       
% % % 
%  c=ncread('LE_4DEnVar_state_20180702_0000.nc','c');   % read the actual concentration variables from the restart correspond field
% %-----------------
%  figure
%  hold on
% set(gcf,'color','w');
% set(gcf,'position',[x0,y0,width,height]);
% set(gcf,'defaultTextInterpreter','latex');
% 
%  hold on
% set(gcf,'color','w');
% set(gcf,'position',[x0,y0,width,height]);
% set(gcf,'defaultTextInterpreter','latex');
%  
%  sb666=subplot(1,2,1)         % this plot is just to compare the mean actual concentration with the one in the restart file
%  imagesc(c(:,:,1,1));
%  set(sb666,'CLim',[0 20],'Layer','top','XTickLabel',...
%     {'-76.06','-75.88','-75.7','-75.52','-75.34'},'YTickLabel',...
%     {'6.82','6.73','6.64','6.55','6.46','6.37','6.28','6.19','6.1','6.01','5.92'});
% h=colorbar;
% colormap('jet');
% title('Default emission NO$_2$','FontSize',24,'Interpreter','latex')
% xlabel('Longitude ','FontSize',20,'Interpreter','latex')
% ylabel('Latitude ','FontSize',20,'Interpreter','latex')
%  colorbar
%  d(:,:,:,1)=new_no2(:,:,:)*1e9/2;
% sb667= subplot(1,2,2)
%  imagesc(d(:,:,1,1));
%   set(sb667,'CLim',[0 20],'Layer','top','XTickLabel',...
%     {'-76.06','-75.88','-75.7','-75.52','-75.34'},'YTickLabel',...
%     {'6.82','6.73','6.64','6.55','6.46','6.37','6.28','6.19','6.1','6.01','5.92'});
% h=colorbar;
% colormap('jet');
% title('Updated emission NO$_2$','FontSize',24,'Interpreter','latex')
% xlabel('Longitude ','FontSize',20,'Interpreter','latex')
% ylabel('Latitude ','FontSize',20,'Interpreter','latex')
% colorbar
% %------------------
% 
% ncwrite('LE_4DEnVar_state_20180702_0000.nc','c',c);  % Rewrite the matrix in the states to run the model from that point
% 
% % cd /run/media/dirac/Datos/Dropbox/users/arjo/lotos-euros/4DEnVar_LE/LE_2018_ensembles   % Main folder
%    
% 
% % cd /run/media/dirac/Datos/scratch/projects/LOTOS-EUROS/4DEnVar_ANALYSIS/ens-1/run
% 
% % system('rm -rf le.ok')
% 
% cd /run/media/dirac/Datos/Dropbox/users/arjo/lotos-euros/4DEnVar_LE/LE_2018_ensembles   % Main folder
%    
%    
%    fprintf( '   \n                                                     ' );
%    fprintf( '    \n                                                    ' );
%    fprintf( '     \n                                                   ' );
%    [le,ld,c5]=propagate_LOTOSEUROS_analysis(nens,le,ld);
%    
%   
%    muestra=0;   % inicializa este contador para almacenar cada una de las observaciones
%    
% % part to read outputs of the model, concatenate 
%    i=1;
%    cd /run/media/dirac/Datos/scratch/projects/LOTOS-EUROS/4DEnVar_ANALYSIS      % folder of outputs organize info in each ensemble
% 
% s1='ens-';                  
% s2=num2str(i);
% s3='/output';
% s4=strcat(s1,s2,s3);
% cd(s4) 
% system('cdo select,name=no2_column LE_4DEnVar_ANALYSIS_column_2018*.nc no2_column_analy_20180701_20180708.nc');           % Concatenate for columns
% system('cdo select,name=no2 LE_4DEnVar_ANALYSIS_conc-3d_2018*.nc no2_analy_20180701_20180708.nc');                        % Concatenate for states
% 
% no2_column_anal=ncread('no2_column_analy_20180701_20180708.nc','no2_column') ;     % columns   OUTPUTS   (nx,ny,tsim)
% no2_anal=ncread('no2_analy_20180701_20180708.nc','no2') ;                          % concentrations STATES   (nx,ny,nz,tsim)
% 
% cd /run/media/dirac/Datos/Dropbox/users/arjo/lotos-euros/4DEnVar_LE/LE_2018_ensembles           % folder of the main file
% 
% [no2_column_b_analysis] = mat2vect_analysis(no2_column_anal,tsim,nest);                                                % convert to (nest,tsim)
% [no2_state_vector_analysis] = mat3D2vect_analysis(no2_anal,tsim,nest);                                                 % convert to (nest,nz,tsim)
%  
% no2_column_b;                 % Organice in an array the column vectors                   
% no2_state_vector;              % Organice in an array the state vectors
% 
%    
%    
%    
%    
%  for i=St_Ass+1:tsim   
%     if sum(muestreo==i) 
%          muestra=muestra+1;
%           y(:,muestra)=H*EnseMem_column(:,i,end)+sigma*rand(1,1);                       
%           yxk(:,muestra)=H*no2_column_b_analysis(:,i);                % xk corrida anterior media
%           dd(:,muestra)=(y(:,muestra)-yxk(:,muestra));  %Innovation Matriz
%     end
%  end 
%      
%           
%      k=2;
%     norma(1)=tolerance_inner;
%     flaq=0;
%     
%        while norma(k-1)>=tolerance_inner && k<inner && flaq<2
%            sum1=0;
%            sum2=0;
% 
%            for j=1:M
%               sum1=sum1 + squeeze(Yx(:,muestreo(j),:))'*inR*squeeze(Yx(:,muestreo(j),:));     %Sumatory functional cost
%               sum2=sum2+ squeeze(Yx(:,muestreo(j),:))'*inR*dd(:,j);
%            end
%               dx(:,k)= real(pinv(eye(n,1)+sum1))*sum2;               
%             %=Incremental=
% %           xk(:,St_Ass,k)=xk(:,St_Ass,k-1)+squeeze(Xb(:,St_Ass,k-1))*dx(:,k-1);
% %       Santiago hizo esta modificacion
% 
%         xk(:,St_Ass,k)=max(xk(:,St_Ass,k-1)+squeeze(XXb(:,St_Ass,:))*dx(:,k-1),0);
%        disp(['Analysis Step. Porcentage Analysis: ',num2str(k/inner*100),'%'])
%        
%        
%      new_no2=reshape(xk(:,St_Ass,k),11,11,5,1);    %  transform the mean no2 to the matrix form (nx,ny,nz,tsim)
% % % 
% 
% 
% 
% 
% 
% 
% cd /run/media/dirac/Datos/scratch/projects/LOTOS-EUROS/4DEnVarEns1/ens-1/restart;    % go to folder where the restart file is going to consider
% %  pwd
% % % 
%   c=ncread('LE_4DEnVar_state_20180702_0000.nc','c');   % read the actual concentration variables from the restart correspond field
% % %-----------------
% %  figure
% % subplot(1,2,1)                                 % this plot is just to compare the mean actual concentration with the one in the restart file
% % imagesc(c(:,:,1,1));
% % colorbar
% c(:,:,:,1)=new_no2(:,:,:)*1e9;
% % subplot(1,2,2)
% % imagesc(c(:,:,1,1));
% % colorbar
% % %------------------
% % 
%  ncwrite('LE_4DEnVar_state_20180702_0000.nc','c',c);  % Rewrite the matrix in the states to run the model from that point
% % 
%  cd /run/media/dirac/Datos/Dropbox/users/arjo/lotos-euros/4DEnVar_LE/LE_2018_ensembles   % Main folder
% %    
% 
% %    
%    fprintf( '   \n                                                     ' );
%    fprintf( '    \n                                                    ' );
%    fprintf( '     \n                                                   ' );
%    
%    
%    
%   system('rm -rf  /run/media/dirac/Datos/scratch/projects/LOTOS-EUROS/4DEnVar_ANALYSIS')
%    
%    [le,ld,c5]=propagate_LOTOSEUROS_analysis(nens,le,ld);
% %    
% %   
%    muestra=0;   % inicializa este contador para almacenar cada una de las observaciones
% %    
% % % part to read outputs of the model, concatenate 
%    i=1;
%    cd /run/media/dirac/Datos/scratch/projects/LOTOS-EUROS/4DEnVar_ANALYSIS      % folder of outputs organize info in each ensemble
% % 
%  s1='ens-';                  
%  s2=num2str(i);
%  s3='/output';
%  s4=strcat(s1,s2,s3);
%  cd(s4) 
% 
%  system('cdo select,name=no2_column LE_4DEnVar_ANALYSIS_column_2018*.nc no2_column_20180701_20180708.nc');           % Concatenate for columns
%  system('cdo select,name=no2 LE_4DEnVar_ANALYSIS_conc-3d_2018*.nc no2_20180701_20180708.nc');                        % Concatenate for states
% % 
%  no2_column_anal=ncread('no2_column_20180701_20180708.nc','no2_column') ;     % columns   OUTPUTS   (nx,ny,tsim)
%  no2_anal=ncread('no2_20180701_20180708.nc','no2') ;                          % concentrations STATES   (nx,ny,nz,tsim)
% % 
%  cd /run/media/dirac/Datos/Dropbox/users/arjo/lotos-euros/4DEnVar_LE/LE_2018_ensembles           % folder of the main file
% % 
%  [no2_column_b_analysis] = mat2vect_analysis(no2_column_anal,tsim,nest);              % convert to (nest,tsim)
%  [no2_state_vector_analysis] = mat3D2vect_analysis(no2_anal,tsim,nest);     
% % 
%           for i=St_Ass+1:tsim              
%                if sum(muestreo==i) 
%                     muestra=muestra+1;
%                     y(:,muestra)=H*EnseMem_column(:,i,end)+sigma*rand(1,1);
%                     yxk(:,muestra)=H*no2_column_b_analysis(:,i);
%                     dd(:,muestra)=(y(:,muestra)-yxk(:,muestra));  %Innovation Matriz
%                end
%           end
%           norma(k)=norm(squeeze(XXb(:,St_Ass,:))*dx(:,k))
%           if diff(norma(k-1:k))>0 
%               flaq=flaq+1;          
%           end
%           if flaq==2
%              dx(:,k)=dx(:,k-1); 
%              xk(:,k)=xk(:,k-1)xlim([-0.8 0.8]);
%           end
%           k=k+1;
%        end
%        
%       no2_column_a_analysis=no2_column_b_analysis;
%       no2_state_vector_a_analysis=no2_state_vector_analysis;
% %      
% %      
% %==End Inner Loop==
  
%% ======================================  Plot step Column ======================================================
load mean_column.mat
load mean_state.mat
load EnseMem_column.mat
load EnseMem_state.mat



hold on
set(gcf,'color','w');
set(gcf,'position',[x0,y0,width,height]);
set(gcf,'defaultTextInterpreter','latex');



subplot(1,3,1)

nx=5;ny=3;                           % coordinate of the first state to plot  0<nx<=11,0<ny<=11


S.ensemble(:,:,:,kkk)=EnseMem_column;


state1=(ny-1)*gridy+nx;

diffu=10;
%True=plot(squeeze(EnseMem_column(state1,5:end,end)),'-r','LineWidth',1);
%hold on
for i=1:nens-1
 EEE=plot(EnseMem_column(state1,5:end,i),'linewidth',1,'Color',[0 0 0]+0.05*diffu);xlim([0,tsim])
hold on
end
Ensemble_mean=plot(squeeze(mean_column(state1,5:end)),'k','LineWidth',1);
hold on;
%True=plot(squeeze(EnseMem_column(state1,5:end,nens)),'-r','LineWidth',1);
y1=get(gca,'ylim');
iniAss=plot([St_Ass St_Ass],y1,'g--','LineWidth',2);
FinAss=plot([St_Ass+24 St_Ass+24],y1,'g--','LineWidth',2);

%h=legend([Ensemble_mean True iniAss EEE ],'Ensemble mean meanxb','Xtrue','Init/end assimilation','Ensemble member');hold on;
h=legend([Ensemble_mean iniAss EEE ],'Ensemble mean meanxb','Init/end assimilation','Ensemble member');hold on;
set(h,'FontSize',11,'interpreter','latex');grid on;
xlabel('Time [days]','FontSize',Fx,'interpreter','latex');ylabel('Concentration [mlc/cm$^2$]*1e15','FontSize',Fx,'interpreter','latex')
n=tsim;x = linspace(0,n,n)';b=datenum(2016,07,01,0,0,0);m = 0:1/24:n/24;m=datevec(b+m);
xd=datestr(m);xe=x(1:24:end);xd=xd(1:24:end,1:end-14);TimeAxes=gca;set(TimeAxes,'Xtick',xe);set(TimeAxes,'XTickLabel',xd);
mytitle1=['NO$_2$ Column State 1'];
title(mytitle1,'FontSize',Ft,'interpreter','latex')
dim = [.137 .61 .7 .3];
str = ['lat=  ',num2str(lat(ny)),'°';'lon=',num2str(lon(nx)),'°'];
annotation('textbox',dim,'String',str,'FitBoxToText','on','interpreter','latex');
grid on



subplot(1,3,2)

n2x=4;n2y=9;                        % coordinate of the second state to plot  0<nx<=11,0<ny<=11

state2=(n2y-1)*gridy+n2x;
%True=plot(squeeze(EnseMem_column(state2,5:end,end)),'-r','LineWidth',1);

for i=1:nens-1
 EEE=plot(EnseMem_column(state2,5:end,i),'linewidth',1,'Color',[0 0 0]+0.05*diffu);xlim([0,tsim])
hold on
end
Ensemble_mean=plot(squeeze(mean_column(state2,5:end)),'k','LineWidth',1);
hold on;
%True=plot(squeeze(EnseMem_column(state2,5:end,nens)),'-r','LineWidth',1);
y1=get(gca,'ylim');
iniAss=plot([St_Ass St_Ass],y1,'g--','LineWidth',2);
FinAss=plot([St_Ass+24 St_Ass+24],y1,'g--','LineWidth',2);
%h=legend([Ensemble_mean True iniAss EEE ],'Ensemble mean meanxb','Xtrue','Init/end assimilation','Ensemble member');hold on;
h=legend([Ensemble_mean iniAss EEE ],'Ensemble mean meanxb','Init/end assimilation','Ensemble member');hold on;
set(h,'FontSize',11,'interpreter','latex');grid on;
xlabel('Time [days]','FontSize',Fx,'interpreter','latex');ylabel('Concentration [mlc/cm$^2$]*1e15','FontSize',Fx,'interpreter','latex')
n=tsim;x = linspace(0,n,n)';b=datenum(2016,07,01,0,0,0);m = 0:1/24:n/24;m=datevec(b+m);
xd=datestr(m);xe=x(1:24:end);xd=xd(1:24:end,1:end-14);TimeAxes=gca;set(TimeAxes,'Xtick',xe);set(TimeAxes,'XTickLabel',xd);
mytitle1=['NO$_2$ Column State 2'];
title(mytitle1,'FontSize',Ft,'interpreter','latex')
dim = [.418 .61 .7 .3];
str = ['lat=  ',num2str(lat(n2y)),'°';'lon=',num2str(lon(n2x)),'°'];
annotation('textbox',dim,'String',str,'FitBoxToText','on','interpreter','latex');
grid on


subplot(1,3,3)

n3x=8;n3y=6;                        % coordinate of the second state to plot  0<nx<=11,0<ny<=11

state3=(n3y-1)*gridy+n3x;
%True=plot(squeeze(EnseMem_column(state2,5:end,end)),'-r','LineWidth',1);

for i=1:nens-1
 EEE=plot(EnseMem_column(state3,5:end,i),'linewidth',1,'Color',[0 0 0]+0.05*diffu);xlim([0,tsim])
hold on
end
Ensemble_mean=plot(squeeze(mean_column(state3,5:end)),'k','LineWidth',1);
hold on;
%True=plot(squeeze(EnseMem_column(state2,5:end,nens)),'-r','LineWidth',1);
y1=get(gca,'ylim');
iniAss=plot([St_Ass St_Ass],y1,'g--','LineWidth',2);
FinAss=plot([St_Ass+24 St_Ass+24],y1,'g--','LineWidth',2);
%h=legend([Ensemble_mean True iniAss EEE ],'Ensemble mean meanxb','Xtrue','Init/end assimilation','Ensemble member');hold on;
h=legend([Ensemble_mean iniAss EEE ],'Ensemble mean meanxb','Init/end assimilation','Ensemble member');hold on;
set(h,'FontSize',11,'interpreter','latex');grid on;
xlabel('Time [days]','FontSize',Fx,'interpreter','latex');ylabel('Concentration [mlc/cm$^2$]*1e15','FontSize',Fx,'interpreter','latex')
n=tsim;x = linspace(0,n,n)';b=datenum(2016,07,01,0,0,0);m = 0:1/24:n/24;m=datevec(b+m);
xd=datestr(m);xe=x(1:24:end);xd=xd(1:24:end,1:end-14);TimeAxes=gca;set(TimeAxes,'Xtick',xe);set(TimeAxes,'XTickLabel',xd);
mytitle1=['NO$_2$ Column State 3'];
title(mytitle1,'FontSize',Ft,'interpreter','latex')
dim = [.700 .61 .7 .3];
str = ['lat=  ',num2str(lat(n3y)),'°';'lon=',num2str(lon(n3x)),'°'];
annotation('textbox',dim,'String',str,'FitBoxToText','on','interpreter','latex');
grid on