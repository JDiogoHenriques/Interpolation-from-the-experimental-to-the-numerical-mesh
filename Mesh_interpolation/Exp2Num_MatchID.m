%% Input

% J. Henriques 2023

% clear workSpace
clc; close all; clear all;
tic
% ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: %
%User options
specimen = 'Dogbone_0'; % Dogbone_0 ; Dogbone_45 ; Dogbone_90
Interpolation = 'natural'; % linear ; nearest ; natural
Extrapolation = 'nearest'; % linear ; nearest
% ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: %

% ensure root units are pixels and get the size of the screen
set(0,'Units','pixels');
script.scnsize = get(0,'ScreenSize');
% define the size and location of the figures
script.fig_pos = [script.scnsize(3)/5, script.scnsize(4)/5, ...
    script.scnsize(3)*1/2, script.scnsize(4)*1/1.4];

% image info
script.size_font1 = 19; 
script.size_font2 = 22; 
script.size_font3 = 24;
script.nomeF      = 'Times New Roman';
script.imgformat  = '-djpeg';
script.filetype   = '.jpg';
script.filetype2  = '.eps';
script.mat        = '.mat';
script.dat        = '.dat';
script.ext        = '.txt';
script.csv        = '.csv';
script.resol      = '-r300';
script.xls        = '.xls';
script.MatchIDext = '.tif';

% Change default axes fonts.
set(0,'DefaultAxesFontName', script.nomeF)
set(0,'DefaultAxesFontSize', script.size_font2)

% Change default text fonts.
set(0,'DefaultTextFontname', script.nomeF)
set(0,'DefaultTextFontSize', script.size_font2)

% struture variable
Exp = struct;
Exp.specimen = specimen;

% Get the current working directory
current_dir = pwd;

% Search for the experimental and numerical data
switch specimen
        case 'Dogbone_0'
        source_dir = [fileparts(pwd),filesep,'Experimental_Data',filesep,'MatchID',filesep,'0deg']; 
        d = dir([source_dir, '\*.csv']);
        %ID = 'dogBone_0';
        %Centroid data
        Mesh_X = importdata([current_dir,filesep,specimen,filesep,'COORD_1_int_points_all_el_undef_DB-0.csv']); %FEA Mesh X (centroid)
        Mesh_Y = importdata([current_dir,filesep,specimen,filesep,'COORD_2_int_points_all_el_undef_DB-0.csv']); %FEA Mesh Y (centroid)
        %Nodal data
        Mesh_nodal = importdata([current_dir,filesep,specimen,filesep,'DB-0.inp']); %FEA Mesh X (centroid)
        %Defining alignment settings
        Rot_Angle = -0.3; %degree anti clock-wise
        x_trans = -0.65; % translation in the x axis
        y_trans = 0.5; % translation in the y axis
        x_scaling = 1.029; %scaling to match both meshes in x-axis 1.0199
        
        case 'Dogbone_45'
        source_dir = [fileparts(pwd),filesep,'Experimental_Data',filesep,'MatchID',filesep,'45deg']; 
        d = dir([source_dir, '\*.csv']);      
        %ID = 'dogBone_45';
        %Centroid data
        Mesh_X = importdata([current_dir,filesep,specimen,filesep,'COORD_1_int_points_all_el_undef_DB-45.csv']); %FEA Mesh X (centroid)
        Mesh_Y = importdata([current_dir,filesep,specimen,filesep,'COORD_2_int_points_all_el_undef_DB-45.csv']); %FEA Mesh Y (centroid)
        %Nodal data
        Mesh_nodal = importdata([current_dir,filesep,specimen,filesep,'DB-45.inp']); %FEA Mesh X (centroid)
        %Defining alignment settings
        Rot_Angle = -0.2; %degree anti clock-wise
        x_trans = -0.75; % translation in the x axis
        y_trans = -0.68; % translation in the y axis
        x_scaling = 1.0305; %scaling to match both meshes in x-axis
        
        case 'Dogbone_90'
        source_dir = [fileparts(pwd),filesep,'Experimental_Data',filesep,'MatchID',filesep,'90deg']; 
        d = dir([source_dir, '\*.csv']);    
        %ID = 'dogBone_90';
        %Centroid data
        Mesh_X = importdata([current_dir,filesep,specimen,filesep,'COORD_1_int_points_all_el_undef_DB-90.csv']); %FEA Mesh X (centroid)
        Mesh_Y = importdata([current_dir,filesep,specimen,filesep,'COORD_2_int_points_all_el_undef_DB-90.csv']); %FEA Mesh Y (centroid)
        %Nodal data
        Mesh_nodal = importdata([current_dir,filesep,specimen,filesep,'DB-90.inp']); %FEA Mesh X (centroid)
        %Defining alignment settings
        Rot_Angle = -0.3; %degree anti clock-wise
        x_trans = -0.86; % translation in the x axis
        y_trans = -0.5; % translation in the y axis
        x_scaling = 1.03; %scaling to match both meshes in x-axis
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% FEA DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Processing for the numeric data
%Centroid data
Mesh_X = Mesh_X.data(3,1:end).';
Mesh_Y = Mesh_Y.data(3,1:end).';
x_FEA = Mesh_X;
y_FEA = Mesh_Y;

%Nodal data
x_FEA_nodal = Mesh_nodal.data(:,2);
y_FEA_nodal = Mesh_nodal.data(:,3);

%% Import DIC data

% Total number of stages
n_stages = length(d);                
Exp.stages = n_stages;

for i =1:n_stages
Exp.index(i,1) = i;
end

eval(['Exp.x_coord',' = [];'])
eval(['Exp.y_coord',' = [];'])
eval(['Exp.ux',' = [];'])
eval(['Exp.uy',' = [];'])
eval(['Exp.exx',' = [];'])
eval(['Exp.eyy',' = [];'])
eval(['Exp.exy',' = [];'])

%Loading static DIC coord
%Opening files 
fid_s = fopen([source_dir,filesep,'Static',filesep,'imagecam0000_0.tif.csv']);
Static_data = readmatrix([source_dir,filesep,'Static',filesep,'imagecam0000_0.tif.csv'], 'HeaderLines', 1);
fclose(fid_s);

Exp.x_coord.stage0 = Static_data(:,3);
Exp.y_coord.stage0 = Static_data(:,4);

% All Stages

for j = 1:length(Exp.index)
    i = Exp.index(j);
    disp([specimen,' - reading cvs files: ','stage ',num2str(i)])
    filen = [source_dir,filesep,'imagecam',...
         num2str(mod(i,n_stages+1).', '%04d'),'_0.tif.csv'];
    %Opening files 
    fid = fopen(filen);
    matrix_data = readmatrix(filen, 'HeaderLines', 1);
    fclose(fid);
    %Assigning matrixes to struct
    eval(['Exp.x_coord.stage',num2str(j),' = matrix_data(:,3);']) % mm
    eval(['Exp.y_coord.stage',num2str(j),' = matrix_data(:,4);']) % mm
    eval(['Exp.ux.stage',num2str(j),' = matrix_data(:,5);']) % mm
    eval(['Exp.uy.stage',num2str(j),' = matrix_data(:,6).*-1;']) % mm
    eval(['Exp.exx.stage',num2str(j),' = matrix_data(:,7);']) % mm
    eval(['Exp.eyy.stage',num2str(j),' = matrix_data(:,8);']) % mm
    eval(['Exp.exy.stage',num2str(j),' = matrix_data(:,9).*-2;']) % mm
end

%% Processing DIC data %% 
% Shifting the origin to the corner of the specimen
for i = 0:n_stages    
    xdic_aux   = eval(['Exp.x_coord.stage',num2str(i)]);
    xdic_aux = xdic_aux-min(xdic_aux(:));
    eval(['Exp.x_coord.stage',num2str(i),' = xdic_aux;'])
    
    ydic_aux   = eval(['Exp.y_coord.stage',num2str(i)]);
    ydic_aux = -1.*(ydic_aux-max(ydic_aux(:)));
    eval(['Exp.y_coord.stage',num2str(i),' = ydic_aux;'])   
end

% Rotating data to align in the vertical direction
for i = 0:n_stages    
    xdic_aux   = eval(['Exp.x_coord.stage',num2str(i)]);
    ydic_aux   = eval(['Exp.y_coord.stage',num2str(i)]);
    
    xdic_rotaux = xdic_aux.*cosd(Rot_Angle)-ydic_aux.*sind(Rot_Angle);
    eval(['Exp.x_coord.stage',num2str(i),' = xdic_rotaux;'])
    
    ydic_rotaux = ydic_aux.*cosd(Rot_Angle)+xdic_aux.*sind(Rot_Angle);
    eval(['Exp.y_coord.stage',num2str(i),' = ydic_rotaux;'])   
end

% Calculating differences in the coordinate dimensions between exp and num
diff_x = max(Exp.x_coord.stage0) - max(x_FEA_nodal);
diff_y = max(Exp.y_coord.stage0) - max(y_FEA_nodal);

% Centering DIC data
 for i = 0:n_stages    
  	xdic_aux   = eval(['Exp.x_coord.stage',num2str(i)]);
  	xdic_aux = (xdic_aux-diff_x/2+x_trans)*x_scaling;
  	eval(['Exp.x_coord.stage',num2str(i),' = xdic_aux;'])
       
 	ydic_aux   = eval(['Exp.y_coord.stage',num2str(i)]);
 	ydic_aux = ydic_aux-diff_y/2+y_trans;
 	eval(['Exp.y_coord.stage',num2str(i),' = ydic_aux;'])   
 end

%Static time
DIC.x.stage0 = Exp.x_coord.stage0;
DIC.y.stage0 = Exp.y_coord.stage0;

% Removing NaNs and turning matrixes to arrays
for i = 1:n_stages 
    eval([' mask = any(isnan(Exp.exx.stage',num2str(i),'), 2);']) 
    
    aux_X = rmmissing(eval(['Exp.x_coord.stage',num2str(i)]));
    eval(['DIC.x.stage',num2str(i),' = aux_X(~mask, :);'])
    
    aux_Y = rmmissing(eval(['Exp.y_coord.stage',num2str(i)]));
    eval(['DIC.y.stage',num2str(i),' = aux_Y(~mask, :);'])
    
    aux_U = rmmissing(eval(['Exp.ux.stage',num2str(i)]));
    eval(['DIC.ux.stage',num2str(i),' = aux_U(~mask, :);'])
    
    aux_V = rmmissing(eval(['Exp.uy.stage',num2str(i)]));
    eval(['DIC.uy.stage',num2str(i),' = aux_V(~mask, :);']) 
    
    aux_XX = rmmissing(eval(['Exp.exx.stage',num2str(i)]));
    eval(['DIC.epsXX.stage',num2str(i),' = aux_XX;'])

    aux_YY = rmmissing(eval(['Exp.eyy.stage',num2str(i)]));
    eval(['DIC.epsYY.stage',num2str(i),' = aux_YY;'])
    
    aux_XY = rmmissing(eval(['Exp.exy.stage',num2str(i)]));
    eval(['DIC.epsXY.stage',num2str(i),' = aux_XY;']) 
end


%% Interpolation %%
disp('Doing the interpolation...')
for i = 1:n_stages  
    
    %ux
    eval(['ux_DIC2FEM_aux = scatteredInterpolant(DIC.x.stage0,DIC.y.stage0,DIC.ux.stage',num2str(i),');']) 
    ux_DIC2FEM_aux.Method = Interpolation;
    ux_DIC2FEM_aux.ExtrapolationMethod = Extrapolation;
    
    ux_FEM2DIC = ux_DIC2FEM_aux(x_FEA_nodal,y_FEA_nodal);
    
    eval(['DIC2FEM.ux.stage',num2str(i),'= ux_FEM2DIC;'])    
    
    %uy
    eval(['uy_DIC2FEM_aux = scatteredInterpolant(DIC.x.stage0,DIC.y.stage0,DIC.uy.stage',num2str(i),');']) 
    uy_DIC2FEM_aux.Method = Interpolation;
    uy_DIC2FEM_aux.ExtrapolationMethod = Extrapolation;
    
    uy_FEM2DIC = uy_DIC2FEM_aux(x_FEA_nodal,y_FEA_nodal);
    
    eval(['DIC2FEM.uy.stage',num2str(i),'= uy_FEM2DIC;'])    
    
    %epsXX
    eval(['epsXX_DIC2FEM_aux = scatteredInterpolant(DIC.x.stage0,DIC.y.stage0,DIC.epsXX.stage',num2str(i),');']) 
    epsXX_DIC2FEM_aux.Method = Interpolation;
    epsXX_DIC2FEM_aux.ExtrapolationMethod = Extrapolation;
    
    epsXX_FEM2DIC = epsXX_DIC2FEM_aux(x_FEA,y_FEA);
    
    eval(['DIC2FEM.exx.stage',num2str(i),'= epsXX_FEM2DIC;']) 
    
    %epsYY 
    eval(['epsYY_DIC2FEM_aux = scatteredInterpolant(DIC.x.stage0,DIC.y.stage0,DIC.epsYY.stage',num2str(i),');']) 
    epsYY_DIC2FEM_aux.Method = Interpolation;
    epsYY_DIC2FEM_aux.ExtrapolationMethod = Extrapolation;
    
    epsYY_FEM2DIC = epsYY_DIC2FEM_aux(x_FEA,y_FEA);
    
    eval(['DIC2FEM.eyy.stage',num2str(i),'= epsYY_FEM2DIC;']) 
    
    %epsxY  
    eval(['epsXY_DIC2FEM_aux = scatteredInterpolant(DIC.x.stage0,DIC.y.stage0,DIC.epsXY.stage',num2str(i),');']) 
    epsXY_DIC2FEM_aux.Method = Interpolation;
    epsXY_DIC2FEM_aux.ExtrapolationMethod = Extrapolation;
    
    epsXY_FEM2DIC = epsXY_DIC2FEM_aux(x_FEA,y_FEA);
    
    eval(['DIC2FEM.exy.stage',num2str(i),'= epsXY_FEM2DIC;']) 
end


%% Plotting meshes and fields 
disp('Creating the the plots...')
% Comparing meshes from experimental and numerical - nodes
figure('color',[1 1 1]);
plot(x_FEA_nodal,y_FEA_nodal,'r*'); 
hold on;
plot(DIC.x.stage0,DIC.y.stage0,'ko'); 
legend('FEM Mesh','DIC Mesh')
title('Mesh DIC VS FEM - Nodes') 
xlabel('x [mm]')  
ylabel('y [mm]')

% Comparing meshes from experimental and numerical - integration points
figure('color',[1 1 1]);
plot(x_FEA,y_FEA,'r*'); 
hold on;
plot(DIC.x.stage0,DIC.y.stage0,'ko'); 
legend('FEM Mesh','DIC Mesh')
title('Mesh DIC VS FEM - Integration points') 
xlabel('x [mm]')  
ylabel('y [mm]')

% ============================================== %
% ux: displacement components 
% ============================================== %
figure('color', [1 1 1]);        
scatter3(DIC.x.stage0,DIC.y.stage0,DIC.ux.stage137);hold on
plot3(x_FEA_nodal,y_FEA_nodal,DIC2FEM.ux.stage137,'ro','LineWidth',2);
legend('DIC','Interpolated DIC data')
xlabel('x, [mm]');
ylabel('y, [mm]');
zlabel('Ux');

% ============================================== %
% uy: displacement components 
% ============================================== %
figure('color', [1 1 1]);        
scatter3(DIC.x.stage0,DIC.y.stage0,DIC.uy.stage137);hold on
plot3(x_FEA_nodal,y_FEA_nodal,DIC2FEM.uy.stage137,'ro','LineWidth',2);
legend('DIC','Interpolated DIC data')
xlabel('x, [mm]');
ylabel('y, [mm]');
zlabel('Uy');

% ============================================== %
% epsX: strain components 
% ============================================== %
figure('color', [1 1 1]);        
scatter3(DIC.x.stage0,DIC.y.stage0,DIC.epsXX.stage137);hold on
plot3(x_FEA,y_FEA,DIC2FEM.exx.stage137,'ro','LineWidth',2);
legend('DIC','Interpolated DIC data')
xlabel('x, [mm]');
ylabel('y, [mm]');
zlabel('eps XX');

% ============================================== %
% epsY: strain components 
% ============================================== %
figure('color', [1 1 1]);        
scatter3(DIC.x.stage0,DIC.y.stage0,DIC.epsYY.stage137);hold on
plot3(x_FEA,y_FEA,DIC2FEM.eyy.stage137,'ro','LineWidth',2);
legend('DIC','Interpolated DIC data')
xlabel('x, [mm]');
ylabel('y, [mm]');
zlabel('eps YY');

% ============================================== %
% epsS: strain components 
% ============================================== %
figure('color', [1 1 1]);        
scatter3(DIC.x.stage0,DIC.y.stage0,DIC.epsXY.stage137);hold on
plot3(x_FEA,y_FEA,DIC2FEM.exy.stage137,'ro','LineWidth',2);
legend('DIC','Interpolated DIC data')
xlabel('x, [mm]');
ylabel('y, [mm]');
zlabel('eps XY');

%% Displacement shifting for encastre
disp('Performing displacement shifting...')
%indexes of bottom nodes
ind_encastre = find(y_FEA_nodal == min(y_FEA_nodal));
%indexes of top nodes
ind_top = find(y_FEA_nodal == max(y_FEA_nodal));

ux_shift=struct;
uy_shift = struct;

for l=1:n_stages
        eval(['ux_shift_aux(',num2str(l),',1)= mean(DIC2FEM.ux.stage',num2str(l),'(ind_encastre));'])
        eval(['uy_shift_aux(',num2str(l),',1)= mean(DIC2FEM.uy.stage',num2str(l),'(ind_encastre));'])

        eval(['ux_shift.stage',num2str(l),'= mean(DIC2FEM.ux.stage',num2str(l),'(ind_encastre));'])
        eval(['uy_shift.stage',num2str(l),'= mean(DIC2FEM.uy.stage',num2str(l),'(ind_encastre));'])
        eval(['DIC2FEM.ux_shifted.stage',num2str(l),'= DIC2FEM.ux.stage',num2str(l),'-ux_shift.stage',num2str(l),';'])
        eval(['DIC2FEM.uy_shifted.stage',num2str(l),'= DIC2FEM.uy.stage',num2str(l),'-uy_shift.stage',num2str(l),';'])        
        
        eval(['ux_BC(',num2str(l),',1)= mean(DIC2FEM.ux_shifted.stage',num2str(l),'(ind_top));'])
        eval(['uy_BC(',num2str(l),',1)= mean(DIC2FEM.uy_shifted.stage',num2str(l),'(ind_top));'])
end

%Plotting shifted displacements
% ============================================== %
% ux: displacement components 
% ============================================== %
figure('color', [1 1 1]);        
scatter3(DIC.x.stage0,DIC.y.stage0,DIC.ux.stage137);hold on
plot3(x_FEA_nodal,y_FEA_nodal,DIC2FEM.ux_shifted.stage137,'ro','LineWidth',2);
legend('DIC','Interpolated DIC data')
xlabel('x, [mm]');
ylabel('y, [mm]');
zlabel('Ux shifted');

% ============================================== %
% uy: displacement components 
% ============================================== %
figure('color', [1 1 1]);        
scatter3(DIC.x.stage0,DIC.y.stage0,DIC.uy.stage137);hold on
plot3(x_FEA_nodal,y_FEA_nodal,DIC2FEM.uy_shifted.stage137,'ro','LineWidth',2);
legend('DIC','Interpolated DIC data')
xlabel('x, [mm]');
ylabel('y, [mm]');
zlabel('Uy shifted');

% ============================================== %
% displacements shift
% ============================================== %
figure('color', [1 1 1]);        

plot(1:n_stages,-ux_shift_aux,'k','LineWidth',2);hold on
plot(1:n_stages,-uy_shift_aux,'r','LineWidth',2);hold on
legend('Ux','Uy')
xlabel('Timestep');
ylabel('Shifted displacement value [mm]');

% ============================================== %
% Displacements at encastre on last load step
% ============================================== %
figure('color', [1 1 1]);        

plot(x_FEA_nodal(ind_encastre),DIC2FEM.ux_shifted.stage137(ind_encastre),'k','LineWidth',2);hold on
plot(x_FEA_nodal(ind_encastre),DIC2FEM.uy_shifted.stage137(ind_encastre),'r','LineWidth',2);hold on
plot(x_FEA_nodal(ind_encastre),mean(DIC2FEM.ux_shifted.stage137(ind_encastre)).*ones(length(ind_encastre),1),'k--','LineWidth',2);hold on
plot(x_FEA_nodal(ind_encastre),mean(DIC2FEM.uy_shifted.stage137(ind_encastre)).*ones(length(ind_encastre),1),'r--','LineWidth',2);
legend('Ux shifted','Uy shifted','Average Ux shifted','Average Uy shifted')
xlabel('x [mm]');
ylabel('Displacement [mm] (at encastre y=0)');

% ============================================== %
% Displacements at top on last load step
% ============================================== %
figure('color', [1 1 1]);        

plot(x_FEA_nodal(ind_top),DIC2FEM.ux_shifted.stage137(ind_top),'k','LineWidth',2);hold on
plot(x_FEA_nodal(ind_top),DIC2FEM.uy_shifted.stage137(ind_top),'r','LineWidth',2);hold on
plot(x_FEA_nodal(ind_top),mean(DIC2FEM.ux_shifted.stage137(ind_top)).*ones(length(ind_top),1),'k--','LineWidth',2);hold on
plot(x_FEA_nodal(ind_top),mean(DIC2FEM.uy_shifted.stage137(ind_top)).*ones(length(ind_top),1),'r--','LineWidth',2);
legend('Ux shifted','Uy shifted','Average Ux shifted','Average Uy shifted')
xlabel('x [mm]');
ylabel('Displacement [mm] (at top y=80)');

% ============================================== %
% Ux BC
% ============================================== %
figure('color', [1 1 1]);        

plot(1:n_stages,ux_BC,'k','LineWidth',2);hold on
%plot(1:n_stages,uy_BC,'r','LineWidth',2);hold on
%legend('Ux')
xlabel('Timestep');
ylabel('U_x [mm]');

% ============================================== %
% Uy BC
% ============================================== %
figure('color', [1 1 1]);        

plot(1:n_stages,uy_BC,'k','LineWidth',2);hold on
%plot(1:n_stages,uy_BC,'r','LineWidth',2);hold on
%legend('Ux')
xlabel('Timestep');
ylabel('U_y [mm]');

%% Writting results %%%%%
% Create ExpData txt file
disp('Writting the interpolation results...')
% Increments for nodes
for i = 1:length(x_FEA_nodal)
    for j =1:n_stages-1
        inc_nodal (i,1) = 1;
        inc_nodal (i+j*length(x_FEA_nodal),1) = j+1;
    end
end

% Nodes
Nodes = [1:length(x_FEA_nodal)].';
for i = 1:n_stages-1
    Nodes_aux(:,1) = [1:length(x_FEA_nodal)].';
    Nodes = cat(1,Nodes,Nodes_aux);
end

% Increments for elements
for i = 1:length(Mesh_X)
    for j =1:n_stages-1
        inc (i,1) = 1;
        inc (i+j*length(Mesh_X),1) = j+1;
    end
end

% Elements 
Elements = [1:length(Mesh_X)].';
for i = 1:n_stages-1
    Elements_aux(:,1) = [1:length(Mesh_X)].';
    Elements = cat(1,Elements,Elements_aux);
end

% ux 
ux = DIC2FEM.ux.stage1;
for i = 2:n_stages
    eval(['ux_aux(:,1) = DIC2FEM.ux_shifted.stage',num2str(i),';'])   
    ux = cat(1,ux,ux_aux);
end

% uy 
uy = DIC2FEM.uy.stage1;
for i = 2:n_stages
    eval(['uy_aux(:,1) = DIC2FEM.uy_shifted.stage',num2str(i),';'])   
    uy = cat(1,uy,uy_aux);
end

% Exx 
EpsXX = DIC2FEM.exx.stage1;
for i = 2:n_stages
    eval(['EpsXX_aux(:,1) = DIC2FEM.exx.stage',num2str(i),';'])   
    EpsXX = cat(1,EpsXX,EpsXX_aux);
end

% Eyy
EpsYY = DIC2FEM.eyy.stage1;
for i = 2:n_stages
    eval(['EpsYY_aux(:,1) = DIC2FEM.eyy.stage',num2str(i),';'])   
    EpsYY = cat(1,EpsYY,EpsYY_aux);
end

% Exy 
EpsXY = DIC2FEM.exy.stage1;
for i = 2:n_stages
    eval(['EpsXY_aux(:,1) = DIC2FEM.exy.stage',num2str(i),';'])   
    EpsXY = cat(1,EpsXY,EpsXY_aux);
end

% Writing displacement interpolated data
NameFile_nodal = [pwd,filesep,specimen,filesep,'ExpData_displ.inp'];

fid_nodal = fopen(NameFile_nodal,'w');
fprintf(fid_nodal, 'inc,Node,Ux,Uy \n'); % header
% fclose(fid);
ToWrite_from_inc0_nodal = [inc_nodal,Nodes,ux(:,1),uy(:,1)].';
formatSpec_nodal = ['%d,%d,%.6f,%.6f \n'];
fprintf(fid_nodal,formatSpec_nodal, ToWrite_from_inc0_nodal);
status_nodal = fclose(fid_nodal);

% Writing strain interpolated data
NameFile = [pwd,filesep,specimen,filesep,'ExpData_strain.inp'];

fid = fopen(NameFile,'w');
fprintf(fid, 'inc,Element,Epsxx,Epsyy,Epsxy \n'); % header
% fclose(fid);
ToWrite_from_inc0 = [inc,Elements,EpsXX(:,1),EpsYY(:,1),EpsXY(:,1)].';
formatSpec = ['%d,%d,%.6f,%.6f,%.6f \n'];
fprintf(fid,formatSpec, ToWrite_from_inc0);
status = fclose(fid);

% Writing shifted BC
NameFile = [pwd,filesep,specimen,filesep,'BC_top.csv'];
increments = [1:n_stages].';
fid = fopen(NameFile,'w');
fprintf(fid, 'inc,ux,uy \n'); % header
% fclose(fid);
ToWrite_from_inc0 = [increments,ux_BC,uy_BC].';
formatSpec = ['%d,%.6f,%.6f \n'];
fprintf(fid,formatSpec, ToWrite_from_inc0);
status = fclose(fid);

disp('Done')
toc
