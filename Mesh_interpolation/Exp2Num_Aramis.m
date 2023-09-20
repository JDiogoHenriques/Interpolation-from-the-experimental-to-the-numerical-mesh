%% Input

% J. Henriques, M. Conde 2023

% clear workSpace
clc; close all; clear all;

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

% ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: %
specimen = 'Dogbone_0'; % Dogbone_0 ; Dogbone_45 ; Dogbone_90
% ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: %

% struture variable
Exp = struct;
Exp.specimen = specimen;

% Get the current working directory
current_dir = pwd;

% Search for the experimental and numerical data
switch specimen
        case 'Dogbone_0'
        source_dir = [fileparts(pwd),filesep,'Experimental_Data',filesep,'ARAMIS',filesep,'Export_files_0']; 
        d = dir([source_dir, '\*.txt']);
        ID = 'dogBone_0';
        %Centroid data
        Mesh_X = importdata([current_dir,filesep,specimen,filesep,'COORD_1_int_points_all_el_undef_DB-0.csv']); %FEA Mesh X (centroid)
        Mesh_Y = importdata([current_dir,filesep,specimen,filesep,'COORD_2_int_points_all_el_undef_DB-0.csv']); %FEA Mesh Y (centroid)
        %Nodal data
        Mesh_nodal = importdata([current_dir,filesep,specimen,filesep,'DB-0.inp']); %FEA Mesh X (centroid)
        %Defining alignment settings
        Rot_Angle = 1.3; %degree anti clock-wise
        x_trans = 0 ; % translation in the x axis
        y_trans = 3; % translation in the y axis
        
        case 'Dogbone_45'
        source_dir = [fileparts(pwd),filesep,'Experimental_Data',filesep,'ARAMIS',filesep,'Export_files_45']; 
        d = dir([source_dir, '\*.txt']);      
        ID = 'dogBone_45';
        %Centroid data
        Mesh_X = importdata([current_dir,filesep,specimen,filesep,'COORD_1_int_points_all_el_undef_DB-45.csv']); %FEA Mesh X (centroid)
        Mesh_Y = importdata([current_dir,filesep,specimen,filesep,'COORD_2_int_points_all_el_undef_DB-45.csv']); %FEA Mesh Y (centroid)
        %Nodal data
        Mesh_nodal = importdata([current_dir,filesep,specimen,filesep,'DB-45.inp']); %FEA Mesh X (centroid)
        %Defining alignment settings
        Rot_Angle = 1.3; %degree anti clock-wise
        x_trans = 0.3; % translation in the x axis
        y_trans = 0.85; % translation in the y axis
        
        case 'Dogbone_90'
        source_dir = [fileparts(pwd),filesep,'Experimental_Data',filesep,'ARAMIS',filesep,'Export_files_90']; 
        d = dir([source_dir, '\*.txt']);    
        ID = 'dogBone_90';
        %Centroid data
        Mesh_X = importdata([current_dir,filesep,specimen,filesep,'COORD_1_int_points_all_el_undef_DB-90.csv']); %FEA Mesh X (centroid)
        Mesh_Y = importdata([current_dir,filesep,specimen,filesep,'COORD_2_int_points_all_el_undef_DB-90.csv']); %FEA Mesh Y (centroid)
        %Nodal data
        Mesh_nodal = importdata([current_dir,filesep,specimen,filesep,'DB-90.inp']); %FEA Mesh X (centroid)
        %Defining alignment settings
        Rot_Angle = 1.2; %degree anti clock-wise
        x_trans = 0; % translation in the x axis
        y_trans = 0; % translation in the y axis
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% FEA DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Processing for the numeric data
%Centroid data
Mesh_X = Mesh_X.data(2,2:end).';
Mesh_Y = Mesh_Y.data(2,2:end).';
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

% CoordX - All Stages
Mcomp = 'x_coord';
eval(['Exp.',Mcomp,' = [];'])

for j = 1:length(Exp.index)
    i = Exp.index(j);
    disp([specimen,' - reading cvs files: ',Mcomp,'...stage ',num2str(i)])
    filen = [source_dir,filesep,ID,'-Stage-0-',...
         num2str(mod(i,n_stages+1)),'.txt'];
    %Opening files 
    fid = fopen(filen);
    for k = 1:14
        fgetl(fid); % Read and discard the first 6 lines
    end
    data = textscan(fid, '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f', 'CollectOutput', true);
    matrix_data = [data{:}];
    fclose(fid);
    %Assigning matrixes to struct
    eval(['Exp.',Mcomp,'.stage',num2str(j),' = matrix_data(:,3);'])
end


% CoordY - All Stages
Mcomp = 'y_coord';
eval(['Exp.',Mcomp,' = [];'])

for j = 1:length(Exp.index)
    i = Exp.index(j);
    disp([specimen,' - reading cvs files: ',Mcomp,'...stage ',num2str(i)])
    filen = [source_dir,filesep,ID,'-Stage-0-',...
         num2str(mod(i,n_stages+1)),'.txt'];
    %Opening files 
    fid = fopen(filen);
    for k = 1:14
        fgetl(fid); % Read and discard the first 6 lines
    end
    data = textscan(fid, '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f', 'CollectOutput', true);
    matrix_data = [data{:}];
    fclose(fid);
    %Assigning matrixes to struct
    eval(['Exp.',Mcomp,'.stage',num2str(j),' = matrix_data(:,4);'])
end

% Ux - All Stages
Mcomp = 'ux';
eval(['Exp.',Mcomp,' = [];'])

for j = 1:length(Exp.index)
    i = Exp.index(j);
    disp([specimen,' - reading cvs files: ',Mcomp,'...stage ',num2str(i)])
    filen = [source_dir,filesep,ID,'-Stage-0-',...
         num2str(mod(i,n_stages+1)),'.txt'];
    %Opening files 
    fid = fopen(filen);
    for k = 1:14
        fgetl(fid); % Read and discard the first 6 lines
    end
    data = textscan(fid, '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f', 'CollectOutput', true);
    matrix_data = [data{:}];
    fclose(fid);
    %Assigning matrixes to struct
    eval(['Exp.',Mcomp,'.stage',num2str(j),' = matrix_data(:,9);'])
end

% Uy - All Stages
Mcomp = 'uy';
eval(['Exp.',Mcomp,' = [];'])

for j = 1:length(Exp.index)
    i = Exp.index(j);
    disp([specimen,' - reading cvs files: ',Mcomp,'...stage ',num2str(i)])
    filen = [source_dir,filesep,ID,'-Stage-0-',...
         num2str(mod(i,n_stages+1)),'.txt'];
    %Opening files 
    fid = fopen(filen);
    for k = 1:14
        fgetl(fid); % Read and discard the first 6 lines
    end
    data = textscan(fid, '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f', 'CollectOutput', true);
    matrix_data = [data{:}];
    fclose(fid);
    %Assigning matrixes to struct
    eval(['Exp.',Mcomp,'.stage',num2str(j),' = matrix_data(:,10);'])
end

% EpsX - All Stages
Mcomp = 'exx';
eval(['Exp.',Mcomp,' = [];'])

for j = 1:length(Exp.index)
    i = Exp.index(j);
    disp([specimen,' - reading cvs files: ',Mcomp,'...stage ',num2str(i)])
    filen = [source_dir,filesep,ID,'-Stage-0-',...
         num2str(mod(i,n_stages+1)),'.txt'];
    %Opening files 
    fid = fopen(filen);
    for k = 1:14
        fgetl(fid); % Read and discard the first 6 lines
    end
    data = textscan(fid, '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f', 'CollectOutput', true);
    matrix_data = [data{:}];
    fclose(fid);
    %Assigning matrixes to struct
    eval(['Exp.',Mcomp,'.stage',num2str(j),' = matrix_data(:,12)/100;'])
end

% EpsY - All Stages
Mcomp = 'eyy';
eval(['Exp.',Mcomp,' = [];'])

for j = 1:length(Exp.index)
    i = Exp.index(j);
    disp([specimen,' - reading cvs files: ',Mcomp,'...stage ',num2str(i)])
    filen = [source_dir,filesep,ID,'-Stage-0-',...
         num2str(mod(i,n_stages+1)),'.txt'];
    %Opening files 
    fid = fopen(filen);
    for k = 1:14
        fgetl(fid); % Read and discard the first 6 lines
    end
    data = textscan(fid, '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f', 'CollectOutput', true);
    matrix_data = [data{:}];
    fclose(fid);
    %Assigning matrixes to struct
    eval(['Exp.',Mcomp,'.stage',num2str(j),' = matrix_data(:,13)/100;'])
end

% EpsXY - All Stages
Mcomp = 'exy';
eval(['Exp.',Mcomp,' = [];'])

for j = 1:length(Exp.index)
    i = Exp.index(j);
    disp([specimen,' - reading cvs files: ',Mcomp,'...stage ',num2str(i)])
    filen = [source_dir,filesep,ID,'-Stage-0-',...
         num2str(mod(i,n_stages+1)),'.txt'];
    %Opening files 
    fid = fopen(filen);
    for k = 1:14
        fgetl(fid); % Read and discard the first 6 lines
    end
    data = textscan(fid, '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f', 'CollectOutput', true);
    matrix_data = [data{:}];
    fclose(fid);
    %Assigning matrixes to struct
    eval(['Exp.',Mcomp,'.stage',num2str(j),' = matrix_data(:,14)*2;'])
end

%% Processing DIC data %% 
% Shifting the origin to the corner of the specimen
for i = 1:n_stages    
    xdic_aux   = eval(['Exp.x_coord.stage',num2str(i)]);
    xdic_aux = xdic_aux-min(xdic_aux(:));
    eval(['Exp.x_coord.stage',num2str(i),' = xdic_aux;'])
    
    ydic_aux   = eval(['Exp.y_coord.stage',num2str(i)]);
    ydic_aux = ydic_aux-min(ydic_aux(:));
    eval(['Exp.y_coord.stage',num2str(i),' = ydic_aux;'])   
end

% Rotating data to align in the vertical direction
for i = 1:n_stages    
    xdic_aux   = eval(['Exp.x_coord.stage',num2str(i)]);
    ydic_aux   = eval(['Exp.y_coord.stage',num2str(i)]);
    
    xdic_rotaux = xdic_aux.*cosd(Rot_Angle)-ydic_aux.*sind(Rot_Angle);
    eval(['Exp.x_coord.stage',num2str(i),' = xdic_rotaux;'])
    
    ydic_rotaux = ydic_aux.*cosd(Rot_Angle)+xdic_aux.*sind(Rot_Angle);
    eval(['Exp.y_coord.stage',num2str(i),' = ydic_rotaux;'])   
end

% Calculating differences in the coordinate dimensions between exp and num
diff_x = max(Exp.x_coord.stage1) - max(x_FEA_nodal);
diff_y = max(Exp.y_coord.stage1) - max(y_FEA_nodal);

% Centering DIC data
for i = 1:n_stages    
	xdic_aux   = eval(['Exp.x_coord.stage',num2str(i)]);
	xdic_aux = xdic_aux-diff_x/2+x_trans;
	eval(['Exp.x_coord.stage',num2str(i),' = xdic_aux;'])
      
	ydic_aux   = eval(['Exp.y_coord.stage',num2str(i)]);
	ydic_aux = ydic_aux-diff_y/2+y_trans;
	eval(['Exp.y_coord.stage',num2str(i),' = ydic_aux;'])   
end

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
for i = 1:n_stages  
    
    %ux
    eval(['ux_DIC2FEM_aux = scatteredInterpolant(DIC.x.stage',num2str(i),',DIC.y.stage',num2str(i),',DIC.ux.stage',num2str(i),');']) 
    ux_DIC2FEM_aux.Method = 'nearest';
    ux_DIC2FEM_aux.ExtrapolationMethod = 'nearest';
    
    ux_FEM2DIC = ux_DIC2FEM_aux(x_FEA_nodal,y_FEA_nodal);
    
    eval(['DIC2FEM.ux.stage',num2str(i),'= ux_FEM2DIC;'])    
    
    %uy
    eval(['uy_DIC2FEM_aux = scatteredInterpolant(DIC.x.stage',num2str(i),',DIC.y.stage',num2str(i),',DIC.uy.stage',num2str(i),');']) 
    uy_DIC2FEM_aux.Method = 'nearest';
    uy_DIC2FEM_aux.ExtrapolationMethod = 'nearest';
    
    uy_FEM2DIC = uy_DIC2FEM_aux(x_FEA_nodal,y_FEA_nodal);
    
    eval(['DIC2FEM.uy.stage',num2str(i),'= uy_FEM2DIC;'])    
    
    %epsXX
    eval(['epsXX_DIC2FEM_aux = scatteredInterpolant(DIC.x.stage',num2str(i),',DIC.y.stage',num2str(i),',DIC.epsXX.stage',num2str(i),');']) 
    epsXX_DIC2FEM_aux.Method = 'nearest';
    epsXX_DIC2FEM_aux.ExtrapolationMethod = 'nearest';
    
    epsXX_FEM2DIC = epsXX_DIC2FEM_aux(x_FEA,y_FEA);
    
    eval(['DIC2FEM.exx.stage',num2str(i),'= epsXX_FEM2DIC;']) 
    
    %epsYY 
    eval(['epsYY_DIC2FEM_aux = scatteredInterpolant(DIC.x.stage',num2str(i),',DIC.y.stage',num2str(i),',DIC.epsYY.stage',num2str(i),');']) 
    epsYY_DIC2FEM_aux.Method = 'nearest';
    epsYY_DIC2FEM_aux.ExtrapolationMethod = 'nearest';
    
    epsYY_FEM2DIC = epsYY_DIC2FEM_aux(x_FEA,y_FEA);
    
    eval(['DIC2FEM.eyy.stage',num2str(i),'= epsYY_FEM2DIC;']) 
    
    %epsxY  
    eval(['epsXY_DIC2FEM_aux = scatteredInterpolant(DIC.x.stage',num2str(i),',DIC.y.stage',num2str(i),',DIC.epsXY.stage',num2str(i),');']) 
    epsXY_DIC2FEM_aux.Method = 'nearest';
    epsXY_DIC2FEM_aux.ExtrapolationMethod = 'nearest';
    
    epsXY_FEM2DIC = epsXY_DIC2FEM_aux(x_FEA,y_FEA);
    
    eval(['DIC2FEM.exy.stage',num2str(i),'= epsXY_FEM2DIC;']) 
end

%% Writting results %%%%%
% Create ExpData txt file

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
    eval(['ux_aux(:,1) = DIC2FEM.ux.stage',num2str(i),';'])   
    ux = cat(1,ux,ux_aux);
end

% uy 
uy = DIC2FEM.uy.stage1;
for i = 2:n_stages
    eval(['uy_aux(:,1) = DIC2FEM.uy.stage',num2str(i),';'])   
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

%% Plotting meshes and fields 
% Comparing meshes from experimental and numerical - integration points
figure('color',[1 1 1]);
plot(x_FEA_nodal,y_FEA_nodal,'r*'); 
hold on;
plot(DIC.x.stage1,DIC.y.stage1,'ko'); 
legend('FEM Mesh','DIC Mesh')
title('Mesh DIC VS FEM - Nodes') 
xlabel('x [mm]')  
ylabel('y [mm]')

% Comparing meshes from experimental and numerical - integration points
figure('color',[1 1 1]);
plot(x_FEA,y_FEA,'r*'); 
hold on;
plot(DIC.x.stage1,DIC.y.stage1,'ko'); 
legend('FEM Mesh','DIC Mesh')
title('Mesh DIC VS FEM - Integration points') 
xlabel('x [mm]')  
ylabel('y [mm]')

% ============================================== %
% ux: displacement components 
% ============================================== %
figure('color', [1 1 1]);        
scatter3(Exp.x_coord.stage146,Exp.y_coord.stage146,Exp.ux.stage146);hold on
plot3(x_FEA_nodal,y_FEA_nodal,DIC2FEM.ux.stage146,'ro','LineWidth',2);
legend('DIC','Interpolated DIC data')
xlabel('x, [mm]');
ylabel('y, [mm]');
zlabel('Ux');

% ============================================== %
% uy: displacement components 
% ============================================== %
figure('color', [1 1 1]);        
scatter3(Exp.x_coord.stage146,Exp.y_coord.stage146,Exp.uy.stage146);hold on
plot3(x_FEA_nodal,y_FEA_nodal,DIC2FEM.uy.stage146,'ro','LineWidth',2);
legend('DIC','Interpolated DIC data')
xlabel('x, [mm]');
ylabel('y, [mm]');
zlabel('Uy');

%%%%% Processing data for plots %%%%%%%%%
tri = delaunay(Exp.x_coord.stage146,Exp.y_coord.stage146);

% ============================================== %
% epsX: strain components 
% ============================================== %
figure('color', [1 1 1]);        
trisurf(tri,Exp.x_coord.stage146,Exp.y_coord.stage146,Exp.exx.stage146);hold on
plot3(x_FEA,y_FEA,DIC2FEM.exx.stage146,'ro','LineWidth',2);
legend('DIC','Interpolated DIC data')
xlabel('x, [mm]');
ylabel('y, [mm]');
zlabel('eps XX');

% ============================================== %
% epsY: strain components 
% ============================================== %
figure('color', [1 1 1]);        
trisurf(tri,Exp.x_coord.stage146,Exp.y_coord.stage146,Exp.eyy.stage146);hold on
plot3(x_FEA,y_FEA,DIC2FEM.eyy.stage146,'ro','LineWidth',2);
legend('DIC','Interpolated DIC data')
xlabel('x, [mm]');
ylabel('y, [mm]');
zlabel('eps YY');

% ============================================== %
% epsS: strain components 
% ============================================== %
figure('color', [1 1 1]);        
trisurf(tri,Exp.x_coord.stage146,Exp.y_coord.stage146,Exp.exy.stage146);hold on
plot3(x_FEA,y_FEA,DIC2FEM.exy.stage146,'ro','LineWidth',2);
legend('DIC','Interpolated DIC data')
xlabel('x, [mm]');
ylabel('y, [mm]');
zlabel('eps XY');

disp('Done')


%% Smoothing of the interpolated data
% 
% a =  griddata(DIC.x.stage146,DIC.y.stage146,DIC.ux.stage146); 
% 
% 
% figure('color', [1 1 1]);        
% scatter3(Exp.x_coord.stage146,Exp.y_coord.stage146,Exp.ux.stage146);hold on
% plot3(x_FEA_nodal,y_FEA_nodal,a,'ro','LineWidth',2);
% legend('DIC','Interpolated DIC data')
% xlabel('x, [mm]');
% ylabel('y, [mm]');
% zlabel('Ux');

