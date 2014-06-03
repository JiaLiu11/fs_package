% read in nodes data and process events. 

% Author: Jia Liu

% History: 
% Jun. 3, 2014  Updated path info, new format of "Epx_initial.dat", missing
%                       vn is assigned to NaN
% Nov. 5, 2013 Only read odd rows of "Epx_initial.dat"
% Aug. 9, 2013  First version.

% Purpose: Calculate event-average for vn, eccn and vn/eccn_init. 
% Prepare to be used by vn_plots_generator.m. 

% Have a clean start
clear all
clc

% Specify the info for running
nodes_list = 1:10;
events_per_node = 40;   %number of events in one node
tau =[1:1:10];     %matching time list


% specify directory structure
rootDir = fileparts(fileparts(pwd())); % root directory is at grand-parent level
iStable_location = fullfile(rootDir, 'fs_package', 'iS/tables');
pt_list = load(fullfile(iStable_location, 'pT_gauss_table.dat'));  %Gaussian points of pt
tau0 = 0.01;           %inital time of free-streaming
particles_list = [211, 321, 2212];  %thermal particle species
particles_name_list = {'Pion+', 'Kaon+', 'Proton'};
orders_vn_eccn_range = 2:3;   %range of anisotropy orders

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nodes_total = length(nodes_list);       %number of nodes
mtimes_total = length(tau);
particles_total = length(particles_list);
events_total = (nodes_total) * events_per_node; 
orders_total = length(orders_vn_eccn_range);
pt_list_total = length(pt_list);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define data structure to store the data for each node.
%   Structure:
%       vn_cell{particle}{order}(event#, matchingTime)
%       vn_pt_cell{particle}{order}{matchingTime}(event#, pt)
%       eccn_cell{order}(event#, matchingTime)
%       eccn_init_cell{order}(event#)
%       dn_dyptdpt2pi_cell{particle}{matchingTime}(event#, pt)
%       dn_dydpt2pi_cell{particle}{matchingTime}(event#, pt)
vn_cell = cell(particles_total, 1); 
vn_real_cell = cell(particles_total, 1); 
vn_img_cell = cell(particles_total, 1); 
angle_psin_cell= cell(particles_total, 1); 
eccn_cell = cell(orders_total, 1);
eccn_init_cell = cell(orders_total, 1);
vn_pt_cell = cell(particles_total, 1);
vn_pt_real_cell = cell(particles_total, 1);
vn_pt_img_cell = cell(particles_total, 1);
dn_ptdpt_cell = cell(particles_total, 1);
dn_dpt_cell = cell(particles_total, 1);
% Store particle multiplicity
particles_mp = cell(particles_total, 1); %thermal particle multiplicity
hydro_run_time = zeros(events_total, mtimes_total);
% Assign space for data structure
for i=1:particles_total
    vn_cell{i} = cell(orders_total, 1);
    vn_real_cell{i} = cell(orders_total, 1);
    vn_img_cell{i} = cell(orders_total, 1);
    angle_psin_cell{i}= cell(orders_total, 1);
    vn_pt_cell{i} = cell(orders_total, 1);
    vn_pt_img_cell{i} = cell(orders_total, 1);
    vn_pt_real_cell{i} = cell(orders_total, 1);
    dn_ptdpt_cell{i}= cell(mtimes_total, 1);
    dn_dpt_cell{i} = cell(mtimes_total, 1);
    particles_mp{i} = zeros(events_total, mtimes_total);
    % Assign space for data matrix
    for j=1:orders_total
        vn_cell{i}{j} = zeros(events_total, mtimes_total);
        angle_psin_cell{i}{j}= zeros(events_total, mtimes_total);
        vn_real_cell{i}{j} = zeros(events_total, mtimes_total);
        vn_img_cell{i}{j} = zeros(events_total, mtimes_total);
        vn_pt_cell{i}{j} = cell(mtimes_total, 1);
        vn_pt_real_cell{i}{j} = cell(mtimes_total, 1);
        vn_pt_img_cell{i}{j} = cell(mtimes_total, 1);
        for k=1:mtimes_total
            vn_pt_cell{i}{j}{k} = zeros(events_total, pt_list_total);
            vn_pt_real_cell{i}{j}{k} = zeros(events_total, pt_list_total);
            vn_pt_img_cell{i}{j}{k} = zeros(events_total, pt_list_total);
            dn_ptdpt_cell{i}{k} = zeros(events_total, pt_list_total);
            dn_dpt_cell{i}{k} = zeros(events_total, pt_list_total);
        end
    end    
end
for i=1:orders_total
    eccn_cell{i} = zeros(events_total, mtimes_total);
    eccn_init_cell{i} = zeros(events_total, 1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Read in ecc_init from fs
events_count = 1;
for i =1: nodes_total
    eccn_init_node_folder = strcat(rootDir,'/dataBase/node', num2str(nodes_list(i)));
    eccn_init_file = fullfile(eccn_init_node_folder, 'Epx_initial.dat');
    eccn_init_raw = load(eccn_init_file);
    for j = 1:orders_total
        eccn_init_cell{j}(events_count: events_count + events_per_node -1 , 1)...
            = eccn_init_raw(1: events_per_node ,2*orders_vn_eccn_range(j)+3);
    end
    events_count = events_count + events_per_node;
end
disp('Initial eccentricity has been read in!');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Loop over the nodes
events_count = 1;
for i=1:nodes_total
    % Specify node folders
    node_folder = strcat(rootDir, '/dataBase/node', num2str(nodes_list(i)));
    
    % Loop over events
    for j=1:events_per_node
        event_folder = strcat(node_folder, '/', 'event_', num2str(j));
        % Loop over matching times
        for k=1:mtimes_total
            % Loop over orders to get eccentricity data
            eccn_file = strcat(event_folder, '/', num2str(tau(k)), '/', ...
                'ecc-init-r_power-2.dat');
            eccn_raw = load(eccn_file);
            % Read off hydro run time
            hydro_run_aniso = load(strcat(event_folder, '/', num2str(tau(k)), '/', ...
                'anisotropy.dat'));
            hydro_run_time(events_count, k) = hydro_run_aniso(end, 1);
            % Read off eccentricities
            for n = 1: orders_total
                eccn_cell{n}(events_count, k) ...
                    = eccn_raw(orders_vn_eccn_range(n), end-1);
            end
            % Loop over particle type
        	for m=1:particles_total
                vn_files = strcat(event_folder, '/', num2str(tau(k)), '/', ...
                    'thermal_', num2str(particles_list(m)), ...
                    '_integrated_vndata.dat');
                if exist(vn_files, 'file')
                    vn_raw = load(vn_files);
                else
                    vn_raw = NaN(10,6);
                end
                vn_pt_file = strcat(event_folder, '/', num2str(tau(k)), '/', ...
                    'thermal_', num2str(particles_list(m)), ...
                    '_vndata.dat');
                if exist(vn_pt_file, 'file')
                    vn_pt_raw = load(vn_pt_file);
                else
                    vn_pt_raw = NaN(length(pt_list),30);
                end
                % Loop over orders to get vn data
                for n = 1:orders_total
                    vn_cell{m}{n}(events_count, k)...
                        = vn_raw(orders_vn_eccn_range(n)+1, end);
                    angle_psin_cell{m}{n}(events_count, k)...
                        = atan2(vn_raw(orders_vn_eccn_range(n)+1, end-1), ...
                                 vn_raw(orders_vn_eccn_range(n)+1, end-2));
                    vn_real_cell{m}{n}(events_count, k)...
                        = vn_raw(orders_vn_eccn_range(n)+1, end-2);
                    vn_img_cell{m}{n}(events_count, k)...
                        = vn_raw(orders_vn_eccn_range(n)+1, end-1);
                    vn_pt_cell{m}{n}{k}(events_count, :)  ...
                        = vn_pt_raw(:, 3+3*orders_vn_eccn_range(n));  %read in vn(pt) data
                    vn_pt_real_cell{m}{n}{k}(events_count, :)  ...
                        = vn_pt_raw(:, 1+3*orders_vn_eccn_range(n));  %read in vn(pt) data
                    vn_pt_img_cell{m}{n}{k}(events_count, :)  ...
                        = vn_pt_raw(:, 2+3*orders_vn_eccn_range(n));  %read in vn(pt) data
                    dn_ptdpt_cell{m}{k}(events_count, :)  ...
                        = vn_pt_raw(:, 3);  %read in dn spectral data
                    dn_dpt_cell{m}{k}(events_count, :)  ...
                        = vn_pt_raw(:, 3).*vn_pt_raw(:, 1);  %pt* dn/dyptdpt spectral data
                    %debug
                end
                % Get thermal particle multiplicity
                particles_mp{m}(events_count, k) ...
                    = vn_raw(1, 2);
            end            
        end
        % disp(['Event count: ', num2str(events_count)]);
        events_count = events_count +1;        
    end        
end

% Sanity check
events_count = events_count -1;
if events_count~=events_total
    disp(['Events count ', num2str(), ', it does not equal to events expected: ', num2str(events_total)]);
    quit
end

disp('Data readin success!');
disp(['Total events: ', num2str(events_count)]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Processing data to prepare tables for plot
% Calculate mean and standard deviation of vn, eccn and eccn_init
vn_mean = cell(particles_total, 1);
vn_std = cell(particles_total, 1);
vn_pt_mean = cell(particles_total, 1);
vn_pt_std = cell(particles_total, 1);

angle_psin_mean = cell(particles_total, 1);
angle_psin_std = cell(particles_total, 1);

vn_pt_real_mean = cell(particles_total, 1);
vn_pt_real_std = cell(particles_total, 1);
vn_pt_img_mean = cell(particles_total, 1);
vn_pt_img_std = cell(particles_total, 1);
dn_spectral_mean = cell(particles_total, 1);
dn_spectral_std = cell(particles_total, 1);
dn_dpt_spectral_mean = cell(particles_total, 1);
dn_dpt_spectral_std = cell(particles_total, 1);

eccn_mean = cell(orders_total, 1);
eccn_std = cell(orders_total, 1);
eccn_init_mean = cell(orders_total, 1);
eccn_init_std = cell(orders_total, 1);
% Ratio tables for <vn/eccn> and <vn/eccn_init>
vn_eccn_ratio_mean = cell(orders_total, 1);
vn_eccn_ratio_std = cell(orders_total, 1);
vn_eccn_init_ratio_mean = cell(orders_total, 1);
vn_eccn_init_ratio_std = cell(orders_total, 1);

for i = 1:particles_total
    vn_mean{i} = cell(orders_total, 1);
    vn_std{i} = cell(orders_total, 1);
    angle_psin_mean{i}= cell(orders_total, 1);
    angle_psin_std{i}= cell(orders_total, 1);
    vn_pt_mean{i} = cell(orders_total, 1);
    vn_pt_std{i} = cell(orders_total, 1);
    vn_pt_real_mean{i} = cell(orders_total, 1);
    vn_pt_real_std{i} =cell(orders_total, 1);
    
    vn_pt_img_mean{i} = cell(orders_total, 1);
    vn_pt_img_std{i} = cell(orders_total, 1);
    for j=1:orders_total
        vn_pt_mean{i}{j} = cell(mtimes_total,1);
        vn_pt_std{i}{j} = cell(mtimes_total,1);
        
        dn_spectral_mean{i} = cell(mtimes_total,1);
        dn_spectral_std{i} = cell(mtimes_total,1);
        dn_dpt_spectral_mean{i} = cell(mtimes_total,1);
        dn_dpt_spectral_std{i} = cell(mtimes_total,1);
        
        vn_pt_real_mean{i}{j} = cell(mtimes_total,1);
        vn_pt_real_std{i}{j} =cell(mtimes_total,1);
        vn_pt_img_mean{i}{j} = cell(mtimes_total,1);
        vn_pt_img_std{i}{j} = cell(mtimes_total,1);
        for k = 1:mtimes_total
            vn_pt_mean{i}{j}{k} = zeros(pt_list_total, 1);
            vn_pt_std{i}{j}{k} = zeros(pt_list_total, 1);
            
            dn_spectral_mean{i}{k} = zeros(pt_list_total, 1);
            dn_spectral_std{i}{k} = zeros(pt_list_total, 1);            
            dn_dpt_spectral_mean{i}{k} = zeros(pt_list_total, 1);
            dn_dpt_spectral_std{i}{k} = zeros(pt_list_total, 1);   
            
            vn_pt_real_mean{i}{j}{k} = zeros(pt_list_total, 1);
            vn_pt_real_std{i}{j}{k} =zeros(pt_list_total, 1);
            vn_pt_img_mean{i}{j}{k} = zeros(pt_list_total, 1);
            vn_pt_img_std{i}{j}{k} = zeros(pt_list_total, 1);
        end
    end
    vn_eccn_ratio_mean{i} = cell(orders_total, 1);
    vn_eccn_ratio_std{i} = cell(orders_total, 1);
    vn_eccn_init_ratio_mean{i} = cell(orders_total, 1);
    vn_eccn_init_ratio_std{i} = cell(orders_total, 1);
end

% Calculate ratio tables
for i=1:orders_total
    % Calculate eccentricity tables
    eccn_mean{i} = mean(eccn_cell{i}, 1);
    eccn_std{i} = std(eccn_cell{i}, 0, 1);
    eccn_init_mean{i} = mean(eccn_init_cell{i}, 1);
    eccn_init_std{i} = std(eccn_init_cell{i}, 0, 1);
    % Calculate vn_eccn ratio tables
    for j=1:particles_total
         vn_mean{j}{i} = mean(vn_cell{j}{i}, 1);
         vn_std{j}{i} = std(vn_cell{j}{i}, 0, 1);
         angle_psin_mean{j}{i} = mean(angle_psin_cell{j}{i}, 1);
         angle_psin_std{j}{i} = std(angle_psin_cell{j}{i}, 1);
         vn_eccn_ratio_temp = vn_cell{j}{i}./eccn_cell{i};
         vn_eccn_ratio_mean{j}{i} = mean(vn_eccn_ratio_temp, 1);
         vn_eccn_ratio_std{j}{i} = std(vn_eccn_ratio_temp, 0, 1);
         % Calculate <vn/eccn_init>
         vn_eccn_init_ratio_temp = zeros(events_total, mtimes_total);
         for k=1:mtimes_total
             vn_eccn_init_ratio_temp(:, k) = vn_cell{j}{i}(:, k)./eccn_init_cell{i};
             vn_pt_mean{j}{i}{k} = mean(vn_pt_cell{j}{i}{k}, 1);
             vn_pt_std{j}{i}{k} = std(vn_pt_cell{j}{i}{k}, 0, 1);
             vn_pt_real_mean{j}{i}{k} = mean(vn_pt_real_cell{j}{i}{k}, 1);
             vn_pt_real_std{j}{i}{k} = std(vn_pt_real_cell{j}{i}{k}, 0, 1);
             vn_pt_img_mean{j}{i}{k} = mean(vn_pt_img_cell{j}{i}{k}, 1);
             vn_pt_img_std{j}{i}{k} = std(vn_pt_img_cell{j}{i}{k}, 0, 1);   
             dn_spectral_mean{j}{k} = mean(dn_ptdpt_cell{j}{k}, 1);
             dn_spectral_std{j}{k} = std(dn_ptdpt_cell{j}{k}, 0, 1); 
             dn_dpt_spectral_mean{j}{k} = mean(dn_dpt_cell{j}{k}, 1);
             dn_dpt_spectral_std{j}{k} = std(dn_dpt_cell{j}{k}, 0, 1);
         end
         vn_eccn_init_ratio_mean{j}{i} = mean(vn_eccn_init_ratio_temp, 1);
         vn_eccn_init_ratio_std{j}{i} = std(vn_eccn_init_ratio_temp, 0, 1);
    end
end


% Process total momentum anisotropy e2p
e2pdataFile = strcat('e2p_data_',num2str(events_total),'events.mat');
if exist(e2pdataFile, 'file')
  load(e2pdataFile, 'e2p_cell');
  e2p_totalOrderLength = size(e2p_cell, 1);
    e2p_mean=cell(e2p_totalOrderLength, 1);
    e2p_std=cell(e2p_totalOrderLength, 1);

    for i=1:e2p_totalOrderLength
        e2p_mean{i}=zeros(mtimes_total, 1);
        e2p_std{i}=zeros(mtimes_total, 1);
    end

    for i=1:e2p_totalOrderLength
        e2p_mean{i}=mean(e2p_cell{i}, 1);
        e2p_std{i}=std(e2p_cell{i}, 0, 1);
    end
    disp('e2p file has been processed!');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save results
file_name = strcat('fs_package_', num2str(events_total), 'events.mat');
save(file_name);
disp(['Data has been saved to file: ', file_name]);

