% Calculate anisotropies from azimuthally energy distribution dEdydphip
% tables from various events.
% Path: should be put in the parent folder of main/.

% Author:   Jia Liu, liu.2053@osu.edu

% History: 
% Jan. 27, 2014 First version.

% Have a clean start
clear all
tic 
% Specify the info for running
nodes_list = 1:10;
events_list = 1:40;   %number of events in one node
tau =1:1:10;     %matching time list
tau0 = 0.01;           %inital time of free-streaming
order_list = 2:3;  % anisotropy order list

nodes_total = length(nodes_list);       %number of nodes
mtimes_total = length(tau);
events_total = nodes_total* length(events_list); 
orders_total = length(order_list);

backup_length = 40;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define data structure to store the data for each node.
%   Structure:
% wn_cell{order}(event, matching time)
% wn_th_cell{order}(event, matching time)
% wn_fo_cell{order}(event, matching time)
wn_cell = cell(orders_total, 1);
wn_psin_cell = cell(orders_total, 1);
wn_th_numerator_cell = cell(orders_total, 1);   
wn_th_denominator_tbl = zeros(events_total, mtimes_total);
wn_fo_numerator_cell = cell(orders_total, 1);
wn_fo_denominator_tbl = zeros(events_total, mtimes_total);
% pre-allocate space
for i=1:orders_total
wn_cell{i} = zeros(events_total, mtimes_total);
wn_psin_cell{i} = zeros(events_total, mtimes_total);
wn_th_numerator_cell{i} = zeros(events_total, mtimes_total);
wn_fo_numerator_cell{i} = zeros(events_total, mtimes_total);
end

% prepare data file names
rootDir = fileparts(fileparts(pwd())); %grand-parent folder of the current directory
tableDir = fullfile(rootDir, 'fs_package','tables');
dataBaseDir = fullfile(rootDir, 'dataBase');
eventDataDir = fullfile(dataBaseDir, 'node%d', 'event_%d', '%g') ;
eccnInitDir = fullfile(dataBaseDir, 'node%d');

% read in phip Gaussian Table:
phip_th = load(fullfile(tableDir,'phip_gauss_table.dat'));  %100 points Gaussian table, (phip, w(phip))
phip_fo = load(fullfile(tableDir, 'phi_gauss_table.dat'));  %48 points Gaussian table

% Loop over nodes
events_count = 1;
for i=1:nodes_total
    % load in initial eccentricity file
    eccn_init_file = fullfile(sprintf(eccnInitDir, nodes_list(i)), 'Epx_initial.dat');
    eccn_init_data = load(eccn_init_file);
    eccn_init_data = eccn_init_data(1:length(events_list), :);   % suppose this node hasn't been finished.
    % Loop over events per node
    for j=1:length(events_list);
        for k=1:mtimes_total
            eventDataDir_now = sprintf(eventDataDir, nodes_list(i), events_list(j), tau(k));
            dEdydphip_th_file =  fullfile(eventDataDir_now, 'dEdydphipThermal.dat'); 
            dEdydphip_fo_file = fullfile(eventDataDir_now, 'dEdydphipFO.dat');
            % safty check: both file must exists
            if(exist(dEdydphip_th_file, 'file')==0)
                disp(['File does not exist: ', 'dEdydphipThermal.dat']);
                return
            elseif(exist(dEdydphip_fo_file, 'file')==0)
                 %disp(['File does not exist: ', 'dEdydphipFO.dat']);
                 dEdydphip_fo_data = zeros(length(phip_fo),1);  
            else
                dEdydphip_fo_data = load(dEdydphip_fo_file);               
            end
            dEdydphip_th_data = load(dEdydphip_th_file);
            
            % denominators of wn, independent of order
            wn_th_denominator_tbl(events_count, k) = sum(dEdydphip_th_data.*phip_th(:,2)); %\int dE/dydphip dphip=\sum dEdydphip*w(phip)
            wn_fo_denominator_tbl(events_count, k) = sum(dEdydphip_fo_data.*phip_fo(:,2)); 
            % Loop over all orders
            for iorder=1:orders_total
                %eccn_phin = eccn_init_data(j, 2+2*(order_list(iorder)+1)); % applies to new 'Epx_initial.dat' format.
                eccn_phin = eccn_init_data(j, 3+order_list(iorder));  % only applies to old 'Epx_initial.dat' and only ecc2 and ecc3 are there.
                eccn_phin = 0 ;%debug
                
                wn_th_numerator_real = sum(dEdydphip_th_data.*cos(order_list(iorder).*phip_th(:,1)-eccn_phin).*phip_th(:,2));
                wn_th_numerator_img = sum(dEdydphip_th_data.*sin(order_list(iorder).*phip_th(:,1)-eccn_phin).*phip_th(:,2));
                
                wn_fo_numerator_real = sum(dEdydphip_fo_data.*cos(order_list(iorder).*phip_fo(:,1)-eccn_phin).*phip_fo(:,2));
                wn_fo_numerator_img = sum(dEdydphip_fo_data.*sin(order_list(iorder).*phip_fo(:,1)-eccn_phin).*phip_fo(:,2));
                
                wn_cell{iorder}(events_count, k) = sqrt((wn_th_numerator_real+wn_fo_numerator_real)^2 ...
                    +(wn_th_numerator_img+wn_fo_numerator_img)^2)./abs(wn_th_denominator_tbl(events_count, k) ...
                        +wn_fo_denominator_tbl(events_count, k));
                wn_psin_cell{iorder}(events_count, k) = atan2(wn_th_numerator_img+wn_fo_numerator_img, ...
                    wn_th_numerator_real+wn_fo_numerator_real)/order_list(iorder);
                % save the intermediate data
                wn_th_numerator_cell{iorder}(events_count, k) = sqrt(wn_th_numerator_real^2 ...
                    +wn_th_numerator_img^2);
                wn_fo_numerator_cell{iorder}(events_count, k) = sqrt(wn_fo_numerator_real^2 ...
                    +wn_fo_numerator_img^2);
            end %<-> order
        end  %<-> matching time
        events_count = events_count +1;
        % save temperare files
        if(~mod(events_count, backup_length))
            tempFileName = 'wn_backup.mat';
            save(tempFileName);
        end
    end %<-> events per node
end %<-> nodes

% post processing: doing statistics

% clean temp backup file
tempFileName='wn_backup.mat';
if(exist(tempFileName, 'file'))
    delete(tempFileName);
end

% save all data to file
saveFileName =sprintf('wn_data_%devents.mat', events_total);
save(saveFileName)   % test mode, save all variables
% closure
disp(['All calculations finish! Save to file: ', saveFileName]);
toc
