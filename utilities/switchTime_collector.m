% Collect hydro run data for testing how the hydro quantities changes for
% different matching times combine with switchTime_plotter.m to visulize data

% Revise history: 
%       Aug.20,2014   collect real and imaginary part of epsilon_p and
%       epsilon_p'
%       Jun.04, 2014   path updates and integrated to fs_package
%       Nov.12.2013   clean temporary variables and Pi00_scaled_*,
%            Pi1122_scaled_*, Pi11_*, Pi22_*, before saving to the data file.
%       Oct.1. 2013     Add two new outputs from "avg_points.dat"
%       Aug.30. 2013  first version

clear all
clc

% Specify the info for running
nodes_list = 1:10;   
events_per_node = 40;   %number of events in one node
tau =[0.4:0.2:1.8];    %matching time list
tau0 = 0.01;           %inital time of free-streaming
finding_accuracy = 1e-8;  %relate to matchine error of comparing two float numbers

% specify directory structure
rootDir = fileparts(fileparts(pwd())); % root directory is at grand-parent level

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nodes_total = length(nodes_list);       %number of nodes
mtimes_total = length(tau);
events_total = (nodes_total) * events_per_node; 

% define data structure:
% Time-TT0, EpsP, EdAvg, Vaver, Pi00Avg, Pi01Avg, Pi02Avg, 
%    Pi11Avg, Pi12Avg, Pi22Avg, Pi33Avg
% (1) since the hydro run time list is different for different t0, making
%        evo_time_cell{k} dynamical
% (2) momentum anisotropy <Txx-Tyy>/<Txx+Tyy>: 
%        EpsP{matching time}(events_num, hydro time)
% (3) radial flow v, or energy density ed
%       radialV{matching time}(events_num, hydro time)
% (4) shear tensor: pi00, Pi01, Pi02, Pi11, Pi12, Pi22, Pi33
% (5) scaled Pi00 and Pi11+Pi22
evo_time_cell = cell(mtimes_total, 1); 
TEpsP_cell = cell(mtimes_total, 1); %total momentum anisotropy
EpsP_cell = cell(mtimes_total, 1);
TEpsP_full_cell = cell(mtimes_total, 1); %total momentum anisotropy: counting imaginary part
EpsP_full_cell = cell(mtimes_total, 1);
EdAvg_cell = cell(mtimes_total, 1);
RadialFlow_cell = cell(mtimes_total, 1);
Pi00_cell = cell(mtimes_total, 1);
Pi01_cell = cell(mtimes_total, 1);
Pi02_cell = cell(mtimes_total, 1);
Pi11_cell = cell(mtimes_total, 1);
Pi12_cell = cell(mtimes_total, 1);
Pi22_cell = cell(mtimes_total, 1);
Pi33_cell = cell(mtimes_total, 1);
TEpsP_full_weighted_cell = cell(mtimes_total, 1);
TEpsP_inside_cell = cell(mtimes_total, 1);

EpsP_mag_cell = cell(mtimes_total, 1);
EpsP_angle_cell = cell(mtimes_total, 1);
TEpsP_mag_cell = cell(mtimes_total, 1);
TEpsP_angle_cell = cell(mtimes_total, 1);

for k=1:mtimes_total
    if tau(k)<0.59
        evo_time_cell{k} = [tau(k):0.002:0.588, 0.59:0.02:16.0];
    else 
       evo_time_cell{k} = tau(k):0.02:16.0; 
    end
    % assign spacee for EpsP at one matching time
    evo_time_total = length(evo_time_cell{k});
    EpsP_cell{k} = zeros(events_total, evo_time_total);
    TEpsP_cell{k} = zeros(events_total, evo_time_total);
    EpsP_full_cell{k} = zeros(events_total, evo_time_total);
    TEpsP_full_cell{k} = zeros(events_total, evo_time_total);    
    RadialFlow_cell{k} = zeros(events_total, evo_time_total);
    EdAvg_cell{k} = zeros(events_total, evo_time_total);
    Pi00_cell{k}  = zeros(events_total, evo_time_total);
    Pi01_cell{k}  = zeros(events_total, evo_time_total);
    Pi02_cell{k}  = zeros(events_total, evo_time_total);
    Pi11_cell{k}  = zeros(events_total, evo_time_total);
    Pi12_cell{k}  = zeros(events_total, evo_time_total);
    Pi22_cell{k}  = zeros(events_total, evo_time_total);
    Pi33_cell{k}  = zeros(events_total, evo_time_total);
    TEpsP_full_weighted_cell{k} = zeros(events_total, evo_time_total);
    TEpsP_inside_cell{k} = zeros(events_total, evo_time_total);
    
    EpsP_mag_cell{k} = zeros(events_total, evo_time_total);
    EpsP_angle_cell{k} = zeros(events_total, evo_time_total);
    TEpsP_mag_cell{k} = zeros(events_total, evo_time_total);
    TEpsP_angle_cell{k} = zeros(events_total, evo_time_total);
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read data from the database

events_count = 1;   
for i=1:nodes_total
    % Specify node folders
    node_folder = strcat(rootDir, '/dataBase/node', num2str(nodes_list(i)));
    % Loop over events
    for j=1:events_per_node
        event_folder = strcat(node_folder, '/', 'event_', num2str(j));
        % Loop over matching times
        for k=1:mtimes_total
            anisotrop_fileName = strcat(event_folder, '/',...
                num2str(tau(k)),'/', 'avg_points_corrected.dat'); %'avg_points.dat' is still in the folder for backup
            anisotrop_raw = load(anisotrop_fileName);
            
            TEpsP_fileName= strcat(event_folder, '/',...
                num2str(tau(k)),'/', 'anisotropy.dat');
            TEpsP_raw = load(TEpsP_fileName);
            % pre-processing: change the first column of anisotropy file
            % from tau-tau0 to tau
            anisotrop_raw(:,1) = anisotrop_raw(:,1) + tau(k);
            TEpsP_raw(:,1) = TEpsP_raw(:,1) + tau(k);
            if max(anisotrop_raw(:, 1)) > max(evo_time_cell{k})
                disp(['Estimate of maximum evolution time is wrong: estimate time ',...
                    num2str(max(evo_time_cell{k})), ...
                    ',   real final time', num2str(max(anisotrop_raw(:, 1)))]);
                quit force;
            end
            % calculate the position in the cell structure to store the
            % data
            time_min = anisotrop_raw(1,1);
            time_max = anisotrop_raw(end,1);
            time_min_idx = find(abs(evo_time_cell{k}-time_min)<finding_accuracy);
            time_max_idx = find(abs(evo_time_cell{k}-time_max)<finding_accuracy);
            % disp([num2str(time_min_idx), '  ',num2str(time_max_idx), '  ', num2str(length(anisotrop_raw))])
            % %debug

            % start recording from the second row since the 1st and 2nd row
            % are repetitive.            
            EpsP_cell{k}(events_count, time_min_idx:time_max_idx) = anisotrop_raw(2:end, 2);
            TEpsP_cell{k}(events_count, time_min_idx:time_max_idx) = TEpsP_raw(2:end, 5);
            EpsP_full_cell{k}(events_count, time_min_idx:time_max_idx) = anisotrop_raw(2:end, 12);
            TEpsP_full_cell{k}(events_count, time_min_idx:time_max_idx) = anisotrop_raw(2:end, 13);
            EdAvg_cell{k}(events_count, time_min_idx:time_max_idx) = anisotrop_raw(2:end, 3);
            RadialFlow_cell{k}(events_count, time_min_idx:time_max_idx) = anisotrop_raw(2:end, 4);
            Pi00_cell{k}(events_count, time_min_idx:time_max_idx) = anisotrop_raw(2:end, 5);
            Pi01_cell{k}(events_count, time_min_idx:time_max_idx) = anisotrop_raw(2:end, 6);
            Pi02_cell{k}(events_count, time_min_idx:time_max_idx) = anisotrop_raw(2:end, 7);
            Pi11_cell{k}(events_count, time_min_idx:time_max_idx) = anisotrop_raw(2:end, 8);
            Pi12_cell{k}(events_count, time_min_idx:time_max_idx) = anisotrop_raw(2:end, 9);
            Pi22_cell{k}(events_count, time_min_idx:time_max_idx) = anisotrop_raw(2:end, 10);
            Pi33_cell{k}(events_count, time_min_idx:time_max_idx) = anisotrop_raw(2:end, 11);
            TEpsP_full_weighted_cell{k}(events_count, time_min_idx:time_max_idx) = anisotrop_raw(2:end, 14);
            TEpsP_inside_cell{k}(events_count, time_min_idx:time_max_idx) = anisotrop_raw(2:end, 15); 
            
            EpsP_mag_cell{k}(events_count, time_min_idx:time_max_idx) = ...
                sqrt(anisotrop_raw(2:end, 16).^2 +anisotrop_raw(2:end, 17).^2); 
            EpsP_angle_cell{k}(events_count, time_min_idx:time_max_idx) = ...
                atan2( anisotrop_raw(2:end, 17), anisotrop_raw(2:end, 16)); 
            TEpsP_mag_cell{k}(events_count, time_min_idx:time_max_idx) = ...
                sqrt(anisotrop_raw(2:end, 18).^2 +anisotrop_raw(2:end, 19).^2); 
            TEpsP_angle_cell{k}(events_count, time_min_idx:time_max_idx)=  ...
                atan2( anisotrop_raw(2:end, 19), anisotrop_raw(2:end, 18)); 
        end
        events_count = events_count + 1;
    end
    disp([node_folder, ' has been recorded!'])
end

disp('Inital momentum anisotropy data has been read in!');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Post-processing
% find mean and standard deviation of momentum anisotropy
% since there are many zeros in the tail of each event, need to extrac
% non-zero elements of each column at first
EpsP_mean = cell(mtimes_total, 1);
EpsP_std = cell(mtimes_total, 1);
TEpsP_mean = cell(mtimes_total, 1);
TEpsP_std = cell(mtimes_total, 1);
EpsP_full_mean = cell(mtimes_total, 1);
EpsP_full_std = cell(mtimes_total, 1);
TEpsP_full_mean = cell(mtimes_total, 1);
TEpsP_full_std = cell(mtimes_total, 1);
RadialFlow_mean = cell(mtimes_total, 1);
RadialFlow_std = cell(mtimes_total, 1);
EdAvg_mean = cell(mtimes_total, 1);
EdAvg_std = cell(mtimes_total, 1);
Pi00_mean = cell(mtimes_total, 1);
Pi01_mean = cell(mtimes_total, 1);
Pi02_mean = cell(mtimes_total, 1);
Pi11_mean = cell(mtimes_total, 1);
Pi12_mean = cell(mtimes_total, 1);
Pi22_mean = cell(mtimes_total, 1);
Pi33_mean = cell(mtimes_total, 1);
TEpsP_full_weighted_mean = cell(mtimes_total, 1);
TEpsP_inside_mean = cell(mtimes_total, 1);

EpsP_mag_mean = cell(mtimes_total, 1);
TEpsP_mag_mean = cell(mtimes_total, 1);


Pi00_std = cell(mtimes_total, 1);
Pi01_std = cell(mtimes_total, 1);
Pi02_std = cell(mtimes_total, 1);
Pi11_std = cell(mtimes_total, 1);
Pi12_std = cell(mtimes_total, 1);
Pi22_std = cell(mtimes_total, 1);
Pi33_std = cell(mtimes_total, 1);
TEpsP_full_weighted_std = cell(mtimes_total, 1);
TEpsP_inside_std = cell(mtimes_total, 1);

EpsP_mag_std = cell(mtimes_total, 1);
TEpsP_mag_std = cell(mtimes_total, 1);

for i=1:mtimes_total
    evo_time_total = length(evo_time_cell{i});
    EpsP_mean{i} = zeros(1, evo_time_total);
    EpsP_std{i} = zeros(1, evo_time_total);
    TEpsP_mean{i} = zeros(1, evo_time_total);
    TEpsP_std{i} = zeros(1, evo_time_total);   
    
    EpsP_full_mean{i} = zeros(1, evo_time_total);
    EpsP_full_std{i} = zeros(1, evo_time_total);
    TEpsP_full_mean{i} = zeros(1, evo_time_total);
    TEpsP_full_std{i} = zeros(1, evo_time_total);     
    
    RadialFlow_mean{i} = zeros(1, evo_time_total);
    RadialFlow_std{i} = zeros(1, evo_time_total);
    EdAvg_mean{i} = zeros(1, evo_time_total);
    EdAvg_std{i} = zeros(1, evo_time_total);
    
    Pi00_mean{i} = zeros(1, evo_time_total);
    Pi01_mean{i} = zeros(1, evo_time_total);
    Pi02_mean{i} = zeros(1, evo_time_total);
    Pi11_mean{i} = zeros(1, evo_time_total);
    Pi12_mean{i} = zeros(1, evo_time_total);
    Pi22_mean{i} = zeros(1, evo_time_total);
    Pi33_mean{i} = zeros(1, evo_time_total);
    TEpsP_full_weighted_mean{i} = zeros(1, evo_time_total);
    TEpsP_inside_mean{i} = zeros(1, evo_time_total);
    EpsP_mag_mean{i} = zeros(1, evo_time_total);
    TEpsP_mag_mean{i} = zeros(1, evo_time_total);
    
    Pi00_std{i} = zeros(1, evo_time_total);
    Pi01_std{i} = zeros(1, evo_time_total);
    Pi02_std{i} = zeros(1, evo_time_total);
    Pi11_std{i} = zeros(1, evo_time_total);
    Pi12_std{i} = zeros(1, evo_time_total);
    Pi22_std{i} = zeros(1, evo_time_total);
    Pi33_std{i} = zeros(1, evo_time_total);
    TEpsP_full_weighted_std{i} = zeros(1, evo_time_total);
    TEpsP_inside_std{i} = zeros(1, evo_time_total);
    EpsP_mag_std{i} = zeros(1, evo_time_total);
    TEpsP_mag_std{i} = zeros(1, evo_time_total);    
    
    for j=1:evo_time_total
        % extract non-zero elements of momentum anisotropy
        epsp_nonzero=nonzeros(EpsP_cell{i}(:, j));
        if isempty(epsp_nonzero)  %all events end from this time
            break
        end
        EpsP_mean{i}(j) = mean(epsp_nonzero, 1);
        EpsP_std{i}(j) = std(epsp_nonzero, 0, 1);
        
        % extract non-zero elements of total momentum anisotropy
        tepsp_nonzero=nonzeros(TEpsP_cell{i}(:, j));
        if isempty(tepsp_nonzero)  %all events end from this time
            break
        end
        TEpsP_mean{i}(j) = mean(tepsp_nonzero, 1);
        TEpsP_std{i}(j) = std(tepsp_nonzero, 0, 1);

         % extract non-zero elements of momentum anisotropy
        epsp_full_nonzero=nonzeros(EpsP_full_cell{i}(:, j));
        if isempty(epsp_full_nonzero)  %all events end from this time
            break
        end
        EpsP_full_mean{i}(j) = mean(epsp_full_nonzero, 1);
        EpsP_full_std{i}(j) = std(epsp_full_nonzero, 0, 1);
 
        % extract non-zero elements of total momentum anisotropy
        tepsp_full_nonzero=nonzeros(TEpsP_full_cell{i}(:, j));
        if isempty(tepsp_full_nonzero)  %all events end from this time
            break
        end
        TEpsP_full_mean{i}(j) = mean(tepsp_full_nonzero, 1);
        TEpsP_full_std{i}(j) = std(tepsp_full_nonzero, 0, 1);        
        
         % extract non-zero elements of energy density
        ed_nonzero=nonzeros(EdAvg_cell{i}(:, j));
        if isempty(ed_nonzero)  %all events end from this time
            break
        end
        EdAvg_mean{i}(j) = mean(ed_nonzero, 1);
        EdAvg_std{i}(j) = std(ed_nonzero, 0, 1); 

        % extract non-zero elements of radial flow
        radialflow_nonzero=nonzeros(RadialFlow_cell{i}(:, j));
        if isempty(radialflow_nonzero)  %all events end from this time
            break
        end
        RadialFlow_mean{i}(j) = mean(radialflow_nonzero, 1);
        RadialFlow_std{i}(j) = std(radialflow_nonzero, 0, 1); 
        
         % extract non-zero elements of pi00
        pi00_nonzero=nonzeros(Pi00_cell{i}(:, j));
        if isempty(pi00_nonzero)  %all events end from this time
            break
        end
        Pi00_mean{i}(j) = mean(pi00_nonzero, 1);
        Pi00_std{i}(j) = std(pi00_nonzero, 0, 1); 

         % extract non-zero elements of pi00
        pi01_nonzero=nonzeros(Pi01_cell{i}(:, j));
        if isempty(pi01_nonzero)  %all events end from this time
            break
        end
        Pi01_mean{i}(j) = mean(pi01_nonzero, 1);
        Pi01_std{i}(j) = std(pi01_nonzero, 0, 1);   
        
         % extract non-zero elements of pi00
        pi02_nonzero=nonzeros(Pi02_cell{i}(:, j));
        if isempty(pi02_nonzero)  %all events end from this time
            break
        end
        Pi02_mean{i}(j) = mean(pi02_nonzero, 1);
        Pi02_std{i}(j) = std(pi02_nonzero, 0, 1);        
        
          % extract non-zero elements of pi00
        pi11_nonzero=nonzeros(Pi11_cell{i}(:, j));
        if isempty(pi11_nonzero)  %all events end from this time
            break
        end
        Pi11_mean{i}(j) = mean(pi11_nonzero, 1);
        Pi11_std{i}(j) = std(pi11_nonzero, 0, 1);       
 
         % extract non-zero elements of pi00
        pi12_nonzero=nonzeros(Pi12_cell{i}(:, j));
        if isempty(pi12_nonzero)  %all events end from this time
            break
        end
        Pi12_mean{i}(j) = mean(pi12_nonzero, 1);
        Pi12_std{i}(j) = std(pi12_nonzero, 0, 1);   

         % extract non-zero elements of pi00
        pi22_nonzero=nonzeros(Pi22_cell{i}(:, j));
        if isempty(pi22_nonzero)  %all events end from this time
            break
        end
        Pi22_mean{i}(j) = mean(pi22_nonzero, 1);
        Pi22_std{i}(j) = std(pi22_nonzero, 0, 1); 

         % extract non-zero elements of pi00
        pi33_nonzero=nonzeros(Pi33_cell{i}(:, j));
        if isempty(pi33_nonzero)  %all events end from this time
            break
        end
        Pi33_mean{i}(j) = mean(pi33_nonzero, 1);
        Pi33_std{i}(j) = std(pi33_nonzero, 0, 1);     
        
        % extract non-zero elements of scaled pi00
        pi00_scaled_nonzero=nonzeros(TEpsP_full_weighted_cell{i}(:, j));
        if isempty(pi00_scaled_nonzero)  %all events end from this time
            break
        end
        TEpsP_full_weighted_mean{i}(j) = mean(pi00_scaled_nonzero, 1);
        TEpsP_full_weighted_std{i}(j) = std(pi00_scaled_nonzero, 0, 1);         
        
        % extract non-zero elements of scaled pi11+pi22
        pi1122_scaled_nonzero=nonzeros(TEpsP_inside_cell{i}(:, j));
        if isempty(pi1122_scaled_nonzero)  %all events end from this time
            break
        end
        TEpsP_inside_mean{i}(j) = mean(pi1122_scaled_nonzero, 1);
        TEpsP_inside_std{i}(j) = std(pi1122_scaled_nonzero, 0, 1);          
 
        
        % extract non-zero elements of momentum anisotropy
        epsp_mag_nonzero=nonzeros(EpsP_mag_cell{i}(:, j));
        if isempty(epsp_mag_nonzero)  %all events end from this time
            break
        end
        EpsP_mag_mean{i}(j) = mean(epsp_mag_nonzero, 1);
        EpsP_mag_std{i}(j) = std(epsp_mag_nonzero, 0, 1);              

        
        % extract non-zero elements of total momentum anisotropy
        tepsp_mag_nonzero=nonzeros(TEpsP_mag_cell{i}(:, j));
        if isempty(tepsp_mag_nonzero)  %all events end from this time
            break
        end
        TEpsP_mag_mean{i}(j) = mean(tepsp_mag_nonzero, 1);
        TEpsP_mag_std{i}(j) = std(tepsp_mag_nonzero, 0, 1);            
    end
end
disp('Post processing complete!');

%clear epsp_nonzero  ed_nonzero radialflow_nonzero pi00_nonzero pi01_nonzero ...
%       pi02_nonzero pi11_nonzero pi12_nonzero pi22_nonzero pi33_nonzero;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clean results before saving
clear *_nonzero   % clean temp variables
clear Pi11_* Pi22_* 

% Save results
file_name = strcat('st_data_', num2str(events_total), 'events.mat');
save(file_name);
disp(['Data has been saved to file: ', file_name]);

