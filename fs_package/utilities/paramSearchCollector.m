% Name: paramSearchCollector.m

% Purpose: collect, process and make table for 3-parameter search
% 1. collect stable particle spectral and calculate mean pT;
% 2. collect stable particle v2ch, v3ch and make information  table.

clear all

% run parameters
node_list= 1:10;
event_list = 1:100;
event_id = 99;
particle_collection = {'pion_p', 'proton'};

% define folder structure
rootDir = pwd();
event_folder_pattern = fullfile(rootDir, 'dataBase','node%d', 'run_%d');
event_folder_v3_pattern = fullfile(rootDir, 'dataBase2', 'node%d', 'run_%d'); %file contains v3

% pre-allocate space
resultTable = zeros(length(node_list)*length(event_list), ...
    6+length(particle_collection)); % event_id, tau_s, eta/s, Tdec, V2ch, V3ch, pion meanPT, proton meanPT

% pT gaussian points, (0.3, 3) GeV/c
pT_new = [0.33401871,0.37228685,0.44019231,0.53615897, ...
             0.65793937,0.80267954,0.96698706,1.14701076, ...
             1.33853107,1.53705894,1.73794106,1.93646893, ...
             2.12798924,2.30801294,2.47232046,2.61706063, ...
             2.73884103,2.83480769,2.90271315,2.94098129];
pT_weight = [2.31183844e-02,5.32893766e-02,8.22570634e-02, ...
              1.09300723e-01,1.33783282e-01,1.55130323e-01, ...
              1.72841338e-01,1.86501143e-01,1.95789545e-01, ...
              2.00488821e-01,2.00488821e-01,1.95789545e-01, ...
              1.86501143e-01,1.72841338e-01,1.55130323e-01, ...
              1.33783282e-01,1.09300723e-01,8.22570634e-02, ...
              5.32893766e-02,2.31183844e-02];


% loop over all parameters
fprintf('---- start to collect data ----\n');
idx=1; % indicator of current line
for inode = 1:length(node_list)  
    for ievent = 1:length(event_list)       
        % process one event
        node_now = node_list(inode);
        event_now= event_list(ievent);
        event_folder_now = sprintf(event_folder_pattern, node_now, event_now);   
        event_folder_v3_now = sprintf(event_folder_v3_pattern, node_now, event_now);

        % record run paramters
        params_now = dlmread(fullfile(event_folder_now, 'params.dat'), ...
            '',1,0);
        taus_now = params_now(1); etas_now = params_now(2); tdec_now = params_now(3);
        
        % record all charged particle v2 and v3
        charged_v2_file = fullfile(event_folder_now, 'Charged_ptcut0510_eta_integrated_vndata.dat');
        charged_v2_data=load(charged_v2_file);
        v2ch_now = charged_v2_data(3, end);

        charged_v3_file = fullfile(event_folder_v3_now, 'Charged_ptcut0510_eta_integrated_vndata.dat');
        charged_v3_data=load(charged_v3_file);            
        v3ch_now = charged_v3_data(4, end);

        % calculate mean pT for specified particles
        meanpT_now = zeros(1, length(particle_collection));
        for iparticle = 1:length(particle_collection)
            spectra_file = fullfile(event_folder_now, ...
                sprintf('%s_vndata.dat', particle_collection{iparticle}));
            spectra_data=load(spectra_file);
            pT_array = spectra_data(:,1);
            dN_array = spectra_data(:,3);%dN/(2pi pT dpT)
            % do interpolation on pT Gaussian points
            dN_interped  = exp(interp1(pT_array, log(dN_array), pT_new, 'linear'));
            dN_interped  = dN_interped.*2.*pi.*pT_new;
            % find mean pT
            totalpT_now  = sum(pT_new.*pT_weight.*dN_interped);
            totalNum_now=sum(pT_weight.*dN_interped);                
            meanpT_now(iparticle) = totalpT_now/totalNum_now;
        end % per stable particle

        resultTable(idx, :) = [event_id, taus_now, etas_now, tdec_now, v2ch_now, v3ch_now, meanpT_now]; 
        idx = idx+1; % next line
    end % per run
    fprintf('%d / %d finished!\n', inode, length(node_list));
end % per node

% save table
filename = '3paramSearchResult.dat';
dlmwrite(filename, resultTable, 'delimiter', '\t', 'precision', '%14.8E');
fprintf('---------finished!---------\n');
