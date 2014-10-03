% Name: paramSearchCollector.m

% Purpose: collect, process and make table for 3-parameter search
% 1. collect stable particle spectral and calculate mean pT;
% 2. collect stable particle v2ch, v3ch and make information  table.

clear all

% run parameters
taus_list = 0.4:0.2:2.0;
etas_list = 0.08:0.04:0.4;
tdec_list = 100:10:160;
event_id = 99;
particle_collection = {'pion_p', 'proton'};

% define folder structure
rootDir = pwd();
event_folder_pattern = fullfile(rootDir, 'dataBase', 'taus_%g', 'etas_%g', 'tdec_%g');
% pre-allocate space
resultTable = zeros(length(taus_list)*length(etas_list)*length(tdec_list), ...
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
for itaus = 1:length(taus_list)
    taus_now  = taus_list(itaus);    
    
    for ietas = 1:length(etas_list)
        etas_now = etas_list(ietas);
        
        for itdec = 1:length(tdec_list)
            tdec_now = tdec_list(itdec);
            
            % process one event
            event_folder_now = sprintf(event_folder_pattern, taus_now, etas_now, tdec_now);
            
            % record all charged particle v2 and v3
            charged_vn_file = fullfile(event_folder_now, 'Charged_integrated_vndata.dat');
            charged_vn_data=load(charged_vn_file);
            v2ch_now = charged_vn_data(3, end);
            v3ch_now = charged_vn_data(4, end);
            
            % calculate mean pT for specified particles
            meanpT_now = zeros(1, length(particle_collection));
            for iparticle = 1:length(particle_collection)
                spectra_file = fullfile(event_folder_now, ...
                    sprintf('%s_vndata.dat', particle_collection{iparticle}));
                spectra_data=load(spectra_file);
                pT_array = spectra_data(:,1);
                dN_array = spectra_data(:,3);%dN/(2pi pT dpT)
                % do interpolation on pT Gaussian points
                dN_forInterp = dN_array.*2.*pi.*pT_array; %dN/dpT
                dN_interped  = exp(interp1(pT_array, log(dN_forInterp), pT_new, 'linear'));
                % find mean pT
                totalpT_now  = sum(pT_new.*pT_weight.*dN_interped);
                totalNum_now=sum(pT_weight.*dN_interped);                
                meanpT_now(iparticle) = totalpT_now/totalNum_now;
            end % per stable particle
            
            resultTable(idx, :) = [event_id, taus_now, etas_now, tdec_now, v2ch_now, v3ch_now, meanpT_now]; 
            idx = idx+1; % next line
        end % per tdec
    end % per eta/s
    fprintf('%d / %d finished!\n', itaus, length(taus_list));
end % per tau_s

% save table
filename = '3paramSearchResult.dat';
dlmwrite(filename, resultTable, 'delimiter', '\t', 'precision', 5);
fprintf('---------finished!---------\n');