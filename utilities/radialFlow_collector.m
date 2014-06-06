% Calculate radial flow on the freeze-out surface.
% v_aver = u^\mu d^3\sigma_\mu*v/u^\mu d^3\sigma_\mu

% Author: Jia Liu
 
clear all
disp('radialFlow_collector: begins running---->')
% run parameters
events_list = 1:40;   %number of events in one node
nodes_list = 1:10;
nodes_total=length(nodes_list);
events_total = length(events_list)*nodes_total; 
backup_length = 100;
tau=1:1:10;

% define directory structure
rootDir = pwd();
dataBaseDir = fullfile(rootDir,'nodes_data');
eventDataDir = fullfile(dataBaseDir, 'node%d', 'event_%d', '%g') ;
                                    %sprintf(eventDataDir, node_num, event_num, matchingTime_dir)

% pre-allocate space
vaver_tbl = zeros(events_total, length(tau));

% Loop over all the nodes, events and matching times
idx = 1;
for i=1:length(nodes_list)
    for k=1:length(events_list)
        for j=1:length(tau)            
            %processing one event and one matching time:
            %  load in files and reconstruct Txx and Tyy
            eventDataDir_now = sprintf(eventDataDir, nodes_list(i), events_list(k), tau(j));
            decData_file = fullfile(eventDataDir_now, 'decdat2.dat');
            decData_fo = load(decData_file);

            % format: Time,DA0,DA1,DA2,VZCM,
            %              VRCM, Ed*HbarC,BN, Temp*HbarC,BAMU,
            %              SMU, PDec2,  CPi33*HbarC, CPi00*HbarC,CPi01*HbarC,
            %              CPi02*HbarC, CPi11*HbarC,CPi12*HbarC,CPi22*HbarC,  CPPI
            if(isempty(decData_fo))
               vaver_tbl(idx, j) = 0;
                continue
            else
                gamma_fo = 1./sqrt(1-decData_fo(:,5).^2 - decData_fo(:,6).^2);
                ux_fo = decData_fo(:,5).*gamma_fo;
                uy_fo = decData_fo(:,6).*gamma_fo;   
                vr_fo = sqrt(decData_fo(:,5).^2+decData_fo(:,6).^2);
                udsigma = decData_fo(:,1).*(gamma_fo.*decData_fo(:,2)+ux_fo.*decData_fo(:,3)...
                    +uy_fo.*decData_fo(:,4));%./gamma_fo;   %u^\mu d^3\sigma_\mu=\tau_f*(u0*DA0+u1*DA1+u2*DA2) 
                vaver_tbl(idx, j) = sum(udsigma.*vr_fo)/sum(udsigma);
            end %if decdat2.dat is empty
        end  % per matching time
    disp([' event ', num2str(idx),...
    ' processing completes!']);
    idx = idx+1;
    end %per event
end
% disp(['node', num2str(nodes_list(i)), ' completes!'])

disp(['All calculation finished!']);

filename=sprintf('radialFlowFO_%devents.dat', events_total);
dlmwrite(filename, vaver_tbl, 'delimiter', '\t', 'precision', '%12.6f')  % test mode, save all variables
