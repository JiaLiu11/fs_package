%collect mean pt data
%meanPT file contains two different range of switching times: 2:1:10 and
%0.4:0.2:1.8. 40 events in each node. Modify this file to adjust to this
%change.
clear all

rootDir = fileparts(fileparts(pwd()));% grand-parent folder is the root folder
dataBaseDir = fullfile(rootDir, 'dataBase');
node_list = 1:1:10;
event_list = 1:1:40;
events_per_node = length(event_list);
events_total = length(node_list)*length(event_list);

matchingTime_list1 = 2:1:10;
matchingTime_list2 = 0.4:0.2:1.8;
matchingTime_total1 = length(matchingTime_list1);
matchingTime_total2 = length(matchingTime_list2);

%pre-allocate space
meanpt1 = zeros(events_total, matchingTime_total1); % for mt_time1
meanpt2 = zeros(events_total, matchingTime_total2); % for mt_time2
k1=1; k2=1;

%loop over nodes, events
for i=1:length(node_list)
    meanptFile = fullfile(dataBaseDir, sprintf('node%d', node_list(i)),'meanPT.dat');
    meanptData = load(meanptFile);
    %split data according to two matching time range
    meanptData1 = meanptData(1:events_per_node*matchingTime_total1,:);
    meanptData2 = meanptData(events_per_node*matchingTime_total1+1:end,:);
    %assign data to tables
    for j=1:matchingTime_total1:(events_per_node*matchingTime_total1)
        meanpt1(k1,:) = meanptData1(j:j+matchingTime_total1-1,end);
        k1=k1+1;
    end
    for j=1:matchingTime_total2:(events_per_node*matchingTime_total2)
        meanpt2(k2,:) = meanptData2(j:j+matchingTime_total2-1,end);
        k2=k2+1;
    end
end

%combine two tables
meanpt = cat(2, meanpt2, meanpt1);
%save to file
filename = sprintf('meanPT_%devents.dat',events_total);
dlmwrite(filename, meanpt, 'precision','%10.6f','delimiter','\t')
disp(['All complete! Data saved to ', filename])
