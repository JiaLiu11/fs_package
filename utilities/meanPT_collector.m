%collect mean pt data

clear all

rootDir = fileparts(fileparts(pwd()));
dataBaseDir = fullfile(rootDir, 'dataBase');
node_list = 1:1:10;
event_list = 1:1:40;
matchingTime_list = 1:1:10;
events_total = length(node_list)*length(event_list);
matchingTime_total = length(matchingTime_list);

%pre-allocate space
meanpt = zeros(events_total, matchingTime_total);
k=1;

% begin looping
for i=1:length(node_list)
    meanptFile = fullfile(dataBaseDir, sprintf('node%d', node_list(i)),'meanPT.dat');
    meanptData = load(meanptFile);
    for j=1:matchingTime_total:(length(event_list)*matchingTime_total)
        meanpt(k,:) = meanptData(j:j+9,end);
        k=k+1;
    end
end

%save to file
filename = sprintf('meanPT_%devents.dat',events_total);
dlmwrite(filename, meanpt, 'precision','%10.6f','delimiter','\t')
disp(['All complete! Data saved to ', filename])