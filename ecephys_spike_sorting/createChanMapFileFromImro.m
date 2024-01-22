% main data
directory = 'Z:\WT_Sequences\2024_winter\Raw_Data\Neural_Traces\IMRO';
name = 'Payne_HPC';

%%
% read imro file
imroData = readmatrix(fullfile(directory, [name, '.imro']), 'FileType','text','Delimiter',')(','OutputType','string');
channel_list = zeros(length(imroData)-1,1);
shank = zeros(length(imroData)-1,1);
bank = zeros(length(imroData)-1,1);
for i = 1:length(imroData)-1
    c = double(split(imroData(i+1),' '));
    if i == length(imroData)-1
        c_end = split(imroData(i+1),' ');
        c_end = c_end{end};
        channel_list(i) = str2double(c_end(1:end-1)) + 1;
    else
        channel_list(i) = c(5) + 1;
    end
    shank(i) = c(2) + 1;
    bank(i) = c(3);
end
kcoords = shank;

%%
% create a channel map
pos = [32, 0]; % possible x-positions on shank
Nchannels = 384; % if there is a sync channel, otherwise 384
connected = true(Nchannels, 1);
connected(192, 1) = false; % exclude reference and sync channel
chanMap = [1:Nchannels]';
chanMap0ind = chanMap - 1;
xcoords = zeros(Nchannels,1);
ycoords = zeros(Nchannels,1);
for i = 1:length(channel_list)
    ycoords(i,1) = round((channel_list(i)+bank(i)*(Nchannels-1))/2) * 15;
    xcoords(i,1) = pos(mod(channel_list(i),2)+1) + (shank(i)-1)*250;
end

fs = 30000; % sampling frequency

%%
% save data in a mat file
save(fullfile(directory, [name, '.mat']), 'chanMap','connected', 'xcoords', 'ycoords', 'kcoords', 'chanMap0ind', 'name')
