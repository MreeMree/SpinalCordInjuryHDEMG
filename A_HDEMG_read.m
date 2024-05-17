%% RPM study

%% 1a: select directory

clear all;

switch 4
    case 1 % if working from PC
        filePath = '\\its-rds.bham.ac.uk\rdsprojects\c\chious-arm-cycling-study\ISRT Subacute SCI\';
        addpath(filePath)
    case 2 % if working from Mac
        filePath = 'Volumes/chious-arm-cycling-study/ISRT Subacute SCI/';
        cd(filePath)
    case 3 % if working from BEAR portal (Linux)
        filePath = '/rds/projects/c/chious-arm-cycling-study/ISRT Subacute SCI';
        addpath(filePath)
    case 4
        filePath = input('Enter file path for participant HDEMG data: ',"s");
        addpath(filePath); cd(filePath);
end

%% 1b: select task to be analysed

switch 1
    case 1
        task = 'FR';
    case 2
        task = 'LR';
    case 3
        task = 'SS';
    case 4
        task = 'SSF';
    case 5
        task = 'RSF';
    case 6
        task = 'CMEP';
end

%% 1c: indicate participants and time points (skip if looking at pilot data)

%participants = {input("Enter participant number (e.g. P03): ","s")}; % set participant no
cd(filePath)

timePoints = {'01','02','03','04'};
timePoint = input('Enter time point (1,2,3 or 4): ');

file = dir(['*' task '*.otb+']); 
if size(file,1) > 1
    dir(['*' task '*.otb+'])
    fileChoice = input('Choose file, enter number: ');
else
    fileChoice = 1;
end

disp(file(fileChoice).name)

participant = {input('Enter participant ID (e.g. JP01): ',"s")};

disp(['Participant: ' participant{:}])
disp(['Time point: ' timePoints{timePoint}])

% checks participant selected matches with file folder
fileCheck = NaN(1,length(file(fileChoice).folder)-4);
for i = 1:length(file(fileChoice).folder)-3
    fileCheck(i) = strcmp(participant,file(fileChoice).folder(i:i+3));
    if sum(fileCheck) < 1
        disp('Warning: check if correct participant file.')
    else
    end
end

% checks timepoint selected matches with file folder
fileCheck = NaN(1,length(file(fileChoice).folder)-2);
for i = 1:length(file(fileChoice).folder)-2
    fileCheck(i) = strcmp(['_' timePoints{timePoint}],file(fileChoice).folder(i:i+2));
    if sum(fileCheck) < 1
        disp('Warning: check time point of file.')
    else
    end
end

%% 1d: untar

filePath=[ pwd '/' file(fileChoice).name(1:end-5) ];
untar(file(fileChoice).name,filePath);
cd(filePath)

%% 1e: look at monopolar signals

% parameters
xmlFiles = dir('*.xml');
paramInfo = readtable(xmlFiles(1).name); % info file should always be first because filename is in numbers
param = struct;
param.noChannels = max(paramInfo.ChannelStartIndexAttribute)+1; % finds number of channels from .xml file
param.gain = paramInfo.GainAttribute(1); % obtains gain
param.highPass = paramInfo.HighPassFilterAttribute; % obtains high pass filter for EMG channels
param.lowPass = paramInfo.LowPassFilterAttribute; % obtains low pass filter for EMG channels
param.sampleFreq = 2000;

% open sig file
sigFile = dir('*.sig')
sig = fopen(sigFile.name);
sig = fread(sig,[param.noChannels inf],'short');

% get seconds
seconds = (1:length(sig))/param.sampleFreq;

% CHECK TRIGGER, CHANNEL 65?
%HDsEMG_trigger=sig(65,:)-median(sig(65,:)); figure; plot(seconds,HDsEMG_trigger); title('Does this look like the trigger channel?');

% CHECK TRIGGER, CHANNEL 66?
%HDsEMG_trigger=sig(64,:)-median(sig(64,:)); figure; plot(seconds,HDsEMG_trigger); ('Trigger channel');

% filter
[B A]=butter(4,[30 500]/1000);
[Bacc Aacc]=butter(2,[.01 2]/75);

emg=filtfilt(B,A,sig(1:64,:)')';
emg = emg*0.286; %µV
% (can get dynamic range, no. bits, gain from xml)
% dynamic range = 5

%% 1f: process IMU data

% From the SQ+ manual:
% The first four additional channels (channels 65, 66, 67 and 68) are the 
% data relating to the IMU (Inertial Measurement Unit) and corresponding 
% respectively to the W, X, Y and Z quaternions derived from the 3 
% integrated sensors: accelerometer, gyroscope and magnetometer.

orientChans = sig(67:70,:);
subtitle = [participant{:} '_' timePoints{timePoint}];

angleOrder = 'xyz';
angles_3D = quat2eul(orientChans',angleOrder);
angles = unwrap(rad2deg(angles_3D));
T4 = angles;
%T4(isnan(T4)) = 0;
[Bacc,Aacc] = butter(2,[.01 2]/1000);
T4 = filtfilt(Bacc,Aacc,T4);
T4 = T4';

figure; 
ax1 = subplot(211); hold on; 
for i=1:4; plot(ax1,seconds,orientChans(i,:)); end
ax2 = subplot(212); hold on; 
for i=1:3; plot(ax2,seconds,T4(i,:)); end
legend(ax2,'X','Y','Z'); xlabel('Seconds'); title(['T4 - ' subtitle]);
hold off

%% 2a: plot all channels in same arrangement as electrodes

yMax = 700;

figure;
for i=1:64
    ax = subplot(5,13,i+1);
    plot(ax,emg(i,:))
    title(num2str(i))
    switch 2
        case 1 % off
        case 2 % on
            ylim([0-yMax yMax])
    end
end

%% OPTIONAL: check PRD of each channel vs mean of the rest of the channels
% monopolar
% if PRDs are all really low, e.g. below 30, don't bother marking as bad.


clear PRDall

switch 1
    case 1
        emgData = emgInt;
    case 2
        emgData = emg;
end

for i = 1:size(emgData,1)
    first = emgData(i,:);
    second = mean(emgData(setdiff(1:height(emgData),i),:));

    diffSum = ((second - first).^2)/(second.^2);
    PRDall(i) = sqrt(diffSum)*100; % percent residual difference
end

y = PRDall;

figure; bar(y); 
hold on; plot([1 64],[mean(y(y~=0))+2*std(y(y~=0)) mean(y(y~=0))+2*std(y(y~=0))])
title('PRD of each channel and mean of all other channels')
mean2sd = mean(y)+2*std(y);

disp('Channel monopolars with PRD > mean+2*sd = '); yInd = 1:length(y);
disp(yInd(y>mean2sd))
disp(['Mean PRD: ' num2str(mean(y))])
disp(['Median PRD: ' num2str(median(y))])

%% 2b: User Interface to observe and select bad channels (monopolar)

scaleFactor = 900;

%tabGroup = {'a1','a2','a3','a4'};
axesStruct = struct;

h = figure;
tg = uitabgroup(h);
%t1 = uitab(tg,"Title","Select bad channels");
t1 = uitab(tg,"Title","1-16"); axesStruct(1).ax = uiaxes(t1);
t2 = uitab(tg,"Title","17-32"); axesStruct(2).ax = uiaxes(t2);
t3 = uitab(tg,"Title","33-48"); axesStruct(3).ax = uiaxes(t3);
t4 = uitab(tg,"Title","49-64"); axesStruct(4).ax = uiaxes(t4);

cyc = 0; % used to cycle through all channels

for ind=1:4
    i = (1:16)+cyc;
    plot(axesStruct(ind).ax,seconds,emg(i,:)+scaleFactor*i'); hold on;
    ylim(axesStruct(ind).ax,[min(scaleFactor*i)-scaleFactor max(scaleFactor*i)+scaleFactor]); hold on;
    leg = cell(1,16); leg2 = (-15:0)+ind*16; for idx=1:16; leg{idx}=num2str(leg2(idx)); end
    legend(axesStruct(ind).ax,leg);
    hold off;
    cyc = cyc + 16;
end

%disp("hit return to bring up channel selection"); pause

% fig = uifigure;
% bgroup = uibuttongroup(fig,'Title','Select bad channels','Scrollable','on');
% % bg1 = uibutton(bgroup,"state","Text",'Channel 1');
% % bg2 = uibutton(bgroup,"state","Text",'Channel 2');
% % bg3 = uibutton(bgroup,"state","Text",'Channel 3');
% 
% bg = struct;
% for i=1:64
%     bg(i).button = uibutton(bgroup,"state","Text",['Channel ' num2str(i)]);
%     bg(i).button.Position = [100 i*30 100 22];
% end

%bg = uibuttongroup(t1,"Title","Ch 1");


%tb = uitogglebutton(t1,"Text","Ch 1");

badChannels = input("Enter bad channel numbers in square brackets, e.g. [16,18,42]: ");
%for 1:64

%% 2c: plot a selection of channels (OPTIONAL)

for i = 1:length(badChannels)
    chanToPlot = badChannels(i);
    %chanToPlot = [16];
    
        figure;
        for i=chanToPlot
            plot(seconds,emg(i,:)+mean(range(emg,2))*i); hold on;
            % the mean(range()) bit spaces the lines out appropriately
        end
        leg = cell(1,length(chanToPlot)); 
        for idx=1:length(chanToPlot); leg{idx}=num2str(chanToPlot(idx)); end
        legend(leg)
end

%% INTERPOLATE – select adjacent channels
% interpolate here if bad channels at this point, messes up differentials
% e.g. if dropping 3, interpolation with all surrounding channels

% Electrode arrangement for RPM study and ACEIII
% 12 25 38 51 64
% 11 24 37 50 63
% 10 23 36 49 62
%  9 22 35 48 61
%  8 21 34 47 60
%  7 20 33 46 59
%  6 19 32 45 58
%  5 18 31 44 57
%  4 17 30 43 56
%  3 16 29 42 55
%  2 15 28 41 54
%  1 14 27 40 53
%    13 26 39 52

% example
%badChannels = [28; 43];
%adjacentChannels = [42 44; 27 29];

% if bad channel is 16, 32, 48, don't interpolate as channel gets removed
% anyway when calculating differentials
% NOT TRUE - 17, 33, 49

%adjacentChannels = input("Enter the two adjacent channels to the bad channel: ");

adjacentChanStruct = struct;
for i = 1:length(badChannels)
    adjacentChanStruct(i).chans = input(['Enter the adjacent channels to interpolate channel ' num2str(badChannels(i)) ': ']);
end

%% ACTUALLY INTERPOLATE version 1:
% only works for 2 adjacent channels for each one

emgInt = emg;

for i = 1:length(badChannels)
    adjacentChannels = adjacentChanStruct(i).chans;
    if isempty(badChannels) == 0
        emgInt(badChannels(i),:) = NaN; % clears bad channel if not already NaN    
        % method 1
        signalOne = emgInt(adjacentChannels(1),:);
        signalTwo = emgInt(adjacentChannels(2),:);
        array1 = (1:length(signalOne))*2-1;
        array2 = (1:length(signalTwo))*2;
        array3 = (array1+0.5);
        x = [ array1 array2 ];
        v = [ signalOne signalTwo ];
        xq = [ array3 ];
        emgInt(badChannels(i),:) = interp1(x,v,xq);
        figure; plot(emgInt(badChannels(i),:)); 
        title(['Reconstructed channel (' num2str(badChannels(i)) ')'] )
    else
    end
end

%% ACTUALLY INTERPOLATE version 2:

emgInt = emg;

row = length(emg);
[XX,YY] = ndgrid(1:13,1:5);
electrodeGrid = reshape(0:64,13,5);
electrodeGrid = flip(electrodeGrid);
vq = NaN(1,row); % value array to store newly generated interpolated vals


for index = 1:length(adjacentChanStruct)
    
    switch 1
        case 1 
            adjacentChannels = adjacentChanStruct(index).chans;
        case 2 % interpolate with all other channels
            adjacentChannels = setdiff(1:64,badChannels);
    end

    col = length(adjacentChannels);
    
    x = NaN(1,col);
    y = NaN(1,col);
    v = NaN(1,col);
    for k=1:length(x); x(k) = XX(electrodeGrid == adjacentChannels(k)); end
    for k=1:length(y); y(k) = YY(electrodeGrid == adjacentChannels(k)); end
    
    xq = XX(electrodeGrid == badChannels(index));
    yq = YY(electrodeGrid == badChannels(index));

    for j = 1:row
        
        for k=1:length(v); v(k) = emgInt(adjacentChannels(k),j); end
        
        %vq(j) = interp2(x,y,v,xq,yq);
        vq(j) = griddata(x,y,v,xq,yq,"v4");
        
    end

    emgInt(badChannels(index),:) = vq;

    figure; plot(emgInt(badChannels(index),:)); 
    title(['Reconstructed channel (' num2str(badChannels(index)) ')'] )

end                

% x = x.*ones(row,col);
% y = y.*ones(row,col);
% v = [emgInt(adjacentChannels(1),:);emgInt(adjacentChannels(2),:)]';
% xq = xq.*ones(row,col);
% yq = yq.*ones(row,col);

F = griddedInterpolant(emgInt(1,:)',1:row);
F = griddedInterpolant(emgInt(1,:)',1:row);

%% ACTUALLY INTERPOLATE version 3:

emgInt = emg;
electrodeGrid = flip(reshape(0:64,[13 5]));

for i = 1:length(badChannels)
    emgInt(badChannels(i),:) = interpolateHDEMG(emgInt,badChannels(i),adjacentChanStruct(i).chans,electrodeGrid);
end

%% CHECK INTERPOLATION (OPTIONAL)
% check against nearby channels

nearbyChannels = input(['Enter a nearby channel for each bad ' ...
    'channel ('  num2str(badChannels)  '): ']);

% choose a channel that has NOT been used for interpolation
% a channel can be used twice, e.g. 63 and 64 both bad channels, use 46.

% 12 25 38 51 64
% 11 24 37 50 63
% 10 23 36 49 62
%  9 22 35 48 61
%  8 21 34 47 60
%  7 20 33 46 59
%  6 19 32 45 58
%  5 18 31 44 57
%  4 17 30 43 56
%  3 16 29 42 55
%  2 15 28 41 54
%  1 14 27 40 53
%    13 26 39 52

for i = 1:length(badChannels)
    channelChoice = i;
    
    interpolated = emgInt(badChannels(channelChoice),:);
    comparison = emgInt(nearbyChannels(channelChoice),:);
    
    diffSum = ((comparison - interpolated).^2)/(comparison.^2);
    PRD = sqrt(diffSum)*100; % percent residual difference
    
    figure;
    ax=subplot(141); bar(ax,PRD); title('PRD'); ylabel('%')
    ax=subplot(142); plot(ax,interpolated); title('Reconstructed – Linear horizontal');
    ax=subplot(143); plot(ax,comparison); title(['Comparison signal: ' num2str(nearbyChannels(channelChoice))]);
    ax=subplot(144); plot(ax,emg(badChannels(channelChoice),:)); title(['Original bad signal: ' num2str(badChannels(channelChoice))]);
end

%% 3a: calculate differentials

switch 2
    case 1 % calculate differentials and then interpolate if needed
        vert_mat = diff(emg);
        figure;
        for i=1:63
            ax = subplot(4,16,i);
            plot(ax,vert_mat(i,:))
            title(num2str(i))
        end
        vert_mat([16 32 48],:)=[]; % get rid of last electrode of each column
    case 2 % interpolate first
        vert_mat = diff(emgInt);
        vert_mat([12 25 38 51],:)=[]; % get rid of last electrode of each column
        %vert_mat([12 25 38 51],:)=[]; % get rid of last electrode of each column
end

figure;
for i=1:59
    ax = subplot(5,12,i+1);
    plot(ax,vert_mat(i,:))
    title(num2str(i))
    ylim([-100 100])
end

%% OPTIONAL: check PRD of each channel vs mean of the rest of the channels
% differential
% this is not that useful a section because by calculating the
% differentials, the channels are all naturally very different from one
% another now, so PRD values very high. outliers may still be noted and
% kept an eye on in following scripts, but no need to alter at this point.

clear PRDall

for i = 1:size(vert_mat,1)
    first = vert_mat(i,:);
    second = mean(vert_mat(setdiff(1:height(vert_mat),i),:));

    diffSum = ((second - first).^2)/(second.^2);
    PRDall(i) = sqrt(diffSum)*100; % percent residual difference
end

y = PRDall;

figure; bar(y); 
hold on; plot([1 60],[mean(y(y~=0))+2*std(y(y~=0)) mean(y(y~=0))+2*std(y(y~=0))])
title('PRD of each channel and mean of all other channels')
mean2sd = mean(y)+2*std(y);

disp('Channel differentials with PRD > mean+2*sd = '); yInd = 1:length(y);
disp(yInd(y>mean2sd))
disp(['Mean PRD: ' num2str(mean(y))])
disp(['Median PRD: ' num2str(median(y))])

% RMS

%% 3b i: save processed emg to .mat – select directory

monopolar = emgInt;
differential = vert_mat;
name = file(fileChoice).name;
%trigger = HDsEMG_trigger;
sampleFreq = param.sampleFreq;
IMU = T4;

%targetName = ['HDEMG_' baseline(i).participant '_' timePoints{timePoint}];
targetName = ['HDEMG_' task '_' participant{:} '_' timePoints{timePoint} '.mat'];

disp('Select folder to write file to then run next section.')

%% 3b ii: save processed emg to .mat – write file

response = input('Do you want to save in current folder? (y/n): ',"s");

if strcmp(response,'y')
    save(targetName,"differential","monopolar","name","seconds","sampleFreq","IMU")
    disp('DONE. Onto next script.')
else
    disp('Pick folder and run section again')
end

%% 3c: plot again for visual inspection (OPTIONAL)

scaleFactor = 100;

%tabGroup = {'a1','a2','a3','a4'};
axesStruct = struct;

h = figure;
tg = uitabgroup(h);
%t1 = uitab(tg,"Title","Select bad channels");
t1 = uitab(tg,"Title","1-15"); axesStruct(1).ax = uiaxes(t1);
t2 = uitab(tg,"Title","16-30"); axesStruct(2).ax = uiaxes(t2);
t3 = uitab(tg,"Title","31-45"); axesStruct(3).ax = uiaxes(t3);
t4 = uitab(tg,"Title","46-59"); axesStruct(4).ax = uiaxes(t4);

cyc = 0; % used to cycle through all channels

vert_mat(60,:) = 1;
for ind=1:4
    i = (1:15)+cyc;
    plot(axesStruct(ind).ax,seconds,vert_mat(i,:)+scaleFactor*i'); hold on;
    ylim(axesStruct(ind).ax,[min(scaleFactor*i)-scaleFactor max(scaleFactor*i)+scaleFactor]); hold on;
    leg = cell(1,15); leg2 = (-14:0)+ind*15; for idx=1:15; leg{idx}=num2str(leg2(idx)); end
    legend(axesStruct(ind).ax,leg);
    hold off;
    cyc = cyc + 15;
end
vert_mat(60,:) = [];

%% 3d: plot a selection of channels (OPTIONAL)

chanToPlot = [16,40];

    figure;
    for i=chanToPlot
        plot(seconds,vert_mat(i,:)+mean(range(vert_mat,2))*i); hold on;
        % the mean(range()) bit spaces the lines out appropriately
    end
    leg = cell(1,length(chanToPlot)); 
    for idx=1:length(chanToPlot); leg{idx}=num2str(chanToPlot(idx)); end
    legend(leg)

%% 3f: EXCLUDE bad channels

% 2 options 
% – interpolate bad channels (perhaps better for spatial distribution, and 
% if lots of bad channels), RMSE
% – just exclude the channels (if just 2-3 channels)

%% 4: write to .mat

