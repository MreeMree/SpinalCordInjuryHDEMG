%% RPM study

%% 1a: select directory

clear all

switch 4
    case 1 % if working from PC
        filePath = '\\its-rds.bham.ac.uk\rdsprojects\c\chious-arm-cycling-study\ISRT Subacute SCI\';
        addpath([filePath 'Scripts']); addpath([filePath 'Processed'])
        option = 1;
    case 2 % if working from Mac
        filePath = 'Volumes/chious-arm-cycling-study/ISRT Subacute SCI/';
        addpath([filePath 'Scripts']); addpath([filePath 'Processed'])
        option = 2;
    case 3 % if working from BEAR portal (Linux)
        filePath = '/rds/projects/c/chious-arm-cycling-study/ISRT Subacute SCI';
        addpath([filePath 'Scripts']); addpath([filePath 'Processed'])
        option = 3;
    case 4 % custom file path
        filePath = input("Enter file path for Processed data: ","s");
        addpath(filePath);
        filePath = input("Enter file path for Scripts: ","s");
        addpath(filePath);
        %filePath = input("Enter file path for Delsys data for chosen participant: ","s");
        %addpath(filePath);
        %option = 4;
end

%% 1b: select task, participant and time point

switch 1
    case 1
        task = 'FR';
    case 2
        task = 'LR';
end

timePoints = {'01','02','03','04'};
participants = {input('Enter participant ID (e.g. JP01): ',"s")}; % select participant, one at a time
p = 1;

switch input('Enter timepoint (1-4): ') % select time point ( baseline / pre / post / washout )
    case 1
        timePoint = 1; % baseline
        filename = ['HDEMG_' task '_' participants{p} '_' timePoints{timePoint}  '.mat'];
        load(filename)
    case 2
        timePoint = 2; % pre-test
        filename = ['HDEMG_' task '_' participants{p} '_' timePoints{timePoint}  '.mat'];
        load(filename)
    case 3
        timePoint = 3;
        filename = ['HDEMG_' task '_' participants{p} '_' timePoints{timePoint}  '.mat'];
        load(filename)
    case 4
        timePoint = 4;
        filename = ['HDEMG_' task '_' participants{p} '_' timePoints{timePoint}  '.mat'];
        load(filename)
end

%HDsEMG_trigger = trigger;

% change HDsEMG_trigger to trigger for neatness

%% quickly process deltoid EMG for extra guidance for trial selection
% 
% deltoidEMG = delsys.AD_EMG4(delsys.X_s__9>=trigStartDelsysSeconds & delsys.X_s__9<=trigStopDelsysSeconds);
% deltoidSeconds = delsys.X_s__9(delsys.X_s__9>=trigStartDelsysSeconds & delsys.X_s__9<=trigStopDelsysSeconds);
% 
% % filter deltoidEMG
% [B A]=butter(4,[30 500]/1000);
% [Bacc Aacc]=butter(2,[.01 2]/75);
% 
% deltoidEMG=filtfilt(B,A,deltoidEMG')';
% deltoidEMG = deltoidEMG*0.286; %µV

%% 3a: plot IMU data with HDEMG
% 
% x_min = trigStartDelsysSeconds;
% x_max = max(secondsT4Crop);

subtitle = [participants{:} '_' timePoints{timePoint}];

yMax = 30;
channels = [1,9,5,25,19];
f = figure; set(f,'WindowStyle','docked')
%1
ax = subplot(311);
plot(ax,seconds,IMU(1,:)); hold on; %xlim([trigStartDelsysSeconds max(secondsT4Crop)]);
legend('X'); xlabel('Seconds'); title(['T4 - ' subtitle]); ylabel('Degrees')
%2
ax = subplot(312);
plot(ax,seconds,differential(channels,:)); hold on; %xlim([trigStartSeconds trigStopSeconds]); 
legend(num2str(channels)); title('Selection of HDEMG channels'); xlabel('Seconds')
%3
emgArray = 1:length(differential); % to easily index segmentTimes
ax = subplot(313);
plot(ax,emgArray,mean(differential,1)); hold on; %xlim([emgCropArray(1) emgCropArray(end)]);
ylim([0-yMax yMax])
legend(num2str(channels)); title('Mean of al EMG channels, custom y limits'); xlabel('Data point')
%4
% yMax_delt = 0.0001;
% ax = subplot(414);
% plot(ax,deltoidSeconds,deltoidEMG); ylim([0-yMax_delt yMax_delt]);
% xlim([min(deltoidSeconds) max(deltoidSeconds)])
% title('Anterior deltoid EMG (Delsys)'); xlabel('Seconds'); 


% define start and end point in seconds of 3 best segments

% start end
% start end
% start end

% old, in seconds: segmentTimes = [6.91875 11.07; 14.7487 18.4073; 21.924 25.245]; % these need to be for 

% plot EMG with highlight from segmentTimes

% plot Delsys big to select trial segments

%emgArray = 1:length(emg); % to easily index segmentTimes
%delsysCropArray = 1:length(T4CropX); % to easily index segmentTimes
h = figure; plot(emgArray,IMU(1,:)); set(h,'WindowStyle','docked')

%% 3b: select start and end of trials

%dxa = [27245,51616,67977,103574,132736,147932,276630,293860,307895]; % JP01_01
dxa = [103574,132736,147932,276630,293860,307895,409291,432609,442291]; % JP01_01
dxa = [12887,21865,33690,48182,55573,69233,75620,88041,97449]; % JP01_02
% first one JP didn't keep hand up and leant down so selected other trials

% dxa short for Delsys X annotated
dxa = input("Enter x values for the start, middle and end" + ...
    " of each of the three selected trials in square brackets: ");
segmentFull = [dxa(1),dxa(3),dxa(4),dxa(6),dxa(7),dxa(9)];
segmentConcentric = [dxa(2),dxa(3),dxa(5),dxa(6),dxa(8),dxa(9)];
segmentEccentric = [dxa(1),dxa(2),dxa(4),dxa(5),dxa(7),dxa(8)];

% % 1: segment EMG for whole movement 
% segmentDelsysPortion = segmentFull/delsysCropArray(end);
% segmentEMG = segmentDelsysPortion*emgCropArray(end);
% 
% % 2: segment EMG for eccentric portion of movement
% segmentDelsysPortion = segmentEccentric/delsysCropArray(end);
% segmentEMGEccentric = segmentDelsysPortion*emgCropArray(end);
% 
% % 3: segment EMG for concentric portion of movement
% segmentDelsysPortion = segmentConcentric/delsysCropArray(end);
% segmentEMGConcentric = segmentDelsysPortion*emgCropArray(end);


%% 3c: plot again with highlights 

yMax = 30;

switch 1
    case 1
        m = segmentFull; l = segmentFull;
    case 2
        m = segmentEccentric; l = segmentEccentric;
    case 3
        m = segmentConcentric; l = segmentConcentric;
end

channels = [5,7];
cs = struct; cs.Col = [0.2 0.2 0.4]; cs.Opa = 0.3; cs.Edg = 0; %colourSettings struct, colour, opacity, edge opacity
figure; 
%1
ax = subplot(211);
plot(ax,IMU(1,:)); hold on; %xlim([delsysCropArray(1) delsysCropArray(end)]);
for i=[1 3 5] % index start of each trial
    p=patch([l(i),l(i),l(i+1),l(i+1)],[ax.YLim(1),ax.YLim(2),ax.YLim(2),ax.YLim(1)],cs.Col); hold on;
    p.FaceAlpha = cs.Opa;
    p.EdgeAlpha = cs.Edg; 
end
legend('X'); xlabel('Datapoint'); title([participants{:} ' T4 - ' task ' | Time point: ' timePoints{timePoint}]); ylabel('Degrees')
%2
ax = subplot(212);
plot(ax,mean(differential,1)); hold on; %xlim([emgCropArray(1) emgCropArray(end)]);
for i=[1 3 5] % index start of each trial
    p=patch([m(i),m(i),m(i+1),m(i+1)],[ax.YLim(1),ax.YLim(2),ax.YLim(2),ax.YLim(1)],cs.Col); hold on;
    p.FaceAlpha = cs.Opa;
    p.EdgeAlpha = cs.Edg;
end
ylim([0-yMax yMax])
legend('Mean EMG signal from 59 differential signals'); xlabel('Datapoint'); 
title([participants{:} ' Erector Spinae - ' task ' | Time point: ' timePoints{timePoint}]); ylabel('EMG (µV)')

% %3
% ax = subplot(313);
% plot(ax,deltoidSeconds,deltoidEMG); %ylim([0-yMax_delt yMax_delt]);
% xlim([min(deltoidSeconds) max(deltoidSeconds)])
% title('Anterior deltoid EMG (Delsys)'); xlabel('Seconds');


%% 4: calculate and store trunk displacement for selected trials

trunkDisplacement = nan(1,3);
endTrunkAngles = nan(1,3);
T4_X = IMU(1,:);

for i=1:3
    trunkDisplacement(i) = T4_X( dxa(-1+3*i) ) - T4_X( dxa(-2+3*i) );
    endTrunkAngles(i) = T4_X( dxa(-1+3*i) );
end

meanTrunkDisplacement = abs(mean(trunkDisplacement));
disp(['Mean angular trunk displacement: ' num2str(meanTrunkDisplacement)])


%% 5d: 5c alternative
% 
% channels = [1,2,3];
% cs = struct; cs.Col = [0.2 0.2 0.4]; cs.Opa = 0.3; cs.Edg = 0; %colourSettings struct, colour, opacity, edge opacity
% figure; 
% %1
% ax = subplot(211);
% plot(ax,delsysCropArray,T4CropX); hold on; xlim([delsysCropArray(1) delsysCropArray(end)]);
% for i=[1 3 5] % index start of each trial
%     p=patch([segmentDelsys(i),segmentDelsys(i),segmentDelsys(i+1),segmentDelsys(i+1)],[ax.YLim(1),ax.YLim(2),ax.YLim(2),ax.YLim(1)],cs.Col); hold on;
%     p.FaceAlpha = cs.Opa;
%     p.EdgeAlpha = cs.Edg; 
% end
% legend('X'); xlabel('Datapoint'); title(['T4 - ' task]); ylabel(y_label_name)
% %2
% ax = subplot(212);
% plot(ax,emgCropArray,emgCrop(channels,:)); hold on; xlim([emgCropArray(1) emgCropArray(end)]);
% for i=[1 3 5] % index start of each trial
%     p=patch([segmentEMG(i),segmentEMG(i),segmentEMG(i+1),segmentEMG(i+1)],[ax.YLim(1),ax.YLim(2),ax.YLim(2),ax.YLim(1)],cs.Col); hold on;
%     p.FaceAlpha = cs.Opa;
%     p.EdgeAlpha = cs.Edg;
% end
% legend(num2str(channels)); xlabel('Datapoint'); title(['Erector Spinae - ' task]); ylabel('EMG (µV)')


%% 6a: save segmented data for trials

% EMG data from 3 best trials
% for i=1:3
% emgSegments

hdemg = struct;

for i = 1:3
    j = i*2;
    hdemg(i).full = differential(:,segmentFull(j-1):segmentFull(j));
    hdemg(i).eccentric = differential(:,segmentEccentric(j-1):segmentEccentric(j));
    hdemg(i).concentric = differential(:,segmentConcentric(j-1):segmentConcentric(j));
    %hdemg(i).emgCropArray = emgCropArray;
end

%% 6b: calculate RMS for each channel

% rmsEMG = nan(height(emgCrop),3);
% rmsEMG_ecc = nan(height(emgCrop),3);
% rmsEMG_con = nan(height(emgCrop),3);
% 
% for i = 1:3
%     %for j = 1:height(emgCrop) % number of channels left, 60 after differential calc 
%     rmsEMG(:,i) = rms(hdemg(i).full,2,"omitnan");
%     rmsEMG_ecc(:,i) = rms(hdemg(i).eccentric,2,"omitnan");
%     rmsEMG_con(:,i) = rms(hdemg(i).concentric,2,"omitnan");
% end
% 
% % channel 15 has nan values

%% 6c: write to .mat

filePath = input('Enter file path of folder for file to be saved: ',"s");
cd(filePath)
targetName = [participants{:} '_' timePoints{timePoint} '_' task '.mat'];
save(targetName,"differential","name","seconds",...
    "sampleFreq","hdemg","segmentFull","segmentEccentric","segmentConcentric",...
    "dxa","trunkDisplacement","meanTrunkDisplacement","endTrunkAngles");

disp('DONE. Onto the next script.')

%% old section
% redo seconds for delsys

figure; plot(emgCropArray,emgCrop(channels,:)); hold on; 
legend(num2str(channels))

% segmentTimes = [39507 52143; 68443 79870; 99152 114713]; % these need to be for 

% emgCrop1 = emgCrop(:,segmentTimes(1,1):segmentTimes(1,2));
% emgOut1 = mean(emgCrop(:,segmentTimes(1,1):segmentTimes(1,2)),2)';


% %% 7: for now, quick plots
% 
% electrodeGrid = ndgrid(1:15,1:4);
% electrodeGrid(:,2) = electrodeGrid(:,2) + 15;
% electrodeGrid(:,3) = electrodeGrid(:,3) + 15*2;
% electrodeGrid(:,4) = electrodeGrid(:,4) + 15*3;
% 
% % electrodeGrid(:,2) = flip(electrodeGrid(:,2));
% % electrodeGrid(:,4) = flip(electrodeGrid(:,4));
% 
% meanRMS = mean(rmsEMG,2);
% 
% mappedRMS = electrodeGrid;
% for i = 1:60
%     mappedRMS(mappedRMS == i) = meanRMS(i);
% end
% 
% figure; imagesc(mappedRMS); title([task ' ' timePoints{timePoint}]); colorbar;
