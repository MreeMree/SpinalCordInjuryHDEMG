%% RPM

%% select directory

clear all;

switch 4
    case 1 % if working from PC
        filePath = '\\its-rds.bham.ac.uk\rdsprojects\c\chious-arm-cycling-study\ISRT Subacute SCI\';
        cd([filePath 'Processed']);
    case 2 % if working from Mac
        %filePath = 'Volumes/chious-arm-cycling-study/ISRT Subacute SCI/';
        filePath = '/Users/joshuakearney/Documents/PostdocSCI/Data/0703Analysis';
        %cd([filePath 'Processed']);
        cd(filePath);
    case 3 % if working from BEAR portal (Linux)
        filePath = '/rds/projects/c/chious-arm-cycling-study/ISRT Subacute SCI';
        cd([filePath 'Processed']);
    case 4 % custom file path
        filePath = input("Enter file path of processed files: ","s");
        cd(filePath)
end

%% select task

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
        task = 'CMEP';
end

%% load data for selected participants

timePoints = {'01','02','03','04'};
%participants = input("Enter participant numbers, e.g. {'P01','P04'}: "); %participantNo = [1];
participants = {input("Enter participant number, e.g. JP01: ","s")}; %participantNo = [1];
clear baseline; clear pre; clear post; clear washout;

for i = 1:length(participants) % participants
        
        timePoint = 1;
        filename = [participants{i} '_' timePoints{timePoint} '_' task '.mat'];
        if exist(filename) == 2 
        baseline(i) = load(filename);
        baseline(i).participant = participants(i);
        else
        end

        timePoint = 2;
        filename = [participants{i} '_' timePoints{timePoint} '_' task '.mat'];
        if exist(filename) == 2 
        pre(i) = load(filename);
        pre(i).participant = participants(i);
        else
        end

        timePoint = 3;
        filename = [participants{i} '_' timePoints{timePoint} '_' task '.mat'];
        if exist(filename) == 2 
        post(i) = load(filename);
        post(i).participant = participants(i);
        else
        end

        timePoint = 4;
        filename = [participants{i} '_' timePoints{timePoint} '_' task '.mat'];
        if exist(filename) == 2 
        washout(i) = load(filename);
        washout(i).participant = participants(i);
        else
        end

end

%% calculate raw RMS: select time point for same participant

switch input("Enter time point (1-4): ")
    case 1 
        struct = baseline;
        tp = "Baseline"; timePoint = 1;
        disp('Looking at baseline data')
    case 2 
        struct = pre;
        tp = "Pre-test"; timePoint = 2;
        disp('Looking at pre-test data')
    case 3 
        struct = post;
        tp = "Post-test"; timePoint = 3;
        disp('Looking at post-test data')
    case 4 
        struct = followUp;
        tp = "Follow up"; timePoint = 4;
        disp('Looking at follow-up data')
end

rmsEMG = nan(height(struct.differential),3); %,length(participantNo));
rmsEMG_ecc = nan(height(struct.differential),3);
rmsEMG_con = nan(height(struct.differential),3);

switch 2
    case 1
        methodRMS = "block";
    case 2
        methodRMS = "rolling";
end

for ind = 1:length(participants)

    if methodRMS == "block"
        for i = 1:3
            %for j = 1:height(emgCrop) % number of channels left, 60 after differential calc 
            rmsEMG(:,i) = rms(struct.hdemg(i).full,2,"omitnan");
            rmsEMG_ecc(:,i) = rms(struct.hdemg(i).eccentric,2,"omitnan");
            rmsEMG_con(:,i) = rms(struct.hdemg(i).concentric,2,"omitnan");
        end
    elseif methodRMS == "rolling"
        for i = 1:3
            tempArray1 = nan(height(struct.differential),floor(length(struct.hdemg(i).full)/500));
            tempArray2 = nan(height(struct.differential),floor(length(struct.hdemg(i).eccentric)/500));
            tempArray3 = nan(height(struct.differential),floor(length(struct.hdemg(i).concentric)/500));
            for j = 1:size(tempArray1,2)
                tempArray1(:,j) = rms(struct.hdemg(i).full(:,(1:500)+500*(j-1)),2,"omitnan");
            end
            for j = 1:size(tempArray2,2)
                tempArray2(:,j) = rms(struct.hdemg(i).eccentric(:,(1:500)+500*(j-1)),2,"omitnan");
            end
            for j = 1:size(tempArray3,2)
                tempArray3(:,j) = rms(struct.hdemg(i).concentric(:,(1:500)+500*(j-1)),2,"omitnan");
            end
            rmsEMG(:,i) = mean(tempArray1,2);
            rmsEMG_ecc(:,i) = mean(tempArray2,2);
            rmsEMG_con(:,i) = mean(tempArray3,2);
        end
    end

end

% % monopolar RMS
% for ind = 1:length(participants)
% 
%     for i = 1:3
%         %for j = 1:height(emgCrop) % number of channels left, 60 after differential calc 
%         rmsEMG_m(:,i) = rms(baseline.monopolar,2,"omitnan");
%         rmsEMG_ecc_m(:,i) = rms(baseline.hdemg(i).eccentric,2,"omitnan");
%         rmsEMG_con_m(:,i) = rms(baseline.hdemg(i).concentric,2,"omitnan");
%     end
% 
% end

%% topoplot

% electrodeGrid = ndgrid(1:15,1:4);
% electrodeGrid(:,2) = electrodeGrid(:,2) + 15;
% electrodeGrid(:,3) = electrodeGrid(:,3) + 15*2;
% electrodeGrid(:,4) = electrodeGrid(:,4) + 15*3;

switch 1
    case 1
        %meanRMS = mean(rmsEMG,2);
        meanRMS = mean(rmsEMG(:,:),2); % whole movement
        titleMovement = 'Reach and return'; disp(titleMovement)
    case 2
        meanRMS = mean(rmsEMG_con,2); % return (concentric)
        titleMovement = 'Return'; disp(titleMovement)
    case 3
        meanRMS = mean(rmsEMG_ecc,2); % reach (eccentric)
        titleMovement = 'Reach only'; disp(titleMovement)
    case 4 % just playing looking at whole movement without trials
        meanRMS = mean(rms(struct.monoCrop,2,"omitnan"),2);
    case 5
        meanRMS = mean(rms(struct.differential,2,"omitnan"),2);
end

meanRMS(60) = NaN;
m = round(size(meanRMS,1)/5);
electrodeGrid = reshape(1:size(meanRMS,1),[m,5]);
%electrodeGrid = flip(electrodeGrid);

% electrodeGrid(:,2) = flip(electrodeGrid(:,2));
% electrodeGrid(:,4) = flip(electrodeGrid(:,4));

mappedRMS = electrodeGrid-1;
for i = 1:size(meanRMS,1)
    mappedRMS(mappedRMS == i) = meanRMS(i);
end

mappedRMS(1) = NaN;
%mappedRMS = flip(mappedRMS)
figure; imagesc(flip(mappedRMS)); title([task ' ' participants{1} ': ' titleMovement]); colorbar;

clear m;

%% OPTIONAL: check PRD of each channel vs mean of the rest of the channels
% differential

% y = struct.hdemg(1).full(:,:);
% PRD = NaN(1,size(y,1));
% 
% % for i = 1:length(y)
% %     first = y(i);
% %     second = mean(y(setdiff(1:length(y),i)));
% % 
% %     diffSum = ((second - first).^2)/(second.^2);
% %     PRD(i) = sqrt(diffSum)*100; % percent residual difference
% % end
% 
% for i = 1:size(y,1)
%     first = y(i,:);
%     second = mean(y(setdiff(1:height(y),i),:));
% 
%     diffSum = ((second - first).^2)/(second.^2);
%     PRD(i) = sqrt(diffSum)*100; % percent residual difference
% end
% 
% y = PRD;
% figure; bar(y); 
% hold on; plot([1 60],[mean(y(y~=0))+2*std(y(y~=0)) mean(y(y~=0))+2*std(y(y~=0))])
% title('PRD of each channel and mean of all other channels')
% mean2sd = mean(y)+2*std(y);
% 
% disp('Channel differentials with PRD > mean+2*sd = '); yInd = 1:length(y);
% disp(yInd(y>mean2sd))
% disp(['Mean PRD: ' num2str(mean(y))])
% disp(['Median PRD: ' num2str(median(y))])
% 
% clear y

%% ROUND 1: get rid of any outlier channels

y = mappedRMS;
y = reshape(y,[1 60]);
% y = flip(mappedRMS);
y = y(2:60);
[~,~,~,stats] = ttest(y);

figure; bar(y); xticks(1:60); hold on;
plot([1 60],[mean(y)+2*std(y) mean(y)+2*std(y)])
%plot([1 60],[stats.tstat*2 stats.tstat*2])
plot([1 60],[mean(y)+stats.tstat*2 mean(y)+stats.tstat*2])

exclude = input("Enter channels to exclude: ");
y(exclude) = 0;

y = flip(y);
y(60) = NaN;
y = flip(y);
mappedRMS_inc = reshape(y,[12 5]);

disp('Channel differentials with PRD > mean+2*sd = '); yInd = 1:length(y);
mean2sd = mean(y(~isnan(y)))+2*std(y(~isnan(y))); disp(yInd(y>mean2sd))

figure; imagesc(flip(mappedRMS_inc)); title([task ' ' participants{1} ': ' titleMovement] ); colorbar;

%% ROUND 2: any further channels to remove? (optional, discouraged)

y = mappedRMS_inc;
y = reshape(y,[1 60]);
% y = flip(mappedRMS);
%y = y(2:60);
[~,~,~,stats] = ttest(y);

figure; bar(0:59,y); xticks(1:60); hold on;
plot([0 59],[mean(y(y>0))+2*std(y(y>0)) mean(y(y>0))+2*std(y(y>0))])
%plot([1 60],[stats.tstat*2 stats.tstat*2])
plot([0 59],[mean(y(y>0))+stats.tstat*2 mean(y(y>0))+stats.tstat*2])

exclude = input("Enter channels to exclude: ");
y(exclude+1) = 0;

mappedRMS_inc = reshape(y,[12 5]);

figure; imagesc(flip(mappedRMS_inc)); title([task ' ' participants{1} ': ' titleMovement] ); colorbar;

%% plot monopolar trial by trial



%% Active channels

%mappedRMS(1,4) = 0;

% keep array flipped now as will not need to index channel by number:
mappedRMS_inc_RWU = flip(mappedRMS_inc); % RWU = right way up

% Calculate active channels
[yActive,xActive]= find( mappedRMS_inc_RWU > 0.7*max(max(mappedRMS_inc_RWU)) ); % coordinates of active channels
IQR_yActive = iqr(unique(yActive,'first'));%numel(YactchanR); % number of active channels

% Calculate the centroid coordinates
x = reshape(1:60,[12 5]); % creating a support matrix to find centroid coordinates
            
for i = 1:numel(yActive)
    x(yActive(i),xActive(i)) = 0; % the active channels in the support matrix 'x' are equal to 0
end

% Calculating the global RMS amplitue for active channels
RMSActive = mean( mappedRMS_inc_RWU(x == 0) );

%% calculate x,y barycentre of active channels

mappedActiveRMS = mappedRMS_inc_RWU;
mappedActiveRMS(find(x)) = 0; % The function 'find' gives the non-zero elements in the matrix.

% How to compute the weighted centroid:
% (i) multiplying the number of each row/col of active channels by
% their corresponding RMS amplitude; (ii) summing up the products; 
% (iii) dividing the resulting value by the sum of RMS amplitude 
% values sum of RMS amplitude values        
% Normalizing with the sum of RMS amplitude 
mappedActiveRMS = mappedActiveRMS/sum(mappedActiveRMS(:));
% Finding centroid coordinates
[m,n] = size(mappedActiveRMS);
[I,J] = ndgrid(1:m,1:n);
centroidCoord(1:2) =[dot(J(:),mappedActiveRMS(:)), dot(I(:),mappedActiveRMS(:))];

%% heat map

figure;
imagesc(mappedRMS_inc_RWU);
hold on;
plot(xActive,yActive,'o','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',6); % plotting the active channels        
plot(centroidCoord(1),centroidCoord(2),'o','MarkerEdgeColor','k','MarkerFaceColor','w',...
        'MarkerSize',14,'Linewidth',2);
colormap('jet'); colorbar;  % changing the colormap 'jet' 'gray'  
title([ struct.participant{1} ' ' tp ' ' task ' - ' titleMovement])

%% write to .csv

channels = (1:60)';
RMS = reshape(mappedRMS_inc,[60,1]);
activeChanRMS = reshape(mappedActiveRMS,[60,1]);

tab = table(channels,RMS,activeChanRMS);

write(tab,[struct.participant{1} '_' timePoints{timePoint} '_' task '_ActiveChannels.csv'])

%% save key variables for quicker looking next time

save([struct.participant{1} '_' timePoints{timePoint} '_' task '_RMS.mat'],...
    'RMSActive','IQR_yActive',"mappedActiveRMS","mappedRMS_inc","centroidCoord")

disp('DONE. Onto the next dataset.')

% END