%% interpolate bad HDEMG channels with nearby channels
% uses linear method but weights adjacents channels according to their
% location on the electrode grid

% Joshua Kearney May 2024

% electrode arrays this was written for are OTBio Elettronica ones, where
% electrodes are 8 mm apart, meaning ~11 mm apart diagonally.

% interpolates one bad channel at a time

function intChan = interpolateHDEMG(emg,badChannel,adjacentChannels,electrodeGrid)

intChan = NaN(1,length(emg));

% 1: weight each adjacent channel before finding the mean
distances = NaN(size(adjacentChannels));
weights = NaN(size(adjacentChannels));

[Y,X] = meshgrid(1:5,1:13);
rowBad = X(electrodeGrid == badChannel);
colBad = Y(electrodeGrid == badChannel);

for i = 1:length(adjacentChannels)
    
    row = X(electrodeGrid == adjacentChannels(i));
    col = Y(electrodeGrid == adjacentChannels(i));
    a = abs(rowBad-row);
    b = abs(colBad-col);

    distances(i) = sqrt(a^2+b^2)*8;

end

weights = distances/sum(distances); % 8 mm shortest distance between electrodes

% 2: cycle through each data point and interpolate (find mean)

for i = 1:length(emg)
    intChan(i) =  sum(emg(adjacentChannels,i).*weights');
end

end