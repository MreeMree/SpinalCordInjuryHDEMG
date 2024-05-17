%% interpolate bad HDEMG channels with nearby channels
% uses linear method but weights adjacents channels according to their
% location on the electrode grid

% Joshua Kearney May 2024

% electrode arrays this was written for are OTBio Elettronica ones, where
% electrodes are 8 mm apart, meaning ~11 mm apart diagonally.

% interpolates one bad channel at a time

% INPUTS:
% emg = multi-channel array of emg data, e.g. 64(channels)x258909(datapoints)
% badChannel = a single number indicating which channel you wish to replace
% adjacentChannels = a 1D array, series of numbers indicating channels adjacent to one being replaced
% electrodeGrid = channels in same arrangement as electrode array.

% electrodeGrid can easily be created with reshape function, e.g.:
% electrodeGrid = reshape(0:64,13,5)

% OUTPUT:
% intChan = 1D array same length as emg data, repping new interpolated channel

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
