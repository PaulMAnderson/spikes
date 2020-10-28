

function [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spikeTimes, eventTimes, window, psthBinSize)
% function [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spikeTimes, eventTimes, window, psthBinSize)
%
% Fast computation of psth and spike counts in a window relative to some
% events. Also returns rasters you can plot. 
%
% Notes on inputs:
% - eventTimes is nEvents x 1
% - window is length 2, e.g. [-0.1 0.3] for a window from -0.1 to +0.3 sec
% relative to events. 
%
% Notes on outputs:
% - psth can be plotted with plot(bins, psth);
% - rasters can be plotted with plot(rasterX, rasterY);
% - spikeCounts is nEvents x 1, where each entry is the number of spikes
% that occurred within the window around that event. 

spikeTimes = spikeTimes(:);
eventTimes = sort(eventTimes(:));

% first we'll subselect spikes that are between the end and beginning of
% all the ranges - in some cases this will be useless but often instead
% reduces spike count by a lot.
spikeTimes = spikeTimes(spikeTimes>min(eventTimes+window(1)) & spikeTimes<max(eventTimes+window(2)));

[binnedArray, bins] = timestampsToBinned(spikeTimes, eventTimes, psthBinSize, window);

spikeCounts = sum(binnedArray,2);
psth = mean(binnedArray./psthBinSize); % normalize to Hz

[tr,b] = find(binnedArray);
[rasterX,yy] = rasterize(bins(b));
rasterY = yy+reshape(repmat(tr',3,1),1,length(tr)*3); % yy is of the form [0 1 NaN 0 1 NaN...] so just need to add trial number to everything


function [binArray, binCenters] = timestampsToBinned(timeStamps, referencePoints, binSize, window)
% [binArray, binCenters] = timestampsToBinned(timeStamps, referencePoints, binSize,
% window)
%
% Returns an array of binned spike counts. If you use a large enough
% binSize, it may well be possible that there is more than one spike in a
% bin, i.e. the value of some bin may be >1. 

binBorders = window(1):binSize:window(2);
numBins = length(binBorders)-1;

if isempty(referencePoints)
    binArray = [];
    binCenters = binBorders(1:end-1)+binSize/2; % not sure if this is the same as what you get below?
    return;
end

binArray = zeros(length(referencePoints), numBins);

if isempty(timeStamps)
    binCenters = binBorders(1:end-1)+binSize/2; % not sure if this is the same as what you get below?
    return;
end

for r = 1:length(referencePoints)
    [n,binCenters] = histdiff(timeStamps, referencePoints(r), binBorders);
    binArray(r,:) = n;
end
