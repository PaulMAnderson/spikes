

function myPSTHViewer(spikeTimes, clu, eventTimes, window, trGroups, sortIdx)
% function psthViewer(spikeTimes, clu, eventTimes, window, trGroups)
%
% Controls:
% - c: dialog box to pick a new cluster ID number
% - t: toggle showing psth traces for each grouping variable or just the
% overall. If showing just overall, raster is sorted chronologically. If
% showing by grouping variable, raster is sorted by that variable.
% - r: select a new range within which to count spikes for the tuning curve
% -
%
% TODO:
% - if plotting all traces, color raster and sort
% - don't replot traces just change X,Y data for speed and for keeping zoom
% levels
% - indicate on tuning curve which color is which
% - add ability to switch from graded to distinguishable color scheme, also
% cyclical
% - add support for a second grouping variable - in that case, should use a
% different graded color scheme for each, overlay the traces in the tuning
% curve view, and also provide a 2-d image of tuning curve
% - add support for plot labels (traces in psth, x-axis in tuning curve)

fprintf(1, 'Controls:\n')
fprintf(1, '- left/right arrow: select previous/next cluster\n')
fprintf(1, '- up/down arrow: change smoothing of psth curves\n')
fprintf(1, '- c: dialog box to pick a new cluster ID number\n')
fprintf(1, '- t: toggle showing psth traces for each grouping variable or just the\n')
fprintf(1, 'overall. If showing just overall, raster is sorted chronologically. If\n')
fprintf(1, 'showing by grouping variable, raster is sorted by that variable.\n')
fprintf(1, '- r: select a new range within which to count spikes for the tuning curve\n')

% Set parameters, need to update so these can be passed into the function
params.histBinSize = 25; % in msec, width of histogram bins
params.smoothSize = 15; % in msec, stdev of gaussian smoothing filter
params.clusterIndex = 1;
params.rasterScale = floor(numel(eventTimes)/100); % height of ticks
params.window = window;
params.showAllTraces = false;
params.showErrorShading = false;
params.startRange = window(1);
params.stopRange = window(2);
params.rasterBinSize = 0.001;

figData.spikeTimes = spikeTimes;
figData.clu = clu;

if size(eventTimes,2) > size(eventTimes,1)
    eventTimes = eventTimes';
end

if nargin < 6
    if ~issorted(eventTimes)
        [~, ii] = sort(eventTimes(:,1));
        eventTimes = eventTimes(ii,:);
        trGroups = trGroups(ii);
    end
else
    sortIdx = sortIdx(:);
    eventTimes = eventTimes(sortIdx,:);
    trGroups = trGroups(sortIdx);
end

figData.eventTimes = eventTimes;
figData.trGroups = trGroups(:);
figData.clusterIDs = unique(clu);
figData.trGroupLabels = unique(figData.trGroups);
figData.nGroups = length(figData.trGroupLabels);
figData.plotAxes = [];

figData.params = params;

f = figure;
set(f, 'UserData', figData);
set(f, 'KeyPressFcn', @(f,k)psthViewerCallback(f, k));

psthViewerPlot(f)
end

function psthViewerPlot(f)
% fprintf(1,'plot with fig %d\n', get(f,'Number'));
figData = get(f,'UserData');

nChildren = length(f.Children);
while nChildren > 0 
    delete(f.Children(1));
    nChildren = length(f.Children);
end

% pick the right spikes
spikeTimes = figData.spikeTimes(figData.clu==figData.clusterIDs(figData.params.clusterIndex));

% compute everything
[rasters, spikeCounts, ba, bins] = ...
    rastersAndBins(spikeTimes, figData.eventTimes, figData.params.window, figData.params.rasterBinSize);

trGroupLabels = figData.trGroupLabels;
nGroups = figData.nGroups;
inclRange = bins>figData.params.startRange & bins<=figData.params.stopRange;
spikeCounts = nansum(ba(:,inclRange),2)./(figData.params.stopRange-figData.params.startRange);

% construct data for rasters
if figData.params.showAllTraces  
    rasterColors = figData.trGroups;
else
    rasterColors = ones(size(figData.trGroups));
end

% Calculate bin size
totalTimeInms = diff(figData.params.window) * 1000;
nBinsfor1ms   = length(bins) ./ totalTimeInms;
widthHistBins = nBinsfor1ms * figData.params.histBinSize;
nBins         = floor(length(bins) / widthHistBins);


% Make plots - using gramm
% Use Gramms stat functions to bin the acutal rasters
grammObj(1,1) = gramm('x',rasters,'color',rasterColors);
grammObj(1,1).stat_bin('nbins',nBins,'geom','line','normalization','countdensity');
% grammObj(1,1).stat_smooth('npoints',100);
% grammObj(1,1).stat_smooth('method','eilers','lambda','auto','npoints',length(bins)/100,'geom','line');
% grammObj(1,1).stat_smooth;

grammObj(2,1) = gramm('x',rasters,'color',rasterColors);
grammObj(2,1).geom_raster(); %'geom','line'

if ~figData.params.showAllTraces
    grammObj(1,1).set_color_options('map',[0.3 0.3 0.3],'n_color',1,'n_lightness',1);
    grammObj(2,1).set_color_options('map',[0.3 0.3 0.3],'n_color',1,'n_lightness',1);
end

grammObj(1,1).set_names('x','Time(s)','y','Firing Rate (Hz)','color','Trial Type');
grammObj(2,1).set_names('x','Time(s)','y','Trial #','color','Trial Type');
grammObj.set_title(['Cluster #' num2str(figData.params.clusterIndex)]);
grammObj.draw();

figData.grammObj = grammObj;

set(f,'UserData',figData)

end

function psthViewerCallback(f, keydata)

% fprintf('callback on %d with source %d\n', f.Number, keydata.Source.Number);


updateOtherFigs = false;

myData = get(f, 'UserData');
% myData.params

switch keydata.Key
    case 'rightarrow' % increment cluster index
        
        myData.params.clusterIndex = myData.params.clusterIndex+1;
        if myData.params.clusterIndex>length(myData.clusterIDs)
            myData.params.clusterIndex=1;
        end
        updateOtherFigs = true;
        
    case 'leftarrow' % decrement cluster index
        
        myData.params.clusterIndex = myData.params.clusterIndex-1;
        if myData.params.clusterIndex<1
            myData.params.clusterIndex=length(myData.clusterIDs);
        end
        updateOtherFigs = true;
        
    case 'uparrow' % increase smoothing
        myData.params.histBinSize = myData.params.histBinSize*1.2;
        
    case 'downarrow' % decrease smoothing
        myData.params.histBinSize = myData.params.histBinSize/1.2;
        
    case 'e' % whether to show standard error as shading
        myData.params.showErrorShading = ~myData.params.showErrorShading;
        
    case 't' % whether to plot the psth trace for each condition or just the overall one
        myData.params.showAllTraces = ~myData.params.showAllTraces;
        updateOtherFigs = true;
        
    case 'r'
        ax = subplot(3,1,1); title('click start and stop of range')
        %         [startRange,~] = ginput();
        %         [stopRange,~] = ginput();
        waitforbuttonpress;
        q = get(ax, 'CurrentPoint');
        myData.params.startRange = q(1,1);
        waitforbuttonpress;
        q = get(ax, 'CurrentPoint');
        myData.params.stopRange = q(1,1);
        if myData.params.stopRange<myData.params.startRange
            tmp = myData.params.startRange;
            myData.params.startRange = myData.params.stopRange;
            myData.params.stopRange = tmp;
        end
        
    case 'c'
        newC = inputdlg('cluster ID?');
        ind = find(myData.clusterIDs==str2num(newC{1}),1);
        if ~isempty(ind)
            myData.params.clusterIndex = ind;
        end
        
        updateOtherFigs = true;
        
end

set(f, 'UserData', myData);

% plot with new settings
psthViewerPlot(f)

if updateOtherFigs && f==keydata.Source.Number && isfield(myData, 'otherFigs')
    % checking that the current figure matches the source number prevents
    % doing this when called *not* as the original fig
    setOtherFigsClusterIndex(myData, myData.params.clusterIndex)
    plotOtherFigs(f)
end

end


function setOtherFigsClusterIndex(myData, cInd)

for thatFig = myData.otherFigs
    thatData = get(thatFig, 'UserData');
    thatData.params.clusterIndex = cInd;
    set(thatFig, 'UserData', thatData);
    
end

end


function plotOtherFigs(f)
myData = get(f, 'UserData');
for thatFig = myData.otherFigs
%     thatFigFcn = get(thatFig, 'KeyPressFcn');
    figs = get(0, 'Children');
    figNumsCell = get(figs, 'Number');
    figNums = [figNumsCell{:}];
    thatFigObj = figs(figNums==thatFig);
    psthViewerPlot(thatFigObj);
end
figure(f) % return focus here
end



function [rasters, spikeCounts, binnedArray, bins] = ...
    rastersAndBins(spikeTimes, eventTimes, window, psthBinSize)

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

    % first we'll subselect spikes that are between the end and beginning of
    % all the ranges - in some cases this will be useless but often instead
    % reduces spike count by a lot.
    spikeTimes = spikeTimes(spikeTimes>min(eventTimes(:,1)+window(1)) & spikeTimes<max(eventTimes(:,1)+window(2)));

    [binnedArray, bins, rasters] = timestampsToBinned(spikeTimes, eventTimes, psthBinSize, window);
    spikeCounts = nansum(binnedArray,2);
end % End of function rastersAndBins


function [binArray, binCenters, rasters] = timestampsToBinned(timeStamps, referencePoints, binSize, window)
    % [binArray, binCenters] = timestampsToBinned(timeStamps, referencePoints, binSize,
    % window)
    %
    % Returns an array of binned spike counts. If you use a large enough
    % binSize, it may well be possible that there is more than one spike in a
    % bin, i.e. the value of some bin may be >1. 

    binBorders = window(1):binSize:window(2);
    numBins = length(binBorders)-1;

    binArray = zeros(length(referencePoints), numBins);
    rasters  = cell(1,length(referencePoints));

    if size(referencePoints,2) == 2
        % in this case there is a reference point and a second point to
        % exclude values from
        for r = 1:length(referencePoints)
            [n,binCenters] = histdiff(timeStamps, referencePoints(r,1), binBorders);
            % zero out bins that fall before/after the second reference
            refWindow = referencePoints(r,2) - referencePoints(r,1);
            binMatch  = dsearchn(binBorders',refWindow);
            if refWindow < 0
                n(1:binMatch-1) = 0;
                rasters{r} = binCenters(logical(n));
               %  n(1:binMatch-1) = nan;
            else
                n(binMatch+1:end) = 0;
                rasters{r} = binCenters(logical(n));
               %  n(binMatch+1:end) = nan;
            end
            binArray(r,:) = n;
        end
        
    else
        % The standard single reference point
        for r = 1:length(referencePoints)
             [n,binCenters] = histdiff(timeStamps, referencePoints(r), binBorders);
             rasters{r} = binCenters(logical(n));
             binArray(r,:) = n;  
        end
    end
    
end % End function timestamps to be binned
