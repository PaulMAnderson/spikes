% Inputs/outputs: mostly self explanatory 
% localizedSpikesOnly (false by default) - if true, only spikes with no discrepancy between depth and site are returned. 
function  [spikes, templates] = spikeTemplateDepths(spikeStruct)
% [spikeTimes, spikeAmps, spikeDepths, spikeSites]

%% Need to check that the principal components exist
pcFeat = squeeze(spikeStruct.pcFeat(:,1,:)); % take first PC only
pcFeat(pcFeat<0) = 0; % some entries are negative, but we don't really want to push the CoM away from there.

%% compute center of mass of these features

% which channels for each spike?
spikeFeatInd = spikeStruct.pcFeatInd(spikeStruct.spikeTemplates+1,:);

% Ycoords of those channels?
spikeFeatYcoords = spikeStruct.ycoords(spikeFeatInd+1); % 2D matrix of size #spikes x 32?
% Xcoords of those channels?
spikeFeatXcoords = spikeStruct.xcoords(spikeFeatInd+1); % 2D matrix of size #spikes x 32

% center of mass is sum(coords.*features)/sum(features)
spikeDepths = sum(spikeFeatYcoords.*pcFeat.^2,2)./sum(pcFeat.^2,2);
spikeWidths = sum(spikeFeatXcoords.*pcFeat.^2,2)./sum(pcFeat.^2,2);

% unwhiten all the templates
tempsUnW = zeros(size(spikeStruct.temps));
for t = 1:size(spikeStruct.temps,1)
    tempsUnW(t,:,:) = squeeze(spikeStruct.temps(t,:,:))*spikeStruct.winv;
end

% The amplitude on each channel is the positive peak minus the negative
tempChanAmps = squeeze(max(tempsUnW,[],2))-squeeze(min(tempsUnW,[],2));

% The template amplitude is the amplitude of its largest channel (but see
% below for true tempAmps)
[tempAmpsUnscaled, maxChan] = max(tempChanAmps,[],2);


% need to zero-out the potentially-many low values on distant channels ...
threshVals = tempAmpsUnscaled*0.3; 
tempChanAmps(bsxfun(@lt, tempChanAmps, threshVals)) = 0;

% ... in order to compute the depth as a center of mass
templateDepths = sum(bsxfun(@times,tempChanAmps,spikeStruct.ycoords'),2)./sum(tempChanAmps,2);
% Also calculate X-Pos
templateWidths = sum(bsxfun(@times,tempChanAmps,spikeStruct.xcoords'),2)./sum(tempChanAmps,2);

% assign all spikes the amplitude of their template multiplied by their
% scaling amplitudes (templates are zero-indexed)
spikeAmps = tempAmpsUnscaled(spikeStruct.spikeTemplates+1).*spikeStruct.tempScalingAmps;

% take the average of all spike amps to get actual template amps (since
% tempScalingAmps are equal mean for all templates)
ta = clusterAverage(spikeStruct.spikeTemplates+1, spikeAmps);
tids = unique(spikeStruct.spikeTemplates);
tempAmps(tids+1) = ta; % because ta only has entries for templates that had at least one spike
tempAmps = tempAmps'; % for consistency, make first dimension template number

% Get channel with largest amplitude, take that as the waveform
[~,tempChannel] = max(max(abs(spikeStruct.temps),[],2),[],3);
templates_max = nan(size(spikeStruct.temps,1),size(spikeStruct.temps,2));
spikeChannel = tempChannel(spikeStruct.spikeTemplates+1);
for curr_template = 1:size(spikeStruct.temps,1)
    templates_max(curr_template,:) = ...
        spikeStruct.temps(curr_template,:,tempChannel(curr_template));
end
waveforms = templates_max;



% Get trough-to-peak time for each template
[~,waveform_trough] = min(templates_max,[],2);
[~,templateDuration] = arrayfun(@(x) ...
    max(templates_max(x,waveform_trough(x):end),[],2), ...
    transpose(1:size(templates_max,1)));

%% Place data in structures

spikes.Amplitudes = spikeAmps;
spikes.Depths     = spikeDepths;
spikes.Widths     = spikeWidths;
spikes.Cluster    = spikeStruct.clu;
spikes.Times      = spikeStruct.st;

templates.Depths     = templateDepths;
templates.Widths     = templateWidths;
templates.Amplitudes = tempAmps;
templates.Channel    = tempChannel;
templates.Duration   = templateDuration;
templates.Waveforms  = waveforms;





