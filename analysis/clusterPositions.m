% function to extract position information on kilosort clusters 
function  cluster = clusterPositions(spikeStruct)

% unwhiten all the templates
tempsUnW = zeros(size(spikeStruct.temps));
for t = 1:size(spikeStruct.temps,1)
    tempsUnW(t,:,:) = squeeze(spikeStruct.temps(t,:,:))*spikeStruct.winv;
end

% The amplitude on each channel is the positive peak minus the negative
tempChanAmps = squeeze(max(tempsUnW,[],2))-squeeze(min(tempsUnW,[],2));

% The template amplitude is the amplitude of its largest channel (but see
% below for true tempAmps)
[tempAmpsUnscaled, templateChannel] = max(tempChanAmps,[],2);


% need to zero-out the potentially-many low values on distant channels ...
threshVals = tempAmpsUnscaled*0.3; 
tempChanAmps(bsxfun(@lt, tempChanAmps, threshVals)) = 0;

% ... in order to compute the depth as a center of mass
templateDepths = sum(bsxfun(@times,tempChanAmps,spikeStruct.ycoords'),2)./sum(tempChanAmps,2);
% Also calculate X-Pos
templateSpan = sum(bsxfun(@times,tempChanAmps,spikeStruct.xcoords'),2)./sum(tempChanAmps,2);

% assign all spikes the amplitude of their template multiplied by their
% scaling amplitudes (templates are zero-indexed)
spikeAmps = tempAmpsUnscaled(spikeStruct.spikeTemplates+1).*spikeStruct.tempScalingAmps;

% take the average of all spike amps to get actual template amps (since
% tempScalingAmps are equal mean for all templates)
templateAmps = nan(size(spikeStruct.temps,1),1);
ta = clusterAverage(spikeStruct.spikeTemplates+1, spikeAmps);
tids = unique(spikeStruct.spikeTemplates);
templateAmps(tids+1) = ta; % because ta only has entries for templates that had at least one spike

% Get the waveform from the channel with the largest amplitude
% [~,tempChannel] = max(max(abs(spikeStruct.temps),[],2),[],3);
waveforms = nan(size(spikeStruct.temps,1),size(spikeStruct.temps,2));
for curr_template = 1:size(spikeStruct.temps,1)
    waveforms(curr_template,:) = ...
        spikeStruct.temps(curr_template,:,templateChannel(curr_template));
end

% Get trough-to-peak time for each template
[~,waveform_trough] = min(waveforms,[],2);
[~,templateDuration] = arrayfun(@(x) ...
    max(waveforms(x,waveform_trough(x):end),[],2), ...
    transpose(1:size(waveforms,1)));

%% Place data in structure
cluster.Depth     = templateDepths;
cluster.Span      = templateSpan;
cluster.Amplitude = templateAmps;
cluster.Channel   = templateChannel;
cluster.Duration  = templateDuration;
cluster.Waveform  = waveforms;
