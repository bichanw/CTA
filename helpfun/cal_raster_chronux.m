function Raster = cal_raster_chronux(SpikeTime,EventTime,TOI)
% input spike time and event time stamps, get raster plot in TOI

	% rotate spiketime and eventtime
	if size(SpikeTime,1) == 1 SpikeTime = SpikeTime'; end
	if size(EventTime,1) == 1 EventTime = EventTime'; end

	% locate all events
	NTrials = length(EventTime);
	Raster  = struct;

	% calculate time from spikes to events
	seDiff = int64(SpikeTime) - int64(EventTime)';

	% spikes within TOI
	Spikes_inTOI = (seDiff <= TOI(2)) & (seDiff > TOI(1));

	for iTrial = 1:NTrials
	    
	    % skip trial without spike in TOI
	    if sum(Spikes_inTOI(:,iTrial))==0  
	    	continue;  
	    end

		% keep all spike times
		tmp = seDiff(:,iTrial);
		Raster(iTrial).times = double(tmp(Spikes_inTOI(:,iTrial)))/1000; % ms to s, for chronux
	    
	end

	% add up remaining trials without spikes
	for iTrial = (numel(Raster)+1):NTrials
		Raster(iTrial).times = [];
	end

end
