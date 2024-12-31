function [resp,resp_err,RR,raster] = cal_psth(spiketime,eventtime,varargin)
% calculate psth
% Input: spiketime, eventtime in ms
% 	     tp, kernel_width in sec
	p = inputParser;
	addParameter(p,'tp',-0.2:0.01:0.7);
	addParameter(p,'kernel_width',0.02);
	parse(p,varargin{:});
	tp = p.Results.tp;


	% NTrials * 1 cell, each with the spike time relative to event onset
	% take extra 100 ms before and after so the psth won't look cut off
	raster = cal_raster_chronux(spiketime,eventtime,[tp(1) tp(end)]*1000 + [-100 100]);

	% calculate Gaussian response
	if isempty(fieldnames(raster))
	    % empty location (no spikes)
	    resp     = zeros(size(tp));
	    resp_err = resp;
	    RR = zeros(numel(tp),numel(eventtime));
	else
	    % calculate response if there're spikes
	    [resp,~,resp_err,RR] = psth(raster,p.Results.kernel_width,'n',[tp(1) tp(end)],1,tp);
	end
end