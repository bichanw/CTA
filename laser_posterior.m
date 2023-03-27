function laser_posterior(data,ops,prefix)

% initiation
if nargin < 3
	prefix = '';
end

% count spikes with larger temporal precision
ops.tp = getOr(ops,'tp',[0 1]);
ops.posterior_t_edges = data.video(1):0.1:data.video(end);
[spk_count,ops] = classifier.count_spk.time_course(data,ops);


% calculate posterior
if ~isfield(ops,'Mdl')
	[~,ops] = ops.classifier.train(data,ops);
end
Posterior = ops.classifier.predict(ops.Mdl,spk_count',ops);


% time locked to laser onset
if isstruct(data.laser)
	% in develop
	error('need to add onsets / offsets conversion');
else
	% convert time to relative time

	signal = Posterior;
	eventtime = data.laser(:,1);
	raw_t = ops.posterior_t;
	toi = [-1 3];

	% flip signal if needed
	if size(signal,1)==numel(raw_t) && size(signal,2)~=numel(raw_t)
		signal = signal';
	end

	% writine a for loop to go through right now
	% to prevent using too much memory
	t = [];
	resp = [];
	for ii = 1:numel(eventtime)
		relative_t = raw_t - eventtime(ii);

		% extract data
		ind  = (relative_t>=toi(1)) & (relative_t<=toi(2));
		t    = [t, relative_t(ind)];
		resp = [resp, signal(:,ind)];

	end
	[R,~,~,tbin] = running_average(t,resp,0.15,0.05);
	[V,~,~,tbin] = running_average(t,resp,0.15,0.05,[],@std); 
	plot(tbin,R(:,1:2));


end