function scatter_event(data,ax)

% change front / familiar color
if isfield(data,'port_is_water')

	% front port is water
	if data.port_is_water(1)
		Colors = [0 0 0; 1 0 0];
	% front port is flavor
	else
		Colors = [1 0 0; 0 0 0];
	end

% by default use front as novel
else
	Colors = [1 0 0; 0 0 0];
end

n_cells = numel(data.spikes);
% Y = repmat([0 n_cells+1 n_cells+1 0]',1,numel(data.rewards.all.front));  % [0.1 0.7]

fill(ax,(data.rewards.all.front+[0.1 0.1 0.7 0.7])',repmat([0 n_cells+1 n_cells+1 0]',1,numel(data.rewards.all.front)),Colors(1,:),'FaceAlpha',0.1,'EdgeColor','none');
fill(ax,(data.rewards.all.rear+[0.1 0.1 0.7 0.7])',repmat([0 n_cells+1 n_cells+1 0]',1,numel(data.rewards.all.rear)),Colors(2,:),'FaceAlpha',0.1,'EdgeColor','none');

% laser
if isfield(data,'laser')
	fill(ax,(data.laser.onsets+[0 0 0.5 0.5])',repmat([0 n_cells+1 n_cells+1 0]',1,numel(data.laser.onsets)),[0 0 1],'FaceAlpha',0.1,'EdgeColor','none');
end
