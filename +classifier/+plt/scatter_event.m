function h = scatter_event(data,ax,y)
% Input: y - location on y axis 

	% initiation 
	if nargin < 3
		y = 0;
	end

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

	% actual reward
	h(1,1) = my_scatter(data.rewards.all.front,y,ax,'.','MarkerEdgeColor',Colors(1,:));
	h(2,1) = my_scatter(data.rewards.all.rear,y,ax,'.','MarkerEdgeColor',Colors(2,:));

	% cues
	h(1,2) = my_scatter(data.cues.all.front,y,ax,'*','MarkerEdgeColor',Colors(1,:));
	h(2,2) = my_scatter(data.cues.all.rear,y,ax,'*','MarkerEdgeColor',Colors(2,:));

	% port entries
	h(1,3) = my_scatter(data.entries.unrewarded.front,y,ax,'^','MarkerEdgeColor',Colors(1,:));
	h(2,3) = my_scatter(data.entries.unrewarded.rear,y,ax,'^','MarkerEdgeColor',Colors(2,:));

	% laser
	if isfield(data,'laser')
		h(1,4) = my_scatter(data.laser.onsets,y,ax,'.','MarkerEdgeColor',[0 0 1]);
	end

end