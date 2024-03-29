function h = scatter_event(data,ax,y)
% Input: y - location on y axis 

	% initiation 
	if nargin < 3
		y = 0;
	end

	% change front / familiar color
	Colors = getOr(data,'port_color',[0.9412    0.2118    0.2745;   0.0667    0.5020    0.9451]);
	% if isfield(data,'port_is_water')

	% 	% front port is water
	% 	if data.port_is_water(1)
	% 		Colors = [0 0 0; 1 0 0];
	% 	% front port is flavor
	% 	else
	% 		Colors = [1 0 0; 0 0 0];
	% 	end

	% % by default use front as novel
	% else
	% 	Colors = [1 0 0; 0 0 0];
	% end

	% actual reward
	h(1,1) = my_scatter(data.rewards.all.front,y,ax,'.','MarkerEdgeColor',Colors(1,:));
	h(2,1) = my_scatter(data.rewards.all.rear,y,ax,'.','MarkerEdgeColor',Colors(2,:));

	% cues
	h(1,2) = my_scatter(data.cues.all.front,y,ax,'*','MarkerEdgeColor',Colors(1,:));
	h(2,2) = my_scatter(data.cues.all.rear,y,ax,'*','MarkerEdgeColor',Colors(2,:));

	% port entries
	if isfield(data.entries,'unrewarded')
		tmp = data.entries.unrewarded;
	else
		% trim down data for plotting
		tmp = data.entries.all;
		tmp.front = tmp.front(:,1);
		tmp.rear  = tmp.rear(:,1);
	end
	h(1,3) = my_scatter(tmp.front,y,ax,'^','MarkerEdgeColor',Colors(1,:));
	h(2,3) = my_scatter(tmp.rear,y,ax,'^','MarkerEdgeColor',Colors(2,:));

	% laser
	if isfield(data,'laser')
		if isstruct(data.laser)
			h(1,4) = my_scatter(data.laser.onsets,y,ax,'.','MarkerEdgeColor',[0 0 1]);
		else
			% h(1,4) = my_scatter(data.laser(:,1),y,ax,'.','MarkerEdgeColor',[0 0 1]);
			% plot from onset to offset
			plot(ax,data.laser',repmat(y,size(data.laser')),'Color',[54 161 86]/255);
		end
	end

end