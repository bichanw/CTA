function h = drinking_event(data,ax,y)
% Input: y - location on y axis 
% Only plot actual drinking for figure

	% initiation 
	if nargin < 3
		y = 0;
	end

	% change front / familiar color
	Colors = getOr(data,'port_color',[0.9412    0.2118    0.2745;   0.0667    0.5020    0.9451]);

	% actual reward
	h(1) = my_scatter(data.rewards.all.front,y,ax,'v','MarkerEdgeColor',Colors(1,:));
	h(2) = my_scatter(data.rewards.all.rear,y,ax,'v','MarkerEdgeColor',Colors(2,:));


	% laser
	tmp = plot(ax,data.laser',repmat(y,size(data.laser')),'Color',[54 161 86]/255);
	h(3) = tmp(1); % only label one of them

end