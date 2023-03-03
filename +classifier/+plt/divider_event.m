function divider_event(data,ax)

% last reward if applied
if isfield(data,'rewards')
	t_last_reward = max([data.rewards.all.rear(end) data.rewards.all.front(end)]);
	plot(ax,[1 1]*t_last_reward,ax.YLim,'k-','LineWidth',0.7);
	h = text(ax,t_last_reward, ax.YLim(1), sprintf('last\nreward'),'FontSize',ax.FontSize-1,'horizontalalignment', 'center', 'verticalalignment', 'top'); 
end

% need to add for opto stimulation

% licl if applied
if isfield(data,'licl') 
	plot(ax,[1 1]*data.licl,ax.YLim,'k-','LineWidth',0.7);
	text(ax,data.licl, ax.YLim(1), 'licl','FontSize',ax.FontSize-1,'horizontalalignment', 'center', 'verticalalignment', 'top'); 
end
% cgrp if applied
if isfield(data,'laser') 
	plot(ax,[1 1]*data.laser(1,1),ax.YLim,'k-','LineWidth',0.7);
	text(ax,data.laser(1,1), ax.YLim(1), sprintf('first\nstim'),'FontSize',ax.FontSize-1,'horizontalalignment', 'center', 'verticalalignment', 'top'); 
end
