function h = plot_multiple_lines(data,ax,varargin)
% h = plot_multiple_lines(data,ax,varargin)
% plot the average lines while keeping individual traces
% rows are individiual samples
% data: 		Nsample * Ntime 
% parameters:	base_color
% 				x
p = inputParser;
addParameter(p,'base_color',[0.24,0.15,0.66]);
addParameter(p,'x',nan);
addParameter(p,'to_plt',{'shade','mean'});
parse(p,varargin{:});
Color = p.Results.base_color;
to_plt = p.Results.to_plt;

%% initiation
% axes
if nargin < 2 || isempty(ax)
	ax = gca; hold on;
end

% x data
if isnan(p.Results.x)
	x = (1:size(data,2)) - 1;
else
	x = p.Results.x;
end


%% calculate and plot
% plot std bounds
M = mean(data,1,'omitnan');
V = std(data,[],1,'omitnan') / sqrt(size(data,1));

% plot multiple lines as data, with rows as 
for ii = 1:numel(to_plt)
	switch to_plt{ii}
	case 'sep'
		if size(data,2) > 1
			h.h_sep = plot(ax,x,data','-','Marker','none','Color',[Color 0.1],'LineWidth',1);
		else
			h.h_sep = my_scatter(x,data,ax,'.','MarkerEdgeColor',Color,'MarkerEdgeAlpha',0.1);
		end
	case 'shade'
		h.shade = fill(ax,[x flip(x)],[M+V, flip(M-V)],Color,'FaceAlpha',0.1,'EdgeColor','none');
	case 'err'
		h.err = errorbar(ax,x,M,V,'.-','Color',Color,'LineWidth',2);
	case 'mean'
		h.h_m   = plot(ax,x,M,'Color',Color,'LineWidth',2);
	end
end
% h.h_v   = plot(ax,x,[M+V; M-V]','Color',Color,'LineStyle','--','LineWidth',1);

% h.h_sep = plot(ax,x,data','-','Marker','none','Color',[Color 0.1],'LineWidth',1);
% % h.shade = fill(ax,[x flip(x)],[M+V, flip(M-V)],Color,'FaceAlpha',0.1,'EdgeColor','none');
% h.err = errorbar(ax,x,M,V,'.-','Color',Color,'LineWidth',2);
% h.h_m   = plot(ax,x,M,'Color',Color,'LineWidth',2);
% % h.h_v   = plot(ax,x,[M+V; M-V]','Color',Color,'LineStyle','--','LineWidth',1);

end
