f = 'mat/280_230214_assignments.jld2';
final_assignments = h5read(f,'/final_assignments');
final_events      = h5read(f,'/final_events');


% assign spikes to event type
[~,Locb]   = ismember(final_assignments,final_events.assignment_id);

tmp = find(Locb~=0);
final_type = final_assignments;
final_type(tmp) = final_events.seq_type(Locb(tmp));


% load spikes
load('mat/julia.mat');

% plot spikes
ax = np;

Colors = [0 0 0;1 0 0;0 0 1];

% plot spikes
cell_oi  = unique(clu);
% final_type = ones(size(final_type)); % plot all black to see raw raster
event_oi = unique(final_type);
for ineuron = 1:numel(cell_oi)
	for ievent = 1:numel(event_oi)
		h = my_scatter(t(clu==cell_oi(ineuron)&final_type==event_oi(ievent)),ineuron,ax,...
						3,'MarkerEdgeColor','none','MarkerFaceColor',Colors(ievent,:),'MarkerFaceAlpha',0.1);
	end
end

% plot other events
% ?
data = load_data.all(datetime(2023,2,14),'280');
h = classifier.plt.scatter_event(data,ax,numel(cell_oi)+1);

for t_start = 300:100:900
	set(ax,'XLim',[0 100] + t_start);
	export_fig(sprintf('results/ppseq_%d.png',t_start),'-r300');
end


