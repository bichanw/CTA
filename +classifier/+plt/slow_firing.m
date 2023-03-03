function slow_firing(data,ops,prefix,ax)

if nargin < 3
	prefix = '';
end

% no training, zscore and calculate firing; count spikes with non-overlapping 1 sec window
ops.tp = getOr(ops,'tp',[0 1]);
% save temporary file for skipping counting spikes
posterior_t_edges = getOr(ops,'posterior_t_edges',data.video(1):0.5:data.video(end));
[spk_count,ops] = classifier.count_spk.time_course(data,ops);


% calculate posterior
if ~isfield(ops,'Mdl')
	[~,ops] = ops.classifier.train(data,ops);
end
Posterior = ops.classifier.predict(ops.Mdl,spk_count',ops);


% brute-force z-score firing rate
spk_count_z = (spk_count - mean(spk_count,2))./ std(spk_count,[],2);
% average more
bin_width = 60;
[slow_change,~,slow_edges,t] = running_average(ops.posterior_t,spk_count_z,bin_width,60);

switch 1
case 1 % order them by novel vs fam results
	cell_group_plt = ops.novel_vs_fam;
case 2 % plot firing rate of the non-zero coef
	cell_group_plt = classifier.select_cells.non_zero_coef(data,ops);
end
[~,ind] = ismember(cell_group_plt.ordered_id(1:cell_group_plt.ordered_div(3)),ops.spk_count_id);
slow_change = slow_change(:,ind);

% count posterior peaks
[~,locs] = arrayfun(@(ii) findpeaks(Posterior(:,ii),'MinPeakHeight',0.3), 1:2, 'UniformOutput',false);
[~,npeaks] = running_average(ops.posterior_t(locs{1}),[],bin_width,[],slow_edges);
[~,npeaks(:,end+1)] = running_average(ops.posterior_t(locs{2}),[],bin_width,[],slow_edges);


% plot
ax = np(2,1); 
% slow firing
	clear h;
	Colors = getOr(data,'port_color',[1 0 0;0 0 0]);
	h    = plot_multiple_lines(slow_change(:,1:cell_group_plt.ordered_div(2))',ax(1),'x',t,'base_color',Colors(1,:));
	h(2) = plot_multiple_lines(slow_change(:,(cell_group_plt.ordered_div(2)+1):cell_group_plt.ordered_div(3))',ax(1),'x',t,'base_color',Colors(2,:));
	classifier.plt.divider_event(data,ax(1));
	legend(ax(1),[h(1).h_m h(2).h_m],cell_group_plt.cell_cat_name,'box','off');
	ylabel(ax(1),'z-scored firing');
% slow posterior
	plot(ax(2),t,npeaks(:,1),'Color',[1 0 0]);
	plot(ax(2),t,npeaks(:,2),'Color',[0 0 0]);
	xlabel(ax(2),'time (s)'); 
	ylabel(ax(2),'# peaks');

arrayfun(@(h) set(h,'XLim',t([1 end])), ax);
set(gcf,'Position',[0 0 400 300]);
export_fig(sprintf('results/%sslow_firing_%s_%s.pdf',prefix,data.subject,datestr(data.session,'YYmmdd')));


