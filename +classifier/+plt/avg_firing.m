function avg_firing(data,ops)

% rerun cell seletion
if ~isfield(ops,'novel_vs_fam')
	ops.novel_vs_fam.n_sig = Inf;
	ops = classifier.select_cells.novel_vs_fam(data,ops);
end

% parameter
tbin = data.video(1):30:data.video(end);
bin_width = 60;
% average slow dynamic
spk_rate = nan(numel(ops.novel_vs_fam.ordered_id),numel(tbin));
for ii = 1:numel(ops.novel_vs_fam.ordered_id)
	[~,R,~,t] = running_average(data.spikes{ops.novel_vs_fam.ordered_id(ii)},[],bin_width,[],tbin);
	spk_rate(ii,:) = R / bin_width;
end

% plot average firing
ax = np(4,1);
for ii = 1:numel(ax)
	plot_multiple_lines(spk_rate(1:ops.novel_vs_fam.ordered_id((ops.novel_vs_fam.ordered_div(ii)+1):ops.novel_vs_fam.ordered_div(ii+1)),:),ax(ii),'x',t);
	classifier.plt.divider_event(data,ax(ii));
end
% figure setting
arrayfun(@(h) set(h,'YLim',[0 10]), ax);


export_fig(sprintf('results/avr_firing_%s_%s.pdf',data.subject,datestr(data.session,'YYmmdd')));

end