function psth(data,ops)

% initiation
ops.tp = getOr(ops,'tp',[0 1]);

% select cells
ops = classifier.select_cells.novel_vs_fam(data,ops);

% parameters
events_oi = {data.rewards.all.front,data.rewards.all.rear};
toi = -0.2:0.01:3;
cell_oi = ops.novel_vs_fam.ordered_id((ops.novel_vs_fam.ordered_div(1)+1):ops.novel_vs_fam.ordered_div(4));

% count spikes for z-score
ops.posterior_t_edges = min([data.cues.all.front(1) data.cues.all.rear(1)]):max([data.rewards.all.front(end) data.rewards.all.rear(end)]); % first cue to last reward
[spk_count,ops] = classifier.count_spk.time_course(data,ops);

% calculate psth
for jj = 1:numel(events_oi)
	resp = nan(numel(cell_oi),numel(toi));
	for ii = 1:numel(cell_oi)
		[resp(ii,:),resp_err,RR,raster] = cal_psth(data.spikes{cell_oi(ii)}*1000,events_oi{jj}*1000,'tp',toi,'kernel_width',0.1);
	end
	% z score by precue spike count
	resp_zscored{jj} = (resp - mean(spk_count(cell_oi,:),2)) ./ var(spk_count(cell_oi,:),[],2);
	resp_zscored{jj}(resp_zscored{jj}==Inf) = 0;
end



% plotting
ax = np(1,2);
% plot psth
imagesc(ax(1),toi,1:numel(cell_oi),resp_zscored{1});
imagesc(ax(2),toi,1:numel(cell_oi),resp_zscored{2});
% plot divider
arrayfun(@(h) plot(h,toi([1 end])',repmat(ops.novel_vs_fam.ordered_div(2:end-1)+0.5,2,1),'-','LineWidth',0.7,'Color',[0.3 0.3 0.3]),ax);
arrayfun(@(h) plot(h,[0;0],[0.5 numel(cell_oi)+0.5],'--','LineWidth',0.5,'Color',[0.3 0.3 0.3]),ax);
% other figure setting
colormap(flip(cbrewer2('RdBu')));
arrayfun(@(h) set(h,'CLim',[-1 1]*max(abs(h.CLim)),'YDir','reverse','YLim',[0.5 numel(cell_oi)+0.5]),ax);
set(ax(1),'YTick',(ops.novel_vs_fam.ordered_div(1:3)+ops.novel_vs_fam.ordered_div(2:4))/2,'YTickLabel',{'front higher','back higher','no difference'});
set(ax(2),'YTick',[]);
arrayfun(@(h) xlabel(h,'time (s)'), ax);
title(ax(1),'front reward'); title(ax(2),'back reward');

% position
ax(1).Position = [0.2409    0.3328    0.2153    0.4997];
ax(2).Position = [0.5    0.3328    0.2153    0.4997];

export_fig(sprintf('results/psth_%s_%s.pdf',data.subject,datestr(data.session,'YYmmdd')));

end