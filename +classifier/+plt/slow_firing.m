function [ops,to_save] = slow_firing(data,ops,prefix,ax)

if nargin < 3
	prefix = '';
end

% count spikes with non-overlapping 1 sec window
ops.tp = getOr(ops,'tp',[0 1]);
posterior_t_edges = getOr(ops,'posterior_t_edges',data.video(1):0.5:data.video(end));
[spk_count,ops] = classifier.count_spk.time_course(data,ops);


% calculate posterior
if ~isfield(ops,'Mdl')
	[~,ops] = ops.classifier.train(data,ops);
end
spk_count = spk_count(ismember(ops.spk_count_id,ops.decoder_id),:);
Posterior = ops.classifier.predict(ops.Mdl,spk_count',ops);


% brute-force z-score firing rate
spk_count_z = (spk_count - mean(spk_count,2))./ std(spk_count,[],2);
% average more
bin_width = 60;
[slow_change,~,slow_edges,t] = running_average(ops.posterior_t,spk_count_z,bin_width,bin_width);


% select cells to plot for firing rate
ops.plot = getOr(ops,'plot',struct());
ops.plot.slow_firing_cell = getOr(ops.plot,'slow_firing_cell','novel_vs_fam');
switch ops.plot.slow_firing_cell
	case 'novel_vs_fam' % order them by novel vs fam results
		ops.novel_vs_fam.n_sig = Inf;
		[~,cell_group_plt] = classifier.select_cells.novel_vs_fam(data,ops);
	case 'non_zero_coef' % plot firing rate of the non-zero coef
		cell_group_plt = classifier.select_cells.non_zero_coef(data,ops);
	case 'by_coef'
		[~,cell_group_plt] = classifier.select_cells.by_coef(data,ops);
end
[~,ind] = ismember(cell_group_plt.ordered_id(1:cell_group_plt.ordered_div(3)),ops.decoder_id);
slow_change = slow_change(:,ind(ind~=0));
% save info
prefix = [ops.plot.slow_firing_cell '_' prefix];


% count posterior peaks
[~,locs] = arrayfun(@(ii) findpeaks(Posterior(:,ii),'MinPeakHeight',0.5), 1:2, 'UniformOutput',false);
[~,npeaks] = running_average(ops.posterior_t(locs{1}),[],bin_width,[],slow_edges);
[~,npeaks(:,end+1)] = running_average(ops.posterior_t(locs{2}),[],bin_width,[],slow_edges);



% sanity check
	% ax = np; 
	% plot(ops.posterior_t,Posterior(:,1));
	% my_scatter(ops.posterior_t(locs{1}),1.25,ax);

% plot
% change label
if isfield(data,'port_is_water') && sum(data.port_is_water)
	tmp = {'novel','fam'};
	cell_group_plt.cell_cat_name(1:2) = tmp(data.port_is_water+1);
end


% output to save
if nargout > 1
	% select time points within period of interest
	% not rebinning spike count, which is much less efficient; current margin of error small anyway
	% t_2_plot = [10 30 45]; % how many data points to take
	t_2_plot = [15 30 45]; % how many data points to take
	t_partition = [common_t.last_reward(data) common_t.first_laser(data)];
	bin_width = 60; step_size = 30;

	% better just to recount the peaks, right now npeaks is smaller
	for ii = 1:2
		to_save.npeaks(ii,:) = {hist_overlap(ops.posterior_t(locs{ii}),bin_width,[(t_partition(1)-t_2_plot(1)*60):step_size:(t_partition(1)-bin_width)]),...
					  		hist_overlap(ops.posterior_t(locs{ii}),bin_width,[(t_partition(1)):step_size:(t_partition(1)+t_2_plot(2)*60-bin_width)]),...
					  		hist_overlap(ops.posterior_t(locs{ii}),bin_width,[(t_partition(2)):step_size:(t_partition(2)+t_2_plot(3)*60-bin_width)])};
	end

	% also save all peaks time points after CGRP starts
	for ii = 1:2
		to_save.t_peak{ii} = ops.posterior_t(locs{ii}) - common_t.first_laser(data);
		to_save.t_peak{ii} = to_save.t_peak{ii}(to_save.t_peak{ii} > 0);
	end

	% switch novel to first 
	if data.port_is_water(1)
		to_save.npeaks = to_save.npeaks([2 1],:);
		to_save.t_peak = to_save.t_peak([2 1]);
	end

	% save parameters
	to_save.t_2_plot = t_2_plot;	
	to_save.bin_width = bin_width;
	to_save.step_size = step_size;

	% plot for examination
	% [ax,r,c] = np(2,3);
	% for ii = 1:2
	% 	for jj = 1:3
	% 		plot(ax(sub2ind([c r],jj,ii)),to_save.npeaks{ii,jj});
	% 	end
	% end
	% ef;

else

	% average posterior (don't need this for final figure so putting it here)
	[avg_post] = running_average(ops.posterior_t,Posterior,bin_width,[],slow_edges);

	% do no save, just plot
	ax = np(3,1); 
	% slow firing
		clear h;
		Colors = getOr(data,'port_color',[1 0 0;0 0 0]);
		h    = plot_multiple_lines(slow_change(:,1:cell_group_plt.ordered_div(2))',ax(1),'x',t,'base_color',Colors(1,:));
		h(2) = plot_multiple_lines(slow_change(:,(cell_group_plt.ordered_div(2)+1):cell_group_plt.ordered_div(3))',ax(1),'x',t,'base_color',Colors(2,:));
		classifier.plt.divider_event(data,ax(1));
		legend(ax(1),[h(1).h_m h(2).h_m],cell_group_plt.cell_cat_name,'box','off');
		ylabel(ax(1),'z-scored firing');
		
	% posterior peaks
		plot(ax(2),t,npeaks(:,1),'Color',Colors(1,:));
		plot(ax(2),t,npeaks(:,2),'Color',Colors(2,:));
		xlabel(ax(2),'time (s)'); 
		ylabel(ax(2),'# peak');
	% slow posterior
		plot(ax(3),t,avg_post(:,1),'Color',Colors(1,:));
		plot(ax(3),t,avg_post(:,2),'Color',Colors(2,:));
		xlabel(ax(3),'time (s)'); 
		ylabel(ax(3),'avg post');

	arrayfun(@(h) set(h,'XLim',t([1 end])), ax);
	set(gcf,'Position',[0 0 400 400]);
	title(ax(1),{sprintf('%s %s',data.subject,datestr(data.session,'YYmmdd')),...
				sprintf('n = %d',size(ops.Mdl.coef,1))});
	export_fig(sprintf('results/%sslow_firing_%s_%s.pdf',prefix,data.subject,datestr(data.session,'YYmmdd')));
	
end

