function ops = posterior_raster(data,Mdl,ops,prefix)
% plot postieror traces

% initiation
if nargin < 4
	prefix = ''; % figure prefix
end



% posterior
switch getOr(ops,'posterior_method',1)
case 1 % calculate posterior

	% set classifier
	Classifier = getOr(ops,'classifier',classifier.mnr());

	% retrain a model if there's no input
	if isempty(Mdl)
		[Mdl, ops] = Classifier.train(data,ops);
	end

	% count spikes time course
	[count,ops] = classifier.count_spk.time_course(data,ops);

	% calculate posterior
	Posterior = Classifier.predict(Mdl,count',ops);

case 2 % load from Chris
	load(['/jukebox/witten/Chris/matlab/cz/neuropixels-cta/calca_' data.subject '_decoders.mat']);
	Posterior = [decoders.posteriors.Front_Novel;decoders.posteriors.Rear_Water;decoders.posteriors.Before_cue]';
	ops.posterior_t = (1:size(Posterior,1)) * 0.05;
end



% order cells for plotting raster
switch 2
case 1 % sort by ranksum
	if ~isfield(ops,'novel_vs_fam') || ~isfield(ops.novel_vs_fam,'ordered_id')
		ops = classifier.select_cells.novel_vs_fam(data,ops);
	end
case 2 % sort by coefficient

	% ops = classifier.select_cells.by_coef(data,ops);
	ops.novel_vs_fam = classifier.select_cells.non_zero_coef(data,ops);
	% % front vs back
	% d_coef = Mdl.coef(:,1) - Mdl.coef(:,2);
	% [~,I] = sort(d_coef,'descend');
	% ind   = find(~ops.exclude_id);

	% % % reward vs control
	% % d_coef = Mdl.coef(:,1) + Mdl.coef(:,2) - Mdl.coef(:,3)*2;
	% % [~,I] = sort(abs(d_coef),'descend');

	% % replace results
	% ops.novel_vs_fam.ordered_id    = [ind(I(1:15)) ind(I(end-14:end))];
	% ops.novel_vs_fam.ordered_div   = [0 15 30];
	% ops.novel_vs_fam.cell_cat_name = {'front higher','back higher'};
end
spk_2_plt = data.spikes(ops.novel_vs_fam.ordered_id);
ytick_loc = [-2.5 (ops.novel_vs_fam.ordered_div(1:end-1)+ops.novel_vs_fam.ordered_div(2:end))/2];


% plot initiation
t_step  = 100;
t_start = min(ops.posterior_t):t_step:max(ops.posterior_t);
max_ax = 1;
for ibatch = 1:ceil(numel(t_start)/max_ax)
	ax = np(max_ax,1);

	for ii = 1:numel(ax)
		
		% add raster below
		toi = [0 t_step]+t_start(ii)+(ibatch-1)*t_step*ii;
		tmp = cellfun(@(x) x(find(x>toi(1) & x<toi(2))) , spk_2_plt,'UniformOutput',false);
		plt.raster_smooth(tmp,toi,ax(ii),'kernel_width',200); 

		% line for posterior probability
		h = plot(ax(ii),ops.posterior_t,Posterior.*-5); h(2-data.port_is_water(2)).Color = [1 0 0]; h(2-data.port_is_water(1)).Color = [0 0 0]; 
		if numel(h)>2 
			h(3).Color(4) = 0.5;h(3).LineWidth = 0.7;
		end


		% add division for different cells
		plot(ax(ii),toi',repmat(ops.novel_vs_fam.ordered_div(2:end-1)+0.5,2,1),'-','LineWidth',0.7,'Color',[0.3 0.3 0.3]);


		% licl if applied
		if isfield(data,'licl') && (data.licl>=toi(1)&&data.licl<=toi(2))
			plot([1 1]*data.licl,[0 numel(spk_2_plt)+1],'k-','LineWidth',0.7);
			text(ax,data.licl, numel(spk_2_plt)+1, 'licl','FontSize',ax.FontSize-1,'horizontalalignment', 'center', 'verticalalignment', 'top'); 
		end

		% events
		h = classifier.plt.scatter_event(data,ax(ii),-5.5);
		legend([h(:,1)' h(1,2) h(1,3)],{'front reward','back','front cue','front entry'},'Location','northeastoutside');
	end

	% figure setting
	% change x axis to progress in time
	set(gcf,'Position',[0 0 800 250*max_ax]);
	arrayfun(@(i) set(ax(i),'XLim',[0 t_step]+t_start(i)+(ibatch-1)*t_step*max_ax,...
							'YLim',[-6 numel(spk_2_plt)+0.5],...
				            'YTick',ytick_loc,'YTickLabel',['posterior',ops.novel_vs_fam.cell_cat_name]), 1:numel(ax))

	% save figure
	export_fig(sprintf('results/%sposterior_raster_%s_%s_%d.pdf',prefix,data.subject,datestr(data.session,'yymmdd'),ibatch));

end
append_script(sprintf('results/%sposterior_raster_%s_%s',prefix,data.subject,datestr(data.session,'yymmdd')));




%% average posterior
% [Post_avg,~,~,t] = running_average(ops.posterior_t,Posterior>0.1,60,30);

% % plotting
% ax = np;
% plot(ax,t,Post_avg(:,1),'Color',[1 0 0]);
% plot(ax,t,Post_avg(:,2),'Color',[0 0 0]);
% classifier.plt.divider_event(data,ax);

% set(ax,'XLim',t([1 end]),'XTick',round(t([1 end])));
% xlabel('time (s)');

% export_fig(sprintf('results/%savg_post_%s_%s.png',prefix,data.subject,datestr(data.session,'YYmmdd')),'-m3');


end