function ops = raster(data,ops)
ops = classifier.select_cells.novel_vs_fam(data,ops);

% initiation
ytick_loc = (ops.novel_vs_fam.ordered_div(1:end-1)+ops.novel_vs_fam.ordered_div(2:end))/2;
spk_2_plt = data.spikes(ops.novel_vs_fam.ordered_id);
T_STARTS = colon_right(0,100,max(cellfun(@(x) max(x),spk_2_plt))-100);
ax = np;
for it = 1:numel(T_STARTS)
	t_start = T_STARTS(it);
	toi = [0 100] + t_start;
	fig_name = sprintf('results/raster_%s_%s_%d.pdf',data.subject,datestr(data.session,'yymmdd'),t_start/100+1);

	if exist(fig_name,'file') continue; end

	% raster
	tmp = cellfun(@(x) x(find(x>toi(1) & x<toi(2))) , spk_2_plt,'UniformOutput',false);
	cla(ax);
	
	% plot raster
	plt.raster_smooth(tmp,toi,ax,'kernel_width',200); 
	% plt.raster2(tmp,[],[],ax,0,'k-','LineWidth',0.2); 

	% draw color background
	% classifier.plt.color_event(data,ax);

	% event legend at top
	h = classifier.plt.scatter_event(data,ax);
	% licl if applied
	if isfield(data,'licl') && (data.licl>=toi(1)&&data.licl<=toi(2))
		plot([1 1]*data.licl,[0 numel(spk_2_plt)+1],'k-','LineWidth',0.7);
		text(ax,data.licl, numel(spk_2_plt)+1, 'licl','FontSize',ax.FontSize-1,'horizontalalignment', 'center', 'verticalalignment', 'top'); 
	end
	
	% add division for different cells
	plot(toi',repmat(ops.novel_vs_fam.ordered_div(2:end-1)+0.5,2,1),'-','LineWidth',0.7,'Color',[0.3 0.3 0.3]);

	% other figure setting
	legend([h(:,1)' h(1,2) h(1,3)],{'front reward','back','front cue','front entry'},'Location','northeastoutside');
	set(gcf,'Position',[0 0 diff(toi)*7 numel(spk_2_plt)*5]);
	set(gca,'YLim',[0 numel(spk_2_plt)+1],'YDir','reverse','XLim',toi,'YTick',ytick_loc,'YTickLabel',ops.novel_vs_fam.cell_cat_name);

	export_fig(fig_name);
end

append_script(sprintf('results/raster_%s_%s',data.subject,datestr(data.session,'yymmdd')));
