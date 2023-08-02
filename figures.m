%% paper figure
% general figure setting
Colors = [228 45 38;   55 135 192; 54 161 86] / 255; % novel, familiar, CGRP



% raster
	h_posterior = 3;
	Axes = np(1,2);

	% add raster below
	Tois = [670 720; 7970 8010]; % 620 720
	for ii = 1:2
		% convolve raw spikes and average
		ax = Axes(ii);
		toi = Tois(ii,:);
		tmp = cellfun(@(x) x(find(x>toi(1) & x<toi(2))) , spk_2_plt,'UniformOutput',false);
		switch 'smooth_colored'
		case 'smooth_bw'
			% same color
			plt.raster_smooth(tmp,toi,ax,'kernel_width',300); 
		case 'smooth_colored'
			ops.novel_vs_fam.ordered_color = [data.port_color; 0 0 0];
			% color raster by cell category
			for jj = 1:(numel(ops.novel_vs_fam.ordered_div)-1)
				plt.raster_smooth(tmp((ops.novel_vs_fam.ordered_div(jj)+1):(ops.novel_vs_fam.ordered_div(jj+1))),...
								toi,ax,'dy',ops.novel_vs_fam.ordered_div(jj),'kernel_width',300,'base_color',ops.novel_vs_fam.ordered_color(jj,:,:)); 
			end
		case 'lines'
			ops.novel_vs_fam.ordered_color = [data.port_color; 0 0 0];
			for jj = 1:(numel(ops.novel_vs_fam.ordered_div)-1)
				plt.raster2(tmp((ops.novel_vs_fam.ordered_div(jj)+1):(ops.novel_vs_fam.ordered_div(jj+1))),...
								[],[],ax,ops.novel_vs_fam.ordered_div(jj),'Color',ops.novel_vs_fam.ordered_color(jj,:),'LineWidth',0.1); 
			end
		end

		% line for posterior probability
		h = plot(ax,ops.posterior_t,Posterior(:,1:2) .* -h_posterior); 
		arrayfun(@(ii) set(h(ii),'Color',data.port_color(ii,:)), 1:2);

		% add division for different cells
		plot(ax,toi',repmat(ops.novel_vs_fam.ordered_div(2:end-1)+0.5,2,1),'-','LineWidth',0.7,'Color',[0.3 0.3 0.3]);


		% licl if applied
		if isfield(data,'licl') && (data.licl>=toi(1)&&data.licl<=toi(2))
			plot([1 1]*data.licl,[0 numel(spk_2_plt)+1],'k-','LineWidth',0.7);
			text(ax,data.licl, numel(spk_2_plt)+1, 'licl','FontSize',ax.FontSize-1,'horizontalalignment', 'center', 'verticalalignment', 'top'); 
		end

		% events
		h = classifier.plt.drinking_event(data,ax,-h_posterior * 1.2);

	end

	sprintf('\\color[rgb]{%f,%f,%f}Novel',data.port_color(1,:))

	% figure setting
	% change x axis to progress in time
	legend(Axes(2),h,{'novel','water','CGRP'},'Location','northeastoutside');
	set(gcf,'Position',[0 0 600 250]);
	set(Axes(1),'XLim',Tois(1,:),'XTick',[],...
			'YLim',[-h_posterior-1 numel(spk_2_plt)+0.5],'YDir','reverse',...
			'YTick',ytick_loc,'YTickLabel',{'posterior',sprintf('\\color[rgb]{%f,%f,%f}Novel',data.port_color(1,:)),sprintf('\\color[rgb]{%f,%f,%f}Water',data.port_color(2,:)),'Non-selective'},'YTickLabelRotation',90);
	set(Axes(2),'XLim',Tois(2,:),'XTick',[],...
			'YLim',[-h_posterior-1 numel(spk_2_plt)+0.5],...
			'YTick',[],'YDir','reverse');
	Axes(2).YAxis.Visible = 'off';		

	Axes(1).Position = [0.1300    0.1100    0.3347    0.8150];
	Axes(2).Position = [0.49    0.1100    0.3347    0.8150];
	return
% averaged laser posterior
	load('figures/laser_posterior.mat');
	% resp - 6 sessions * 3 time periods * 201 time points * 3 conditions
	% tbin - time stamps for 201 time points
	% ax = np;
	% plot(tbin,squeeze(mean(resp(:,1,:,1),1)));ef;

	% actual plot
	ax = np(2,1);
	% 280
	arrayfun(@(ii) plot(ax(1),tbin,squeeze(resp(2,3,:,ii)),'Color',Colors(ii,:)), 1:2);
	% average
	arrayfun(@(ii) plot_multiple_lines(squeeze(resp(:,3,:,ii)),ax(2),'x',tbin,'base_color',Colors(ii,:)), 1:2);

	% figure setting
	set(ax(1),'XLim',tbin([1 end]),'XTick',[],'YLim',[0 0.8],'YTick',[0 0.8]);
	set(ax(2),'XLim',tbin([1 end]),'XTick',[0 1 2 3],'XTickLabel',{'onset','1','2','offset'});
	xlabel(ax(2),'Time (s)');

	ax(1).Position = [0.1485    0.5    0.7565    0.2630];
	ax(2).Position = [0.1485    0.1882    0.7565    0.2630];
	set(gcf,'Position',[0 0 150 250]);

	% inspection
		% ax = np(3,1);
		% for iplot = 1:3
		% 	plot_multiple_lines(squeeze(resp(:,iplot,:,1)),ax(iplot),'x',tbin,'base_color',Colors(1,:));
		% 	plot_multiple_lines(squeeze(resp(:,iplot,:,2)),ax(iplot),'x',tbin,'base_color',Colors(2,:));
		% end
		% align_ax(ax,true,true);
		% % temporary figure setting (might want to select 1 plot in the future)
		% xlabel(ax(3),'Time (s)');
		% ylabel(ax(2),'Posterior probability');

	
% averaged peak number
	load('figures/decoder_peaks.mat');

	% concatenate into 2 traces
	clear tmp
	for ii = 1:2
		for jj = 1:3
			tmp{ii,jj} = catcell(arrayfun(@(kk) to_save(kk).npeaks{ii,jj}, 1:6,'uni',0),2);
		end
	end
	npeaks_all = catcell(tmp(1,:),1);
	npeaks_all(:,:,2) = catcell(tmp(2,:),1);
	l_period = [0 cumsum(cellfun(@(x) size(x,1),tmp(1,:)))] * to_save(1).step_size / 60;

	% plot 1 trace
	h_top = 12.5; % top of the trace
	ax = np;
	% calculate time
	t = ((1:size(npeaks_all,1))-1) * to_save(1).step_size / 60;

	% plot peak number
	arrayfun(@(ii) plot_multiple_lines(squeeze(npeaks_all(:,:,ii))',ax,'x',t,'base_color',Colors(ii,:)), 1:2);
	% time period
	color_period = [0 0 0; 0.5 0.5 0.5; Colors(3,:,:)];
	arrayfun(@(ii) plot(ax,l_period([1:2]+ii-1) + [1 -1]*0.5,[1 1]*(h_top+0.5),'Color',color_period(ii,:),'LineWidth',2), 1:3);
	label_period = {'Drinking','Delay','CGRP stim'};
	h = arrayfun(@(ii) text(ax,mean(l_period([1:2]+ii-1))+0.5,h_top+0.6,label_period{ii},'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',ax.FontSize-1),1:3);


	% figure setting
	set(ax,'XLim',t([1 end]),'YLim',[0 (h_top+0.5)]);
	set(gcf,'Position',[0 0 400 225]);
	xlabel('Time (min)'); ylabel('# of (re)-activations');
	ef;


	% inspect as 3 subplots
	% [ax,r,c] = np(1,3);
	% for ii = 1:2
	% 	for jj = 1:3
	% 		tmp = catcell(arrayfun(@(kk) to_save(kk).npeaks{ii,jj}, 1:6,'uni',0),2);
	% 		plot_multiple_lines(tmp',ax(jj),'base_color',Colors(ii,:));
	% 	end
	% end
	% align_ax(ax,false,true);
	% ef;

return
%% K99 figures

% chris k99, figure 2c raster
	data = load_data.all(datetime(2022,12,4),'002');
	ops = struct('tp',[0 1]);
	ops  = classifier.select_cells.novel_vs_fam(data,ops);
	cell_oi = ops.novel_vs_fam.ordered_id((ops.novel_vs_fam.ordered_div(1)+1):ops.novel_vs_fam.ordered_div(3));
	ind_front = (ops.novel_vs_fam.ordered_div(1)+1):ops.novel_vs_fam.ordered_div(2);
	ind_rear  = (ops.novel_vs_fam.ordered_div(2)+1):ops.novel_vs_fam.ordered_div(3);
	% retrain decoder
	ops.exclude_id = ~ismember(1:numel(data.spikes),cell_oi);
	ops.mnr = struct('penalty','l2','lambda',1,'zscore',false);
	[Mdl,ops] = classifier.mnr.train(data,ops);

	% figure 2c psth
		% count spikes for z-score
		ops.posterior_t_edges = min([data.cues.all.front(1) data.cues.all.rear(1)]):max([data.rewards.all.front(end) data.rewards.all.rear(end)]); % first cue to last reward
		[spk_count,ops] = classifier.count_spk.time_course(data,ops);

		% parameters
		events_oi = {data.rewards.all.front,data.rewards.all.rear};
		toi = -1:0.01:5;
		cell_oi = ops.novel_vs_fam.ordered_id((ops.novel_vs_fam.ordered_div(1)+1):ops.novel_vs_fam.ordered_div(3));
	
		% calculate psth
		m = mean(spk_count(cell_oi,:),2);
		v = var(spk_count(cell_oi,:),[],2);
		for jj = 1:numel(events_oi)
			resp = nan(numel(cell_oi),numel(toi));
			for ii = 1:numel(cell_oi)
				[resp(ii,:),resp_err,RR,raster] = cal_psth(data.spikes{cell_oi(ii)}*1000,events_oi{jj}*1000,'tp',toi,'kernel_width',0.1);
			end

			% z score by overall
			resp_zscored{jj} = (resp - m) ./ v;
			resp_zscored{jj}(resp_zscored{jj}==Inf) = 0;
		end

		% figure setting
		ind_front = (ops.novel_vs_fam.ordered_div(1)+1):ops.novel_vs_fam.ordered_div(2);
		ind_rear  = (ops.novel_vs_fam.ordered_div(2)+1):ops.novel_vs_fam.ordered_div(3);

		
		ax = np(2,2);
		imagesc(ax(1),toi,1:numel(cell_oi),resp_zscored{1}(ind_front,:)); % novel resp, novel cell
		imagesc(ax(3),toi,1:numel(cell_oi),resp_zscored{1}(ind_rear,:)); % novel resp, water cell
		imagesc(ax(2),toi,1:numel(cell_oi),resp_zscored{2}(ind_front,:)); % novel resp, novel cell
		imagesc(ax(4),toi,1:numel(cell_oi),resp_zscored{2}(ind_rear,:)); % novel resp, novel cell

		arrayfun(@(h) set(h,'CLim',[-1 1]*2,'YDir','reverse','Visible','off'),ax);
		colormap(flip(cbrewer2('RdBu')));
		ef;

		% need average posterior
			events_oi = {data.rewards.all.front,data.rewards.all.rear};
			jj = 1;
			toi = -2:0.1:6;

			ax = np(1,2);
			for jj = 1:2
				Posterior = nan(numel(events_oi{jj}),numel(toi),3);
				for ii = 1:numel(events_oi{jj})
					ops.posterior_t_edges = toi + events_oi{jj}(ii);
					[count,ops] = classifier.count_spk.time_course(data,ops);
					Posterior(ii,:,:) = classifier.mnr.predict(Mdl,count',ops);
				end

				avg_post = squeeze(mean(Posterior,1));
				fig_line_shade(ax(jj),toi,avg_post(:,1),[1 0 0]); 
				fig_line_shade(ax(jj),toi,avg_post(:,2),[0 0 0]);
			end
			arrayfun(@(h) set(h,'Visible','off','XLim',[-1 5]),ax);
			efig;
			ef;


	% figure 2c raster
		% toi = [-1 5]+ 479;
		% toi = [-1 5]+ 1596.3;
		% toi = [-1 5]+ 969.5;

		% for t = [479, 1596.3, 969.5]
		for t = [2740,4281]
			ax = fig_raster(data,ops,t,cell_oi,ind_front,ind_rear);
			export_fig(sprintf('%.1f.pdf',t));
		end
		tmp = cellfun(@(x) x(find(x>toi(1) & x<toi(2))) , data.spikes(cell_oi),'UniformOutput',false);
		ax = np(2,1);
		plt.raster_smooth(tmp(ind_front),toi,ax(1),'kernel_width',100); 
		plt.raster_smooth(tmp(ind_rear),toi,ax(2),'kernel_width',100); 
		set(gcf,'Position',[0 0 250 400]);
		arrayfun(@(h) set(h,'Visible','off','XLim',toi),ax);
		colorbar off;
		ef;

	% figure 2c posterior
		toi = [479 969.5 1596.3];
		Colors = [1 0 0;0 0 0];

		ax  = np(1,3);
		for ii = 1:3
			% count spikes time course
			ops.posterior_t_edges = (-2:0.2:6)+toi(ii);
			[count,ops] = classifier.count_spk.time_course(data,ops);
			Posterior = classifier.mnr.predict(Mdl,count',ops);

			for jj = 1:2
				h = plot(ax(ii),ops.posterior_t,Posterior(:,jj),'Color',Colors(jj,:),'LineWidth',0.7);
				area(ax(ii),ops.posterior_t,Posterior(:,jj),'FaceColor',Colors(jj,:),'EdgeColor','none','FaceAlpha',0.1);
			end
		end


		arrayfun(@(ii) set(ax(ii),'XLim',[-1 5]+toi(ii),'Visible','off'), 1:3);
		ef;

	% figure 2c reactivation, just use the new function

		

		ax = fig_raster(data,ops,1966,cell_oi,ind_front,ind_rear);
		ax = fig_raster(data,ops,t,cell_oi,ind_front,ind_rear);
		ax = fig_raster(data,ops,3230,cell_oi,ind_front,ind_rear);


	% load trained model and corresponding parameters
	load('figures/k99_2c.mat');

	% define time points to plot
	ops.posterior_t_edges = 3780:0.5:3810; % left edges for spike counting windows

	% count spk
	classifier.plt.posterior_raster(data,ops);
	[spk_count,ops] = classifier.count_spk.time_course(data,ops);

