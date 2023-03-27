

%% Code to produce required figures

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

