%% paper figure
% general figure setting
Colors = [228 45 38;   55 135 192; 54 161 86] / 255; % novel, familiar, CGRP


% extended figure, log likelihood results
	load('mat/230923.mat');
	clear M V;
	for jj = 1:numel(to_save)
		M(jj,:) = arrayfun(@(ii) mean(to_save(jj).cv_results(ii).LL(:,2)), 1:numel(to_save(jj).cv_results));
		% V = arrayfun(@(ii) std(to_save(jj).cv_results(ii).LL(:,2)), 1:numel(to_save(jj).cv_results));
		% errorbar(ax,1:numel(to_save(jj).cv_results),M,V,'o-','Color',[0 0 0],'LineWidth',1);
		% errorbar(ax,-4:4,M,V,'o-','LineWidth',1);
	end
	V = std(M,[],1);
	M = mean(M,1);
	ax = np;
	errorbar(ax,-4:4,M,V,'o-k','LineWidth',1);
	set(ax,'FontSize',8,'XLim',[-4.5 4.5],'XTick',[-4 0 4],'XTickLabel',{'10^{-4}','10^0','10^{4}'});
	% ylim([-120 0]);
	xlabel(ax,'\lambda','FontSize',11);
	ylabel(ax,'test LL','FontSize',11);
	set(gcf,'Position',[0 0 150 100]);ef;
	ef;
% extended figure, confusion matrix
	label = [];
	for jj = 1:numel(to_save)
		label = [label; to_save(jj).cv_results(5).true_label, to_save(jj).cv_results(5).post_label];
	end
	N = histcounts2(label(:,1),label(:,2));
	N = N ./ sum(N,2);


	ax = np;
	% need to use histcounts 2 and calculate probability
	imagesc(ax,1:3,1:3,N,[0 1]);
	cmap = cbrewer2('Reds'); cmap(1,:) = 1; colormap(cmap);
	c = colorbar; c.Label.String = 'classification frequency';
	set(ax,'FontSize',8,'XLim',[0.5 3.5],'YLim',[0.5 3.5],...
				'XTick',[1:3],'XTickLabel',{'novel','water','baseline'},...
				'YTick',[1:3],'YTickLabel',{'novel','water','baseline'},'YDir','reverse');
	xlabel(ax,'predicted','FontSize',11,'FontWeight','bold');
	ylabel(ax,'true','FontSize',11,'FontWeight','bold');
	set(gcf,'Position',[0 0 150 110]);
	ef;
% extended figure, posterior locked to drinking
	load('mat/230925_3s.mat');
	tbin = tbin + 0.5; 
	close all; clear ax;
	ax(1) = subplot(2,2,1,'NextPlot','add','FontSize',8);
	ax(2) = subplot(2,2,2,'NextPlot','add','FontSize',8);
	% move by 2/bin (if we're still using 1 s window for the binning)
	% 280
	for jj = 1:2
		for ii = 1:2
			M = squeeze(resp(2,jj,:,ii));
			V = squeeze(resp_err(2,jj,:,ii));
			fill(ax(jj),[tbin flip(tbin)],[M+V; flip(M-V)],Colors(ii,:),'FaceAlpha',0.1,'EdgeColor','none');
			plot(ax(jj),tbin,M,'Color',Colors(ii,:),'LineWidth',2);
		end
	end
	% done with this section

	% averaged across all 180 drinking events
	ax(3) = subplot(2,2,3,'NextPlot','add','FontSize',8);
	ax(4) = subplot(2,2,4,'NextPlot','add','FontSize',8);
	cla(ax(3)); cla(ax(4));
	for ievent = 1:2
		for iPost = 1:2
			t = cat(2,raw_data(:,ievent,:).t) + 0.5;
			sig_oi = catcell(arrayfun(@(s) s.sig_oi(iPost,:),raw_data(:,ievent,:),'uni',0),2);
			[M,~,~,tbin2,V] = running_average(t(:),sig_oi(:),0.15,0.15);
			fill(ax(ievent+2),[tbin2 flip(tbin2)],[M+V; flip(M-V)],Colors(iPost,:),'FaceAlpha',0.1,'EdgeColor','none');
			plot(ax(ievent+2),tbin2,M,'Color',Colors(iPost,:),'LineWidth',2);
		end
	end

	
	% figure setting
	arrayfun(@(h) set(h,'XLim',[-5 10]),ax);
	arrayfun(@(h) set(h,'XTick',[],'YLim',[0 1],'YTick',[0 0.5 1],'YTickLabel',{'0','50','100'}),ax);
	arrayfun(@(h) set(h,'YTick',[]),ax([2 4]));
	ylabel(ax(1),'Flavor decoder output (%)','FontSize',11);
	arrayfun(@(h) set(h,'XTick',-5:5:15), ax([3 4]));
	xlabel(ax(3),'Time (s)','FontSize',11);
	title(ax(1),'Drinking novel','FontSize',11);
	title(ax(2),'Drinking water','FontSize',11);
	h = text(ax(1),7,0.5,'novel','Color',Colors(1,:),'FontSize',8);
	h = text(ax(1),7,0.3,'water','Color',Colors(2,:),'FontSize',8);

	ylabel(ax(2),'Example animal','Rotation',0);
	ylabel(ax(4),'All animals','Rotation',0);



	set(gcf,'Position',[0 0 330 200]); ef;

	ef;



% examine firing rate for included neurons
	% FR = cal_FR(data.spikes(ops.novel_vs_fam.ordered_id));
	% ax = np; plot(FR);
	% xline(ax,ops.novel_vs_fam.ordered_div(2:end-1),'k--');
	% xticks([]);
	% xlabel('neurons');
	% ylabel('Firing rate (Hz)');


	% edges = 0:0.5:30;
	% ax = np;
	% for ii = 1:3
	% 	if ii < 3
	% 		c = Colors(ii,:);
	% 	else
	% 		c = [0 0 0];
	% 	end
	% 	histogram(ax,cal_FR(data.spikes((ops.novel_vs_fam.ordered_id(ii)+1):ops.novel_vs_fam.ordered_id(ii+1))),'FaceColor',c,...
	% 				'BinEdges',edges,'Normalization','probability');
	% end
	% export_fig tmp.png -r300;

% fig 4h: raster
	load('figures/280_old.mat');
	h_posterior = 3;
	Axes = np(1,2);

	% add raster below
	% Tois = [673.5 723.5; 7970 8010]; 
	Tois = [673.5 753.5; 7970 8040]; 
	for ii = 1:2
		% convolve raw spikes and average
		ax = Axes(ii);
		ax.FontSize = 14;
		toi = Tois(ii,:);
		spks_in_toi = cellfun(@(x) x(find(x>toi(1) & x<toi(2))) , spk_2_plt,'UniformOutput',false);
		FR  = cellfun(@(x) numel(x), spks_in_toi); % spike count for sorting
		switch 'smooth_colored'
		case 'smooth_bw'
			% same color
			plt.raster_smooth(spks_in_toi,toi,ax,'kernel_width',300); 
		case 'smooth_colored'
			ops.novel_vs_fam.ordered_color = [data.port_color; 0 0 0];
			% color raster by cell category
			for jj = 1:(numel(ops.novel_vs_fam.ordered_div)-1)
				plt.raster_smooth(spks_in_toi((ops.novel_vs_fam.ordered_div(jj)+1):(ops.novel_vs_fam.ordered_div(jj+1))),...
								toi,ax,'dy',ops.novel_vs_fam.ordered_div(jj),'kernel_width',80,'base_color',ops.novel_vs_fam.ordered_color(jj,:,:),'bin_width',100); 
			end
		case 'lines'
			ops.novel_vs_fam.ordered_color = [data.port_color; 0 0 0];
			for jj = 1:(numel(ops.novel_vs_fam.ordered_div)-1)
				% plt.raster2(spks_in_toi((ops.novel_vs_fam.ordered_div(jj)+1):(ops.novel_vs_fam.ordered_div(jj+1))),...
				% 				[],[],ax,ops.novel_vs_fam.ordered_div(jj),'MarkerEdgeColor',ops.novel_vs_fam.ordered_color(jj,:),'MarkerFaceColor',ops.novel_vs_fam.ordered_color(jj,:),'LineWidth',0.3); 
				tmp = (ops.novel_vs_fam.ordered_div(jj)+1):(ops.novel_vs_fam.ordered_div(jj+1));
				[~,I] = sort(FR(tmp));
				plt.raster2(spks_in_toi(tmp(I)),...
								[],[],ax,ops.novel_vs_fam.ordered_div(jj),'MarkerEdgeColor',ops.novel_vs_fam.ordered_color(jj,:),'MarkerFaceColor',ops.novel_vs_fam.ordered_color(jj,:),'LineWidth',0.3); 
			end
		end

		% line for posterior probability
		h = plot(ax,ops.posterior_t,Posterior(:,1:2) .* -h_posterior,'LineWidth',2); 
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


	% figure setting
	% change x axis to progress in time
	title(Axes(1),'Drinking Period'); title(Axes(2),'CGRP');
	legend(Axes(2),h,{'novel','water','CGRP'},'Location','northeastoutside');
	set(Axes(1),'XLim',Tois(1,:),'XTick',[],...
			'YLim',[-h_posterior-1 numel(spk_2_plt)+0.5],'YDir','reverse',...
			'YTick',ytick_loc,'YTickLabel',{'decoder',sprintf('\\color[rgb]{%f,%f,%f}Novel',data.port_color(1,:)),sprintf('\\color[rgb]{%f,%f,%f}Water',data.port_color(2,:)),'Non-selective'},'YTickLabelRotation',90);
	set(Axes(2),'XLim',Tois(2,:),'XTick',[],...
			'YLim',[-h_posterior-1 numel(spk_2_plt)+0.5],...
			'YTick',[],'YDir','reverse');
	Axes(2).YAxis.Visible = 'off';		


	ax_width = 0.3347;
	ax_l = 0.13;
	h1 = annotation('line',ax_l + [0 ax_width/16],[0.19 0.19]);
	h2 = annotation('textbox',[ax_l+0.0008 0.2-0.11 ax_width/10+0.05 0.1],'String','5 s','FontSize',14,'FitBoxToText','off','LineStyle','none');
	% set(gcf,'Position',[0 0 1375 556]);
	set(gcf,'Position',[0 0 1097 556]);
	Axes(1).Position = [ax_l    0.200    ax_width    0.7];
	Axes(2).Position = [0.49    0.200    ax_width    0.7];

	% export_fig('fig4h.pdf');
	export_fig('tmp.pdf');
	return

% tmp figure averaged CGRP
	load('figures/0.15.mat');
	% actual plot
	close all; clear ax;
	ax(1) = subplot(2,1,1,'NextPlot','add','FontSize',11);
	% move by 2/bin (if we're still using 1 s window for the binning)
	tbin = tbin + 0.5; 
	% 280
	M = squeeze(resp(2,3,:,:));
	V = squeeze(resp_err(2,3,:,:));
	arrayfun(@(ii) fill(ax,[tbin flip(tbin)],[M(:,ii)+V(:,ii); flip(M(:,ii)-V(:,ii))],Colors(ii,:),'FaceAlpha',0.1,'EdgeColor','none'), 1:2);
	arrayfun(@(ii) plot(ax,tbin,M(:,ii),'Color',Colors(ii,:),'LineWidth',2), 1:2);
	% mark laser
	y = 0.83;
	h = plot(ax,[0 3],[1 1]*y,'LineWidth',2.5,'Color',Colors(3,:));

	% averaged across all 100 activations
	ax(2) = subplot(2,1,2,'NextPlot','add','FontSize',11);
	t = cat(1,raw_data(:).t);
	sig_oi = catcell(arrayfun(@(ii) raw_data(ii).sig_oi(1,:),1:600,'uni',0));
	[M,~,~,tbin2,V] = running_average(t(:),sig_oi(:),0.15,0.15);
	sig_oi = catcell(arrayfun(@(ii) raw_data(ii).sig_oi(2,:),1:600,'uni',0));
	[M(:,end+1),~,~,tbin2,V(:,end+1)] = running_average(t(:),sig_oi(:),0.15,0.15);
	arrayfun(@(ii) fill(ax(2),[tbin2 flip(tbin2)],[M(:,ii)+V(:,ii); flip(M(:,ii)-V(:,ii))],Colors(ii,:),'FaceAlpha',0.1,'EdgeColor','none'), 1:2);
	arrayfun(@(ii) plot(ax(2),tbin2,M(:,ii),'Color',Colors(ii,:),'LineWidth',2), 1:2);
	% mark laser
	y = 0.83;
	h = plot(ax(1),[0 3],[1 1]*y,'LineWidth',2.5,'Color',Colors(3,:));
	h = plot(ax(2),[0 3],[1 1]*0.3,'LineWidth',2.5,'Color',Colors(3,:));

	
	% ax = np; plot(tbin2,R);ef;

	% figure setting
	arrayfun(@(h) set(h,'XLim',[tbin(1) 3.5]),ax);
	set(ax(1),'XTick',[],'YLim',[0 y],'YTick',[0 0.8]);
	set(ax(2),'XTick',[],'YLim',[0 0.3],'YTick',[0 0.3]);
	ylabel(ax(1),'Example subject');
	ylabel(ax(2),{'Average','across subjects'});
	xlabel(ax(2),'Time (s)');
	set(ax(2),'XTick',[0 1 2 3],'XTickLabel',{'onset','1','2','offset'});
	title(ax(1),{'Average decoder probability','time locked to CGRP stim'},'FontWeight','Bold');
	set(gcf,'Position',[0 0 300/2 300]);

	ef;


% fig 4i: 

	
	clear ax;
	% example subject reactivation
		load('figures/0.15.mat');
		% actual plot
		close all;
		ax(1) = subplot(2,2,1,'NextPlot','add','FontSize',11);
		% move by 2/bin (if we're still using 1 s window for the binning)
		tbin = tbin + 0.5; 
		% 280
		M = squeeze(resp(2,3,:,:));
		V = squeeze(resp_err(2,3,:,:));
		arrayfun(@(ii) fill(ax,[tbin flip(tbin)],[M(:,ii)+V(:,ii); flip(M(:,ii)-V(:,ii))],Colors(ii,:),'FaceAlpha',0.1,'EdgeColor','none'), 1:2);
		arrayfun(@(ii) plot(ax,tbin,M(:,ii),'Color',Colors(ii,:),'LineWidth',2), 1:2);
		% mark laser
		y = 0.83;
		h = plot(ax,[0 3],[1 1]*y,'LineWidth',2.5,'Color',Colors(3,:));


		% figure setting
		set(ax(1),'XLim',[tbin(1) 4],'XTick',[],'YLim',[0 y],'YTick',[0 0.8]);
		ylabel(ax(1),'decoder probability');
		xlabel(ax(1),'Time (s)');
		set(ax,'XLim',tbin([1 end]),'XTick',[0 1 2 3],'XTickLabel',{'onset','1','2','offset'});
		title(ax,{'Average decoder probability','time locked to CGRP stim'},'FontWeight','Normal');

	% time course of peaks
		load('mat/decoders_peaks_bin60_step30.mat');
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
		% ax(2) = subplot(2,2,[3 4],'NextPlot','add','FontSize',14);
		ax(2) = np;
		% calculate time
		t = ((1:size(npeaks_all,1))-1) * to_save(1).step_size / 60;

		% plot peak number
		arrayfun(@(ii) plot_multiple_lines(squeeze(npeaks_all(:,:,ii))',ax(2),'x',t,'base_color',Colors(ii,:)), 1:2);
		% time period
		color_period = [0 0 0; 0.5 0.5 0.5; Colors(3,:,:)];
		arrayfun(@(ii) plot(ax(2),l_period([1:2]+ii-1) + [1 -1]*0.5,[1 1]*(h_top+0.5),'Color',color_period(ii,:),'LineWidth',3), 1:3);
		label_period = {'Drinking','Delay','CGRP stim'};
		h = arrayfun(@(ii) text(ax(2),mean(l_period([1:2]+ii-1))+0.5,h_top+0.6,label_period{ii},'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',ax(2).FontSize-1),1:3);


		% figure setting
		set(ax(2),'XLim',t([1 end]),'YLim',[0 (h_top+1.5)],'XTick',0:15:90);
		xlabel('Time (min)'); ylabel('# of (re)-activations');
		set(gcf,'Position',[0 0 400 300]);
		% h1 = annotation('line',ax_l + [0 ax_width/16],[0.19 0.19]);

		% title(ax(2),'Number of reactivations across different periods','FontWeight','Normal');

	% reactivation comparison
		% take first 10 min
		n_sub = numel(to_save);
		toi = 10*60;
		npeaks_in_toi = nan(n_sub,2);
		for ii = 1:n_sub
			npeaks_in_toi(ii,:) = cellfun(@(x) sum(x < toi), to_save(ii).t_peak);
		end

		% ax(3) = subplot(2,2,2,'NextPlot','add','FontSize',11);
		ax = np;
		scatter(ax,npeaks_in_toi(:,1),npeaks_in_toi(:,2),[],Colors(3,:),'filled');
		% xlabel(ax,'# of novel reactivations');
		% ylabel(ax,'# of water reactivations');
		xlabel(ax,sprintf('\\color[rgb]{%f,%f,%f}Novel',Colors(1,:)))
		ylabel(ax,sprintf('\\color[rgb]{%f,%f,%f}Water',Colors(2,:)))
		% title(ax,{'Number of reactivations','within first 10 min of stimulation'},'FontWeight','Normal');
		title(ax,'reactivation #')
		equalize_plot(ax);
		set(gca,'XTick',[0 200],'YTick',[0 200]);
		ax.Position = [0.3 0.4 0.6 0.5];
		set(gcf,'Position',[0 0 125 100]);
		ef;



	set(gcf,'Position',[0 0 556 556]);
	ax_b = 0.62;
	ax_h = 0.24;
	ax_w = 0.3;
	ax(1).Position = [0.1300    ax_b    ax_w    ax_h];
	ax(2).Position = [0.1300    0.1100    0.7750    0.3412];
	ax(3).Position = [0.5703    ax_b    ax_w    ax_h];
	ef;

	return


% averaged laser posterior
	% load('figures/laser_posterior.mat');
	load('figures/0.15.mat');
	% resp - 6 sessions * 3 time periods * 201 time points * 3 conditions
	% tbin - time stamps for 201 time points
	% ax = np;
	% plot(tbin,squeeze(mean(resp(:,1,:,1),1)));ef;

	% actual plot
	% ax = np(1,2);
	ax = subplot(2,2,1,'NextPlot','add','FontSize',14);
	% 280
	arrayfun(@(ii) plot(ax(1),tbin,squeeze(resp(2,3,:,ii)),'Color',Colors(ii,:)), 1:2);
	% average
	% arrayfun(@(ii) plot_multiple_lines(squeeze(resp(:,3,:,ii)),ax(2),'x',tbin,'base_color',Colors(ii,:)), 1:2);

	% plot laser


	% figure setting
	set(ax(1),'XLim',tbin([1 end]),'XTick',[],'YLim',[0 0.8],'YTick',[0 0.8]);
	ylabel(ax(1),'decoder reactivation');
	xlabel(ax(1),'Time (s)');
	set(ax(2),'XLim',tbin([1 end]),'XTick',[0 1 2 3],'XTickLabel',{'onset','1','2','offset'});


	% ax(1).Position = [0.1485    0.5    0.7565    0.2630];
	% ax(2).Position = [0.1485    0.1882    0.7565    0.2630];
	% set(gcf,'Position',[0 0 150 250]);

	% rmanova
	load fisheriris;
	t = table(species,meas(:,1),meas(:,2),meas(:,3),meas(:,4),...
			'VariableNames',{'species','meas1','meas2','meas3','meas4'});
	Meas = table([1 2 3 4]','VariableNames',{'Measurements'});
	rm = fitrm(t,'meas1-meas4~species','WithinDesign',Meas);
	ranovatbl = ranova(rm)

	% how about adding bars for first 100 stim posterior
		[nsub, nstim] = size(raw_data);
		mean_reactivation = nan(nsub,nstim,3);
		for isub = 1:nsub
			for istim = 1:nstim
				ind = (raw_data(isub,istim).t>0 & raw_data(isub,istim).t<3);
				mean_reactivation(isub,istim,:) = mean(raw_data(isub,istim).sig_oi(:,ind),2);
			end
		end

		% some inspection plot
		[nsub, nstim] = size(raw_data);
		ax = np(nsub);
		arrayfun(@(ii) histogram(ax(ii),squeeze(mean_reactivation(ii,:,1))), 1:nsub);
		arrayfun(@(ii) histogram(ax(ii),squeeze(mean_reactivation(ii,:,2))), 1:nsub);

		ax = np;
		histogram(ax,reshape(squeeze(mean_reactivation(:,:,1)),[],1));
		histogram(ax,reshape(squeeze(mean_reactivation(:,:,2)),[],1));

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
	% load('figures/decoder_peaks.mat');
	load('mat/decoders_peaks_bin60_step30.mat');

	% raw count for inspection
		tbin = 0:0.5:30;
		n_sub = numel(to_save);
		[ax,r,c] = np(n_sub);

		% individual time courses by spike count
			% for ii = 1:n_sub
			% 	tmp = histcounts(to_save(ii).t_peak{1}/60,tbin);
			% 	tmp(end+1,:) = histcounts(to_save(ii).t_peak{2}/60,tbin);
			% 	plot(ax(ii),tbin(1:end-1),tmp');
			% end
			% ef;

		% take first 5 min?
		n_sub = numel(to_save);
		toi = 10*60;
		npeaks_in_toi = nan(n_sub,2);
		for ii = 1:n_sub
			npeaks_in_toi(ii,:) = cellfun(@(x) sum(x < toi), to_save(ii).t_peak);
		end

		fprintf('sign rank test: %.3f\n',signrank(npeaks_in_toi(:,1),npeaks_in_toi(:,2)));
		fprintf('t-test: %d\n',ttest(npeaks_in_toi(:,1),npeaks_in_toi(:,2)));
		
		ax = np;
		% arrayfun(@(jj) scatter(ax,ones(1,6)*jj,npeaks_in_toi(:,jj),[],Colors(jj,:),'filled'), 1:2);
		% plot(ax,1:2,npeaks_in_toi','Color',[0 0 0 0.5],'LineWidth',1);
		scatter(ax,npeaks_in_toi(:,1),npeaks_in_toi(:,2));
		equalize_plot(ax);
		ef;



	% bar plot
		% count peaks
		npeaks_sum = nan(numel(to_save),size(to_save(1).npeaks,1),size(to_save(1).npeaks,2)); % 6 sessions * 2 novel/fam * 3 time periods 
		for ii = 1:numel(to_save)
			npeaks_sum(ii,:,:) = cellfun(@(x) sum(x) / numel(x) / (60/to_save(1).step_size), to_save(ii).npeaks);
		end
		% change CGRP period to first 20 time points
		n = 10;
		for ii = 1:numel(to_save)
			npeaks_sum(ii,:,3) = cellfun(@(x) sum(x(1:n)) / n / (60/to_save(1).step_size), to_save(ii).npeaks(:,3));
			% npeaks_sum(ii,:,3)
		end

		
		% plot 
		ind_sub = 1:6;
		ax = np(1,3);
		for ii = 1:3
			arrayfun(@(jj) scatter(ax(ii),ones(1,numel(ind_sub))*jj,squeeze(npeaks_sum(ind_sub,jj,ii)),[],Colors(jj,:),'filled'), 1:2);
			plot(ax(ii),1:2,squeeze(npeaks_sum(ind_sub,:,ii))','Color',[0 0 0 0.5],'LineWidth',1);
		end

		for ii = 1:3
			[h(ii),p_ttest(ii)] = ttest(squeeze(npeaks_sum(ind_sub,1,ii)),squeeze(npeaks_sum(ind_sub,2,ii)));
			p_signrank(ii) = signrank(squeeze(npeaks_sum(ind_sub,1,ii)),squeeze(npeaks_sum(ind_sub,2,ii)));
		end
		p_ttest
		p_signrank



	% line plot
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
		ax = subplot(2,2,[3 4],'NextPlot','add','FontSize',14);
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
		set(gcf,'Position',[0 0 556 556]);
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

