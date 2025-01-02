%% paper figure
% general figure setting
addpath('helpfun');
addpath('pkgs/export_fig');
addpath(genpath('pkgs/cbrewer2'));
Colors = [228 45 38;   55 135 192; 54 161 86] / 255; % novel, familiar, CGRP


% new fig extended fig 8: control vs. CGRP ablation
	load('figures/240815.mat');
	Colors = [228 45 38;   55 135 192; 54 161 86] / 255; % novel, familiar, CGRP
	clear tmp
	for ii = 1:2
		for jj = 1:3
			tmp{ii,jj} = catcell(arrayfun(@(kk) to_save(kk).npeaks{ii,jj}, 1:8,'uni',0),2);
		end
	end
	
	npeaks_all = catcell(tmp(1,:),1);
	npeaks_all(:,:,2) = catcell(tmp(2,:),1);
	l_period = [0 cumsum(cellfun(@(x) size(x,1),tmp(1,:)))] * to_save(1).step_size / 60;
	

	% plot 1 trace
	h_top = 17.5; % top of the trace
	% ax(2) = subplot(2,2,[3 4],'NextPlot','add','FontSize',14);
	ax = np(2,1);
	% calculate time
	t = ((1:size(npeaks_all,1))-1) * to_save(1).step_size / 60;

	% plot peak number
	arrayfun(@(ii) plot_multiple_lines(squeeze(npeaks_all(:,1:4,ii))',ax(1),'x',t,'base_color',Colors(ii,:)), 1:2);
	arrayfun(@(ii) plot_multiple_lines(squeeze(npeaks_all(:,5:8,ii))',ax(2),'x',t,'base_color',Colors(ii,:)), 1:2);
	% time period
	color_period = [0 0 0; 0.5 0.5 0.5; Colors(3,:,:)];
	arrayfun(@(ii) plot(ax(1),l_period([1:2]+ii-1) + [1 -1]*0.5,[1 1]*(h_top+0.5),'Color',color_period(ii,:),'LineWidth',3), 1:3);
	label_period = {'Drinking','Delay','LiCl'};
	h = arrayfun(@(ii) text(ax(1),mean(l_period([1:2]+ii-1))+0.5,h_top+0.6,label_period{ii},'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',ax(2).FontSize-1),1:3);


	% figure setting
	set(ax,'XLim',t([1 end]),'YLim',[0 (h_top+1.5)],'XTick',0:15:90);
	xlabel('Time (min)'); ylabel('# of (re)-activations');
	title(ax(1),'control'); title(ax(2),'ablation');
	set(gcf,'Position',[0 0 400 400]);
	export_fig('tmp.pdf');
	% h1 = annotation('line',ax_l + [0 ax_width/16],[0.19 0.19]);


% fig 4g: raster
	load('figures/280_old.mat');

	ops.novel_vs_fam = classifier.select_cells.by_chris(data,ops);
	spk_2_plt = data.spikes(ops.novel_vs_fam.ordered_id);
	ytick_loc = [-2.5 (ops.novel_vs_fam.ordered_div(1:end-1)+ops.novel_vs_fam.ordered_div(2:end))/2];

	h_posterior = 3;
	ops.posterior_t = ops.posterior_t + diff(ops.tp)/2; % make it causal
	Axes = np(1,2);

	% add raster below
	% Tois = [673.5 723.5; 7970 8010]; % CGRP
	Tois = [673.5 753.5; 7970 8040];  % drinking 
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
		arrayfun(@(ii) fill(ax,[ops.posterior_t flip(ops.posterior_t)],[zeros(size(Posterior(:,1))); flip(Posterior(:,ii).* -h_posterior)],Colors(ii,:),'FaceAlpha',0.1,'EdgeColor','none'), 1:2);
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


% extended figure 8d, posterior locked to drinking
	load('figures/231003.mat');
	close all; clear ax;
	ax(1) = subplot(2,2,1,'NextPlot','add','FontSize',8);
	ax(2) = subplot(2,2,2,'NextPlot','add','FontSize',8);
	% move by 2/bin (if we're still using 1 s window for the binning)
	% 280
	for jj = 1:2
		for ii = 1:2
			M = squeeze(resp(2,jj,:,ii));
			V = squeeze(resp_err(2,jj,:,ii));
			% fill(ax(jj),[tbin flip(tbin)],[M+V; flip(M-V)],Colors(ii,:),'FaceAlpha',0.1,'EdgeColor','none');
			fill(ax(jj),[tbin flip(tbin)],[M; zeros(size(M))],Colors(ii,:),'FaceAlpha',0.1,'EdgeColor','none');
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
			t = cat(2,raw_data(:,ievent,:).t);
			sig_oi = catcell(arrayfun(@(s) s.sig_oi(iPost,:),raw_data(:,ievent,:),'uni',0),2);
			[M,~,~,tbin2,V] = running_average(t(:),sig_oi(:),0.15,0.15);
			tbin2 = tbin2 + 0.15/2;
			% fill(ax(ievent+2),[tbin2 flip(tbin2)],[M+V; flip(M-V)],Colors(iPost,:),'FaceAlpha',0.1,'EdgeColor','none');
			fill(ax(ievent+2),[tbin2 flip(tbin2)],[M; zeros(size(M))],Colors(iPost,:),'FaceAlpha',0.1,'EdgeColor','none');
			plot(ax(ievent+2),tbin2,M,'Color',Colors(iPost,:),'LineWidth',0.7);
		end
	end
	ax(3).XAxis.Visible = 'off';
	ax(4).XAxis.Visible = 'off';
	% ax(3).YAxis.Visible = 'off';

	
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



	set(gcf,'Position',[0 0 330 200]); 

	export_fig('tmp.pdf');
	return



% fig 4h: averaged CGRP
	% load('figures/0.15.mat'); % v0: 0.15 s binning
	load('figures/231003_cgrp.mat'); % v1: use right time to label averaing traces as well
	% actual plot
	% tbin = tbin + 0.5; 
	close all; clear ax;
	XLim = [-0.5 3.5];
	[~,I1] = min(abs(tbin+0.5)); [~,I2] = min(abs(tbin-3.5)); ind_XLim = [I1 I2];
	% ind_XLim = [15 175];
	ax(1) = subplot(2,1,1,'NextPlot','add','FontSize',11);
	% move by 2/bin (if we're still using 1 s window for the binning)
	% 280
	M = squeeze(resp(2,3,:,:));
	V = squeeze(resp_err(2,3,:,:));
	arrayfun(@(ii) plot(ax,tbin,M(:,ii),'Color',Colors(ii,:),'LineWidth',1), 1:2);
	arrayfun(@(ii) plot(ax,[XLim; XLim],[0 0; M(ind_XLim(1),ii) M(ind_XLim(2),ii)],'Color',Colors(ii,:),'LineWidth',1), 1:2);
	% add shade
	% arrayfun(@(ii) fill(ax,[tbin flip(tbin)],[M(:,ii)+V(:,ii); flip(M(:,ii)-V(:,ii))],Colors(ii,:),'FaceAlpha',0.1,'EdgeColor','none'), 1:2);
	fill(ax,[tbin flip(tbin)],[M(:,1); flip(M(:,2))],Colors(1,:),'FaceAlpha',0.1,'EdgeColor','none');
	fill(ax,[tbin flip(tbin)],[M(:,2); zeros(size(M(:,2)))],Colors(2,:),'FaceAlpha',0.1,'EdgeColor','none');
	% mark laser
	y = 0.83;
	h = plot(ax,[0 3],[1 1]*y,'LineWidth',2.5,'Color',Colors(3,:));

	% averaged across all subjcts
	ax(2) = subplot(2,1,2,'NextPlot','add','FontSize',11);
	bin_window2 = 0.15;
	t = cat(1,raw_data(:).t);
	sig_oi = catcell(arrayfun(@(ii) raw_data(ii).sig_oi(1,:),1:600,'uni',0));
	[M,~,~,tbin2,V] = running_average(t(:),sig_oi(:),0.15,0.15);
	sig_oi = catcell(arrayfun(@(ii) raw_data(ii).sig_oi(2,:),1:600,'uni',0));
	[M(:,end+1),~,~,tbin2,V(:,end+1)] = running_average(t(:),sig_oi(:),bin_window2,bin_window2);
	% [~,I1] = min(abs(tbin2+0.5)); [~,I2] = min(abs(tbin2-3.5)); ind_XLim = [I1 I2];
	arrayfun(@(ii) plot(ax(2),tbin2+bin_window2/2,M(:,ii),'Color',Colors(ii,:),'LineWidth',1), 1:2);
	% arrayfun(@(ii) plot(ax(2),[XLim; XLim],[0 0; M(ind_XLim(1),ii) M(ind_XLim(2),ii)],'Color',Colors(ii,:),'LineWidth',1), 1:2);
	% add shade, either use error or just fill the bottom section
	% arrayfun(@(ii) fill(ax(2),[tbin2 flip(tbin2)],[M(:,ii)+V(:,ii); flip(M(:,ii)-V(:,ii))],Colors(ii,:),'FaceAlpha',0.1,'EdgeColor','none'), 1:2);
	fill(ax(2),[tbin2 flip(tbin2)]+bin_window2/2,[M(:,1); flip(M(:,2))],Colors(1,:),'FaceAlpha',0.1,'EdgeColor','none');
	fill(ax(2),[tbin2 flip(tbin2)]+bin_window2/2,[M(:,2); zeros(size(M(:,2)))],Colors(2,:),'FaceAlpha',0.1,'EdgeColor','none');
	% mark laser
	y = 0.8;
	h = plot(ax(1),[0 3],[1 1]*y,'LineWidth',2.5,'Color',Colors(3,:));
	h = plot(ax(2),[0 3],[1 1]*0.3,'LineWidth',2.5,'Color',Colors(3,:));

	
	% ax = np; plot(tbin2,R);export_fig('tmp.pdf')

	% figure setting
	ax(1).XAxis.Visible = 'off'; ax(2).XAxis.Visible = 'off'; 
	arrayfun(@(h) set(h,'XLim',XLim),ax);
	set(ax(1),'XTick',[],'YLim',[0 y],'YTick',[0 0.8]);
	set(ax(2),'XTick',[0 1 2 3],'YLim',[0 0.4],'YTick',[0 0.4]);
	ylabel(ax(1),'Example subject');
	ylabel(ax(2),{'Average','across subjects'});
	xlabel(ax(2),'Time (s)');
	set(ax(2),'XTick',[0 1 2 3],'XTickLabel',{'onset','1','2','offset'});
	title(ax(1),{'Average decoder probability','time locked to CGRP stim'},'FontWeight','Bold');
	set(gcf,'Position',[0 0 300/2 300]);

	set(gcf,'Color','white');
	export_fig('tmp.pdf');
	return
	


% extended figure 4c, log likelihood results
	load('figures/230923.mat');
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
	set(gcf,'Position',[0 0 150 100]);
	export_fig('tmp.pdf'); 

% extended figure 8e, confusion matrix
	load('figures/230923.mat');
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
	export_fig('tmp.pdf');







	
% fig 4h, top example subject reactivation
	clear ax;
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

% figure 4i: time course of peaks
	load('figures/decoders_peaks_bin60_step30.mat');  % main figure version when counting peaks
	% load('mat/230928.mat'); % another version when counting the start of 0.5 periods
	% load('mat/240815.mat'); % ablation experiment, not in paper
	Colors = [228 45 38;   55 135 192; 54 161 86] / 255; % novel, familiar, CGRP
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


	
