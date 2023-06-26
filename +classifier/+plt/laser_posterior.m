function ops = laser_posterior(data,ops,prefix)

% initiation
if nargin < 3
	prefix = '';
end

bin_width = 0.15;
toi       = -1:0.025:4;

% count spikes with larger temporal precision
ops.tp = getOr(ops,'tp',[0 1]);
% ops.posterior_t_edges = data.video(1):0.5:data.video(end);
ops.posterior_t_edges = data.video(1):0.1:data.video(end);
[spk_count,ops] = classifier.count_spk.time_course(data,ops);


% calculate posterior
if ~isfield(ops,'Mdl')
	[~,ops] = ops.classifier.train(data,ops);
end
Posterior = ops.classifier.predict(ops.Mdl,spk_count',ops);

% time locked to laser onset
if isstruct(data.laser)
	% in development
	error('need to add onsets / offsets conversion');
else
	% onset and offset
	% [resp,resp_err,tbin] = psth_time_series(Posterior,data.laser(:,1),ops.posterior_t,'bin_width',bin_width,'toi',toi);
	[resp(1,:,:),resp_err(1,:,:),tbin] = psth_time_series(Posterior,data.laser(1:100,1),ops.posterior_t,'bin_width',bin_width,'toi',toi);
	[resp(2,:,:),resp_err(2,:,:),tbin] = psth_time_series(Posterior,data.laser(end-99:end,1),ops.posterior_t,'bin_width',bin_width,'toi',toi);
	[resp(3,:,:),resp_err(3,:,:),tbin] = psth_time_series(Posterior,data.laser(:,1),ops.posterior_t,'bin_width',bin_width,'toi',toi);
end

Colors = getOr(data,'port_color',[1 0 0; 0 0 0]);
% ax = np;
% for ii=1:2
% 	h(ii).shade = fill(ax,[tbin flip(tbin)],[squeeze(resp(:,ii))+squeeze(resp_err(:,ii)); flip(squeeze(resp(:,ii))-squeeze(resp_err(:,ii)))]',Colors(ii,:),'FaceAlpha',0.1,'EdgeColor','none');
% 	h(ii).h_m   = plot(ax,tbin,squeeze(resp(:,ii)),'Color',Colors(ii,:),'LineWidth',1.5);
% end
ax = np(3,1);
for iplot = 1:size(resp,1)
	for ii=1:2
		h(ii).shade = fill(ax(iplot),[tbin flip(tbin)],[squeeze(resp(iplot,:,ii))+squeeze(resp_err(iplot,:,ii)), flip(squeeze(resp(iplot,:,ii))-squeeze(resp_err(iplot,:,ii)))]',Colors(ii,:),'FaceAlpha',0.1,'EdgeColor','none');
		h(ii).h_m   = plot(ax(iplot),tbin,squeeze(resp(iplot,:,ii)),'Color',Colors(ii,:),'LineWidth',1.5);
	end
end
align_ax(ax,true,true);

set(ax(3),'XLim',tbin([1 end]),'XTick',[0 1.5 3]);
ax(3).XTickLabel(ax(3).XTick==0) = {'onset'}; 
ax(3).XTickLabel(ax(3).XTick==3) = {'offset'};


tmp = {'novel','fam'}; legend(ax(1),[h(1).h_m h(2).h_m],tmp(data.port_is_water+1),'box','off');
xlabel(ax(3),'time (s)'); ylabel(ax(3),'average posteior');


title(ax(1),{['calca\_' data.subject],sprintf('step: %d ms',bin_width*1e3),'first 100 lasers'});
title(ax(2),'last 100 lasers');
title(ax(3),'all lasers');
set(gcf,'Position',[0 0 150 400]);
export_fig(sprintf('results/%slaserlock_%s_%s.pdf',prefix,data.subject,datestr(data.session,'YYmmdd')));

