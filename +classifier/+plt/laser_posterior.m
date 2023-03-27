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
	% in develop
	error('need to add onsets / offsets conversion');
else
	
	% onset and offset
	[resp,resp_err,tbin] = psth_time_series(Posterior,data.laser(:,1),ops.posterior_t,'bin_width',bin_width,'toi',toi);

end

ax = np;
Colors = getOr(data,'port_color',[1 0 0; 0 0 0]);
for ii=1:2
	h(ii).shade = fill(ax,[tbin flip(tbin)],[squeeze(resp(:,ii))+squeeze(resp_err(:,ii)); flip(squeeze(resp(:,ii))-squeeze(resp_err(:,ii)))]',Colors(ii,:),'FaceAlpha',0.1,'EdgeColor','none');
	h(ii).h_m   = plot(ax,tbin,squeeze(resp(:,ii)),'Color',Colors(ii,:),'LineWidth',1.5);
end


set(ax,'XLim',tbin([1 end]),'XTick',[0 1.5 3]);
ax.XTickLabel(ax.XTick==0) = {'onset'}; 
ax.XTickLabel(ax.XTick==3) = {'offset'};


tmp = {'novel','fam'}; legend(ax,[h(1).h_m h(2).h_m],tmp(data.port_is_water+1),'box','off');
xlabel(ax(1),'time (s)'); ylabel(ax(1),'average posteior');
title(ax(1),{['calca\_' data.subject],sprintf('step: %d ms',bin_width*1e3)});
set(gcf,'Position',[0 0 200 150]);
export_fig(sprintf('results/%slaserlock_%s_%s.pdf',prefix,data.subject,datestr(data.session,'YYmmdd')));

