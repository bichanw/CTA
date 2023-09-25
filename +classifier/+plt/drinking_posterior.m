function [ops,resp,resp_err,tbin,raw_data] = drinking_posterior(data,ops,prefix)

% initiation
if nargin < 3
	prefix = '';
end

bin_width = 0.15;
toi       = -5.5:0.05:10;

% count spikes with larger temporal precision
ops.tp = getOr(ops,'tp',[0 1]);
% ops.posterior_t_edges = data.video(1):0.5:data.video(end);
ops.posterior_t_edges = getOr(ops,'posterior_t_edges',data.video(1):0.1:data.video(end));


% calculate posterior
if ~isfield(ops,'Mdl')
	[~,ops] = ops.classifier.train(data,ops);
end
[spk_count,ops] = classifier.count_spk.time_course(data,ops); % count spike after training classifier since the training might remove more cells
Posterior = ops.classifier.predict(ops.Mdl,spk_count',ops);

% time locked to laser onset
% onset and offset (plot all 3)
clear resp resp_err;
[resp(1,:,:),resp_err(1,:,:),tbin,raw_data] = psth_time_series(Posterior,data.rewards.all.front,ops.posterior_t,'bin_width',bin_width,'toi',toi);
[resp(2,:,:),resp_err(2,:,:),tbin,raw_data(end+1,:)] = psth_time_series(Posterior,data.rewards.all.rear,ops.posterior_t,'bin_width',bin_width,'toi',toi);
% [resp(3,:,:),resp_err(3,:,:),tbin] = psth_time_series(Posterior,data.laser(:,1),ops.posterior_t,'bin_width',bin_width,'toi',toi);

% flip if front is water (so first is novel)
if data.port_is_water(1)
	resp = resp(data.port_is_water+1,:,[data.port_is_water+1 3]);
	resp_err = resp_err(data.port_is_water+1,:,[data.port_is_water+1 3]);
	raw_data = raw_data(data.port_is_water+1,:);
	for ii = 1:size(raw_data,1)
		for jj = 1:size(raw_data,2)
			raw_data(ii,jj).sig_oi = raw_data(ii,jj).sig_oi([data.port_is_water+1 3],:);
		end
	end
end


% Colors = [228 45 38; 55 135 192] / 255;
% close all; clear ax;
% ax(1) = subplot(1,2,1,'NextPlot','add','FontSize',11);
% ax(2) = subplot(1,2,2,'NextPlot','add','FontSize',11);
% % move by 2/bin (if we're still using 1 s window for the binning)
% % 280
% for jj = 1:2
% 	for ii = 1:2
% 		M = squeeze(resp(jj,:,ii))';
% 		V = squeeze(resp_err(jj,:,ii))';
% 		fill(ax(jj),[tbin flip(tbin)]+0.5,[M+V; flip(M-V)],Colors(ii,:),'FaceAlpha',0.1,'EdgeColor','none');
% 		plot(ax(jj),tbin+0.5,M,'Color',Colors(ii,:),'LineWidth',2);
% 	end
% end

% export_fig('tmp.pdf');

