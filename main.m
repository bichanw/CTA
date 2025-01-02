clearvars -except sessions_oi; clc; close all; addpath('helpfun');



sessions; % initialize sessions to run
for iSession = sessions_oi
% for iSession = 2


	% load data
	data = load_data.all(Sessions.session(iSession),Sessions.subject{iSession});
	
	% parameter section
	params;


	% ----- session level plot ----- %
	for lambda = 1
		ops.mnr.lambda = lambda;
		prefix = sprintf('amp%d_%dcat_%s_%s_%.1e_',ops.amplitude_cutoff,numel(ops.amplitude_cat),ops.zscore_method,ops.mnr.penalty,ops.mnr.lambda);
		% prefix = 'testremove_';

		% classifier.plt.posterior_raster(data,[],ops,prefix); % plot posterior
		[ops,to_save(iSession)] = classifier.plt.slow_firing(data,ops,prefix);
		% ops.posterior_t_edges = data.video(1):0.15:data.video(end); [ops,resp(iSession,:,:,:),resp_err(iSession,:,:,:),tbin,raw_data(iSession,:)] = classifier.plt.laser_posterior(data,ops,prefix);
		% ops.posterior_t_edges = data.video(1):0.15:data.video(end); [ops,resp(iSession,:,:,:),resp_err(iSession,:,:,:),tbin,raw_data(iSession,:,:)] = classifier.plt.drinking_posterior(data,ops,prefix);
		% classifier.plt.examine_coef(data,ops,prefix);
		% classifier.plt.brain_region(data,ops,prefix); 
		% return
	end
	% return
	fprintf('%d session, front is %s\n',iSession,data.ports.front);
	% to_save(iSession).cv_results = classifier.plt.cv(data,ops,prefix);
	% return
end

% save processing setting
saveops(ops);
% return

if exist('to_save')
	% save(sprintf('mat/decoders_peaks_bin%d_step%d.mat',to_save(1).bin_width,to_save(1).step_size),'to_save');
	save(sprintf('mat/%s.mat',datestr(now,'YYmmdd')),'to_save');
end
if exist('resp')
	save(sprintf('mat/%s.mat',datestr(now,'YYmmdd')),'resp','resp_err','tbin','raw_data');
	% save('figures/0.15.mat','resp','resp_err','tbin','raw_data');
end

% append_script(['results/' prefix],true);
% append_script(['results/novel_vs_fam_' prefix],true);
% fprintf('Finish running %s\n',prefix);
return



