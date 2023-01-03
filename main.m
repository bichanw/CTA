% ----- individual cell plotting ----- %
% parameters
plt_name = 'FR_amp';
extension = 'png';


sessions; % initialize sessions to run
for iSession = 1:numel(Sessions.subject)

	% load data
	data = load_data.all(Sessions.session(iSession),Sessions.subject{iSession});

	% plotting initiation
	cells_2_plt = 1:numel(data.spikes); 
	N = numel(cells_2_plt);
	N_in_batch = 100;
	for ibatch = 1:ceil(N/N_in_batch)
		[ax,r,c] = np(N_in_batch);

		% plot cell by cell
		for ii = 1:min([N_in_batch N-(ibatch-1)*N_in_batch])
			eval(['plt.' plt_name '(data,ii+(ibatch-1)*N_in_batch,ax(ii));'])
		end

		% figure setting
		eval(['set_fig.' plt_name '(ax(sub2ind([c r],1,r)),data);'])
		
		% save
		export_fig(sprintf('results/%s_%s_%s_%d.%s',plt_name,datestr(data.session,'YYmmdd'),data.subject,ibatch,extension));

	end
	% combine plots
	if strcmp(extension,'pdf') append_script(sprintf('results/%s_%s_%s',plt_name,datestr(data.session,'YYmmdd'),data.subject)); end


	% % each cell as a pdf?
	% cells_2_plt = 1:numel(data.spikes); 
	% for ii = cells_2_plt
	% 	ax = np;
	% 	plt.FR_amp(data,ii,ax);
	% 	set_fig.FR_amp(data,ax);
	% 	return
	% end


end


return


% classifier pipeline
data = load_data.all(datetime(2022,12,4),'001');

% classifier setting
ops.exclude_id = false(size(data.spikes)); % include all cells


[spk_count,ops] = classifier.count_spk(data,struct());

[Mdl, ops] = classifier.nb.train(spk_count,ops);


% count spikes for all
toi = data.video(1):0.5:data.video(end);
count = NaN(sum(~ops.exclude_id),numel(toi)); t = [];

cell_oi = find(~ops(1).exclude_id);
for ii = 1:numel(cell_oi)
	[~,count(ii,:),~,t] = running_average(data.spikes{cell_oi(ii)},[],diff(ops.tp),[],toi); 
end


% basic plotting

[label,Posterior,Cost] = predict(Mdl,count');
% Posterior = classifier.posterior_bw(Mdl{j},count');

t_step  = 200;
t_start = min(t):t_step:max(t);

ax = np(numel(t_start),1);
set(gcf,'Position',[0 0 1000 150*numel(t_start)]);

clear s;
mk_sz = 12;
for i = 1:numel(ax)
	% line for posterior probability
	h = plot(ax(i),t,Posterior); h(1).Color = [1 0 0]; h(2).Color = [0 0 0 ]; h(3).Color(4) = 0.5;h(3).LineWidth = 0.7;
	% events
	classifier.plt.scatter_event(data,ax(i));
	
end	
arrayfun(@(i) set(ax(i),'XLim',[0 t_step]+t_start(i),'YLim',[0 1.5],'YTick',[0 0.5 1]), 1:numel(ax))

export_fig(sprintf('calca%s_%s.pdf',data.subject,datestr(data.session,'yymmdd')));


return


% examine if amplitude change
	% first 5 min vs last 5 min?
	t_2_cmp = 5 * 60; % 300 s

	% number of spikes available?
	N = numel(data.spikes);
	n_spk = [arrayfun(@(ii) sum(data.spikes{ii}<data.video(1)+t_2_cmp), 1:N)',...
			 arrayfun(@(ii) sum(data.spikes{ii}>=data.video(end)-t_2_cmp), 1:N)'];

	% send report
	fprintf('Session %s_%s\n',datestr(data.session,'YYmmdd'),data.subject);
	fprintf('Start video frame %.1f s, end video frame %.1f s, time of interest %d sec.\n',data.video(1),data.video(end),t_2_cmp);
	arrayfun(@(ii) fprintf('Cell %d available spikes: %d, %d\n',ii,n_spk(ii,1),n_spk(ii,2)), 1:numel(data.spikes));





