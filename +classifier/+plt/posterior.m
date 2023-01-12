function ops = posterior(data,Mdl,ops)
% plot postieror traces

% set classifier
Classifier = getOr(ops,'classifier',classifier.nb());

% retrain a model if there's no input
if isempty(Mdl)
	[Mdl, ops] = Classifier.train(data,ops);
end

% count spikes time course
[count,ops] = classifier.count_spk.time_course(data,ops);

% calculate posterior
Posterior = Classifier.predict(Mdl,count');


% plot initiation
t_step  = 200;
t_start = min(ops.posterior_t):t_step:max(ops.posterior_t);
max_ax = 10;
for ibatch = 1:ceil(numel(t_start)/max_ax)
	ax = np(max_ax,1);
	set(gcf,'Position',[0 0 1000 1500]);

	for i = 1:numel(ax)
		% line for posterior probability
		h = plot(ax(i),ops.posterior_t,Posterior); h(2-data.port_is_water(2)).Color = [1 0 0]; h(2-data.port_is_water(1)).Color = [0 0 0]; h(3).Color(4) = 0.5;h(3).LineWidth = 0.7;
		% events
		classifier.plt.scatter_event(data,ax(i),1.25);
	end

	% change x axis to progress in time
	arrayfun(@(i) set(ax(i),'XLim',[0 t_step]+t_start(i)+(ibatch-1)*t_step*max_ax,'YLim',[0 1.5],'YTick',[0 0.5 1]), 1:numel(ax))
	% figure setting
	legend(ax(1),{'front','back'});

	% save figure
	export_fig(sprintf('results/posterior_%s_%s_%d.pdf',data.subject,datestr(data.session,'yymmdd'),ibatch));

end
append_script(sprintf('results/posterior_%s_%s',data.subject,datestr(data.session,'yymmdd')));


end