function ops = reactivation(data,Mdl,ops,prefix)
% plot postieror traces

% initiation
if nargin < 4
	prefix = ''; % figure prefix
end

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
% average posterior
[Post_avg,~,~,t] = running_average(ops.posterior_t,Posterior,60,30);


% plotting
ax = np;
plot(ax,t,Post_avg(:,1),'Color',[1 0 0]);
plot(ax,t,Post_avg(:,2),'Color',[0 0 0]);
classifier.plt.divider_event(data,ax);

set(ax,'XLim',t([1 end]),'XTick',round(t([1 end])));
xlabel('time (s)');

export_fig(sprintf('results/%savg_post_%s_%s.pdf',prefix,data.subject,datestr(data.session,'YYmmdd')));

end