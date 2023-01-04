function ops = dprime(data,ops)
	% count spikes
	[count,ops] = classifier.count_spk.events(data,ops);

	% calculate dp
	dp_front = my_dp(count{3}',count{1}');
	dp_rear = my_dp(count{3}',count{2}');

	% plot
	ax = np(1,2); 
	% raw dp
	scatter(ax(1),dp_front,dp_rear);
	M = max([ax(1).XLim(2) ax(1).YLim(2)]); plot(ax(1),[0 M],[0 M],'k--');
	xlabel(ax(1),'dp front');ylabel(ax(1),'dp rear');
	% absolute dp
	scatter(ax(2),abs(dp_front),abs(dp_rear));
	M = max([ax(2).XLim(2) ax(2).YLim(2)]); plot(ax(2),[0 M],[0 M],'k--');
	xlabel(ax(2),'abs dp front');ylabel(ax(2),'abs dp rear');

	% export figure
	export_fig(sprintf('results/dprime1s_%s_%s.pdf',data.subject,datestr(data.session,'YYmmdd')));

end