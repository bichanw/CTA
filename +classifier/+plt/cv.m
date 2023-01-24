function cv(data,ops,prefix)
% cross validation to justify result

% initiation
if nargin < 3
	prefix = '';
end

% prepare data
[X,Y,ops] = classifier.data_2_XY(data,ops);
n_repeat = 10;
n_sample = numel(Y);

% value of interest
predicted = nan(n_repeat,numel(Y));
cost  = nan(n_repeat,1);
valid = nan(n_repeat,1);

% run cross validation
for ii = 1:n_repeat
	[cost(ii),cv] = classifier.CV_bw(X,Y,10,getOr(ops,'classifier',classifier.nb()));
	predicted(ii,:) = cv.post_label;
	valid(ii) = sum(~isnan(cv.post_label)) / n_sample;
end
predicted(isnan(predicted)) = 0;


% plot label
ax = np(1,2);
imagesc(ax(1),1,1,Y');
imagesc(ax(1),1,3,predicted);
xlabel(ax(1),'trial'); 
set(ax(1),'YTick',[1 2+n_repeat/2],'YTickLabel',{'true','predicted'});


% plot cv results
% always valid
if prod(valid==n_sample)
	histogram(ax(2),cost);
	xlabel('misclassification rate');
% invalid cv
else
	histogram2(ax(2),valid,cost,'DisplayStyle','tile','ShowEmptyBins','on','EdgeColor','none','XBinEdges',linspace(0,1,ceil(n_sample/10)),'YBinEdges',linspace(0,1,n_sample));
	% histogram2(ax(2),valid,cost,'DisplayStyle','tile','ShowEmptyBins','on','EdgeColor','none');
	xlabel(ax(2),'% valid');
	ylabel(ax(2),{'misclassification','rate'});
end
cmap = colormap('gray');
colormap(ax(2),flip(cmap));
c = colorbar(ax(2));


% figure setting
if any(predicted==0)
	colormap(ax(1),[1 1 1;1 0 0; 0 0 0; 0.8 0.8 0.3]); % need to be after calling colormap('gray')
else
	colormap(ax(1),[1 0 0; 0 0 0; 0.8 0.8 0.3]); % need to be after calling colormap('gray')
end
export_fig(sprintf('results/%scv_%s_%s.png',prefix,data.subject,datestr(data.session,'YYmmdd')),'-m3');

end