function cv(data,ops,prefix)
% cross validation to justify result

% initiation
if nargin < 3
	prefix = '';
end

% prepare data
[X,Y,ops] = classifier.data_2_XY(data,ops);
n_repeat = 10;
n_fold   = 10;
n_sample = numel(Y);
Lambdas = logspace(-4,4,9);

% value of interest
predicted = nan(n_repeat,numel(Y));
valid = nan(n_repeat,1);
LL = nan(numel(Lambdas),n_fold,2); % LL for train / test

% run cross validation
for jj = 1:numel(Lambdas)

	% set lambda
	ops.mnr.lambda = Lambdas(jj);

	% LL for 10 fold
	[~,cv] = classifier.CV_bw(X,Y,n_fold,ops);
	LL(jj,:,:) = cv.LL;


	% for ii = 1:n_repeat
		% [cost(ii),cv] = classifier.CV_bw(X,Y,10,ops);
		% predicted(ii,:) = cv.post_label;
		% valid(ii) = sum(~isnan(cv.post_label)) / n_sample;
	% end
end




% ----plot LL
ax = np(1,2);
Colors = cbrewer2('Accent',2);
for jj = 1:numel(Lambdas)
	my_scatter(Lambdas(jj),squeeze(LL(jj,:,1)),ax(1),'MarkerEdgeColor',Colors(1,:));
	my_scatter(Lambdas(jj),squeeze(LL(jj,:,2)),ax(2),'MarkerEdgeColor',Colors(2,:));
end

xlabel(ax(2),'\lambda','Interpreter','tex');
ylabel(ax(1),'LL');
title(ax(1),'train'); title(ax(2),'test');
% legend(ax(2),{'train','test'},'Location','eastoutside','box','off');
arrayfun(@(h) set(h,'XScale','log'), ax);
set(gcf,'Position',[0 0 450 200]);



% ----original, plot label
% predicted(isnan(predicted)) = 0;
% % plot label
% ax = np(1,2);
% imagesc(ax(1),1,1,Y');
% imagesc(ax(1),1,3,predicted);
% xlabel(ax(1),'trial'); 
% set(ax(1),'YTick',[1 2+n_repeat/2],'YTickLabel',{'true','predicted'});


% % plot cv results
% % always valid
% if prod(valid==n_sample)
% 	histogram(ax(2),cost);
% 	xlabel('misclassification rate');
% % invalid cv
% else
% 	histogram2(ax(2),valid,cost,'DisplayStyle','tile','ShowEmptyBins','on','EdgeColor','none','XBinEdges',linspace(0,1,ceil(n_sample/10)),'YBinEdges',linspace(0,1,n_sample));
% 	% histogram2(ax(2),valid,cost,'DisplayStyle','tile','ShowEmptyBins','on','EdgeColor','none');
% 	xlabel(ax(2),'% valid');
% 	ylabel(ax(2),{'misclassification','rate'});
% end
% cmap = colormap('gray');
% colormap(ax(2),flip(cmap));
% c = colorbar(ax(2));


% % figure setting
% if any(predicted==0)
% 	colormap(ax(1),[1 1 1;1 0 0; 0 0 0; 0.8 0.8 0.3]); % need to be after calling colormap('gray')
% else
% 	colormap(ax(1),[1 0 0; 0 0 0; 0.8 0.8 0.3]); % need to be after calling colormap('gray')
% end
export_fig(sprintf('results/%scv_%s_%s.pdf',prefix,data.subject,datestr(data.session,'YYmmdd')));

end