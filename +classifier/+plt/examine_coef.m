function examine_coef(data,ops,prefix)
if nargin < 3
	prefix = '';
end

if ~isfield(ops,'Mdl')
	[~,ops] = ops.classifier.train(data,ops);
end

% initiation
pairs = [1 2;1 3;2 3];
n = size(pairs,1);
Colors = [1 0 0;0 0 0;1 1 0.2];

% color scatter plots based on ranksum categories
if isfield(ops,'novel_vs_fam')
	ind = {ismember(find(~ops.exclude_id),ops.novel_vs_fam.ordered_id((ops.novel_vs_fam.ordered_div(1)+1):ops.novel_vs_fam.ordered_div(2))),...
		   ismember(find(~ops.exclude_id),ops.novel_vs_fam.ordered_id((ops.novel_vs_fam.ordered_div(2)+1):ops.novel_vs_fam.ordered_div(3))),...
		   ~ismember(find(~ops.exclude_id),ops.novel_vs_fam.ordered_id((ops.novel_vs_fam.ordered_div(1)+1):ops.novel_vs_fam.ordered_div(3)))};
else
	ind = {1:size(ops.Mdl.coef,1)};
end

% plotting
% scatter
[ax,r,c] = np(2,n);
for ii = 1:n
	for jj = 1:numel(ind)
		% skip if not data
		if sum(ind{jj})>0
			scatter(ax(sub2ind([c r],ii,1)),ops.Mdl.coef(ind{jj},pairs(ii,1)),ops.Mdl.coef(ind{jj},pairs(ii,2)),[],Colors(jj,:));
		end
	end
	% figure setting
	xlabel(ax(sub2ind([c r],ii,1)),ops.events{pairs(ii,1)});
	ylabel(ax(sub2ind([c r],ii,1)),ops.events{pairs(ii,2)});
end
% histogram
m = min(ops.Mdl.coef(:));
M = max(ops.Mdl.coef(:));
for ii = 1:n
	histogram2(ax(sub2ind([c r],ii,2)),ops.Mdl.coef(:,pairs(ii,1)),ops.Mdl.coef(:,pairs(ii,2)),...
				'XBinEdges',linspace(m,M,25),'YBinEdges',linspace(m,M,25),...
				'EdgeColor','none','DisplayStyle','tile','ShowEmptyBins','on');
	
	% label
	xlabel(ax(sub2ind([c r],ii,2)),ops.events{pairs(ii,1)});
	ylabel(ax(sub2ind([c r],ii,2)),ops.events{pairs(ii,2)});
	
	% other axes settings
	ax(sub2ind([c r],ii,2)).ColorScale = 'log';
	axis(ax(sub2ind([c r],ii,2)),'square');
	colorbar(ax(sub2ind([c r],ii,2)));
end
colormap(cbrewer2('Reds'));


% figure settings
title(ax(sub2ind([c r],2,1)),'decoder coefficient');
l = legend(ax(sub2ind([c r],3,1)),{'front','back','other'},'box','off'); title(l,'by ranksum');

% set xlim ylim
arrayfun(@(h) set(h,'XLim',[m M],'YLim',[m M]),ax);
% arrayfun(@(h) set(h,'XLim',[min([h.XLim(1) h.YLim(1)]) max([h.XLim(2) h.YLim(2)])],'YLim',[min([h.XLim(1) h.YLim(1)]) max([h.XLim(2) h.YLim(2)])]),ax);

set(gcf,'Position',[0 0 505 350]);
export_fig(sprintf('results/%scoef_%s_%s.pdf',prefix,data.subject,datestr(data.session,'YYmmdd')));


end