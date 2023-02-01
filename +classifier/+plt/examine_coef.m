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
Colors = [1 0 0;0 0 0;0.2 1 0.2];

if isfield(ops,'novel_vs_fam')
	ind = {ismember(find(~ops.exclude_id),ops.novel_vs_fam.ordered_id((ops.novel_vs_fam.ordered_div(1)+1):ops.novel_vs_fam.ordered_div(2))),...
		   ismember(find(~ops.exclude_id),ops.novel_vs_fam.ordered_id((ops.novel_vs_fam.ordered_div(2)+1):ops.novel_vs_fam.ordered_div(3))),...
		   ~ismember(find(~ops.exclude_id),ops.novel_vs_fam.ordered_id((ops.novel_vs_fam.ordered_div(1)+1):ops.novel_vs_fam.ordered_div(3)))};
else
	ind = {1:size(ops.Mdl.coef,1)};
end

% plotting
ax = np(n);
for ii = 1:n
	for jj = 1:numel(ind)
		% skip if not data
		if sum(ind{jj})>0
			scatter(ax(ii),ops.Mdl.coef(ind{jj},pairs(ii,1)),ops.Mdl.coef(ind{jj},pairs(ii,2)),[],Colors(jj,:));
		end
	end
	xlabel(ax(ii),ops.events{pairs(ii,1)});
	ylabel(ax(ii),ops.events{pairs(ii,2)});
end

title(ax(2),'decoder coefficient');
arrayfun(@(h) set(h,'XLim',[min([h.XLim(1) h.YLim(1)]) max([h.XLim(2) h.YLim(2)])],'YLim',[min([h.XLim(1) h.YLim(1)]) max([h.XLim(2) h.YLim(2)])]),ax);
set(gcf,'Position',[0 0 505 170]);
export_fig(sprintf('results/%scoef_%s_%s.pdf',prefix,data.subject,datestr(data.session,'YYmmdd')));


end