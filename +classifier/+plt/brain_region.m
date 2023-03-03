function brain_region(data,ops,prefix)
if nargin < 3
	prefix = '';
end

ax = np(4,1);

Regions = unique(data.brain_region);
% histogram(ax(1),categorical(data.brain_region(:)),Regions);
histogram(ax(1),categorical(data.brain_region(unique(sub2ind(size(data.brain_region),data.sh+1,data.ch+1)))),Regions);
title(ax(1),'channels with units');
ylabel(ax(1),'# channels');

histogram(ax(2),categorical(data.brain_region(sub2ind(size(data.brain_region),data.sh+1,data.ch+1))),Regions);
title(ax(2),'unit distribution');
ylabel(ax(2),'# units');


% plot non-zero coef based on mnr
% only applies for l1
if isfield(ops,'mnr') && strcmp(ops.mnr.penalty,'l1')

	% retrain model if necessary
	if ~isfield(ops,'Mdl')
		[~,ops] = ops.classifier.train(data,ops);
	end
	ind = find(~ops.exclude_id);

	% front non-zero
	title(ax(3),'front units');
	histogram(ax(3),categorical(data.brain_region(sub2ind(size(data.brain_region),data.sh(ind(ops.Mdl.coef(:,1)~=0))+1,data.ch(ind(ops.Mdl.coef(:,1)~=0))+1))),Regions);

	% rear non-zero
	title(ax(4),'rear units');
	histogram(ax(4),categorical(data.brain_region(sub2ind(size(data.brain_region),data.sh(ind(ops.Mdl.coef(:,2)~=0))+1,data.ch(ind(ops.Mdl.coef(:,2)~=0))+1))),Regions);
end


set(gcf,'Position',[0 0 400 600]);
export_fig(sprintf('results/%sBrainRegion_%s_%s.pdf',prefix,data.subject,datestr(data.session,'YYmmdd')));
