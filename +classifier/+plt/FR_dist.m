FR = cellfun(@(x) numel(x)/(max(x)-min(x)), data.spikes);

ax = np; histogram(FR); export_fig tmp.pdf;

