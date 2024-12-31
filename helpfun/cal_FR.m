function FR = cal_raster(spikes)
% calculate firing rate based on input spike times
% input: spikes - nneurons * 1 cells, spike times

FR = cellfun(@(x) numel(x) / range(x),spikes);

end
