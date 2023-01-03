ax = np(2,2);
% cue raster
tp = -1:0.01:3;
n_cells = numel(data.spikes);
ax = np(1,2);
eventtime = data.laser.onsets;

% laser onset
ax = np(n_cells);
arrayfun(@(i) raster(data.spikes{i}*1000,eventtime * 1000,ax(i),'tp',tp), 1:n_cells);

ax = np(n_cells);
arrayfun(@(i) plt.psth_line(data.spikes{i}*1000,eventtime*1000,ax(i),'tp',tp,'kernel_width',0.05), 1:n_cells);
% plt.psth_line(data.spikes*1000,data.rewards.all.front*1000,ax(3),'tp',tp,'line_color',[1 0 0],'kernel_width',0.05);

