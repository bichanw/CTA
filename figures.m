%% Code to produce required figures

% chris k99, figure 2c raster
data = load_data.all(datetime(2022,12,4),'001');

% load trained model and corresponding parameters
load('figures/k99_2c.mat');

% define time points to plot
ops.posterior_t_edges = 3780:0.5:3810; % left edges for spike counting windows

% count spk
[spk_count,ops] = classifier.count_spk.time_course(data,ops);