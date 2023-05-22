%% user input
mouse = 'calca_280';
session = '2023-02-14';
%% load data
close all; clearvars -except mouse session; clc
addpath(genpath('Z:\Chris\matlab\numpy-matlab\'));
addpath(genpath('Z:\Chris\matlab\heatmap\'));
addpath(genpath('Z:\Chris\matlab\cz\neuropixels-utils\'));
dpath = ['Z:\Chris\data\neuropixels_cta\',mouse,'\',session,'\',mouse,'_',session,'.mat'];
load(dpath,'data')
behavior = data;
binName = 'cz_npxl_g0_t0.imec0.ap.bin';
path = ['Z:\Chris\data\neuropixels_cta\',mouse,'\',session,'\cz_npxl_g0\cz_npxl_g0_imec0\'];
meta = ReadMeta(binName, path);
times = [0:str2double(meta.fileSizeBytes)/(2*str2double(meta.nSavedChans))]./str2double(meta.imSampRate);
myKsDir = ['Z:\Chris\data\neuropixels_cta\',mouse,'\',session,'\catgt_cz_npxl_g0\'];
spike_times = times(readNPY([myKsDir,'spike_times.npy']))';
spike_clusters = readNPY([myKsDir,'spike_clusters.npy']);
metrics = tdfread([myKsDir,'cluster_info.tsv']);
path = ['Z:\Chris\data\neuropixels_cta\',mouse,'\',session,'\cz_npxl_g0\cz_npxl_g0_imec0\'];
good_units = metrics.cluster_id(find(sum(metrics.group=='g',2)|sum(metrics.group=='m',2)));
spike_amps = readNPY([myKsDir,'amplitudes.npy'])*2.3438';
load(['Z:/Chris/data/neuropixels_cta/',mouse,'/channel_locations.mat']);
min_amp = 20; amp = [];
min_fr = 0.1; fr = [];
loc = cell(0,0);
tstart = min([behavior.rewards.all.rear(1) behavior.rewards.all.front(1)])-1;
tend   = behavior.laser(end,2)+1;
spike_clusters(spike_times<tstart | spike_times>tend) = NaN;
load('Z:\Chris\matlab\cz\neuropixels-utils\Allen_CCFv3_Ephys_AST.mat','atlas')
min_isi = 0.0002;
for i = 1:length(good_units)
    spks_idx = find(spike_clusters==good_units(i));
    spks_idx(find(diff(spike_times(spks_idx))<=min_isi)+1) = [];
    amp(i) = median(spike_amps(spks_idx));
    fr(i)  = length(spks_idx)/(tend-tstart);
    loc{i} = channel_locations.(['shank',num2str(metrics.sh(find(metrics.cluster_id==good_units(i))))]).brain_region(metrics.ch(find(metrics.cluster_id==good_units(i)))+1,:);
end
good_units_amygdala = good_units(abs(amp)>=min_amp & fr>=min_fr & (cellfun(@(x) contains(x,'CEA'),loc)|cellfun(@(x) contains(x,'STR'),loc)));
good_units_amygdala = RemoveDuplicateUnits(good_units_amygdala,spike_times,spike_clusters);
good_units = good_units_amygdala;
%good_units_other = good_units(abs(amp)>=min_amp & fr>=min_fr & ~(cellfun(@(x) contains(x,'CEA'),loc)|cellfun(@(x) contains(x,'STR'),loc)));
%good_units_other = RemoveDuplicateUnits(good_units_other,spike_times,spike_clusters);
%good_units = sort([good_units_amygdala;good_units_other]);
for i = 1:length(good_units)
    spks_idx = find(spike_clusters==good_units(i));
    spks_idx(find(diff(spike_times(spks_idx))<=min_isi)+1) = [];
    ephys.spikes{i} = spike_times(spks_idx);
    ephys.cluster_id(i) = good_units(i);
    ephys.ks_label{i} = metrics.KSLabel(find(metrics.cluster_id==good_units(i)),:);
    ephys.phy_label{i} = metrics.group(find(metrics.cluster_id==good_units(i)),:);
    ephys.amplitude(i) = median(spike_amps(spks_idx));
    ephys.firing_rate(i) = length(spks_idx)/(tend-tstart);
    ephys.isi_viol(i) = contamination(ephys.spikes{i},'min_time',tstart,'max_time',tend,'min_isi',min_isi);
    ephys.shank(i) = metrics.sh(find(metrics.cluster_id==good_units(i)));
    ephys.channel(i) = metrics.ch(find(metrics.cluster_id==good_units(i)));
    ephys.region{i} = channel_locations.(['shank',num2str(ephys.shank(i))]).brain_region(ephys.channel(i)+1,:);
    ephys.region{i} = ephys.region{i}(find(~isspace(ephys.region{i})));
    ephys.x(i) = channel_locations.(['shank',num2str(ephys.shank(i))]).x(ephys.channel(i)+1);
    ephys.y(i) = channel_locations.(['shank',num2str(ephys.shank(i))]).y(ephys.channel(i)+1);
    ephys.z(i) = channel_locations.(['shank',num2str(ephys.shank(i))]).z(ephys.channel(i)+1);
    if isequal(ephys.region{i},'STR')
        [ML,AP,DV] = bregma2ccf(ephys.x(i),ephys.y(i),ephys.z(i));
        if atlas(round(ML),round(AP),round(DV)) == 4E5
            ephys.region{i} = 'CEAv';
        elseif atlas(round(ML),round(AP),round(DV)) == 2E5
            ephys.region{i} = 'CEAast';
        end
    end
end
ephys.trajectory = channel_locations;
clearvars -except ephys behavior mouse session
%% generate PSTHs

tstart = 0;
tend = ceil(behavior.laser(end,2))+10;
binsize = 0.010;
edges = tstart:binsize:tend;
centers = edges(1:end-1)+binsize/2;
spikedata = nan(length(ephys.spikes),length(centers));
a1 = round((min([behavior.rewards.all.front(1),behavior.rewards.all.rear(1)])-1)./binsize);
a2 = round((max([behavior.rewards.all.front(end),behavior.rewards.all.rear(end)])+5)./binsize);
mu = [];
sigma = [];
for i = 1:length(ephys.spikes)
    spikedata(i,:) = histcounts(ephys.spikes{i},edges);
    mu(i) = mean(spikedata(i,a1:a2));
    sigma(i) = std(spikedata(i,a1:a2));
end

e = -2:binsize:10;
c = e(1:end-1)+binsize/2;
if isequal(behavior.novelport,'rear')
    event = behavior.rewards.all.rear;
elseif isequal(behavior.novelport,'front')
    event = behavior.rewards.all.front;
end
PSTH.novel = nan(length(ephys.spikes),length(event),length(c));
for i = 1:length(ephys.spikes)
    for j = 1:length(event)
        spks = ephys.spikes{i} - event(j);
        PSTH.novel(i,j,:) = (histcounts(spks,e)-mu(i))./sigma(i);
    end
end

if isequal(behavior.novelport,'rear')
    event = behavior.rewards.all.front;
elseif isequal(behavior.novelport,'front')
    event = behavior.rewards.all.rear;
end
PSTH.familiar = nan(length(ephys.spikes),length(event),length(c));
for i = 1:length(ephys.spikes)
    for j = 1:length(event)
        spks = ephys.spikes{i} - event(j);
        PSTH.familiar(i,j,:) = (histcounts(spks,e)-mu(i))./sigma(i);
    end
end

e = -2:binsize:10;
c = e(1:end-1)+binsize/2;
event = behavior.laser(:,1);
PSTH.CGRP = nan(length(ephys.spikes),length(behavior.laser),length(c));
for i = 1:length(ephys.spikes)
    for j = 1:length(event)
        spks = ephys.spikes{i} - event(j);
        PSTH.CGRP(i,j,:) = (histcounts(spks,e)-mu(i))./sigma(i);
    end
end
%% classify neurons
for i = 1:length(ephys.spikes)
    X = squeeze([PSTH.novel(i,:,:);PSTH.familiar(i,:,:)]);
    if sum(isnan(X(:))) == 0
        decoders.individual.clusterid(i) = ephys.cluster_id(i);
        decoders.individual.region{i} = ephys.region{i};
        p = ranksum(mean(squeeze(PSTH.novel(i,:,201:300)),2),mean(squeeze(PSTH.familiar(i,:,201:300)),2));
        p2 = signrank([mean(squeeze(PSTH.novel(i,:,201:300)),2);mean(squeeze(PSTH.familiar(i,:,201:300)),2)],[mean(squeeze(PSTH.novel(i,:,101:200)),2);mean(squeeze(PSTH.familiar(i,:,101:200)),2)]);
        decoders.individual.pvalue(i) = p;
        decoders.individual.pvalue2(i) = p2;
        decoders.individual.response_novel(i) = mean(mean(squeeze(PSTH.novel(i,:,201:300))))-mean(mean(squeeze(PSTH.novel(i,:,101:200))));
        decoders.individual.response_familiar(i) = mean(mean(squeeze(PSTH.familiar(i,:,201:300))))-mean(mean(squeeze(PSTH.familiar(i,:,101:200))));
        decoders.individual.preference(i) = decoders.individual.response_novel(i)-decoders.individual.response_familiar(i);
        disp(['Finished decoder for neuron #',num2str(i,'%03.f'),', p = ',num2str(decoders.individual.pvalue(i),'%.3f')])
    else
        decoders.individual.preference(i) = 0;
        decoders.individual.performance(i) = 0;
        decoders.individual.response_novel(i) = 0;
        decoders.individual.response_familiar(i) = 0;
        decoders.individual.pvalue(i) = 1;
        decoders.individual.pvalue2(i) = 1;
    end
end
decoders.individual.pvalue(isnan(decoders.individual.pvalue)) = 1;
decoders.individual.pvalue2(isnan(decoders.individual.pvalue2)) = 1;
%% make data structure
PSTHdata = struct;
pval = 0.05;
modulated = find(decoders.individual.pvalue<=pval);
nomod = find(decoders.individual.pvalue>pval & decoders.individual.pvalue2>pval);
nopref = find(decoders.individual.pvalue>pval & decoders.individual.pvalue2<=pval);
novelpref = modulated(decoders.individual.preference(modulated)>0);
familiarpref = modulated(decoders.individual.preference(modulated)<0);

PSTHdata.novelpref.novel = PSTH.novel(novelpref,:,:);
PSTHdata.novelpref.familiar = PSTH.familiar(novelpref,:,:);
PSTHdata.novelpref.CGRP = PSTH.CGRP(novelpref,:,:);
PSTHdata.novelpref.cgrp_period = [];
PSTHdata.novelpref.drinking_period = [];
PSTHdata.novelpref.pvalues = decoders.individual.pvalue(novelpref);
PSTHdata.novelpref.clusterid = decoders.individual.clusterid(novelpref);
PSTHdata.novelpref.region = decoders.individual.region(novelpref);
PSTHdata.novelpref.preference = decoders.individual.preference(novelpref);

PSTHdata.familiarpref.novel = PSTH.novel(familiarpref,:,:);
PSTHdata.familiarpref.familiar = PSTH.familiar(familiarpref,:,:);
PSTHdata.familiarpref.CGRP = PSTH.CGRP(familiarpref,:,:);
PSTHdata.familiarpref.cgrp_period = [];
PSTHdata.familiarpref.drinking_period = [];
PSTHdata.familiarpref.pvalues = decoders.individual.pvalue(familiarpref);
PSTHdata.familiarpref.clusterid = decoders.individual.clusterid(familiarpref);
PSTHdata.familiarpref.region = decoders.individual.region(familiarpref);
PSTHdata.familiarpref.preference = decoders.individual.preference(familiarpref);

PSTHdata.nopref.novel = PSTH.novel(nopref,:,:);
PSTHdata.nopref.familiar = PSTH.familiar(nopref,:,:);
PSTHdata.nopref.CGRP = PSTH.CGRP(nopref,:,:);
PSTHdata.nopref.cgrp_period = [];
PSTHdata.nopref.drinking_period = [];
PSTHdata.nopref.pvalues = decoders.individual.pvalue(nopref);
PSTHdata.nopref.clusterid = decoders.individual.clusterid(nopref);
PSTHdata.nopref.region = decoders.individual.region(nopref);
PSTHdata.nopref.preference = decoders.individual.preference(nopref);

PSTHdata.nomod.novel = PSTH.novel(nomod,:,:);
PSTHdata.nomod.familiar = PSTH.familiar(nomod,:,:);
PSTHdata.nomod.CGRP = PSTH.CGRP(nomod,:,:);
PSTHdata.nomod.cgrp_period = [];
PSTHdata.nomod.drinking_period = [];
PSTHdata.nomod.pvalues = decoders.individual.pvalue(nomod);
PSTHdata.nomod.clusterid = decoders.individual.clusterid(nomod);
PSTHdata.nomod.region = decoders.individual.region(nomod);
PSTHdata.nomod.preference = decoders.individual.preference(nomod);

e = [-30*60:binsize:45*60]+binsize;
c = e(1:end-1)+binsize/2;
PSTH.laser_pd = zeros(length(ephys.spikes),length(c));
mu = []; sigma = [];
for i = 1:length(ephys.spikes)
    event = behavior.laser(1,1);
    t = nan(1,length(c));
    idx = round(event./binsize)+round((e(1:end-1)/binsize));
    t(1,:) = spikedata(i,idx);
    delay_idx1 = 1;
    delay_idx2 = 30*60/binsize;
    mu(i) = mean(t(delay_idx1:delay_idx2));
    sigma(i) = std(t(delay_idx1:delay_idx2));
    PSTH.laser_pd(i,:) = (t - mu(i))/sigma(i);
end
out = [];
for i = 1:length(ephys.spikes)
    y = PSTH.laser_pd(i,:);
    out(i,:) = mean(reshape(y(:),60/binsize,[]));
end

PSTHdata.novelpref.cgrp_period =  out(novelpref,:);
PSTHdata.familiarpref.cgrp_period = out(familiarpref,:);
PSTHdata.nopref.cgrp_period = out(nopref,:);
PSTHdata.nomod.cgrp_period = out(nomod,:);

e = [0:binsize:20*60]+binsize;
c = e(1:end-1)+binsize/2;
PSTH.laser_pd = zeros(length(ephys.spikes),length(c));
for i = 1:length(ephys.spikes)
    event = min([behavior.rewards.all.front(1),behavior.rewards.all.rear(1)]);
    t = nan(1,length(c));
    idx = round(event./binsize)+round((e(1:end-1)/binsize));
    t(1,:) = spikedata(i,idx);
    PSTH.laser_pd(i,:) = (t - mu(i))/sigma(i);
end
out = [];
for i = 1:length(ephys.spikes)
    y = PSTH.laser_pd(i,:);
    out(i,:) = mean(reshape(y(:),60/binsize,[]));
end
PSTHdata.novelpref.drinking_period =  out(novelpref,:);
PSTHdata.familiarpref.drinking_period = out(familiarpref,:);
PSTHdata.nopref.drinking_period = out(nopref,:);
PSTHdata.nomod.drinking_period = out(nomod,:);

save([mouse,'-',session,'.mat'],'PSTHdata','behavior','ephys')
%% plot
mgw = myGaussWin(0.1,1/binsize);
mgw(1:round(numel(mgw)/2)-1) = 0; mgw = mgw./sum(mgw);
cmap = flipud(cbrewer('div','RdBu',1000,'spline')); cmap(cmap<0) = 0;

a = []; b = []; c = [];
for i = 1:length(novelpref)
    a(i,:) = conv(squeeze(mean(PSTHdata.novelpref.novel(i,:,:),2)),mgw,'same');
    b(i,:) = conv(squeeze(mean(PSTHdata.novelpref.familiar(i,:,:),2)),mgw,'same');
    c(i,:) = conv(squeeze(mean(PSTHdata.novelpref.CGRP(i,:,:),2)),mgw,'same');
end
a = a - mean(a(:,101:200),2); b = b - mean(b(:,101:200),2); c = c - mean(c(:,101:200),2);
a = a(:,101:700); b = b(:,101:700); c = c(:,101:600);
plotdata1 = [a nan(length(novelpref),100) b nan(length(novelpref),100) c];
[~,i1] = sort(decoders.individual.response_novel(novelpref),'descend');
plotdata1 = plotdata1(i1,:);

a = []; b = []; c = [];
for i = 1:length(familiarpref)
    a(i,:) = conv(squeeze(mean(PSTHdata.familiarpref.novel(i,:,:),2)),mgw,'same');
    b(i,:) = conv(squeeze(mean(PSTHdata.familiarpref.familiar(i,:,:),2)),mgw,'same');
    c(i,:) = conv(squeeze(mean(PSTHdata.familiarpref.CGRP(i,:,:),2)),mgw,'same');
end
a = a - mean(a(:,101:200),2); b = b - mean(b(:,101:200),2); c = c - mean(c(:,101:200),2);
a = a(:,101:700); b = b(:,101:700); c = c(:,101:600);
plotdata2 = [a nan(length(familiarpref),100) b nan(length(familiarpref),100) c];
[~,i2] = sort(decoders.individual.response_familiar(familiarpref),'descend');
plotdata2 = plotdata2(i2,:);

a = []; b = []; c = [];
for i = 1:length(nopref)
    a(i,:) = conv(squeeze(mean(PSTHdata.nopref.novel(i,:,:),2)),mgw,'same');
    b(i,:) = conv(squeeze(mean(PSTHdata.nopref.familiar(i,:,:),2)),mgw,'same');
    c(i,:) = conv(squeeze(mean(PSTHdata.nopref.CGRP(i,:,:),2)),mgw,'same');
end
a = a - mean(a(:,101:200),2); b = b - mean(b(:,101:200),2); c = c - mean(c(:,101:200),2);
a = a(:,101:700); b = b(:,101:700); c = c(:,101:600);
plotdata3 = [a nan(length(nopref),100) b nan(length(nopref),100) c];
[~,i3] = sort(mean([decoders.individual.response_novel(nopref);decoders.individual.response_familiar(nopref)]),'descend');
plotdata3 = plotdata3(i3,:);

a = []; b = []; c = [];
for i = 1:length(nomod)
    a(i,:) = conv(squeeze(mean(PSTHdata.nomod.novel(i,:,:),2)),mgw,'same');
    b(i,:) = conv(squeeze(mean(PSTHdata.nomod.familiar(i,:,:),2)),mgw,'same');
    c(i,:) = conv(squeeze(mean(PSTHdata.nomod.CGRP(i,:,:),2)),mgw,'same');
end
a = a - mean(a(:,101:200),2); b = b - mean(b(:,101:200),2); c = c - mean(c(:,101:200),2);
a = a(:,101:700); b = b(:,101:700); c = c(:,101:600);
plotdata4 = [a nan(length(nomod),100) b nan(length(nomod),100) c];
[~,i4] = sort(mean([decoders.individual.response_novel(nomod);decoders.individual.response_familiar(nomod)]),'descend');
plotdata4 = plotdata4(i4,:);

plotdata = [plotdata1;nan(2,1900);plotdata2;nan(2,1900);plotdata3;nan(2,1900);plotdata4];

figure('Position', get(0, 'Screensize'))
subplot(1,2,1)
hold on
heatmap(flipud(plotdata),[],[],[],'Colormap',cmap,'ColorLevels',1000,'MaxColorValue',1,'MinColorValue',-1,'NaNColor',[1 1 1]);
plot([101 101]+0.5,[0 length(ephys.spikes)+6]+0.5,'k-','LineWidth',2)
plot([801 801]+0.5,[0 length(ephys.spikes)+6]+0.5,'k-','LineWidth',2)
plot([1501 1501]+0.5,[0 length(ephys.spikes)+6]+0.5,'k-','LineWidth',2)
plot([1801 1801]+0.5,[0 length(ephys.spikes)+6]+0.5,'k-','LineWidth',2)
xticks([300 950 1650])
xticklabels({'Novel','Water','CGRP stim'})
yticks([length(nomod)./2 length(nomod)+length(nopref)./2+2 length(nomod)+length(nopref)+length(familiarpref)./2+4 length(nomod)+length(nopref)+length(familiarpref)+length(novelpref)./2+6])
yticklabels({['Not modulated (',num2str(length(nomod)),')'],['No preference (',num2str(length(nopref)),')'],['Water (',num2str(length(familiarpref)),')'],['Novel (',num2str(length(novelpref)),')']})
ytickangle(90)
xtickangle(0)
xlabel('Time (–1 sec → +5 sec)')
set(gca,'FontSize',16)
h = colorbar;
ylabel(h,'Firing rate (\sigma)');
h.Ticks = -1:.5:1;
h.Box = 'off';
hold off

subplot(1,2,2)
hold on
idx1 = max([behavior.rewards.all.rear(end),behavior.rewards.all.front(end)])/60;
idx2 = behavior.laser(1,1)/60;
idx3 = behavior.laser(end,2)/60;
a = [PSTHdata.novelpref.drinking_period(i1,:),nan(size(PSTHdata.novelpref.cgrp_period(i1,:),1),5),PSTHdata.novelpref.cgrp_period(i1,:)];
b = [PSTHdata.familiarpref.drinking_period(i2,:),nan(size(PSTHdata.familiarpref.cgrp_period(i2,:),1),5),PSTHdata.familiarpref.cgrp_period(i2,:)];
c = [PSTHdata.nopref.drinking_period(i3,:),nan(size(PSTHdata.nopref.cgrp_period(i3,:),1),5),PSTHdata.nopref.cgrp_period(i3,:)];
d = [PSTHdata.nomod.drinking_period(i4,:),nan(size(PSTHdata.nomod.cgrp_period(i4,:),1),5),PSTHdata.nomod.cgrp_period(i4,:)];
plotdata = [a;nan(2,100);b;nan(2,100);c;nan(2,100);d];
heatmap(flipud(plotdata),[],[],[],'Colormap',cmap,'ColorLevels',1000,'MaxColorValue',1,'MinColorValue',-1,'NaNColor',[1 1 1]);
plot([55 55]+0.5,[0 length(ephys.spikes)+6]+0.5,'k-','LineWidth',2)
xticks([10 40 77.5])
xticklabels({'Drinking','Delay','CGRP-stim'})
xlabel(['Time (100 min)'])
set(gca,'FontSize',16)
h = colorbar;
ylabel(h,'Firing rate (\sigma)');
h.Ticks = -1:.5:1;
h.Box = 'off';
hold off

set(gcf,'renderer','Painters')
%saveas(gcf,[dpath(1:end-4),'_drinking'],'png')
%saveas(gcf,[dpath(1:end-4),'_drinking'],'epsc')
%% ReadSGLXData functions

% =========================
% General Utility Functions
% =========================


% =========================================================
% Parse ini file returning a structure whose field names
% are the metadata left-hand-side tags, and whose right-
% hand-side values are MATLAB strings. We remove any
% leading '~' characters from tags because MATLAB uses
% '~' as an operator.
%
% If you're unfamiliar with structures, the benefit
% is that after calling this function you can refer
% to metafile items by name. For example:
%
%   meta.fileCreateTime  // file create date and time
%   meta.nSavedChans     // channels per timepoint
%
% All of the values are MATLAB strings, but you can
% obtain a numeric value using str2double(meta.nSavedChans).
% More complicated parsing of values is demonstrated in the
% utility functions below.
%
function [meta] = ReadMeta(binName, path)

% Create the matching metafile name
[dumPath,name,dumExt] = fileparts(binName);
metaName = strcat(name, '.meta');

% Parse ini file into cell entries C{1}{i} = C{2}{i}
fid = fopen(fullfile(path, metaName), 'r');
% -------------------------------------------------------------
%    Need 'BufSize' adjustment for MATLAB earlier than 2014
%    C = textscan(fid, '%[^=] = %[^\r\n]', 'BufSize', 32768);
C = textscan(fid, '%[^=] = %[^\r\n]');
% -------------------------------------------------------------
fclose(fid);

% New empty struct
meta = struct();

% Convert each cell entry into a struct entry
for i = 1:length(C{1})
    tag = C{1}{i};
    if tag(1) == '~'
        % remake tag excluding first character
        tag = sprintf('%s', tag(2:end));
    end
    meta = setfield(meta, tag, C{2}{i});
end
end % ReadMeta


% =========================================================
% Read nSamp timepoints from the binary file, starting
% at timepoint offset samp0. The returned array has
% dimensions [nChan,nSamp]. Note that nSamp returned
% is the lesser of: {nSamp, timepoints available}.
%
function dataArray = ReadBin(samp0, nSamp, meta, binName, path)

nChan = str2double(meta.nSavedChans);

nFileSamp = str2double(meta.fileSizeBytes) / (2 * nChan);
samp0 = max(samp0, 0);
nSamp = min(nSamp, nFileSamp - samp0);

sizeA = [nChan, nSamp];

fid = fopen(fullfile(path, binName), 'rb');
fseek(fid, samp0 * 2 * nChan, 'bof');
dataArray = fread(fid, sizeA, 'int16=>double');
fclose(fid);
end % ReadBin


% =========================================================
% Return sample rate as double.
%
function srate = SampRate(meta)
if strcmp(meta.typeThis, 'imec')
    srate = str2double(meta.imSampRate);
else
    srate = str2double(meta.niSampRate);
end
end % SampRate


% =========================================================
% Return a multiplicative factor for converting 16-bit
% file data to voltage. This does not take gain into
% account. The full conversion with gain is:
%
%   dataVolts = dataInt * fI2V / gain.
%
% Note that each channel may have its own gain.
%
function fI2V = Int2Volts(meta)
if strcmp(meta.typeThis, 'imec')
    fI2V = str2double(meta.imAiRangeMax) / 512;
else
    fI2V = str2double(meta.niAiRangeMax) / 32768;
end
end % Int2Volts


% =========================================================
% Return array of original channel IDs. As an example,
% suppose we want the imec gain for the ith channel stored
% in the binary data. A gain array can be obtained using
% ChanGainsIM() but we need an original channel index to
% do the look-up. Because you can selectively save channels
% the ith channel in the file isn't necessarily the ith
% acquired channel, so use this function to convert from
% ith stored to original index.
%
% Note: In SpikeGLX channels are 0-based, but MATLAB uses
% 1-based indexing, so we add 1 to the original IDs here.
%
function chans = OriginalChans(meta)
if strcmp(meta.snsSaveChanSubset, 'all')
    chans = (1:str2double(meta.nSavedChans));
else
    chans = str2num(meta.snsSaveChanSubset);
    chans = chans + 1;
end
end % OriginalChans


% =========================================================
% Return counts of each nidq channel type that compose
% the timepoints stored in binary file.
%
function [MN,MA,XA,DW] = ChannelCountsNI(meta)
M = str2num(meta.snsMnMaXaDw);
MN = M(1);
MA = M(2);
XA = M(3);
DW = M(4);
end % ChannelCountsNI


% =========================================================
% Return counts of each imec channel type that compose
% the timepoints stored in binary file.
%
function [AP,LF,SY] = ChannelCountsIM(meta)
M = str2num(meta.snsApLfSy);
AP = M(1);
LF = M(2);
SY = M(3);
end % ChannelCountsIM


% =========================================================
% Return gain for ith channel stored in the nidq file.
%
% ichan is a saved channel index, rather than an original
% (acquired) index.
%
function gain = ChanGainNI(ichan, savedMN, savedMA, meta)
if ichan <= savedMN
    gain = str2double(meta.niMNGain);
elseif ichan <= savedMN + savedMA
    gain = str2double(meta.niMAGain);
else
    gain = 1;
end
end % ChanGainNI


% =========================================================
% Return gain arrays for imec channels.
%
% Index into these with original (acquired) channel IDs.
%
function [APgain,LFgain] = ChanGainsIM(meta)
C = textscan(meta.imroTbl, '(%*s %*s %*s %d %d %*s', ...
    'EndOfLine', ')', 'HeaderLines', 1 );
APgain = double(cell2mat(C(1)));
LFgain = double(cell2mat(C(2)));
end % ChanGainsIM


% =========================================================
% Having acquired a block of raw nidq data using ReadBin(),
% convert values to gain-corrected voltages. The conversion
% is only applied to the saved-channel indices in chanList.
% Remember saved-channel indices are in range [1:nSavedChans].
% The dimensions of the dataArray remain unchanged. ChanList
% examples:
%
%   [1:MN]      % all MN chans (MN from ChannelCountsNI)
%   [2,6,20]    % just these three channels
%
function dataArray = GainCorrectNI(dataArray, chanList, meta)

[MN,MA] = ChannelCountsNI(meta);
fI2V = Int2Volts(meta);

for i = 1:length(chanList)
    j = chanList(i);    % index into timepoint
    conv = fI2V / ChanGainNI(j, MN, MA, meta);
    dataArray(j,:) = dataArray(j,:) * conv;
end
end


% =========================================================
% Having acquired a block of raw imec data using ReadBin(),
% convert values to gain-corrected voltages. The conversion
% is only applied to the saved-channel indices in chanList.
% Remember saved-channel indices are in range [1:nSavedChans].
% The dimensions of the dataArray remain unchanged. ChanList
% examples:
%
%   [1:AP]      % all AP chans (AP from ChannelCountsIM)
%   [2,6,20]    % just these three channels
%
function dataArray = GainCorrectIM(dataArray, chanList, meta)

% Look up gain with acquired channel ID
chans = OriginalChans(meta);
[APgain,LFgain] = ChanGainsIM(meta);
nAP = length(APgain);
nNu = nAP * 2;

% Common conversion factor
fI2V = Int2Volts(meta);

for i = 1:length(chanList)
    j = chanList(i);    % index into timepoint
    k = chans(j);       % acquisition index
    if k <= nAP
        conv = fI2V / APgain(k);
    elseif k <= nNu
        conv = fI2V / LFgain(k - nAP);
    else
        continue;
    end
    dataArray(j,:) = dataArray(j,:) * conv;
end
end