classdef load_data < handle
% all functions to load data


	methods (Static)

		function data = all(session,subject)
			% master script to load all sorts of data on one day
			% session - datetime variable of the day of the recording
			% subject - subject number
			% e.g. data = load_all(datetime(2022,7,15),'237');

			% default example
			if nargin == 0
				session = datetime(2022,7,15);
				subject = '238';
			end

			% main folder where data is stored
			% session = datetime(2022,11,30);
			% subject = '001';
			folder = sprintf('/jukebox/witten/Chris/data/neuropixels_cta/calca_%s/%s',subject,datestr(session,'yyyy-mm-dd'));


			% ----- event data -----
			% main struct
				event_f = [dir(sprintf('%s/data_for_bichan.mat',folder)), dir(sprintf('%s/calca_%s_%s.mat',folder,subject,datestr(session,'yyyy-mm-dd')))];
				if ~isempty(event_f)
					load([event_f.folder filesep event_f.name]);
				else
					data = struct;
				end
			% front / rare - water / grape pair
				port_f = dir(sprintf('%s/reward_location_*.txt',folder));
				if ~isempty(port_f)
					f = fopen([port_f.folder '/' port_f.name],'r');
					C = textscan(f,'%s %s %s\n');
					data.port_is_water = [false false];
					if (strcmpi('Rear',C{1}) && strcmpi('Water',C{3})) || (strcmpi('Front',C{1}) && ~strcmpi('Water',C{3}))
						data.port_is_water(2) = true;
					else
						data.port_is_water(1) = true;
					end
					fclose(f);
				end
			% novel / familiar preference based on 
				if exist([folder '/psths.mat']);
					load([folder '/psths.mat']);
					data.PSTHdata = PSTHdata;
				end


			% ----- spike data -----
			% if session>datetime(2022,11,21)
			% 	data_f = [folder filesep 'catgt_cz_npxl_g0'];
			% 	fprintf(data_f);
			% else
			% 	data_f = [folder '/pykilosortfiles/results'];
			% end
			data_f = load_data.ks_folder(session,subject);
			data = load_data.spike_ks(data_f,data);


			% ----- video data -----
			video_f = dir(sprintf('~/SpikeSorting/Codes/witten/mat/sleap/%s_*%s/*.h5',datestr(session,'yyyymmdd'),subject));
			if ~isempty(video_f)
				f = video_f;
				tracks = h5read([f.folder '/' f.name],'/tracks');
				pos_center = squeeze(tracks(:,2,:));

				% linear interpolation
				nanx = isnan(pos_center(:,1));
				for i = 1:2
					pos_center(nanx,i) = interp1(find(~nanx), pos_center(~nanx,i), find(nanx));
				end
				data.pos = pos_center;
			end


			% save session info
			data.session = session;
			data.subject = subject;

		end


		function folder = ks_folder(session,subject)
			folder = sprintf('/jukebox/witten/Chris/data/neuropixels_cta/calca_%s/%s',subject,datestr(session,'yyyy-mm-dd'));
			if session>datetime(2022,11,21)
				folder = [folder filesep 'catgt_cz_npxl_g0'];
			else
				folder = [folder '/pykilosortfiles/results'];
			end
		end



		% spike data from kilosort
		function data = spike_ks(folder, data, if_good)
			% load kilosort into data struct format
			% folder  - kilosort results folder
			% data    - output struct 
			% if_good - if selecting only good units 

			% if only keeping good units
			if nargin < 3
				if_good = false;
			end

			% load sampling rate
			if exist('pyrunfile')
				SR = double(pyrunfile([folder filesep 'params.py'],'sample_rate'));
			else % default SR
				SR = 3e4;
			end

			% load cluster group
			[data.cids, data.cgs] = readClusterGroups_bw([folder '/cluster_info.tsv']); % unit group
			data.cids(data.cgs==0) = [];
			data.cgs(data.cgs==0) = [];
			% - 0 = noise
			% - 1 = mua
			% - 2 = good
			% - 3 = unsorted

			% load spikes
			spike_times = readNPY([folder '/spike_times.npy']);
			spike_clusters = readNPY([folder '/spike_clusters.npy']);
			% spike_amp = readNPY([folder '/amplitudes.npy']);
			for i = 1:numel(data.cids)
				data.spikes{i} = double(spike_times(spike_clusters==data.cids(i)))/SR;
			end
		end

		% load spike amplitude
		function amp = spike_amp(data)
			% load data information
			data_f = load_data.ks_folder(data.session,data.subject);
			spike_amp = readNPY([data_f '/amplitudes.npy']);
			spike_clusters = readNPY([data_f '/spike_clusters.npy']);

			% assign amplitude
			for i = 1:numel(data.cids)
				amp{i} = double(spike_amp(spike_clusters==data.cids(i)));
			end
		
		end

	
	end


end