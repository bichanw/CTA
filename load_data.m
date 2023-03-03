classdef load_data < handle
% all functions to load data


	methods (Static)


		function brain_region = anatomical(data)

			% load data
			subject = data.subject;
			folder = sprintf('/jukebox/witten/Chris/data/neuropixels_cta/calca_%s',subject);
			for ishank = 1:4
				f = [folder filesep 'channel_locations_shank' num2str(ishank-1) '.csv'];
				
				% when there is data
				if exist(f)
					opts = detectImportOptions(f);
					T = readtable(f,opts);

					% extend region 
					tmp = T.brain_region';
					if numel(tmp) == 383
						tmp(end+1) = tmp(end);
					end
					brain_region(ishank,:) = tmp;
				else
					brain_region(ishank,:) = repmat({nan},1,384);
				end
			end
		end

		function data = all(session,subject,ops)
			% master script to load all sorts of data on one day
			% session - datetime variable of the day of the recording
			% subject - subject number
			% e.g. data = load_all(datetime(2022,7,15),'237');

			% default example
			if nargin == 0
				session = datetime(2022,12,4);
				subject = '001';
			end
			if nargin < 3
				ops = struct();
			end

			% main folder where data is stored
			folder = sprintf('/jukebox/witten/Chris/data/neuropixels_cta/calca_%s/%s',subject,datestr(session,'yyyy-mm-dd'));


			% ----- event data -----
			% main struct
				event_f = [dir(sprintf('%s/data_for_bichan.mat',folder)), dir(sprintf('%s/calca_%s_%s.mat',folder,subject,datestr(session,'yyyy-mm-dd')))];
				if ~isempty(event_f)
					load([event_f.folder filesep event_f.name]);
				else
					data = struct();
				end
			% front / rare - water / grape pair
				if session < datetime(2023,1,1)
					port_f = dir(sprintf('%s/reward_location_*.txt',folder));
					f = fopen([port_f.folder '/' port_f.name],'r');
					C = textscan(f,'%s %s %s\n');
					data.port_is_water = [false false];
					if (strcmpi('Rear',C{1}) && strcmpi('Water',C{3})) || (strcmpi('Front',C{1}) && ~strcmpi('Water',C{3}))
						data.port_is_water(2) = true;
					else
						data.port_is_water(1) = true;
					end
				else
					port_f = dir(sprintf('%s/*_ports.txt',folder));
					f = fopen([port_f.folder '/' port_f.name],'r');
					C = textscan(f,'%s %s %s %s');

					% front novel
					if (strcmpi(C{1},'Novel')&&strcmpi(C{2},'Front,')) || (strcmpi(C{1},'Water')&&strcmpi(C{2},'Rear,'))
						data.port_is_water = [false true];
					% back novel
					else
						data.port_is_water = [true false];
					end
				end
				tmp = [1 0 0; 0 0 0];
				data.port_color = tmp(data.port_is_water+1,:);

				fclose(f);
			% LiCl injection time
				licl_f = dir(sprintf('%s/licl_time*.txt',folder));
				if ~isempty(licl_f)
					f = fopen([licl_f.folder '/' licl_f.name],'r');
					A = fscanf(f,'%d:%d');
					data.licl = A(1)*60 + A(2); % in seconds
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
			data = load_data.spike_ks(data_f,data,ops);



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

			% load amplitude
			% data.amp = load_data.spike_amp(data); 
			data.brain_region = load_data.anatomical(data);

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
		function data = spike_ks(folder, data, ops)
			% load kilosort into data struct format
			% folder  - kilosort results folder
			% data    - output struct 
			% if_good - if selecting only good units 

			% if only keeping good units
			if nargin < 3
				ops = struct();
			end
			if_good = getOr(ops,'if_good',false);

			% load sampling rate
			if exist('pyrunfile')
				SR = double(pyrunfile([folder filesep 'params.py'],'sample_rate'));
			else % default SR
				SR = 3e4;
			end

			% load cluster group
			[data.cids, data.cgs,data.ch, data.sh] = readClusterGroups_bw([folder '/cluster_info.tsv']); % unit group
			data.cids(data.cgs==0) = [];
			data.ch(data.cgs==0)   = [];
			data.sh(data.cgs==0)   = [];
			data.cgs(data.cgs==0)  = [];
			% - 0 = noise
			% - 1 = mua
			% - 2 = good
			% - 3 = unsorted

			% select good units
			if if_good
				data.cids(data.cgs==1) = [];
				data.cgs(data.cgs==1) = [];
				fprintf('good units\n');
			end

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