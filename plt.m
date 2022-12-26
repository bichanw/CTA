classdef plt < handle
% functions making plots


	% by cell
	methods (Static)
		function [ax,ops] = FR(data,ii,ax)
			% bin spike
			ops.bin_width = 100; 
			ops.step_size = 100; % bin and step

			% ----- plotting ----- %
			% plot binned data
			[~,spk_count,~,tbin] = running_average(data.spikes{ii},[],ops.bin_width,ops.step_size);
			plot(ax,tbin,spk_count/ops.bin_width);

			% mark events
			% last reward
			if isfield(data,'rewards')
				t_last_reward = max([data.rewards.all.rear(end) data.rewards.all.front(end)]);
				plot(ax,t_last_reward*[1 1],ax.YLim,'k--'); % last reward
			end

			% figures setting
			% h = text(ax,t_last_reward, ax.YLim(1), sprintf('last\nreward'),'FontSize',ax.FontSize-1,'horizontalalignment', 'center', 'verticalalignment', 'top'); 
			% xlabel(ax,'time (s)');
		end

		function ax = amp(data,ii,ax)
			% parameters
			max_2_plt = 500; % number of maximum spikes displayed
			ind = randperm(numel(data.spikes{ii}),min([max_2_plt numel(data.spikes{ii})]));

			% plot
			scatter(ax,data.spikes{ii}(ind),data.amp{ii}(ind),'MarkerEdgeAlpha',0.1);
		end

		
	end



	methods	(Static)


		function ax = FR_all(data,cells_2_plt)
			% load help function...
			addpath('/mnt/cup/people/bichanw/SpikeSorting/Codes/bucket/helpfun');
			
			% initiation
			if nargin<2
				cells_2_plt = 1:numel(data.spikes);
			end
			N = numel(cells_2_plt);
			t_last_reward = max([data.rewards.all.rear(end) data.rewards.all.front(end)]);

			% bin spike
			bin_width = 100; 
			step_size = 100; % bin and step

			% ----- plotting ----- %
			[ax,r,c] = np(N);
			for ii = 1:N
				% plot binned data
				[~,spk_count,~,tbin] = running_average(data.spikes{cells_2_plt(ii)},[],bin_width,step_size);
				plot(ax(ii),tbin,spk_count);

				% mark events
				% last reward
				if isfield(data,'rewards')
					plot(ax(ii),t_last_reward*[1 1],ax(ii).YLim,'k--'); % last reward
				end
			end

			% figures setting
			ind = sub2ind([c r],1,r);
			h = text(ax(ind),t_last_reward, ax(ind).YLim(1), sprintf('last\nreward'),'FontSize',ax(ind).FontSize,'horizontalalignment', 'center', 'verticalalignment', 'top'); 
			xlabel(ax(ind),'time (s)');
			ef;
		end



	end

	

end