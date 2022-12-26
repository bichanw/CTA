classdef set_fig < handle
% functions setting figure legends, should correspond with plt.m


	% by cell
	methods (Static)
		function ax = FR(data,ax)
			% last reward
			if isfield(data,'rewards')
				t_last_reward = max([data.rewards.all.rear(end) data.rewards.all.front(end)]);
				h = text(ax,t_last_reward, ax.YLim(1), sprintf('last\nreward'),'FontSize',ax.FontSize-1,'horizontalalignment', 'center', 'verticalalignment', 'top'); 
			end

			% figures setting
			xlabel(ax,sprintf('\ntime (s)'));
			ylabel(ax,'FR (Hz)');
		end

		function ax = amp(ax)
			% parameters
			ax.YAxis(2).Label.String = 'spk amplitude';
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