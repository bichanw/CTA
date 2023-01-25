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
			% Licl
			if isfield(data,'licl')
				plot(ax,data.licl*[1 1],ax.YLim,'k--'); % last reward
			end

			% figures setting
			% h = text(ax,t_last_reward, ax.YLim(1), sprintf('last\nreward'),'FontSize',ax.FontSize-1,'horizontalalignment', 'center', 'verticalalignment', 'top'); 
			% xlabel(ax,'time (s)');
		end

		function ax = ampimg(data,ii,ax)
			% ii = 1;
			% ax = np; scatter(ax,data.spikes{ii},data.amp{ii}); ef;

			% ax = np; 
			histogram2(ax,data.spikes{ii},data.amp{ii},'EdgeColor','none','DisplayStyle','tile','ShowEmptyBins','on');
			colorbar(ax);
			% ef;
		
		end

		function ax = amp(data,ii,ax)
			% parameters
			max_2_plt = 500; % number of maximum spikes displayed
			ind = randperm(numel(data.spikes{ii}),min([max_2_plt numel(data.spikes{ii})]));

			% plot
			scatter(ax,data.spikes{ii}(ind),data.amp{ii}(ind),'MarkerEdgeAlpha',0.1);
		end


		
	end


	% combined plot
	methods	(Static)


		function ax = FR_amp(data,ii,ax)
			yyaxis(ax,'right'); plt.amp(data,ii,ax); 
			yyaxis(ax,'left');  plt.FR(data,ii,ax);
		end



	end

	

end