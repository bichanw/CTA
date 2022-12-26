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
			% ax.YAxis(2).Label.String = 'spk amplitude';
			ylabel(ax,'spk amplitude');
		end

		
	end



	methods	(Static)


		function ax = FR_amp(data,ax)
			yyaxis(ax,'left'); set_fig.FR(data,ax);
			yyaxis(ax,'right'); set_fig.amp(ax);
		end



	end

	

end