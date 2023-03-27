classdef common_t < handle

	methods (Static)


		function t = last_reward(data)
			t = max([data.rewards.all.rear(end) data.rewards.all.front(end)]);
		end


		function t = first_laser(data)
			
			if isfield(data,'laser')
				t = data.laser(1,1);
			else
				error('No laser in data');
			end

		end

	end

end