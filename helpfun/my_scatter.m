function h = my_scatter(x,y,ax,varargin)
	if nargin<3 || isempty(ax)
		ax = np;
	end

	% repeat y
	if numel(y)<=1 
		if isempty(y) y = 1; end
		y = y*ones(size(x));

	% repeat x
	else
		if isempty(x) x = 1; end
		x = x*ones(size(y));
		
	end

	h = scatter(ax, x,y,varargin{:});
end