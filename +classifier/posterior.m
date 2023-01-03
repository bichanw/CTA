function [Posterior,t,ops] = posterior(Mdl,count,t,ops)

% load files, by default load 9 cells version
if nargin == 0 || isempty(Mdl)
	load('mat/tmp_mdls.mat');
	Mdl = Mdl{1};
	ops = ops_mdl(1);
end
if nargin < 2 || isempty(count)
	load('mat/run_ave.mat');
end

[~,Posterior,Cost] = predict(Mdl,count');