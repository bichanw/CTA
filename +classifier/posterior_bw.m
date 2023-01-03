function [Posterior,p_ll] = posterior_bw(Mdl,X)

% check input dimension
if size(X,2) ~= size(Mdl.DistributionParameters,2)
	error('Column number of input matrix X is different from the number of neurons.');
end

% calculating likelihood
p_ll = nan([size(X,1) size(Mdl.DistributionParameters)]);
for isample = 1:size(X,1)
	for icell = 1:size(Mdl.DistributionParameters,2)
		% p_ll(isample,:,icell) = arrayfun(@(icond) log(normpdf(X(isample,icell),Mdl.DistributionParameters{icond,icell}(1),Mdl.DistributionParameters{icond,icell}(2))), 1:3);
	
		% try writing in for loop see if it speeds up
		for icond = 1:size(Mdl.DistributionParameters,1)
			p_ll(isample,icond,icell) = log(normpdf(X(isample,icell),Mdl.DistributionParameters{icond,icell}(1),Mdl.DistributionParameters{icond,icell}(2)));
		end

	end
end

% calculate posterior
Posterior = sum(p_ll,3) + log(Mdl.Prior);
Posterior = exp(Posterior) ./ sum(exp(Posterior),2);