classdef mnr < handle
% master class for multinomial logistics classifier


	% by cell
	methods (Static)


		function [Posterior,LL] = predict(Mdl,X,ops,Y_true)

			% initiation
			if nargin < 3
				ops = struct;
			end

			if 0 % matlab
				Posterior = mnrval(Mdl.B,X);
			else % python
				% if scaled data when training model, scale again
				if isfield(ops.mnr,'scale') && ~isfield(ops,'zscore_by_time')
					X = (X - ops.mnr.scale.M)./ ops.mnr.scale.V;
				end

				% predict
				tmp = X * Mdl.coef + Mdl.intercept;
				% stack for binary decision
				if size(Mdl.coef,2) == 1
					tmp = [-tmp, tmp];
				end
				Posterior = exp(tmp) ./ sum(exp(tmp),2);
			end


			% calculate log likelihood
			if nargout > 1 
				% check true label
				if nargin < 4
					error('Provide true label to calculate LL');
				end

				% calculate LL
				LL = sum(log(Posterior(arrayfun(@(ii) sub2ind(size(Posterior),ii,Y_true(ii)), 1:numel(Y_true)))));
			end
		end

		function [Mdl, ops] = train(data,ops)
			% train naive bayes classifier, with preprocessing steps
			% Input: 	data - data struct
			% 				 - or cell wit {X, Y}
			% 			ops - classifier parameters, rescale, if_save, exclude_id

			if nargin < 2
				ops = struct;
			end
			if_cv = getOr(ops,'if_cv',false);

			% convert data to spike count if input is raw data struct
			if isstruct(data)
				[X,Y,ops] = classifier.data_2_XY(data,ops);
			elseif iscell(data)
				X = data{1};
				Y = data{2};
			end
			ops.classifier = classifier.mnr();

			% parameter initation
			ops.exclude_id = getOr(ops,'exclude_id',false(size(X,2),1));
			ops.mnr = getOr(ops,'mnr',struct());
			ops.mnr.penalty = getOr(ops.mnr,'penalty','l2');
			ops.mnr.lambda  = getOr(ops.mnr,'lambda',1);
			ops.mnr.zscore  = getOr(ops.mnr,'zscore',true);

			% zscore for preprocessing
			if ops.mnr.zscore
				% calculate mean and variance
				if ~isfield(ops.mnr,'scale')
					ops.mnr.scale.M = mean(X,1); ops.mnr.scale.V = std(X,[],1);
				end

				% remove nan
				ind = (ops.mnr.scale.V==0);
				X(:,ind) = [];
				ops.mnr.scale.M(ind) = [];
				ops.mnr.scale.V(ind) = [];
				tmp = find(~ops.exclude_id); ops.exclude_id(tmp(ind)) = true;
				ops.decoder_id = find(~ops.exclude_id);

				% zscore
				X = (X - ops.mnr.scale.M)./ ops.mnr.scale.V;
			end


			% train model
			if 0 % using matlab built in function
				[Mdl.B,Mdl.dev,Mdl.stats] = mnrfit(X,Y);
			else % using python
				id = sprintf('%.3f',rand());
				save(sprintf('mat/tmp_%s.mat',id),'X','Y','ops');
				system(['python +classifier/mnr.py ' id]);
				load(['mat/pythonsave_' id '.mat']);
				Mdl = struct('coef',coef,'intercept',intercept);
				system(['rm mat/tmp_' id '.mat mat/pythonsave_' id '.mat']);
			end
			ops.Mdl = Mdl;

			% cross validation by me
			if if_cv
				[cost,cv] = classifier.CV_bw(X,Y,10,ops.classifier);
				fprintf('Cross validation loss: %.2f\n',cost);
			end

			% confusion matrix
			% [isLabels,score] = resubPredict(Mdl);
			% ConfusionMat = confusionchart(Y,isLabels);

			% predict posterior
			% [label,Posterior,Cost] = predict(Mdl,X);
		end


	end
	

end