classdef mnr < handle
% master class for multinomial logistics classifier


	% by cell
	methods (Static)


		function Posterior = predict(Mdl,X)
			Posterior = mnrval(Mdl.B,X);
		end

		function [Mdl, ops] = train(data,ops)
			% train naive bayes classifier, with preprocessing steps
			% Input: 	data - data struct
			% 				 - or cell wit {X, Y}
			% 			ops - classifier parameters, rescale, if_save, exclude_id

			if nargin < 2
				ops = struct;
			end
			if_cv = getOr(ops,'if_cv',true);

			% convert data to spike count if input is raw data struct
			if isstruct(data)
				[X,Y,ops] = classifier.data_2_XY(data,ops);
			elseif iscell(data)
				X = data{1};
				Y = data{2};
			end
			ops.classifier = classifier.mnr();
			[Mdl.B,Mdl.dev,Mdl.stats] = mnrfit(X,Y);

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