classdef nb < handle
% master class for naive bayes classifier


	% by cell
	methods (Static)

		function [Mdl, ops] = train(count,ops)
			% train naive bayes classifier, with preprocessing steps
			% Input: 	count - nclass * 1 cells, with each cell being ntrs * nneurons matrix
			% 			ops - classifier parameters, rescale, if_save, exclude_id


			if nargin < 2
				ops = struct;
			end
			ops.rescale = getOr(ops,'rescale',1);
			if_save = getOr(ops,'if_save',false);


			% include certain cells and rescale, default include all cells
			% include_id = [9 13 28 29 30 34 46 54 67]; ops.exclude_id = getOr(ops,'exclude_id',~ismember(1:size(count{1},1),include_id));
			ops.exclude_id = getOr(ops,'exclude_id',false(1,size(count{1},1)));

			% remove cells with 0 variance
			v = catcell(cellfun(@(x) std(x,[],2), count,'UniformOutput',false));
			ops.exclude_id = ops.exclude_id | (prod(v,2)==0)';
			data_nb = cellfun(@(x) x(~ops.exclude_id,:) / ops.rescale, count,'UniformOutput',false);

			% create table and run classifier
			X = catcell(data_nb,2)'; % 90 cells * 100 trials response
			Y = catcell(arrayfun(@(i) i*ones(size(data_nb{i},2),1),1:numel(data_nb),'UniformOutput',false),1);
			Mdl = fitcnb(X,Y);


			% save file
			if if_save save('mat/tmp_mdl.mat','Mdl','ops'); end

			% cross validation by me
			[cost,cv] = classifier.nb.CV_bw(Mdl.X,Mdl.Y);
			save tmp cost cv
			fprintf('Cross validation loss: %.2f\n',cost);

			% confusion matrix
			% [isLabels,score] = resubPredict(Mdl);
			% ConfusionMat = confusionchart(Y,isLabels);

			% predict posterior
			% [label,Posterior,Cost] = predict(Mdl,X);
		end


		function [cost,cv] = CV_bw(X,Y,n_fold)

			% parameters
			if nargin < 3
				n_fold = 10;
			end

			% partition
			n_sample = size(X,1);
			ind = repmat([1:n_fold],1,n_sample/n_fold);
			ind = ind(randperm(numel(ind)));


			% keep cross-validation information
			cv.ind = ind; % test set label
			cv.post_label = nan(n_sample,1);


			% validate
			n_mis = 0;
			n_valid = 0;
			for ifold = 1:n_fold

				try
					% train
					mdl = fitcnb(X(ind~=ifold,:),Y(ind~=ifold,:));

					% predict
					post = classifier.posterior_bw(mdl,X(ind==ifold,:));
					[~,post_label] = max(post,[],2);

				catch ME
					% if zero variance, continue to the next partition
					save('error.mat','ME');
					continue
				end

				% keep label
				cv.post_label(ind==ifold) = post_label;

				% performance
				n_valid = n_valid + n_sample / n_fold;
				n_mis   = n_mis + sum(post_label~=Y(ind==ifold));
			end

			% no valid solution across any partition
			fprintf('%d valid prediction out of %d\n',sum(~isnan(cv.post_label)), numel(cv.post_label));
			if sum(~isnan(cv.post_label))==0
				fprintf('Warning: no valid solution across any partition for cross-validation\n');
			end

			cost = n_mis / n_valid;

		end
	end
	

end