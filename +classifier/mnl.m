classdef mnl < handle
% master class for multinomial logistics classifier


	% by cell
	methods (Static)

		function [Mdl, ops] = train(data,ops)
			% train naive bayes classifier, with preprocessing steps
			% Input: 	data - data struct
			% 				   alternatively it could be ncategories * 1 cell
			% 			ops - classifier parameters, rescale, if_save, exclude_id

			if nargin < 2
				ops = struct;
			end
			if_save = getOr(ops,'if_save',false);

			% convert data to spike count if input is raw data struct
			if ~isdouble(data)
				[X,Y] = classifier.data_2_XY(data,ops);
			end
			[Mdl.B,Mdl.dev,Mdl.stats] = mnrfit(X,Y);

			% cross validation by me
			[cost,cv] = classifier.nb.CV_bw(Mdl.X,Mdl.Y);
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
			% ind = repmat([1:n_fold],1,n_sample/n_fold);
			% ind = ind(randperm(numel(ind)));

			c = cvpartition(n_sample,'KFold',n_fold);

			% keep cross-validation information
			cv.ind = nan(n_sample,1); % test set label
			cv.post_label = nan(n_sample,1);


			% validate
			n_mis = 0;
			n_valid = 0;
			for ifold = 1:n_fold

				% save index
				cv.ind(c.test(ifold)) = ifold;

				try
					% train
					mdl = fitcnb(X(c.training(ifold),:),Y(c.training(ifold),:));

					% predict
					post = classifier.posterior_bw(mdl,X(c.test(ifold),:));
					[~,post_label] = max(post,[],2);

				catch ME
					% if zero variance, continue to the next partition
					save('error.mat','ME');
					continue
				end

				% keep label
				cv.post_label(c.test(ifold)) = post_label;

				% performance
				n_valid = n_valid + sum(sum(c.test(ifold)));
				n_mis   = n_mis + sum(post_label~=Y(c.test(ifold)));
			end

			% no valid solution across any partition
			fprintf('%d valid prediction out of %d\n',sum(~isnan(cv.post_label)), numel(cv.post_label));
			if sum(~isnan(cv.post_label))==0
				fprintf('Warning: no valid solution across any partition for cross-validation\n');
			end

			cost = n_mis / n_valid;

		end

		function [cost,valid] = individual_cv(data,ops,if_plot)

			% initiation
			if nargin < 3
				if_plot = false;
			end

			% naive bayes cross validation by cv
			f = sprintf('mat/individual_cv_%s_%s.mat',data.subject,datestr(data.session,'YYmmdd'));
			if isempty(dir(f))
				[spk_count,ops] = classifier.count_spk.events(data,ops);

				cost = NaN(numel(data.spikes),1);
				valid = cost;
				for ii = 1:numel(data.spikes)
					tmp = cellfun(@(x) x(ii,:), spk_count,'UniformOutput',false);
					X = catcell(tmp,2)'; % 90 cells * 100 trials response
					Y = catcell(arrayfun(@(i) i*ones(size(tmp{i},2),1),1:numel(tmp),'UniformOutput',false),1);
					[cost(ii),cv(ii)] = classifier.nb.CV_bw(X,Y);
					valid(ii) = sum(~isnan(cv(ii).post_label))/ numel(cv(ii).post_label);
				end

				% save result
				save(f,'cost','valid','cv');
			else
				load(f);
			end

			% plot for inspection
			% represent invalid cell as -1
			if if_plot
				cost_plt = cost;
				cost_plt(isnan(cost_plt)) = -0.1;

				ax = np(1,2);
				histogram(ax(1),cost_plt);
				histogram2(ax(2),valid,cost_plt,'XBinEdges',(0:0.1:1.1)-0.05,'YBinEdges',-0.1:0.05:1,'DisplayStyle','tile','ShowEmptyBins','on','EdgeColor','none');
				
				% figure setting
				c = colorbar(ax(2));
				xlabel(ax(1),'cost'); ylabel(ax(1),'# cells');
				xlabel(ax(2),'% valid sample'); ylabel(ax(2),'cost');
				export_fig(sprintf('results/cv_cells_%s_%s.pdf',data.subject,datestr(data.session,'YYmmdd')));
			end
		end
	end
	

end