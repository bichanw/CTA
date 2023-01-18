classdef nb < handle
% master class for naive bayes classifier


	% by cell
	methods (Static)
		function Posterior = predict(Mdl,X)
			[label,Posterior,Cost] = predict(Mdl,X);
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


			% run classifier
			ops.classifier = classifier.nb();
			Mdl = fitcnb(X,Y);

			% cross validation by me
			if if_cv
				[ops.cv_result.cost,ops.cv_result.cv] = classifier.CV_bw(Mdl.X,Mdl.Y,10,ops.classifier);
				fprintf('Cross validation loss: %.2f\n',ops.cv_result.cost);
			end

			% confusion matrix
			% [isLabels,score] = resubPredict(Mdl);
			% ConfusionMat = confusionchart(Y,isLabels);

			% predict posterior
			% [label,Posterior,Cost] = predict(Mdl,X);
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