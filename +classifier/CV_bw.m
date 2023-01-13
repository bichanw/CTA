function [cost,cv] = CV_bw(X,Y,n_fold,Classifier)
% my script for doing cross validation

	% parameters
	if nargin < 3
		n_fold = 10;
	end
	if nargin < 4
		Classifier = classifier.nb();
	end


	% partition
	n_sample = size(X,1);
	c = cvpartition(n_sample,'KFold',n_fold);

	% keep cross-validation information
	cv.ind = nan(n_sample,1); % test set label
	cv.post_label = nan(n_sample,1);


	% validate
	ops = struct('if_cv',false);
	n_mis = 0;
	n_valid = 0;
	for ifold = 1:n_fold

		% save index
		cv.ind(c.test(ifold)) = ifold;

		try
			% train
			mdl = Classifier.train({X(c.training(ifold),:),Y(c.training(ifold),:)},ops);

			% predict
			post = Classifier.predict(mdl,X(c.test(ifold),:));
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