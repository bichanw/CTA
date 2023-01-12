% do not exclude cell
ops = struct('exclude_id',false(1,numel(data.spikes)));

% run classifier with individual cell and exclude bad ones
[cost,valid] = classifier.nb.individual_cv(data,ops);
ops.exclude_id = ops.exclude_id | ~(cost'<0.6);