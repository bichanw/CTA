ops = struct('if_cv',false,'classifier',classifier.mnr());  % classifier

% multinomial logistic regression settings
ops.mnr.lambda = 1;
ops.mnr.penalty = 'l2';