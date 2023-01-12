function saveops(ops)

ops_save = ops;
save(sprintf('mat/ops/%s.mat',datestr(now,'YYmmdd_hhMM')),'ops_save');