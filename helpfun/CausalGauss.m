function v = CausalGauss(kernel_width)
n_points = ceil(3*kernel_width);
v = 1:n_points;
v = exp(-v.^2 / (2*kernel_width^2)) ./ (kernel_width * sqrt(2*pi));
v = v./sum(v);
end