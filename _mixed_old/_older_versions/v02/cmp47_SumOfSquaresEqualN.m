function [SS_total, SS_treat, SS_error, N, M] = cmp47_SumOfSquaresEqualN(data, trl )

	[N, M] = size(trl);	%[trl num per con, num of cons]
	sz = size(data); %[trl, ...]
	data = reshape(data, [sz(1), prod(sz(2:end))]);
	d0 = zeros([N, M, prod(sz(2:end))]);
	
	for im = 1 : M
		d0(:, im, :) = data(trl(:, im), :);
	end
	
	mt = mean(d0, 1);
	m0 = mean(mt, 2);
	
	SS_total = reshape(sum(sum(bsxfun(@minus, d0, m0).^2, 1), 2), sz(2 : end));
	SS_treat = reshape(N * sum(bsxfun(@minus, mt, m0).^2, 2), sz(2 : end));
	SS_error = reshape(sum(sum(bsxfun(@minus, d0, mt).^2, 1), 2), sz(2 : end));
end

