function w2 = cmp48_OmegaSquaredEqualN(data, trl)

	[SS_total, SS_treat, SS_error, N, M] = cmp47_SumOfSquaresEqualN(data, trl);
	MSE = SS_error /(N * M);
	w2 = (SS_treat - (M - 1) .* MSE) ./ (SS_total + MSE);
end

