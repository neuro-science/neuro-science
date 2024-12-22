function diff = mat_compare_difference_matrix(X, Y)
% % % written on 11/08/2017 by wp, compute new xyz with transform matrix
	if ~all(size(X) == size(Y))
		error('they are not even same in size!');
	else
		NNN = numel(X);
		X = X(:);
		Y = Y(:);
		i1 = isnan(X);
		nn1 = numel(find(i1));
		if nn1
			i2 = isnan(Y);
			if all(i1 == i2)
				rn1 = nn1 * 100 ./ NNN;
				X(i1) = [];
				Y(i1) = [];
				fprintf('\n %d[%0.1f%%] NaNs were removed!\n', nn1, rn1);
			else
				error('they are not even same in nans!');
			end
		end
		d1 = (X - Y).^2;
		d2 = (X.^2 + Y.^2);
		idx0 = d2 < 1e-8;
		if any(idx0) && sum(d1(idx0)) < 1e-8
			d1(idx0) = [];
			d2(idx0) = [];
			n1 = length(find(idx0));
			r1 = 100 * n1 ./ NNN;
			fprintf('\n %d[%0.1f%%] zeros were removed!\n', n1, r1);
		end
		diff = mean(sqrt(d1 ./ d2));
		fprintf('The difference is about %5.1f scale!\n\n\n', log10(diff));
	end

end %end of function