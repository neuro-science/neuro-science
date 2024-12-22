
function [tr, n] = cmp39_select_trials(tr0)
% % % This function select equal number of trials from a trl sets
% % % written in 03/12/14 by wp

	nn = length(tr0);
	ntrs = zeros(nn, 1);
	for k = 1 : nn
		sz = size(tr0{k});
		ntrs(k) = prod(sz);
		if sz(1) == ntrs(k)
			fprintf('Sequence #%d is a vector!\n', k);
		elseif sz(2) == ntrs(k)
			fprintf('Sequence #%d is not a vector, \n converted!\n', k);
			tr0{k} = tr0{k}';
		else
			error('Sequence #%d is not a vector, \n cannot be converted!\n');
		end
	end
	
	[n, I] = min(ntrs);
	tr1 = tr0{I};
	tr = zeros(n, nn);
	clear I;
	
	for k = 1 : nn
		if ntrs(k) == n
			tr(:, k) = tr0{k};
		else
			tmp1 = bsxfun(@minus, tr0{k}, tr1');
			[tmp2, tmp3] = sort(sum(tmp1.^2, 2));
			tr(:, k) = sort(tr0{k}(tmp3(1 : n)));
		end
		clear tmp1 tmp2 tmp3;
	end
end % end of function