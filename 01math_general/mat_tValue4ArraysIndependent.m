	function t = mat_tValue4ArraysIndependent(data, trl, dim)
	% % % 26/03/18	written by wp, force squeezed

		% % % check paras
		if nargin < 3 || isempty(dim) %default last dim
			dim = length(size(data));
		end

		% % % get the classification
		sz = size(data);
		szn = 1:length(sz);
		n = length(trl);
		if sz(dim) ~= n
			error('data and case sizes mismatch!\n');
		end
		[x1, x2] = unique(sort(trl));
		n1 = x2(2) - 1;
		n2 = n - x2(2) + 1;
		
		% % % work for two cons 
		if length(x1) ~= 2
			error('The data were not in two groups!\n');
		else
			szn(dim) = [];
			sz(dim) = [];
			data = reshape(permute(data, [dim, szn]), [n, prod(sz)]);
			d1 = data(1:n1, :);
			d2 = data(x2(2):end, :);
			m1 = mean(d1, 1);
			m2 = mean(d2, 1);
			s1 = var(d1, 0, 1);	%note: here is not std, thus no square needed below
			s2 = var(d2, 0, 1);
			sp = sqrt((s1 * (n1 - 1) + s2 * (n2 - 1)) / (n - 2) * (1./n1 + 1./n2));
			t = (m1 - m2) ./ sp;
		end
		
		% % % reshape
		if numel(sz) < 2
			sz = [sz, 1];
		end
		t = reshape(t, sz);

	end
