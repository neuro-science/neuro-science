function [t2, t1, ttf, te] = cmp51_ed2t(ed, t, N, nf, nt)
% % % 07/12/2014	written by wp
% % % 	ed([v1, v2, f, t])
% % %		t(f, t, v, v)

	%% prepare
	if nargin < 5
		nt = max(ed(:, 4));
	end
	if nargin < 4
		nf = max(ed(:, 3));
	end
	if nargin < 3
		N = max(max(ed(:, 1:2)));
	end
	t1 = zeros(N, 1);
	t2 = zeros(N, N);
	n = size(ed, 1);
	te = zeros(n, 1);
	ttf = zeros(nf, nt);
	
	%% work
	for k = 1 : n
		te(k) = t(ed(k, 3), ed(k, 4), ed(k, 1), ed(k, 2));
		t1(ed(k, 1)) = t1(ed(k, 1)) + t(ed(k, 3), ed(k, 4), ed(k, 1), ed(k, 2));
		t1(ed(k, 2)) = t1(ed(k, 2)) + t(ed(k, 3), ed(k, 4), ed(k, 1), ed(k, 2));
		t2(ed(k, 1), ed(k, 2)) = t2(ed(k, 1), ed(k, 2)) + t(ed(k, 3), ed(k, 4), ed(k, 1), ed(k, 2));
		t2(ed(k, 2), ed(k, 1)) = t2(ed(k, 2), ed(k, 1)) + t(ed(k, 3), ed(k, 4), ed(k, 1), ed(k, 2));
		ttf(ed(k, 3), ed(k, 4)) = ttf(ed(k, 3), ed(k, 4)) + ...
			t(ed(k, 3), ed(k, 4), ed(k, 1), ed(k, 2));
	end
	
	
	
end % end of function
