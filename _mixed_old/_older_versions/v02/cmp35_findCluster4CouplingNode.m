function [c, nc] = cmp35_findCluster4CouplingNode(con)

% % %   17/09/2014: updated by wp fpr better compatibility
	

	sz =size(con);
	for ch1 = 1 : sz(4) - 1
		con(:, :, ch1, ch1 : sz(3)) = 0;
	end
	[x1, x2, x3, x4] = ind2sub(sz, find(con));
	nds = length(x1);
	fg = true(nds, 1);
	ids = [x3, x4, x1, x2];
	nc = 0;
	c = [];
	for ig = 1 : nds
		if fg(ig)
			cid = [];
			[cid, fg] = searchClusterLoopSparse(ids, ig, cid, fg);
			if ~isempty(cid)
				nc = nc + 1;
				c{nc}.ne = length(cid);
				c{nc}.ed = ids(cid, :);
				c{nc}.nd = unique(ids(cid, 1:2));
				c{nc}.nn = length(c{nc}.nd);
				c{nc}.id = sub2ind(sz, x1(cid), x2(cid), x3(cid), x4(cid));
			end
		end
	end
end

function [y, flag] = searchClusterLoopSparse(c4, x, y, flag)

	ng = length(x);
	x2 = [];
	for k = 1 : ng
		if flag(x(k))
			flag(x(k)) = 0;
			y = [y; x(k)];
			s1 = find(sum(bsxfun(@eq, c4, c4(x(k), :)), 2) >= 3);
			s2 = find(sum(bsxfun(@eq, c4, c4(x(k), [2 1 3 4])), 2) >= 3);
			s = unique([s1; s2]);
			s(s==x(k)) = [];
			d = sum(abs(bsxfun(@minus, c4(s, 3:4), c4(x(k), 3:4))), 2) > 1;
			s(d) = [];
			x2 = [x2; s];
		end
	end
	x2 = unique(x2);
	if ~isempty(x2)
		[y, flag] = searchClusterLoopSparse(c4, x2, y, flag);
	end
end
