function [c, nc] = cmp28_findCluster4PowerNeighborFree (binArray)

	sz =size(binArray);
	[x1, x2, x3] = ind2sub(sz, find(binArray));
	nds = length(x1);
	flag = true(nds, 1);
	ids = [x3, x1, x2];
	nc = 0;
	c = [];
	for ig = 1 : nds
		if flag(ig)
			cid = [];
			[cid, flag] = searchClusterLoopSparse(ids, ig, cid, flag);
			if ~isempty(cid)
				nc = nc + 1;
				c{nc}.ne = length(cid);
				c{nc}.ed = ids(cid, :);
				c{nc}.nd = unique(ids(cid, 1));
				c{nc}.nn = length(c{nc}.nd);
				c{nc}.id = sub2ind(sz, x1(cid), x2(cid), x3(cid));
			end
		end
	end
end

function [y, flag] = searchClusterLoopSparse(ids, x, y, flag)
	ng = length(x);
	x2 = [];
	for k = 1 : ng
		if flag(x(k))
			flag(x(k)) = 0;
			y = [y; x(k)];
			s = find(sum(bsxfun(@eq, ids, ids(x(k), :)), 2) >= 2);
			s(s==x(k)) = [];
			d = sum(abs(bsxfun(@minus, ids(s, 2:3), ids(x(k), 2:3))), 2) > 1;
			s(d) = [];
			x2 = [x2; s];
		end
	end
	x2 = unique(x2);
	if ~isempty(x2)
		[y, flag] = searchClusterLoopSparse(ids, x2, y, flag);
	end
end
