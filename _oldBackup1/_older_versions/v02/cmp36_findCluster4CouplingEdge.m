function [c, nc] = cmp36_findCluster4CouplingEdge(con, nb)


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
			[cid, fg] = searchClusterLoopSparseNB(ids, ig, cid, fg, nb);
			if ~isempty(cid)
				nc = nc + 1;
				c{nc}.ne = length(cid);
				c{nc}.ed = ids(cid, :);
				c{nc}.nd = unique(ids(cid, 1:2));
				c{nc}.nn = length(c{nc}.nd);
				c{nc}.id = sub2ind(sz, ids(cid, 3), ids(cid, 4), ids(cid, 1), ids(cid, 2));
			end
		end
	end
end

function [y, flag] = searchClusterLoopSparseNB(c4, x, y, flag, nb)

	ng = length(x);
	x2 = [];
	for k = 1 : ng
		if flag(x(k))
			flag(x(k)) = 0;
			y = [y; x(k)];
			s1 = find(sum(bsxfun(@eq, c4(:, 3:4), c4(x(k), 3:4)), 2) == 2);
			if ~isempty(s1)
				for k1 = 1 : 2
					for k2 = 1 : 2
						s2 = find(ismember(c4(s1, k2), c4(x(k), k1)));
						if ~isempty(s2)
							s3 = find(ismember(c4(s1(s2), 3 - k2), nb{c4(x(k), 3 - k1)}));
							if ~isempty(s3)
								x2 = [x2; s1(s2(s3))];
							end
						end
					end
				end
			end
			clear s1 s2 s3 k1 k2;
			for k1 = 1 : 2
				id = [1, 2, k1 + 2];
				s1 = find(sum(bsxfun(@eq, c4(:, id), c4(x(k), id)), 2) == 3);
				if ~isempty(s1)
					s2 = abs(bsxfun(@minus, c4(s1, 5 - k1), c4(x(k), 5 - k1))) == 1;
					if ~isempty(s2)
						x2 = [x2; s1(s2)];
					end
				end
			end
			clear s1 s2 k1;
		end
	end
	x2 = unique(x2);
	if ~isempty(x2)
		[y, flag] = searchClusterLoopSparseNB(c4, x2, y, flag, nb);
	end
end
