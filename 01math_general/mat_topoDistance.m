function n = mat_topoDistance(nb)
% % % written 07/08/2018
% % % distance matrix among nodes in a topological structure
	nv = length(nb);
	n = zeros(nv) + nan;
	for iv = 1 : nv
		flag = true(nv, 1);
		n(iv, iv) = 0;
		flag(iv) = false;
		n(:, iv) = theTopD(n(:, iv), nb, flag, nb{iv}, 1);
	end
	% % % check contras between cons or sti_vs_bas
end	

function d = theTopD(d, nb, flag, nd, lv)
	nd = intersect(nd, find(flag));
	d(nd) = lv;
	flag(nd) = false;
	if any(flag)
		lv = lv + 1;
		nd = unique(cat(1, nb{nd}));
		d = theTopD(d, nb, flag, nd, lv);
	end
end
