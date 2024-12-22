function nb = meg_neighborhood_redefine(nb0, idx)
	% % % 	written 23/01/2018
	% % % sort input indices to ordered numbers	
	if islogical(idx)
		idx = find(idx);
	else
		idx = unique(sort(idx));
	end
	n = numel(idx);
	% % % remove entries
	nb1 = nb0(idx);
	% % % check each element
	nb = cell(n, 1);
	for k = 1 : n
		[tmp1, tmp2, nb{k}] = intersect(nb1{k}, idx);
	end
	% % % 	clean up
	clear tmp1 tmp2 nb1 n idx k;
		
end
