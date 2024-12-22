function u = meg_neighbor_distance (nb, pr, border)
% % % writen by wp 23/07/2019
% % % nb - neighborhood definition
% % % pr - two columns of node pairs
	if nargin < 3
		border = 0;
	end
	if size(pr, 2) > 2 && size(pr, 1) == 2
		fprintf('\nIt seems the node pairs are in wrong dimensions, reverted!!!!!\n')
	end
	for ie = size(pr, 1) : -1 : 1
		if border && (pr(ie, 1) - border) * (pr(ie, 2) - border) < 0
			u(ie) = nan;
		else
			u(ie) = isNeighbor(nb, pr(ie, 1), pr(ie, 2));
		end
	end
end


function n = isNeighbor(nb, x, y, n)
	if nargin < 4
		n = 0;
	end
	X = unique(cat(1, nb{x}));
	if ismember(y, X)
		n = n + 1;
	else
		n = isNeighbor(nb, X, y, n + 1);
	end
end	