function M = cmp19_triangle_route (tri, sid)
% % % updated 05/01/16 by wp : for more general use
% % % updated 08/08/14 by wp : for more general use
% % % 14/05/14 written by wp

	%% recover the nodes
% % % 	nb = cmp20_tri2nb(tri);
	nb = cmp37_neigborSearch(tri, 't');
	n = size(nb, 1);
	
	%% initialize
	if nargin < 2
		sid = true(n, 1);
	end
	
	src = find(sid);
	m = length(src);
	M = zeros(n, m) + inf;	%not connected, distance is inf
	
	parfor c1 = 1 : m
		% % % initialize search candidates
		iFlag = true(n, 1);
		con = M(:, c1);

		% % % kick out itself
		x = src(c1);
		iFlag(x) = 0;
		con(x) = 0;
		
		% % % para before start
		goFlag = true;
		D = 1;
		
		% % % do it in loop
		while goFlag
% 			tic; %debug only
			[con, x, iFlag, goFlag] = inviteMyNeighbors(x, nb, con, iFlag, D);
% 			nl = length(find(iFlag));
% 			fprintf('level %d done! in  %5.2f seconds, %d nodes left.\n', D, toc, nl);
			D = D + 1;
		end
		M(:, c1) = con;
% 		fprintf('%d\n', c1);
	end
end % end of function 0

function [con, y, iFlag, goFlag] = inviteMyNeighbors(x, nb, con, iFlag, D)
	goFlag = true;
	y = unique(cell2mat(nb(x)));
	id = iFlag(y);
	iFlag(y) = 0;
	y = y(id);
	if ~isempty(y)
		con(y) = D;
	else
		goFlag = false;
	end
end %end of function 1