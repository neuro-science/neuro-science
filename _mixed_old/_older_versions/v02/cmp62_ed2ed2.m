function [spd1, spd2, ed2] = cmp62_ed2ed2(ed, rID, lbs)
% % % written  05/02/16 by wp: edges to edges of rois

	%% prepare
	if nargin < 3
		error('edges rIDs and labels are needed!');
	end

	%% do in look
	ed2 = zeros(size(ed));
	for k = 1 : size(ed, 1)
		[tmp1, tmp2] = ismember([rID(ed(k, 1)), rID(ed(k, 2))], lbs);
		ed2(k, :) = [tmp2, ed(k, 3:4)];
		clear tmp1 tmp2;
	end
	
	[spd1, spd2] = cluster_SP2(ed2, length(lbs));
end

function [spd1, spd2] = cluster_SP2(ed, nvxs)
	tmp = ed(:, 1:2);
	spd2 = zeros(nvxs);
	for k =  1 : size(tmp, 1)
		spd2(tmp(k, 1), tmp(k, 2)) = spd2(tmp(k, 1), tmp(k, 2)) + 1;
		spd2(tmp(k, 2), tmp(k, 1)) = spd2(tmp(k, 2), tmp(k, 1)) + 1;
	end
% 	spd2 = spd2 / 2; %- modified 04/01/2016
	spd1 = sum(spd2, 2);
	clear tmp;
end %end of function
