function vout = plt06_gaussian_topo (data, dThresh, v1, v2)
% % % updated 08/08/14 by wp : topological interpolation
% % % updated 20/6/14 wp free v1 space
% % % written 14/5/14 wp
% % % To interpolate data for surface graph based on topological connection
% % % 3 conditions as below

	%% get distance in various situations
	if nargin < 4	%3rd input is matrix numVert x numSrc
		dist2 = v1 .^ 2;
	elseif size(v2, 2) == 3 %v2 is tri, v1 is index of sources in vertices
		M = cmp19_triangle_route(v2);
		dist2 = M(:, v1) .^ 2;
	else %v2 is distance matrix for all, v1 is index of sources in vertices
		dist2 = v2(:, v1) .^ 2;
	end
	v1 = [];
	
	% % % compute the transform matrix
	alpha2 = (dThresh / sqrt(log(2))) .^ 2;
	w = exp(-dist2 / alpha2);

	% % % obtain data	
	vout = bsxfun(@rdivide, w * data, sum(w, 2));

end %end of function
