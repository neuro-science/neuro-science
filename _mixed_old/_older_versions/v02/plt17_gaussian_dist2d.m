function vout = plt17_gaussian_dist2d (vc2, vc_id, src_id, dThresh, data)
% % % written 04/01/2016 by wp

	% % % check input	
	if nargin > 6 || nargin < 4
		error('input shall be 3 (vertices_2d, 2d indices, source indices)\n or 4(+ distance threshold) or 5 (+data)!');
	elseif nargin == 3
		dThresh = 2;
	end
	
	% % % generate the 2d points
	n = length(vc_id);
	vc = zeros(n, 2) + nan;
	vc(vc_id, :) = vc2; %2d locations of all vertices
	vc0 = vc(src_id, :); %2d locations of sources
	
	% % % compute the transform matrix	
	alpha2 = (dThresh / sqrt(log(2))) .^ 2;
	dist2 = sum(bsxfun(@minus, permute(vc, [1 3 2]), permute(vc0, [3 1 2])).^2, 3);
	w = exp(-dist2 / alpha2);
	w(isnan(w)) = 0;
	
	% % % decide output
	if nargin == 4
		vout = w;
	else
		vout = bsxfun(@rdivide, w * data, sum(w, 2));
	end

end %end of function

    