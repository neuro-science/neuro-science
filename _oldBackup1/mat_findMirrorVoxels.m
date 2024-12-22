function [ids, d1] = mat_findMirrorVoxels (V3)
% % % written 06/07/2017 by wp
% % % V3 is a Nx3 array of x y z coordinates for N points
% % % The function will find out the mirroring points for each point

	% % % get the mirror array
	Vm = V3;
	Vm(:, 1) = -Vm(:, 1);
	
	% % % compute distance
	d = sum(bsxfun(@minus, permute(V3, [1 3 2]), permute(Vm, [3 1 2])).^2, 3);
	
	% % % find the pairs
	[d1, ids] = min(d);
end
