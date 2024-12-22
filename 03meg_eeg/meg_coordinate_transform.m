% % % apply transformation matrix on 3d location data
function out = meg_coordinate_transform(T, in)
	% % % updated on 15/11/2023 by wp, such that the input shall be 3 x n
	% % % written on 11/08/2017 by wp, as mat_matrix_transform_3d, compute new xyz with transform matrix

	%% 1. check inputs and initialize
	if nargin ~= 2 || any(size(T) - [4 4]) || numel(size(in)) ~= 2
		error('We need two inputs: Transformation Matrix(4x4) and input data(3xN or 6xN)');
	end
	
	if size(in, 2) == 3 && size(in, 1) ~= 3 && size(in, 1) ~= 6
		in_loc = [in ones(size(in, 1), 1)]';
		oriFlag = false;
		fprintf('I assume you put it as Nx3, now converted to 3xN !!!\n');
	elseif size(in, 2) == 6 && size(in, 1) ~= 3 && size(in, 1) ~= 6
		in_loc = [in(:, 1:3) ones(size(in, 1), 1)]';
		in_ori = in(:, 4:6)';
		oriFlag = true;
		fprintf('I assume you put it as Nx6, now converted to 6xN !!!\n');
	elseif size(in, 1) == 3
		in_loc = [in; ones(1, size(in, 2))];
		oriFlag = false;
	elseif size(in, 1) == 6
		in_loc = [in(1:3, :); ones(1, size(in, 2))];
		in_ori = in(4:6, :);
		oriFlag = true;
	else
		error('We need two inputs: Transformation Matrix(4x4) and input data(3xN or 6xN)');
	end
		
	%% 2. perform computation
	out_loc = T * in_loc;
	out = out_loc(1:3, :);
	if oriFlag
		out(4:6, :) = T(1:3, 1:3) * in_ori;
	end
		
end %end of function