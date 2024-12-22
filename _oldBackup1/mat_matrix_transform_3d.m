function out = mat_matrix_transform_3d(M, in)
% % % written on 11/08/2017 by wp, compute new xyz with transform matrix

	%% initialize
	if nargin < 2 || nargin > 2 || any(size(M) - [4 4]) || numel(size(in)) ~= 2
		error('We need inputs: M(4x4), in(nx3)');
	end
	
	if size(in, 2) == 3
		out = [in ones(size(in, 1), 1)] * M(1:3, :)';
	else
		fprintf('The data size is not Nx3...');
		if size(in, 1) == 3
			in = in';
			out = [in ones(size(in, 1), 1)] * M(1:3, :)';
			fprintf('   I assume you put it as 3xN, now converted!!!\n');
		elseif size(in, 2) == 6
			out = [in(:, 1:3) ones(size(in, 1), 1)] * M(1:3, :)';
			out(:, 4:6) = in(:, 4:6) * M(1:3, 1:3)';
			fprintf('   I assume the orientations were involved(Nx6), try to work on this!\n');
		elseif size(in, 1) == 6
			in = in';
			out = [in(:, 1:3) ones(size(in, 1), 1)] * M(1:3, :)';
			out(:, 4:6) = in(:, 4:6) * M(1:3, 1:3)';
			fprintf('   I assume you put it as 6xN, try to work on this!!!\n');
		else
			error('We need inputs: M(4x4), in(nx3)');
		end
	end
end %end of function