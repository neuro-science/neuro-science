% Turn those inward orientation to ourtwards, for meg sensors/references, usefull for lead field
% 'pos' is a Nx3 matrix where each row represents x, y, z coordinates of the points
% 'in' is a Nx3 matrix where each row represents the orientation vector of the points
function out = meg_coordinate_ori_outwards(pos, in)
% % % updated on 22/11/2023 by wp, tolerance for other shape
% % % written on 08/11/2023 by wp


	% % % check inputs
	if numel(size(pos)) ~= 2 || nargin > 2 || nargin < 1
		fprinf('We need the iput size to be Nx3, for both position and orientation!\n');
		out = [];
		return;
	elseif nargin == 1
		if size(pos, 2) ~= 6 && size(pos, 1) == 6
			in = pos(4:6, :)';
			pos = pos(1:3, :)';
			fprintf('We accept Nx6 input, assuming last 3 columns indicate orientation and you offer 6xN, inverted!\n');
		elseif size(pos, 2) == 6
			in = pos(:, 4:6);
			pos = pos(:, 1:3);
		else
			out = [];
			fprinf('We need the iput size to be Nx3, for both position and orientation!\n');
			return;
		end		
	else
		if size(pos, 2) ~= 3 && size(pos, 1) == 3
			pos = pos';
			fprintf('We need Nx3 input and you offer 3xN pos, inverted!\n');
		end
		if size(in, 2) ~= 3 && size(in, 1) == 3
			in = in';
			fprintf('We need Nx3 input and you offer 3xN ori, inverted!\n');
		end
	end
		
	% Calculate the center of the sphere (if it's not already known)
	center = mean(pos, 1);

	% Allocate space for the new orientations
	out = in;

	% Loop through each point
	for i = 1 : size(pos, 1)
		 % Vector from center to point
		 center_to_point = pos(i, :) - center;

		 % Normalize the center_to_point vector
		 center_to_point = center_to_point / norm(center_to_point);

		 % Calculate the dot product
		 dot_product = dot(center_to_point, in(i, :));

		 % If dot_product is negative, the orientation is inward. Flip it by multiplying by -1
		 if dot_product < 0
			  out(i, :) = -in(i, :);
		 end
	end

end