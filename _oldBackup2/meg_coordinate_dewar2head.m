% % % This function would compute the transformation matrix based on the fiducials
% % % and (optionally) the new coordinates and orientation of sensors
% % % Input:
% % %		fiducials (3 x 3 array),		3 columns are nasion, lpa and rpa respectively
% % %												3 rows are x, y, z values of their locations
% % %		dewar_pos (3/6 x N array),		N columns are N points in dewar system, 
% % %												3 rows are x, y, z values of their locations
% % %												6 rows are with u, v, w values of their orientations as row 4~6 
% % %		dewar_ori (3 x N array),		N columns are N points in dewar system, 
% % %												3 rows are u, v, w values of their orientations
% % % Output:
% % %		T (4 x 4 matrix),					Transformation matrix
% % %		head_fiducials (3 x 3 array),	3 columns are nasion, lpa and rpa respectively
% % %												3 rows are x, y, z values of their locations in head system
% % %		head_all (3/6 x N array),		N columns are N points in head system, 
% % %												3 rows are x, y, z values of their locations
% % %												6 rows are with u, v, w values of their orientations as row 4~6 

function [T, head_fiducials, head_all] = meg_coordinate_dewar2head (fiducials, dewar_pos, dewar_ori)
% % % written by wp in 07/11/2023


	%% 1. judge the inputs and outputs
	
	% % % 1.1 return with too few or many inputs	
	if nargin < 1 || size(fiducials, 1) ~= 3 || size(fiducials, 2) ~= 3 
		fprintf('We need inputs of fiducial locations in dewar system, the columns shall be: \n');
		fprintf('The nasion, and the left and right preauricular points (LPA and RPA). \n');
		fprintf('The rows shall by x, y, z coordinates. \n');
		return;
	elseif nargin > 3
		fprintf('No more than 3 inputs are accepted!\n');
		return;
	elseif nargout > 3
		fprintf('No more than 3 outputs are accepted!\n');
		return;
	elseif nargout < 1
		fprintf('At least one output!\n');
		return;
	end
	
	% % % 1.2 initiate flags for additional computation
	posFlag = false;
	oriFlag = false;
	fidFlag = false;

	
	% % % 1.3 whether dewar orientations are provided
	if nargin == 3
		if size(dewar_ori, 1) ~= 3  
			fprintf('Orientation shall be 3xN, skipped due to violation!\n');
		elseif size(dewar_pos, 1) ~= 3
			fprintf('Position shall be 3xN with Orientation separately provided, skipped due to violation!\n');
		else
			posFlag = true;
			oriFlag = true;
		end
	end
	
	% % % 1.4 whether dewar positions are provided
	if nargin == 2
		if size(dewar_pos, 1) == 3
			posFlag = true;
		elseif size(dewar_pos, 1) == 6
			oriFlag = true;
			dewar_ori = dewar_pos(4:6, :);
			posFlag = true;
			dewar_pos = dewar_pos(1:3, :);
		else
			fprintf('Position shall be 3xN or 6xN (with orientation), skipped due to violation!\n');
		end
	end
	
	% % % 1.5 output check
	if nargout == 3
		fidFlag = true;
	elseif nargout == 2
		fidFlag = true;
		posFlag = false;
		oriFlag = false;
	end	

	%% 2. computation of transformation matrix T

	% % % Define the coordinates of your fiducials
	nasion = fiducials(:, 1);
	lpa = fiducials(:, 2);
	rpa = fiducials(:, 3);

	% % % Step 1: Find the origin
	origin = (lpa + rpa) / 2;

	% % % Step 2: Create the X-axis
	x_axis = nasion - origin;
	x_axis = x_axis / norm(x_axis);  % Normalize

	% % % Step 3: Start defining the Y-axis
	y_axis_temp = lpa - origin;
	y_axis_temp = y_axis_temp / norm(y_axis_temp);  % Normalize
	
	% % % Step 4: Define the Z-axis
	z_axis = cross(x_axis, y_axis_temp);
	z_axis = z_axis / norm(z_axis);  % Normalize

	% % % Step 5: Correct the Y-axis to ensure orthogonality
	y_axis = cross(z_axis, x_axis);
	y_axis = y_axis / norm(y_axis);  % Normalize

	% % % Step 6: The rotation matrix is formed by the X, Y, and Z axes
	R = [x_axis y_axis z_axis]';

	% % % Step 7: The full transformation matrix including rotation and translation
	T = eye(4);  % Initialize as identity matrix
	T(1:3, 1:3) = R;  % Insert rotation
	T(1:3, 4) = -R * origin;  % Insert translation

	%% 3. Use this transformation matrix T to convert 
	
	% % % 3.1 fiducials
	if fidFlag
		fiducials_homogeneous = [fiducials; ones(1, 3)];
		transformed_fiducials = T * fiducials_homogeneous;
		head_fiducials = transformed_fiducials(1:3, :);
	end
	
	% % % 3.2 head_positions
	if posFlag
		sensor_locations_homogeneous = [dewar_pos; ones(1, size(dewar_pos, 2))];
		transformed_sensor_locations = T * sensor_locations_homogeneous;
		head_all = transformed_sensor_locations(1:3, :);
	end
	
		% % % 3.3 head_orientations
	if oriFlag
		head_all(4:6, :) = T(1:3, 1:3) * dewar_ori;
	end

end