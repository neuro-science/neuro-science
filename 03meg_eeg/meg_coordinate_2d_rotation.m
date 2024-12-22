function loc2d_rotated = meg_coordinate_2d_rotation(loc2d_origin, theta)
	% % % written by wp @11/12/2023
	% % % note that input and output are both Nx2 array (x,y of N points)
	% % % theta pi/4 means 45 degree counter-clockwise	

    % Define the rotation matrix for theta degrees counterclockwise
    R = [cos(theta) -sin(theta); sin(theta) cos(theta)];

    % Apply the rotation to each point
    loc2d_rotated = (R * loc2d_origin')';

    % s2d_rotated now contains the rotated coordinates
end
