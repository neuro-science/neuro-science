% % % This function would get 2d coordinates for topo plots from 3d sensor locations
% % % 
% % % Input:
% % %		loc3d (n x 3 array),		3 columns are x, y, z values of their locations
% % % Output:
% % %		loc2d (4 x 4 matrix),					2 columns are x, y values of their locations
% % % 
% % % closed by wp @11/12/2023 - it does not run well, fall back to guido's function
% % % written by wp @30/11/2023

% 
% function loc2d = meg_coordinate_3d_to_2d (loc3d)
% 	
% 	% Convert Cartesian (x, y, z) to Spherical (r, theta, phi)
% 	[azimuth, elevation, ~] = cart2sph(loc3d(:, 1), loc3d(:, 2), loc3d(:, 3));
% 
% 	% Convert spherical to polar coordinates for 2D plot
% 	theta = pi/2 - elevation; % Theta is the angle from the positive z-axis
% 	phi = azimuth;            % Phi is the azimuth angle from the positive x-axis
% 
% 	% Convert to 2D Cartesian coordinates
% 	[loc2d(:, 1), loc2d(:, 2)] = pol2cart(phi, theta);
% 	
% end
% 
% function loc2d = meg_coordinate_3d_to_2d(loc3d)
%     % Normalize the 3D coordinates (project onto unit sphere)
%     norms = sqrt(sum(loc3d.^2, 2));
%     loc3d_normalized = bsxfun(@rdivide, loc3d, norms);
% 
%     % Convert to spherical coordinates
%     [azimuth, elevation, ~] = cart2sph(loc3d_normalized(:, 1), ...
% 		 loc3d_normalized(:, 2), loc3d_normalized(:, 3));
% 
%     % Azimuthal equidistant projection
%     loc2d(:, 1) = azimuth .* cos(elevation);
%     loc2d(:, 2) = azimuth .* sin(elevation);
% end
% 
% function loc2d = meg_coordinate_3d_to_2d(loc3d)
%     % Project 3D coordinates onto the x-y plane using azimuthal projection
% 
%     % Normalize the 3D coordinates (project onto unit sphere)
%     norms = sqrt(sum(loc3d.^2, 2));
%     loc3d_normalized = bsxfun(@rdivide, loc3d, norms);
% 
%     % Calculate angles
%     azimuth = atan2(loc3d_normalized(:,2), loc3d_normalized(:,1));
%     elevation = asin(loc3d_normalized(:,3));
% 
%     % Azimuthal projection
%     loc2d(:, 1) = cos(elevation) .* cos(azimuth);
%     loc2d(:, 2) = cos(elevation) .* sin(azimuth);
% end
