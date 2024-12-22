function loc3d_new = meg_coordinate_3d_rotation(loc3d_old, angles)
    % Rotate 3D coordinates around the x, y, and z axes.
    % angles: [x_rotation, y_rotation, z_rotation] in radians
    
    % Extract the individual rotation angles
    x_rotation = angles(1);
    y_rotation = angles(2);
    z_rotation = angles(3);

    % Rotation matrix around the x-axis
    Rx = [1 0 0; 
          0 cos(x_rotation) -sin(x_rotation); 
          0 sin(x_rotation) cos(x_rotation)];

    % Rotation matrix around the y-axis
    Ry = [cos(y_rotation) 0 sin(y_rotation); 
          0 1 0; 
          -sin(y_rotation) 0 cos(y_rotation)];

    % Rotation matrix around the z-axis
    Rz = [cos(z_rotation) -sin(z_rotation) 0; 
          sin(z_rotation) cos(z_rotation) 0; 
          0 0 1];

    % Combined rotation matrix (rotation order: z, y, x)
    R = Rx * Ry * Rz;

    % Apply the rotation to the coordinates
    loc3d_new = (R * loc3d_old')';
end
