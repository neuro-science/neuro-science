function M = mat_compute_transform_matrix(Xin, Xout, rigidFlag)
% % % updated on 11/08/2017 by wp, consider head by fiducials
% % % updated on 10/08/2017 by wp, consider rigid transformation
% % % written on 15/05/2017 by wp, compute transform matrix based on two sets of points

	%% initialize
	fidFlag = false;
	if nargin < 1 || nargin > 3
		error('We need input: fiducials(nx3) or Xin(nx3), Xout(nx3) [, rigidfFlag]!');
	elseif nargin == 1
		fidFlag = true;
	elseif nargin == 2
		rigidFlag = true;
	end
	
	if fidFlag
		% % % get the landmarks
		nas = Xin(1, :); lpa = Xin(2, :); rpa = Xin(3, :);
		origin = (lpa + rpa)/2;
		% % % get the directions
		dirx = nas - origin;
		dirz = cross(dirx, lpa - rpa);
		diry = cross(dirz, dirx);
		dirx = dirx / norm(dirx);
		diry = diry / norm(diry);
		dirz = dirz / norm(dirz);	
		% % % rotation
		R = eye(4);
		R(1:3, 1:3) = inv(eye(3) / [dirx; diry; dirz]);
		% % % translation
		T = eye(4);
		T(1:3, 4) = -origin(:);
		% homogeneous transformation matrix
		M = R * T;

	else
		M = zeros(4, 4);
		%% translate	
		% % % 
		m1 = mean(Xin);
		m2 = mean(Xout);
		X1 = bsxfun(@minus, Xin, m1);
		X2 = bsxfun(@minus, Xout, m2);

		%% rotation
		if rigidFlag
			[u, s, v] = svd(X1' * X2);
			m = eye(3);
			m(3, 3) = det(v * u');
			R = v * m * u';
			M(1:3, 1:3) = R;
			M(1:3, 4) = m2' - R * m1';  
		else
			X1(X1 == 0) = 1e-5;
			X2(X2 == 0) = 1e-5;
			tmp = pinv(X1) * X2;
			M(1:3, 1:3) = tmp';
			M(1:3, 4) = m2 - m1 * tmp;
		end
		M(4, 4) = 1;
		clear m1 m2 X1 X2 tmp m R v s u;
	end
end %end of function