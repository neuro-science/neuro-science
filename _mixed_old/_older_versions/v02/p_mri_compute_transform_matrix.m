function M = p_mri_compute_transform_matrix(Xin, Xout, rigidFlag)
% % % updated on 15/05/2017 by wp, consider rigid transformation
% % % written on 15/05/2017 by wp, compute transform matrix based on two sets of points

	%% initialize
	if nargin < 2 || nargin > 3
		error('We need input: Xin(nx3), Xout(nx3), rigidfFlag!');
	elseif nargin == 2
		rigidFlag = true;
	end
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
end %end of function