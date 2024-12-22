function [src_head, M] = fcv09_src_temp2head_nomri(src_temp, fd, zFlag)
% % % 24/11/15 by wp 
% % % src_temp : source locations in template, Nx3
% % % fd			: [center2nasion, lpa-rpa, z];[nas; lpa;rpa] in head;
% % % zFlag		: whether we guess for the z-dimension scaler;
% % % src_head : source locations in head system
% % % M			: transformation matrix, Nx3
	%% 1. verify inputs
	if nargin < 3 || isempty(zFlag)
		zFlag = 0;
	end
	if nargin < 2
		error('Both source in templates and fiducials are needed!');
	end
	if size(src_temp, 2) ~= 3 || length(size(src_temp)) ~= 2
		error('The source locations shall be Nx3 matrix');
	end
	eFlag = all(abs([fd(1, 2); fd(:, 3);fd(2,1) + fd(3,1);fd(2,2) + fd(3,2)]) < 1e-8) & (fd(2, 2) > 0);
	if ~eFlag
		error('The fiducials are not in correct format, be sure it is in head system!');
	end
	z1 = fd(1, 1);z2 = fd(2, 2);z3 = fd(2, 1);
	%% 2. constant definition
	c1 = 11.03676; c2 = 8.8;
	M0 = [0 0.09967 0.0081546 2.9565; -0.1 0 0 0;0 -0.0081546 0.09967 5.47715; 0 0 0 1];
	
	%% 3. computation
	% % % transform in 2-D
	M1 = [z1/c1 z3/c2 0 0; 0 z2/c2 0 0; 0 0 1 0; 0 0 0 1];
	% % % compensate for z dimension?
	if zFlag > 0	
		M1(3, 3) = sqrt((z1 * z2) / (c1 * c2));
	elseif zFlag < 0
		M1(3, 3) = 1;
	elseif ((z1/c1 > 1) && (z2/c2 > 1)) || ((z1/c1 < 1) && (z2/c2 < 1))
		M1(3, 3) = sqrt((z1 * z2) / (c1 * c2));
	else
		M1(3, 3) = 1;
	end
	% % % combination of transforms
	M = M1 * M0;
	% % % transform
	in = src_temp';
	in(4, :) = 1;
	out = M * in;
	out(4, :) = [];
	src_head = out';
	clear in out M1 M0 c1 c2 z1 z2 z3 eFlag;
	
end % end of function