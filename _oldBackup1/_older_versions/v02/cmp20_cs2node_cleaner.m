
function [imc_sum, imc_max] = cmp20_cs2node_cleaner(C, A1, A2)
% % % modified from guido's function <cs2node_mim>
% % % some switches and output text were removed to speed up
% % % output of imaginary coherence with optimization of two orientations
% % % input A(ch, vx, dim) C(ch, ch)

	%% 1. evaluate input & initiate
	% % % get sizes
	[nch, nvx1, ndim]=size(A1);
	[nch2, nvx2, ndim2]=size(A2);
	[nch3, nch4] = size(C);
	% % % check consistance
	if nch ~= nch2 || nch2 ~= nch3 || nch3 ~= nch4
		error('The number of channels do not match!');
	elseif ndim ~= ndim2
		error('Dimension of filters do not match!');
	end
	% % % initiate
	imc_sum = zeros(nvx1, nvx2);
	imc_max = zeros(nvx1, nvx2);
	% % % convert A
	A1 = permute(A1, [1 3 2]);
	A2 = permute(A2, [1 3 2]);
	% % % the regulation factor   
	regu = 0.000001;

	%% 2. prepare data
	% % % voxel level cs 1
	csvoxrealinv = zeros(ndim, ndim, nvx1);
	csvoxvec = zeros(ndim, nch, nvx1);
	parfor iv = 1 : nvx1
		csloc = A1(:, :, iv)' * real(C) * A1(:, :, iv);
		csvoxrealinv(:, :, iv) = inv(csloc + regu * eye(ndim) * mean(diag(csloc)));
		csvoxvec(:, :, iv) = A1(:, :, iv)' * C;
		csloc = [];
	end

	% % % voxel level cs 1
	csvoxrealinv2 = zeros(ndim, ndim, nvx2);
	parfor iv = 1 : nvx2
		csloc = A2(:, :, iv)' * real(C) * A2(:, :, iv);
		csvoxrealinv2(:, :, iv) = inv(csloc + regu * eye(ndim) * mean(diag(csloc)));
		csloc = [];
	end

	%% 3.do it in loop in loop
	for iv = 1 : nvx1
		csvoxvecloc = csvoxvec(:, :, iv);
		caainv = csvoxrealinv(:, :, iv);
		caainvsqrt = sqrtm(csvoxrealinv(:, :, iv));

		parfor ix = 1 : nvx2
			cab = imag(csvoxvecloc * A2(:, :, ix));
			cbbinv = csvoxrealinv2(:, :, ix);
			X = cab * cbbinv * cab';
			imc_sum(iv, ix) = trace(caainv * X);
			Y = caainvsqrt * X * caainvsqrt;
			[u, s, v] = svd(Y);
			imc_max(iv, ix) = sqrt(s(1, 1));
		end
	end


return;
	  

						  