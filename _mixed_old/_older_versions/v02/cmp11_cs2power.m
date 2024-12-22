function [p, A1] = cmp11_cs2power(C, A3)
% % % 15/07/14	written by wp: 
% % % 	C : nChs x nChs;
% % % 	A : nChs x nVxs x 3

	%% 1. check inputs
	% % % spf	
	[nchs, nvxs, ndms] = size(A3); %get size
	if ndms ~= 3
		error('spatial filter shall be nChs x nVXs x 3 in size!');
	else
		A3 = permute(A3, [1, 3, 2]);	%Nx3xM
	end
	% % % C	
	C = real(C);

	%% 2. do computation
	p = zeros(nvxs, 1);
	parfor k = 1 : nvxs;
		c1 = A3(:, :, k)' * C * A3(:, :, k);
		[u, s, v] = svd(c1);
		A1(k, :) = A3(:, :, k) * u(:, 1);
		p(k)=s(1, 1);
	end
		
end %end of function

