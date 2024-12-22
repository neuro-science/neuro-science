function A = cmp31_spf_beamforming(C,  L, sigma)

% % % Input shall be:
%		Cross Spectrum Density: C(nChans, nChans);
%		Leadfield: L(nChans, nGrds, 3), ;
%		A constant for regularization: sigma

	%% headers
	% % % sigma for default
	if nargin < 3
		sigma = 0.05;
	end

	% % % sizes check
	[nChans, tmp2] = size(C);
	if tmp2 ~= nChans
		error('C is not a square matrix!' );
	end
	[tmp1, nGrids, tmp3] = size(L);
	if tmp1 ~=	nChans
		error('incosistent channel number in C and L!');
	elseif tmp3 ~= 3
		error('The leadfield is not 3-D!');
	end
	clear tmp1 tmp2 tmp3;

	% % % initialize filter
	A = zeros([nChans, nGrids, 3]);
% 	A1 = zeros([nChans, nGrids]);
% 	rv = zeros(nGrids, 1);
	
	%% do calculation
	% % % prepare data
	C = real(C);
	Ct = C + sigma * trace(C) / nChans * eye(nChans);

	% % % compute spf in loop
	parfor iG = 1 : nGrids
		L1 = squeeze(L(:, iG, :));
		tmp = (L1' / Ct);
		A0 = (tmp * L1) \ tmp;
		A(:, iG, :) = permute(A0, [2 3 1]);

% 		[u, s, v] = svd(A0 * C * A0');
% 		tmp = A0' * u(:, 1);
% 		rv(iG)=s(1,1);
% 		
% 		if isreal(tmp)
% 			A01 = tmp;
% 		elseif abs(imag(tmp) ./ real(imag(tmp))) < 0.01
% 			A01 = real(tmp);
% 		else
% 			fprintf('======\ncomplex results!\n======\n')
% 			A01 = nan;
% 		end
% 		A1(:, iG) = A01;
	end
	
% 	A1 = A1';
end %end of function

