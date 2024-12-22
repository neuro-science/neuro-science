% % % The code was adopted from MVGC Toolbox 1.0 for my meg data

function [F, G, cfg] = meg_grangerCausality (d, numLag, cfg)
	%% 1. prepare the parameters and data
	if nargin < 2 || isempty(numLag)
		numLag = 10;	% as 400Hz is the default sampling frequency, this is gamma
	end
	
	if nargin < 3
		cfg = [];
	end
	
	if ~isfield(cfg, 'cfg.acMaxLag') || isempty(cfg.acMaxLag)
		cfg.acMaxLag = 0;	% % % default 0 to take minLag
	end
	
	if ~isfield(cfg, 'acDecayTolerance') || isempty(cfg.acDecayTolerance)
		cfg.acDecayTolerance = 1e-8;	% 
	end	
	
	if ~isfield(cfg, 'acMaxRelativeError') || isempty(cfg.acMaxRelativeError)
		cfg.acMaxRelativeError = 1e-8;	% 
	end	
	
	% initialise info struct
	cfg.errID = 0;
	cfg.errMessage = '';
	cfg.warnID = 0;
	cfg.warnMessage = cell(0,1);
	
	[nChans, nPoints, nTrials]   = size(d);     % data size check
	if nPoints < numLag
		fprintf('The designed lag is larger than data size, changed from %d to %d!\n', numLag, nPoints);
		numLag = nPoints;
	end
	d = bsxfun(@minus, d, mean(d(:, :), 2));		%remove mean of the data
	
	%% 2. computation for vector autoregression, LWR (Morf) method
	I = eye(nChans);
	q1n = nChans * (numLag + 1);
	% % % store lags
	dd = zeros(nChans, numLag + 1, nPoints + numLag, nTrials);
	for k = 0 : numLag
		dd(:, k + 1, k + 1 : k + nPoints, :) = d; % k-lagged observations
	end
	
	% % % initialise recursion
	AF = zeros(nChans, q1n); % forward  AR coefficients
	AB = zeros(nChans, q1n); % backward AR coefficients (reversed compared with Morf's treatment)
	
	k  = 1;            % model order is k-1
	kn = k * nChans;
	M  = nTrials * (nPoints - k);
	kf = 1 : kn;         % forward  indices
	kb = q1n - kn + 1 : q1n; % backward indices

	XF = reshape(dd(:, 1:k, k+1:nPoints, :), kn, M);
	XB = reshape(dd(:, 1:k, k:nPoints-1, :), kn, M);

	[CXF, cholp] = chol(XF * XF');
	if cholp; fprintf('CXF\n'); return; end % show-stopper!

	[CXB, cholp] = chol(XB * XB');
	if cholp; fprintf('CXB\n'); return; end % show-stopper! 

	AF(:, kf) = CXF'\I;
	AB(:, kb) = CXB'\I;
	
	% % % loop with recursion	
	while k <= numLag
		EF = AF(:, kf) * reshape(dd(:, 1:k, k+1:nPoints,:), kn, M); % forward  prediction errors
		EB = AB(:, kb) * reshape(dd(:, 1:k, k:nPoints-1,:), kn, M); % backward prediction errors

		[CEF, cholp] = chol(EF * EF');
		if cholp; fprintf('CEF\n'); return; end % show-stopper! 

		[CEB, cholp] = chol(EB * EB');
		if cholp; fprintf('CEB\n'); return; end % show-stopper! 

		R = CEF'\(EF * EB')/CEB; % normalised reflection coefficients

		[RF, cholp] = chol(I - R*R');
		if cholp; fprintf('RF\n'); return; end % show-stopper! 

		[RB, cholp] = chol(I - R'*R);
		if cholp; fprintf('RB\n'); return; end % show-stopper! 

		k  = k + 1;
		kn = k * nChans;
		M  = nTrials * (nPoints - k);
		kf = 1 : kn;
		kb = q1n - kn + 1 : q1n; % backward indices


		AFPREV = AF(:, kf);
		ABPREV = AB(:, kb);

		AF(:, kf) = RF'\(AFPREV - R * ABPREV);
		AB(:, kb) = RB'\(ABPREV - R' * AFPREV);
	end
	clear dd;

	A = reshape(-AF(:, 1 : nChans) \ AF(:,nChans + 1 : end), nChans, nChans, numLag); % so A(:,:,k) is the k-lag coefficients matrix
	E   = AFPREV(:, 1 : nChans) \ EF;   % residuals
	SIG = (E * E') / (M - 1);       % residuals covariance matrix
	[~, cholp] = chol(SIG);
	if cholp; fprintf('SIG\n'); 
		 cfg.errID = 2;
		 cfg.errMessage = 'residuals covariance matrix not positive-definite';
		return; 
	end % show-stopper! 
% 	E   = reshape(E, nChans, nPoints - numLag, nTrials); % put residuals back into per-trial form
	clear AF AB AFPREV ABPREV E EB EF XB XF;

	%% 3. computation for autocovariance sequence
	pn = numLag * nChans;
	pn1 = (numLag - 1) * nChans;
	A = reshape(A, nChans, pn);                                   % coefficients
		
	% construct VAR coefficients for 1-lag problem
	A1 = [A; eye(pn1) zeros(pn1, nChans)];

	% calculate spectral radius
	cfg.rho = max(abs(eig(A1)));
	
	if cfg.rho >= 1
		 cfg.errID = 1;
		 cfg.errMessage = 'unstable VAR (unit root)';
		 return
	end

	% construct residual covariances for 1-lag problem
	SIG1 = [SIG zeros(nChans, pn1); zeros(pn1, nChans) zeros(pn1)];

	% solve the Lyapunov equation for the 1-lag covariance matrix
	try
		G1 = dlyap(A1, SIG1);  % dlyap seems to work better here without balancing, which seems to break positive-definitiveness
		% % %         G1 = lyapslv('D',A1,[],-SIG1); % sometimes. However lyapslv is not an official interface, so this could conceivably break in future.
	catch except
		cfg.errID = 3;
		cfg.errMessage = ['Lyapunov equation solver failed: ' except.message];
		return
	end

	cfg.acRelativeError = norm(A1 * G1 * A1' - G1 + SIG1) / norm(SIG1); % this should be small (see below)

	if cfg.acRelativeError > cfg.acMaxRelativeError
		 cfg.warnID = cfg.warnID+1;
		 cfg.warnMessage{cfg.warnID} = sprintf('large relative error = %g (tolerance = %g)',cfg.acRelativeError, cfg.acMaxRelativeError);
	end

	% estimate number of autocov lags
	cfg.acMinLag = ceil(log(cfg.acDecayTolerance)/log(cfg.rho)); % minimum lags to achieve specified tolerance

	if cfg.acMaxLag < 0  % use exactly -cfg.acMaxLag lags (not encouraged, hence undocumented!)
		cfg.acLag = -cfg.acMaxLag;
	elseif cfg.acMaxLag > 0  % use at most cfg.acMaxLag lags
		cfg.acLag = min(cfg.acMinLag, cfg.acMaxLag);
	else                  % cfg.acMaxLag == 0 - use minimum acceptable lags (recommended)
		cfg.acLag = cfg.acMinLag;
	end

	if cfg.acLag < cfg.acMinLag
		 cfg.warnID = cfg.warnID+1;
		 cfg.warnMessage{cfg.warnID} = sprintf('too few autocovariance lags = %d (minimum = %d)', cfg.acLag, cfg.acMinLag);
	end

	% calculate recursively from 1-lag solution (which supplies up to p-1 lags), from p lags up to q
	assert(cfg.acLag >= numLag, 'number of lags is too small'); % lags must be at least number of VAR lags
	G = cat(3, reshape(G1(1 : nChans, :), nChans, nChans, numLag), zeros(nChans, nChans, cfg.acLag - numLag + 1));   % autocov forward  sequence
	B = [zeros((cfg.acLag - numLag + 1) * nChans, nChans); G1(:, end - nChans + 1 : end)];            % autocov backward sequence
	for k = numLag : cfg.acLag
		 r = cfg.acLag - k + 1;
		 G(:, :, k + 1) = A * B(r * nChans + 1 : r * nChans + pn, :);
		 B((r - 1) * nChans + 1 : r * nChans, :) = G(:, :, k+1);
	end
	clear G1 G0 GB GF B;
	
	%% 4. F statistics
	[~, SIG] = ACV2V(G);
	LSIG = log(diag(SIG));
	F = nan(nChans);
	for j = 1 : nChans
		 % reduced regression
		 jo = [1:j-1 j+1:nChans]; % omit j
		 [~, SIGj] = ACV2V(G(jo, jo, :));
		 LSIGj = log(diag(SIGj));
		 for ii = 1 : nChans - 1
			  i = jo(ii);
			  F(i, j) = LSIGj(ii) - LSIG(i);
		 end
	end	
end


function [AF, SIG] = ACV2V(G)
	% full regression
	[nChans, ~,	q1] = size(G);
	q = q1 - 1;
	qn = q * nChans;

	G0 = G(:,:,1);                                               % covariance
	GF = reshape(G(:, :, 2:end), nChans, qn)';                            % forward  autocov sequence
	GB = reshape(permute(flip(G(:, :, 2:end), 3),[1 3 2]), qn, nChans); % backward autocov sequence

	AF = zeros(nChans, qn); % forward  coefficients
	AB = zeros(nChans, qn); % backward coefficients (reversed compared with Whittle's treatment)

	% initialise recursion

	k = 1;            % model order

	r = q - k;
	kf = 1 : k * nChans;       % forward  indices
	kb = r * nChans + 1 : qn;    % backward indices

	AF(:, kf) = GB(kb, :) / G0;
	AB(:, kb) = GF(kf, :) / G0;

	% and loop

	for k = 2 : q

		 AAF = (GB((r - 1) * nChans + 1 : r * nChans, :) - AF(:, kf) * GB(kb, :)) / (G0 - AB(:, kb) * GB(kb, :)); % DF/VB
		 AAB = (GF((k - 1) * nChans + 1 : k * nChans, :) - AB(:, kb) *GF (kf, :)) / (G0 - AF(:, kf) * GF(kf, :)); % DB/VF

		 AFPREV = AF(:, kf);
		 ABPREV = AB(:, kb);

		 r = q - k;
		 kf = 1 : k * nChans;
		 kb = r * nChans + 1 : qn;

		 AF(:, kf) = [AFPREV - AAF * ABPREV AAF];
		 AB(:, kb) = [AAB ABPREV - AAB * AFPREV];
	end
	
	SIG = G0-AF*GF;
	AF = reshape(AF, nChans, nChans, q);
end


function [f, fres] = ACV2gc(G, fres, useFFT)

	if nargin < 3, useFFT = []; end % force autocov_xform default

	[n, ~, q1] = size(G);
	if nargin < 2 || isempty(fres)
		 fres = q1;
	end

	h = fres + 1;
	f = nan(n, n, h);

	for j = 1:n
		 jo  = [1:j-1 j+1:n]; % omit j
		 joj = [jo j];        % rearrange with omitted j indices first

		 [Aj, SIGj] = ACV2V(G(jo, jo, :));  % reduced regression
		 Gj = autocov_xform(G(joj, joj, :), Aj, SIGj, useFFT); % transform autocov

		 [Ajj, SIGjj] = ACV2V(Gj);        % transformed full regression
		 Hjj = var2trfun(Ajj,fres);               % transfer function

		 for ii=1:n-1;
			  i  = jo(ii);           % i index in omitted j indices
			  io = [1:ii-1 ii+1:n];  % omit i

			  SIGji = SIGjj(io,io)-(SIGjj(io,ii)*SIGjj(ii,io))/SIGjj(ii,ii); % partial covariance
			  Hji = Hjj(ii,io,:);                                            % transfer function
			  Sji = SIGj(ii,ii);                                             % i part of spectrum is flat!

			  LSji = log(Sji);
			  for k = 1:h
					f(i,j,k) = LSji - log(real(Sji-Hji(:,:,k)*SIGji*Hji(:,:,k)'));
			  end
		 end
	end
end

%% below are functions copied from MVGC toolbox (by Lionel Barnett and Anil K. Seth)    

function X = dlyap(A, B, C, E)
	%DLYAP  Solve discrete Lyapunov equations.
	%
	%   X = DLYAP(A,Q) solves the discrete Lyapunov matrix equation:
	%
	%       A*X*A' - X + Q = 0
	%
	%   X = DLYAP(A,B,C) solves the Sylvester equation:
	%
	%       A*X*B - X + C = 0
	%
	%   X = DLYAP(A,Q,[],E) solves the generalized discrete Lyapunov equation:
	%
	%       A*X*A' - E*X*E' + Q = 0
	%
	%   See also DLYAPCHOL, LYAP.

	%	J.N. Little 2-1-86, AFP 7-28-94
	%  Copyright 1986-2018 The MathWorks, Inc.

	narginchk(2,4)
	ni = nargin;
	if ni<3
		C = [];
	end
	if ni<4
		E = [];
	end

	% Validate data 
	try
		[A,B,C,E] = lyapcheckin('dlyap',ni,A,B,C,E);
	catch err
		throw(err)
	end

	% Balance to minimize spectrum distorsions in Hess/Schur/QZ factorizations
	if ni==3
		% Sylvester
		[sA,pA,A] = mscale(A,'fullbal');
		[sB,pB,B] = mscale(B,'fullbal');
		C(pA,pB) = lrscale(C,1./sA,sB);  % TA\C*TB
	else
		[A,~,~,E,s,p] = aebalance(A,[],[],E,'fullbal');
		B(p,p) = lrscale(B,1./s,1./s);  % T\B/T'
	end

	% Solve equation
	if ni==3
		 % Sylvester equation A*X*B-X+C=0
		 X = dsylvester(A, B, -C);
		 X = lrscale(X(pA,pB),sA,1./sB);  % TA*X/TB
	else
		 % Lyapunov equation A*X*A'-X+B=0 or A*X*A'-E*X*E'+B=0
		 if isequal(E, [])
			  X = dlyapunov(A, -B);
		 else
			  X = dgenlyapunov(A, -B, E);
		 end
		 X = lrscale(X(p,p),s,s);  % T*X*T' using T(:,p)=diag(s)
	end

	% Error if X contains Inf or NaN.
	if ~all(isfinite(X(:)))
		 try
			  if ni==3
					ctrlMsgUtils.error('Control:foundation:SingularSylv')
			  else
					ctrlMsgUtils.error('Control:foundation:SingularLyap')
			  end
		 catch err
			  throw(err);
		 end
	end
end


function X = dlyapunov(A, C)
	% Solve simple discrete Lyapunov equation A*X*A' - X = C

	if ishermitian(A)
		 [QA, dA] = eig(A, 'vector');

		 CC = QA'*C*QA;
		 X = CC ./ (dA.*dA' - 1);
		 X = QA*X*QA';

	else
		 % Reduce equation to triangular form
		 flag = 'real';
		 if ~isreal(A) || ~isreal(C)
			  flag = 'complex'; % Need complex Schur form
		 end

		 CC = -C;
		 schurA = matlab.internal.math.isQuasiTriangular(A,flag);
		 if schurA
			  TA = A;
		 else
			  [QA, TA] = schur(A, flag);
			  CC = QA'*CC*QA;
		 end

		 % Solve Lyapunov Equation -TA*X*TA' + X = -QA'*C*QA.
		 X = matlab.internal.math.sylvester_tri(TA, 'I', CC, 'I', -TA, 'transp');

		 % Recover X
		 if ~schurA
			  X = QA*X*QA';
		 end
	end

	if ishermitian(C)
		 X = (X + X') / 2;
	end
end

function X = dsylvester(A, B, C)
	% Solve Sylvester Equation A*X*B - I = C.

	if ishermitian(A) && ishermitian(B)
		 [QA, dA] = eig(A, 'vector');
		 [QB, dB] = eig(B, 'vector');

		 CC = QA'*C*QB;
		 X = CC ./ (dA.*dB' - 1);
		 X = QA*X*QB';

	else
		 % Reduce equation to triangular form
		 flag = 'real';
		 if ~isreal(A) || ~isreal(B) || ~isreal(C)
			  flag = 'complex'; % Need complex Schur form
		 end

		 CC = -C;
		 schurA = matlab.internal.math.isQuasiTriangular(A,flag);
		 if schurA
			  TA = A;
		 else
			  [QA, TA] = schur(A, flag);
			  CC = QA'*CC;
		 end
		 schurB = matlab.internal.math.isQuasiTriangular(B,flag);
		 if schurB
			  TB = B;
		 else
			  [QB, TB] = schur(B, flag);
			  CC = CC*QB;
		 end

		 % Solve Sylvester Equation -TA*X*TB' + X = -QA'*C*QB.
		 X = matlab.internal.math.sylvester_tri(TA, 'I', CC, 'I', -TB, 'notransp');

		 % Recover X
		 if ~schurA
			  X = QA*X;
		 end
		 if ~schurB
			  X = X*QB';
		 end
	end
end

function X = dgenlyapunov(A, C, E)
	% Solve generalized discrete Lyapunov equation A*X*A' - E*X*E' = C

	% Reduce equation to triangular form
	flag = 'real';
	if ~isreal(A) || ~isreal(C) || ~isreal(E)
		 flag = 'complex'; % Need complex Schur form
	end

	CC = C;
	qzAE = matlab.internal.math.isQuasiTriangular(A, E, flag);
	if qzAE
		 TA = A;
		 TE = E;
	else
		 [TA, TE, Q, Z] = qz(A, E, flag);
		 CC = Q*CC*Q';
	end

	% Solve Lyapunov Equation TA*X*TA' - TE*X*TE' = Q*C*Q'.
	X = matlab.internal.math.sylvester_tri(TA, TE, CC, -TE, TA, 'transp');

	% Recover X
	if ~qzAE
		 X = Z*X*Z';
	end

	if ishermitian(C)
		 X = (X + X') / 2;
	end
end

function [A,B,C,E] = lyapcheckin(Caller,ni,A,B,C,E)
	% Validates input arguments to LYAP and DLYAP.

	%   Copyright 1986-2007 The MathWorks, Inc.
	if ~(isnumeric(A) && isnumeric(B) && isnumeric(C) && isnumeric(E))
		ctrlMsgUtils.error('Control:foundation:Lyapunov1',Caller)
	end

	if ni==3
		% Sylvester
		A = double(full(A));  B = double(full(B));  C = double(full(C));
		szA = size(A);  szB = size(B);  szC = size(C);
		if length(szA)>2 || length(szB)>2 || szA(1)~=szA(2) || szB(1)~=szB(2)
			ctrlMsgUtils.error('Control:foundation:Sylvester1',Caller)
		elseif length(szC)>2 || szC(1)~=szA(1) || szC(2)~=szB(1)
			ctrlMsgUtils.error('Control:foundation:Sylvester2',Caller)
		elseif hasInfNaN(A) || hasInfNaN(B) || hasInfNaN(C)
			ctrlMsgUtils.error('Control:foundation:Lyapunov6',Caller)
		end
	else
		% Lyapunov
		A = double(full(A));  B = double(full(B));  E = double(full(E));
		szA = size(A);  szB = size(B);  szE = size(E);
		if length(szA)>2 || length(szB)>2 || szA(1)~=szA(2) || szB(1)~=szB(2) || szA(1)~=szB(1)
			ctrlMsgUtils.error('Control:foundation:Lyapunov4',Caller)
		elseif any(szE) && any(szE~=szA(1))
			ctrlMsgUtils.error('Control:foundation:Lyapunov5',Caller)
		elseif hasInfNaN(A) || hasInfNaN(B) || hasInfNaN(E)
			ctrlMsgUtils.error('Control:foundation:Lyapunov6',Caller)
		end
	end
end


function P = mvgc_cdf(x,X,p,m,N,nx,ny,nz,tstat)
	assert(isvector(x),'evaluation values must be a vector');
	n = length(x);

	assert(isvector(X),'MVGC values must be a vector');
	assert(all(X >= 0),'MVGC values must be non-negative');
	if isscalar(X)
		 X = X*ones(n,1);
	else
		 assert(length(X) == n,'MVGC values must match evaluation values');
	end

	if nargin < 8 || isempty(nz), nz = 0; end % unconditional

	if nargin < 9 || isempty(tstat);
		 ftest = nx == 1; % default: use F-distribution for univariate predictee, chi2 for multivariate
	else
		 switch lower(tstat)
			  case 'f'     % Granger F-test form
					assert(nx == 1,'F-distribution is not appropriate for multivariate predictee');
					ftest = true;
			  case 'chi2'  % Geweke chi2 test form
					ftest = false;
			  otherwise
					error('unknown distribution (must be ''chi2'' or ''F'')');
		 end
	end

	P = zeros(n,1);
	m = N*(m-p);                  % effective number of observations (p-lag autoregression loses p observations per trial)
	if ftest
		 if any(X > 0), fprintf(2,'WARNING (mvgc_cdf): non-central F-distribution is experimental\n'); end
		 d1 = p*ny;                % #{full model parameters} - #{reduced model parameters}
		 d2 = m-p*(1+ny+nz);       % #{observations} - #{full model parameters}
		 mm = d2/d1;
		 for i = 1:n
			  xx = exp(x(i))-1;     % Granger form: (RSS_reduced - RSS_full) / RSS_full
			  if X(i) > 0           % non-central
					XX = exp(X(i))-1; % Granger form
					P(i) = ncfcdf(mm*xx,d1,d2,m*XX); % NOTE: non-centrality parameter factor might reasonably be m, d2 or d2-2
			  else
					P(i) = fcdf(mm*xx,d1,d2);
			  end
		 end
	else
		 d = p*nx*ny;              % note that d does not depend on the number of conditioning variables
		 for i = 1:n
			  if X(i) > 0           % non-central
					P(i) = ncx2cdf(m*x(i),d,m*X(i));
			  else
					P(i) = chi2cdf(m*x(i),d);
			  end
		 end
	end
end
