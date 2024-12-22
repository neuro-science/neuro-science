function [mX, x, y, p, r] = mat_myFitTemporalPrediction (x0, y0, nTrials, pltScale, ctTlrz, mValue, linkFun)
% % % updated 03/01/2021 by wp for temporal prediction
% % % try to be as consistant as Johnathan's method

	%% 1 make data column-wise
	if nargin < 2
		error('You need to give both X and Y for fit!');
	elseif abs(numel(x0) - numel(y0)) > 0.1
		error('X and Y should be of the same size!');
	elseif nargin < 3 || isempty(nTrials)
		nTrials = ones(size(x0));
	end
	
	if nargin < 4 || isempty(pltScale)
		pltScale = 10;
	end
	
	if nargin < 5 || isempty(ctTlrz)
		ctTlrz = 0.001;
	end
	
	if nargin < 6 || isempty(mValue)
		mValue = 0.5;
	end
	
	if nargin < 7 || isempty(linkFun)
% 		linkFun = 'probit';
		linkFun = 'logit';
	end
	
	if size(x0, 1) < size(x0, 2)
		 x0 = x0';
	end
	
	if size(y0, 1) < size(y0, 2)
		 y0 = y0';
	end
	
	if size(nTrials, 1) < size(nTrials, 2)
		 nTrials = nTrials';
	end
	
	%% 2 do fit
	p = glmfit(x0, [y0, nTrials], 'binomial','link',linkFun);	
	
	%% 3 data for plot	
	x = linspace(min(x0), max(x0), numel(x0) * pltScale);
	y = glmval(p, x, linkFun);

	[d, I] = min(abs((y - mValue)));
	mX = x(I);
	while d > ctTlrz && pltScale < 1e4
		pltScale = pltScale * 2;
		fprintf('.');
		xk = linspace(min(x0), max(x0), numel(x0) * pltScale);
		yk = glmval(p, xk, linkFun);
		[d, I] = min(abs((yk - mValue)));
		mX = xk(I);
	end
	
	%% 4 evaluate fit	
	if nargout > 4
		y1 = glmval(p, x0, linkFun);
		m = nanmean(y0);
		d0 = nansum((y0 - m).^2 .* nTrials);
		d1 = nansum((y1 - y0).^2 .* nTrials);
		r = d1 ./ d0;
	end
end
