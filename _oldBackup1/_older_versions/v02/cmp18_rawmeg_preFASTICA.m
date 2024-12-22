function [D, E, infoRatio, nICs] = cmp18_rawmeg_preFASTICA(data, infoRatio, numLim4ICs, stepRatio, rsmpSpace)

% % % 28/07/14	written by wp: work before ICA
% % %	The function was inspired by the original function in FASTICA
% % % toolbox, modified to fit MEG/EEG data, but less flexible

% % % input:	data			(points x channels)
% % %				infoRatio	(the propotion of information shall be explained 0-1)
% % %				numLim4ICs	(The maximum of ICs)
% % %				stepRatio	(When the maximum of ICs cannot be meet, the step to reduce the ratio)
% % %				rsmpSpace	(the resample space)

% % % Output:	D, E			(output of PCA)
% % %				infoRatio	(information explained ratio)
% % %				nICs			(number of components)


	%% 1. check input and set default paras
	% % %	data size check
	sz = size(data);
	if length(sz) > 3
		error('The input data shall be 2-D!');
	elseif sz(1) < sz(2)
		data = data';
		fprintf('Input data may in the wrong shape, now converted, from [%d x %d] to [%d x %d]', ...
			sz(1), sz(2), sz(2), sz(1));
		sz = sz([2, 1]);
	end
	
	% % % resample if needed	
	if nargin > 4 && ~isempty(rsmpSpace)
		data = data(1 : rsmpSpace : end, :);
	end
	
	% % % ICA para defaults
	if nargin < 4 || isempty(stepRatio)
		stepRatio = 0.01;
	end
	if nargin < 3 || isempty(numLim4ICs)
		numLim4ICs = 100;
	end
	if nargin < 2 || isempty(infoRatio)
		infoRatio = 0.99;
	end
	
	%% 2. pca
	% % % Calculate the covariance matrix.
	covarianceMatrix = cov(data, 1);

	% % % Calculate the eigenvalues and eigenvectors of covariance matrix.
	[E, D] = eig (covarianceMatrix);

	% % % Sort the eigenvalues - decending.
	eigenvalues = sort(diag(D), 1, 'descend');
	
	%% 3. select components
	% The rank is determined from the eigenvalues - and not directly by
	% using the function rank - because function rank uses svd, which
	% in some cases gives a higher dimensionality than what can be used
	% with eig later on (eig then gives negative eigenvalues).
	rankTolerance = 1e-7;
	maxLastEig = sum (diag (D) > rankTolerance);
	if maxLastEig == 0,
	  fprintf (['Eigenvalues of the covariance matrix are' ...
			 ' all smaller than tolerance [ %g ].\n' ...
			 'Please make sure that your data matrix contains' ...
			 ' nonzero values.\nIf the values are very small,' ...
			 ' try rescaling the data matrix.\n'], rankTolerance);
	  error ('Unable to continue, aborting.');
	end

	% % %	decide number of components in loop
	nICs = sz(2);
	rExplained = 1;
	while nICs > numLim4ICs
		[nICs, rExplained] = findNumComponents(eigenvalues, infoRatio);
		infoRatio = infoRatio - stepRatio;
	end
	infoRatio = rExplained;
	
	% % % Drop the smaller eigenvalues
	idx = diag(D) > (eigenvalues(nICs) + eigenvalues(nICs + 1)) / 2;

	% % % select subset of the data	
	E = selcol(E, idx);
	D = selcol(selcol(D, idx)', idx);
	
end %end of function

%% sub-function to decide number of ICs
function [nICs, rExplained] = findNumComponents(vector, thresh)

	if isvector(vector) && ~isempty(vector)
		nICs = 1;
		rExplained = 0;
		S = sum(vector);
		while rExplained < thresh
			nICs = nICs + 1;
			rExplained = sum(vector(1 : nICs)) ./ S;
		end
	else
		error('Input data is not a non-empty vector!');
	end
end % end of sub-function 


%% sub-function to select columns, copied from original FASTICA toolbox
function newMatrix = selcol(oldMatrix, maskVector);

	% newMatrix = selcol(oldMatrix, maskVector);
	%
	% Selects the columns of the matrix that marked by one in the given vector.
	% The maskVector is a column vector.

	% 15.3.1998

	if size(maskVector, 1) ~= size(oldMatrix, 2),
	  error ('The mask vector and matrix are of uncompatible size.');
	end

	numTaken = 0;

	for i = 1 : size (maskVector, 1),
	  if maskVector(i, 1) == 1,
		 takingMask(1, numTaken + 1) = i;
		 numTaken = numTaken + 1;
	  end
	end

	newMatrix = oldMatrix(:, takingMask);
end