function [F, p, stats] = mat_rmANOVA2Array(data, FACTNAMES)

% % % wp 06/02/2019 rewrote from rm_anova2 by Aaron Schurger
% % % The purpose is to simplify for n-D arrays
% % % renamed for convenience - wp 30/01/2019
% % % modified for nan values - wp 05/07/2016 dirty change, not accurate
% % % data(F1, F2, SUB, others)
% % % FACTNAMES {'F1', 'F2'}

	%% 1. prepare
	% % % 1.1	check inputs
	if nargin < 2 || isempty(FACTNAMES) || numel(FACTNAMES)~= 2
% 		fprintf('No proper factor names provided, F1 and F2 are used instead!\n');
		FACTNAMES = {'F1', 'F2'};
	end
	% % % 1.2	data check
	sz = size(data);
	if numel(sz) > 4
		data = reshape(data, [sz(1:3), prod(sz(4:end))]);
		sFlag = true;
	elseif numel(sz) < 3
		fprintf('dimension of data is smaller than 3, eixting...');
		return;
	else
		sFlag = false;
	end
	[a, b, n, x] = size(data);
	
	%% 2. computation 	
	AB = nansum(data, 3); %[a b 1 x]
	AS = nansum(data, 2); %[a 1 n x]
	BS = nansum(data, 1); %[1 b n x]
	
	A = nansum(AB, 2); % [a 1 1 x]
	B = nansum(AB, 1); % [1 b 1 x]
	S = nansum(AS, 1); % [1 1 s x]
	T = nansum(A); % [1 1 1 x]

	% degrees of freedom (numbers)
	dfA = a-1;
	dfB = b-1;
	dfAB = (a-1)*(b-1);
	dfS = n-1;
	dfAS = (a-1)*(n-1);
	dfBS = (b-1)*(n-1);
	dfABS = (a-1)*(b-1)*(n-1);

	% bracket terms (expected value)
	expA = nansum(A.^2)./(b*n);	%[1 1 1 x]
	expB = nansum(B.^2)./(a*n);	%[1 1 1 x]
	expAB = nansum(nansum(AB.^2))./n;	%[1 1 1 x]
	expS = nansum(S.^2)./(a*b);	%[1 1 1 x]
	expAS = nansum(nansum(AS.^2))./b;	%[1 1 1 x]
	expBS = nansum(nansum(BS.^2))./a;	%[1 1 1 x]
	expY = nansum(nansum(nansum(data.^2))); 	%[1 1 1 x]
	expT = T.^2 / (a*b*n);	%[1 1 1 x]

	% sums of squares
	ssA = expA - expT;
	ssB = expB - expT;
	ssAB = expAB - expA - expB + expT;
	ssS = expS - expT;
	ssAS = expAS - expA - expS + expT;
	ssBS = expBS - expB - expS + expT;
	ssABS = expY - expAB - expAS - expBS + expA + expB + expS - expT;
	ssTot = expY - expT;

	% mean squares
	msA = ssA / dfA;
	msB = ssB / dfB;
	msAB = ssAB / dfAB;
	msS = ssS / dfS;
	msAS = ssAS / dfAS;
	msBS = ssBS / dfBS;
	msABS = ssABS / dfABS;

	%% 3. results 	
	% f statistic
	fA = msA ./ msAS;
	fB = msB ./ msBS;
	fAB = msAB ./ msABS;

	F = squeeze(cat(1, fA, fB, fAB))';
	if sFlag
		F = reshape(F, [sz(4:end) 3]);
	end
	
	% p values
	if nargout > 1
		pA = 1-fcdf(fA,dfA,dfAS);
		pB = 1-fcdf(fB,dfB,dfBS);
		pAB = 1-fcdf(fAB,dfAB,dfABS);	
		p = squeeze(cat(1, pA, pB, pAB))';
		if sFlag
			p = reshape(p, [sz(4:end) 3]);
		end
	end

	% other return values
	if nargout > 2
		if x == 1
			stats = {'Source','SS','df','MS','F','p';...
						FACTNAMES{1}, ssA, dfA, msA, fA, pA;...
						FACTNAMES{2}, ssB, dfB, msB, fB, pB;...
						[FACTNAMES{1} ' x ' FACTNAMES{2}], ssAB, dfAB, msAB, fAB, pAB;...
						[FACTNAMES{1} ' x Subj'], ssAS, dfAS, msAS, [], [];...
						[FACTNAMES{2} ' x Subj'], ssBS, dfBS, msBS, [], [];...
						[FACTNAMES{1} ' x ' FACTNAMES{2} ' x Subj'], ssABS, dfABS, msABS, [], []};
		else
			fprintf('More dimensions with formula is not convenient, bye!\n')
		end
	end
			
end