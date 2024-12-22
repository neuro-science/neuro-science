	function [F, p]= mat_FThreshFromP(n, dn, p0, delta, tStart)
	% % % 22/06/17 written by wp

		%% 1. check paras
		if nargin < 2
			fprintf('You need input the numerator and denominator!\n')
			return;
		end
		
		if nargin < 3 || isempty(p0)
			p0 = 0.05;
		end
		
		if nargin < 4 || isempty(delta)
			delta = 1e-5;
		end
		
		if nargin < 5 || isempty(tStart)
			FStart = 6.8;
		end
				
		scale = 0.1;
		%% 2. do the recursive loop
		goFlag = 1;
		F = FStart;
		while goFlag
			p = (1 - fcdf(F, n, dn));
			if p - p0 > delta
				F = F * (1 + scale);
			elseif p - p0 < -delta
				F = F * (1 - scale);
			else
				goFlag = 0;
			end
		end
	end
