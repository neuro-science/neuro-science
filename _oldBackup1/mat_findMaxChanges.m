	function [A, As] = mat_findMaxChanges(X, nn)
	% % % 04/07/17 written by wp
		X1 = X([2:end, end]);
		X3 = abs(X1 - X);
		[tmp, As] = sort(X3, 'descend');
		if nargin < 2 || isempty(nn)
			A = As(1);
		else
			id = find(As <= nn);
			A = As(id(1));
		end
		clear X X1 X2 X3;
	end
