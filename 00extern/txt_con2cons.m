% % % turn category texts into all condtions, connecting with symbol
% % % Input in shall be cell arrays, content of which are also cell arrays
function out = txt_con2cons (in, symbol)

	if nargin < 2
		symbol = '-';
	end
	
	n = numel(in);
	for i = n : -1 : 1
		sz(i) = numel(in{i});
	end
	
	out = cell(sz);
	
	for i = prod(sz) : -1 : 1
		[j{1:n}] = ind2sub(sz, i);
		tmp = in{n}{j{n}};
		for k = n-1 : -1 : 1
			tmp = [in{k}{j{k}}, symbol, tmp];
		end
		out{i} = tmp;
	end

end