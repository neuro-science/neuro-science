function out = cmp03_cs2coupling (data, methods)
% % % rewrote 26/06/14 by wp, for newer format
% % % rewrote 13/05/14 by wp, flag1d added to allow 2d output

% % % data format shall be [chan, chan, X, Y, ...]
% % % method could be 'coh', 'imc' or 'c', 'i'
	%% 1. check data size and input
	% % % size etc.
	sz = size(data);
	if sz < 2
		error('Input data shall be at least 2 dimensions!\n');
	elseif sz(1) == sz(2)
		nChs = sz(1);
	else
		error('The first two dimensions shall be channels!\n');
	end
	
	% % % 	default coherence
	if nargin < 2
		methods = [];	%default complex output
	end
	
	%% 2. do coherence computation
	if length(sz) > 2
		N = prod(sz(3:end));
		data = reshape(data, [sz(1), sz(2), N]);
		tmp = zeros([nChs, nChs, N]);
		for k = 1 : N
			tmp(:, :, k) = cs2coh (data(:, :, k));
		end
		c_out = reshape(tmp, sz);
	else
		c_out = cs2coh (data);
	end

	%% 3. output as coherence or imaginary part
	if strcmpi('coh', methods) || strcmpi('c', methods)
		out = abs(c_out);
	elseif strcmpi('imc', methods) || strcmpi('i', methods)
		out = imag(c_out);
	elseif strcmp('both', methods) || strcmpi('b', methods)
		out{1} = abs(c_out);
		out{2} = imag(c_out);
	end
	
end % end of function

	%% 4. subfunction cs2coh
	function y = cs2coh(c)
	
		p1 = diag(c); %power factor 
		p2 = sqrt(p1 * p1');	%extend  to 2d
		y = c ./ p2;	%do it

	end %end of function