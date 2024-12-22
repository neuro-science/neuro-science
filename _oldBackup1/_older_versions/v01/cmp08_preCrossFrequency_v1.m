function d = cmp08_preCrossFrequency(data, fA, fB, type)
% % % 25/01/16	modified by wp(null) 
% % % 08/07/14	written by wp: 
% % % 	data: nTimePoints x ... x ...
% % %		fA, fB: filter parameters
% % %		type: 'p', 'l', 'low' or 'phase' for the low frequency;
% % %				'a', 'h', 'high' or 'amplitude' for the high frequency;

	%% 1. check inputs
	% % % check output type	
	type_char = lower(type(1));
	switch type_char
		case {'p', 'l'}
			cosFlag = 1;
		case {'a', 'h'}
			cosFlag = 0;
		otherwise
			error('unknown types!');
	end
	% % % data format	
	sz = size(data);
	if length(sz) > 2
		data = reshape(data, [sz(1), prod(sz(2:end))]);
	end
	data = double(data);
	
	%% 2. work on the data
	data = FiltFiltM(fB, fA, data, 1);
	data = hilbert(data);
	d = abs(data);
	if cosFlag
		d = real(data) ./ d;
	end
	
	%% 3. return
	d = reshape(d, sz);
		
end %end of function

