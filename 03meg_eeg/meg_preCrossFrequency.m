function d = meg_preCrossFrequency(data, fA, fB, type)
% % % 02/11/19	modified by wp(add source sign consideration for angles) 
% % % 12/01/17	modified by wp(add angles, changed the types) 
% % % 25/01/16	modified by wp(null) 
% % % 08/07/14	written by wp: 
% % % 	data: nTimePoints x ... x ...
% % %		fA, fB: filter parameters
% % %		'legacy'
% % %		type: 'c' or 'cosine' for the cosine of low frequency phase;
% % %				'A' or 'Amplitude' for the high frequency Amplitude;
% % %				'a' or 'angle' for the low frequency angle


	%% 1. check inputs
	sz = size(data);
	if length(sz) > 2
		data = reshape(data, [sz(1), prod(sz(2:end))]);
	end
	data = double(data);
	
	%% 2. work on the data
	data = FiltFiltM(fB, fA, data, 1);
	data = hilbert(data);
	switch type(1)
		case {'c'}
			d = real(data) ./ abs(data);
		case {'A'}
			d = abs(data);
		case {'a'}
			d = angle(data);
	% % % In source reconstruction, when spatial filter was converted from 3D to 1D, 
	% % % the sign is arbitrary during svd,which would cause psuedo double peak
	% % % To solve this problem, the phase was stretched and then averaged.
		case {'s'}	
			d = angle(data .* data);
		otherwise
			error('unknown types!');
	end
	
	%% 3. return
	d = reshape(d, sz);
end %end of function

