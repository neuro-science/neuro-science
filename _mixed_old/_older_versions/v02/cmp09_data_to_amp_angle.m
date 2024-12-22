function [amp, ang] = cmp09_data_to_amp_angle (data, fA, fB)
% % % 13/07/14	written by wp: 
% % % 	data: nTimePoints x ... x ...
% % %		fA, fB: filter parameters

	%% 1. check inputs
	% % % data format	
	sz = size(data);
	if length(sz) > 2
		data = reshape(data, [sz(1), prod(sz(2:end))]);
	end
	data = double(data);
	
	%% 2. work on the data
	data = FiltFiltM(fB, fA, data, 1);
	data = hilbert(data);
	data = reshape(data, sz);
	amp  = abs(data);
	ang = angle(data);
	if length(sz) > 2
		amp = reshape(amp, sz);
		ang = reshape(ang, sz);
	end
		
end %end of function

