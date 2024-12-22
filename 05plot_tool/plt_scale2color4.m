% % % rewritten from plt_scale2color @22/05/2024 by wp 
% 	we now try to write the color array with an alpha channel
function [theColor, alpha] = plt_scale2color4(data, threshold, scale, style)

	%% 1. check input
	if nargin < 1	|| isempty(data)
		fprintf('Default: No inputs, I assume you want a legend!\n');
		data = repmat(-100 : 100, [10 1])';
	end
	
	% % % v shall be a vector, reshape input if necessary
	sz = size(data);
	flagReshape = false;
	if sz(2) == 1 && length(sz) == 2
		v = data;
	elseif sz(1) == 1 && length(sz) == 2
		v = data';
		fprintf('... switching input from row to column vector!\n');
	else
		flagReshape = true;
		v = data(:);
	end
	
	if nargin < 2 || isempty(threshold)
		threshold = 0;
		fprintf('Default: threshold 0%% is set!\n');
	end
	
	if nargin < 3 || isempty(scale)
		scale = max(abs(v));
		fprintf('Default: scaled to data maximum!\n');
	end
	
	if nargin < 4 || isempty(style)
		style = 'm';	%default is my colormap yellow-red-gray-blue-cyan
		fprintf('Default: freesurfer colormap used!\n');
	end
	
	%% 2. data preparation
	% % % 2.1	normalize to [-1 1]
	v = v ./ scale;
	v(v > 1) = 1;
	v(v < -1) = -1;
	
	% % % 2.2 alpha channel
	alpha = abs(v);
	alpha(alpha <= threshold) = 0;
	theColor = zeros(prod(sz), 3);
	
	%% 4. fill colors
	switch lower(style)
		% % % 4.1 my colormap yellow-red-gray-blue-cyan
		case {'mymap', 'm'}
			% % % 4.1.1 base color			[0.781 0.762 0.664];	
			c0 = [0.7 0.7 0.7];
% 			c0 = [0.7 0.7 0.7];
			
			% % % 4.1.2 % positive middle to max:	[1 0 0] to [1 1 0]			
			idx = v > 0.5;	
			theColor(idx, 1) = 1;
			theColor(idx, 2) = (2 * v(idx) - 1);
% 			theColor(idx, 3) = 0;

			% % % 4.1.3 % zero to positive middle:	c0 to [1 0 0]			
			idx = v <=0.5 & v >= 0;
			theColor(idx, :) = bsxfun(@times,	2 * v(idx), [1 0 0] - c0) + c0;
			
			% % % 4.1.4 % zero to negative middle:	c0 to [0 0 1]			
			idx = v >=-0.5 & v < 0;
			theColor(idx, :) = bsxfun(@times,	2 * v(idx), c0 - [0 0 1]) + c0;
			
			% % % 4.1.5 % negative middle to negative max:	[0 0 1]	to [0 1 1]		
			idx = v < -0.5;
% 			theColor(idx, 1) = 0;
			theColor(idx, 2) = (-1 - 2 * v(idx));
			theColor(idx, 3) = 1;

		% % % 4.2 color for freesurfer top, red to blue passing gray		
		case {'freemap', 'f', 'free'}	
			c0 = [0.781 0.762 0.664];
			dp = [1 0 0] - c0;
			dn = c0 - [0 0 1];
			
			idx = v >= 0;
			theColor(idx, :) = bsxfun(@plus,v(idx) * dp, c0);
			
			idx = v < 0;
			theColor(idx, :) = bsxfun(@plus,v(idx) * dn, c0);
				
		% % % 4.3 jet colormap, from matlab
		case {'jet', 'j'}
			idx = v >= 0.75;
				theColor(idx, 1) = 1 - (v(idx) - 0.75) * 2; 
				theColor(idx, 2:3) = 0; 

			idx = v < 0.75 & v >= 0.25;
				theColor(idx, 1) = 1; 
				theColor(idx, 2) = 1 - (v(idx) - 0.25) * 2; 
				theColor(idx, 3) = 0; 

			idx = v < 0.25 & v > -0.25;
				theColor(idx, 1) = (v(idx) + 0.25) * 2; 
				theColor(idx, 2) = 1; 
				theColor(idx, 3) = 1 - (v(idx) + 0.25) * 2; 

			idx = v <= -0.25 & v > -0.75;
				theColor(idx, 1) = 0; 
				theColor(idx, 2) = 1 + (v(idx) + 0.25) * 2; 
				theColor(idx, 3) = 1; 

			idx = v <= -0.75;
				theColor(idx, 1) = 0; 
				theColor(idx, 2) = 0; 
				theColor(idx, 3) = 1 + (v(idx) + 0.75) * 2; 
				
			
		% % % 4.4 hot colormap, from many
		case {'hot', 'h'}
			tmp = find(v < 0, 1);
			if ~isempty(tmp)
				fprintf('Warning: you are using hot scale but the values contain negative!\n');
			end
			
			idx = v < 0;
			sid(idx) = 0;
			
			idx = v <= 3/8;
			theColor(idx, 1) = v(idx) * 8/3;
			theColor(idx, 2) = 0;
			theColor(idx, 3) = 0;
			
			idx = v <= 3/4 & v > 3/8;
			theColor(idx, 1) = 1;
			theColor(idx, 2) = (v(idx) - 3/8) * 8/3;
			theColor(idx, 3) = 0;
			
			idx = v > 3/4;
			theColor(idx, 1) = 1;
			theColor(idx, 2) = 1;
			theColor(idx, 3) = (v(idx) - 3/4) * 4;
			
		otherwise
			error('unkown color map!');
	end
		
	clear idx;

	if flagReshape
		theColor = reshape(theColor, [sz, 3]);
		alpha = reshape(alpha, sz);
	end
end %end of function

