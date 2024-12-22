function [cout, sid, theMap, cMax] = plt_scale2color(d, th, M, style)
% % % updated 15/08/2019 by wp : updated for another color map : 'f' (freesurfer)
% % % updated 27/08/2014 by wp : updated for ND data
% % % updated 05/08/2014 by wp : more general purpose
% % % modified 19/5/14 wp:
% % %		- >= th instead of >
% % % modified 12/5/14 wp:
% % %		- sid will always compute
% % %		- relative threshold
% % %		- colormap selection

% % % modified 9/5/14 wp
% % %		- color: (dark) blue - cyan - green - yellow - red (dark)

	%% 1. check input
	% % % v shall be a vector
	mFlag = 0;
	sz = size(d);
	if sz(2) == 1 && length(sz) == 2
		v = d;
	elseif sz(1) == 1 && length(sz) == 2
		v = d';
		fprintf('Switch from row to column vectors!\n');
	else
% 		warning('Data input shall be Nx1 vectors!');
		mFlag = 1;
		v = d(:);
	end
	
	if nargin < 4 || isempty(style)
% 		style = 'jet';
% 		style = 'm';
% 		style = 'hot';
		style = 'f';
	end
	
	if nargin < 3 || isempty(M)
		M = max(abs(v));
	end
	
	if nargin < 2 || isempty(th)
		th = 0.2;
	end
	
	%% 2. normalize to [-1 1]
	v = v/M;
	v(v > 1) = 1;
	v(v < -1) = -1;
	
	%% 3. check output
	sid = abs(v) >= th;
	v(~sid) = 0;
	cout = zeros(size(v, 1), 3);
	
	%% 4. fill colors
	switch lower(style)
		case {'freemap', 'f', 'free'}
			c0 = [0.7 0.7 0.7];	%base color
% 			c0 = [0.5 0.5 0.5];	%base color
% % % 			c0 = [0.781 0.762 0.664];	%base color
% % % 			cout = repmat(c0, [numel(v), 1]);
% % % 			cout(abs(v) > 0.5, 2) = (2 * abs(v(abs(v) > 0.5)) - 1) * (1 - c0(2)) + c0(2);
% % % 			cout(v > 0, 1) =  2 * v(v > 0) * (1 - c0(1)) + c0(1);
% % % 			cout(v < 0, 3) =  2 * v(v < 0) * (c0(3) - 1) + c0(3);
% % % 			cout(cout > 1) = 1;
% % % 			cout(cout < -1) = -1;
			idx = v > 0.5;
			cout(idx, 1) = 1;
			cout(idx, 2) = (2 * v(idx) - 1);
			cout(idx, 3) = 0;

			idx = v <=0.5 & v >= 0;
			cout(idx, :) = bsxfun(@times,	2 * v(idx), [1 0 0] - c0) + c0;
			
			idx = v >=-0.5 & v < 0;
			cout(idx, :) = bsxfun(@times,	2 * v(idx), c0 - [0 0 1]) + c0;
			
			idx = v < -0.5;
			cout(idx, 1) = 0;
			cout(idx, 2) = (-1 - 2 * v(idx));
			cout(idx, 3) = 1;

% 			cout(abs(v) > 0.5, 2) = (2 * abs(v(abs(v) > 0.5)) - 1);
% 			cout(v > 0, 1) =  2 * v(v > 0);
% 			cout(v < 0, 3) =  -2 * v(v < 0);
% 			cout(cout > 1) = 1;
% 			cout(cout < -1) = -1;
% 			cout = bsxfun(@times,	cout, abs(v)) + bsxfun(@times, 1 - abs(v), [0.781 0.762 0.664]);
		case {'mymap', 'm'}
			c0 = [0.781 0.762 0.664];
			dp = [1 0 0] - c0;
			dn = c0 - [0 0 1];
			
			idx = v >= 0;
			cout(idx, :) = bsxfun(@plus,v(idx) * dp, c0);
			
			idx = v < 0;
			cout(idx, :) = bsxfun(@plus,v(idx) * dn, c0);
				
		case {'jet', 'j'}
			idx = v >= 0.75;
				cout(idx, 1) = 1 - (v(idx) - 0.75) * 2; 
				cout(idx, 2:3) = 0; 

			idx = v < 0.75 & v >= 0.25;
				cout(idx, 1) = 1; 
				cout(idx, 2) = 1 - (v(idx) - 0.25) * 2; 
				cout(idx, 3) = 0; 

			idx = v < 0.25 & v > -0.25;
				cout(idx, 1) = (v(idx) + 0.25) * 2; 
				cout(idx, 2) = 1; 
				cout(idx, 3) = 1 - (v(idx) + 0.25) * 2; 

			idx = v <= -0.25 & v > -0.75;
				cout(idx, 1) = 0; 
				cout(idx, 2) = 1 + (v(idx) + 0.25) * 2; 
				cout(idx, 3) = 1; 

			idx = v <= -0.75;
				cout(idx, 1) = 0; 
				cout(idx, 2) = 0; 
				cout(idx, 3) = 1 + (v(idx) + 0.75) * 2; 
				
			
		case {'hot', 'h'}
			tmp = find(v < 0, 1);
			if ~isempty(tmp)
				fprintf('Warning: you are using hot scale but the values contain negative!\n');
			end
			
			idx = v < 0;
			sid(idx) = 0;
			
			idx = v <= 3/8;
			cout(idx, 1) = v(idx) * 8/3;
			cout(idx, 2) = 0;
			cout(idx, 3) = 0;
			
			idx = v <= 3/4 & v > 3/8;
			cout(idx, 1) = 1;
			cout(idx, 2) = (v(idx) - 3/8) * 8/3;
			cout(idx, 3) = 0;
			
			idx = v > 3/4;
			cout(idx, 1) = 1;
			cout(idx, 2) = 1;
			cout(idx, 3) = (v(idx) - 3/4) * 4;
		otherwise
			error('unkown color map!');
	end
		
	clear idx;
% 	cout(~sid, :) = nan;
	if mFlag
		cout = reshape(cout, [sz, 3]);
		sid = reshape(sid, sz);
	end
	
	if nargout > 2
		cMax = [-M, M]';
		theMap = permute(repmat(plt_scale2color(-100 : 100, th, [], style(1)), [1 1 10]), [1 3 2]);
	end
end %end of function

