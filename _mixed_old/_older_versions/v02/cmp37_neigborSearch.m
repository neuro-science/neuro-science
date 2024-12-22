
function nb = cmp37_neigborSearch(data, method, th)

% % % This function defines the neigborhoods of a node
% % % written in 2013 and updated 26/09/14 by wp

	%% check inputs
	if nargin < 2
		error('You need to specify the connection data, definition method!\n');
	end
	
	%% do for 1 of 3 methods
	switch lower(method)
		% % % topological relation definition
		case {'topo', 't'}
			% % % get nodes and initialize
			src = sort(unique(data(:)));
			ns = length(src);
			if ~src(1)
				if src(end) == ns - 1
					src = [src(2 : end); ns];
					data = data + 1;
				else
					error('Some nodes are left alone!');
				end
			elseif src(end) ~= ns
				error('Some nodes are left alone!');
			end
			nb = cell(ns, 1);
			parfor k = 1 : ns
				id3 = ismember(data, src(k));
				id1 = any(id3, 2);
				tmp = data(id1, :);
				nb{k} = unique(tmp(:));
			end
			
		% % % distance in 3D
		case {'distance3d', 'd3'}
			ns = size(data, 1);
			nb = cell(ns, 1);
			if nargin < 3 || isempty(th)
				th = 17.5;
			end
			x1 = permute(data, [1 3 2]);
			x2 = permute(data, [3 1 2]);
			d = sqrt(sum(bsxfun(@minus, x1, x2).^2, 3));
			for k = 1 : ns
				nb{k} = find(d(k, :) <= th);
			end
			
		% % % distance < th1 in 1D, other 2D shall be identical in th2 cm resolution
		case {'distance1d', 'd1'}
			ns = size(data, 1);
			nb = cell(ns, 1);
			if nargin < 3 || isempty(th) 
				th(1) = 15;
			elseif length(th) < 2
				th(2) = 5;
			end
			x1 = permute(data, [1 3 2]);
			x2 = permute(data, [3 1 2]);
			d = abs(bsxfun(@minus, x1, x2));
			d1 = all(d < th(1), 3);
			d0 = sum(d < th(2), 3) >= 2;
			for k = 1 : ns
				nb{k} = find(d1(k, :) & d0(k, :));
			end
		otherwise
			error('Please specify the method to define neighbors: topo, distance3d or distance1d !');
	end
	
end