function h = plt_scatter_3d(data, varargin)
% % % updated 12/08/2017 - data only

	h = scatter3(data(:, 1), data(:, 2), data(:, 3), varargin{:});
	axis equal;
	
end %end of function