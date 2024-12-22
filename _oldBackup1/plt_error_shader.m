function h = plt_error_shader(dm, de, x0, C, av)
% % % written by wp 04/04/2018 - shaded error bar for mean values
	if nargin < 2
		error('We need the mean and error!');
	else
		dm = dm(:)';
		de = de(:)';
		y = [dm + de, fliplr(dm - de)];
	end
	if nargin < 3 || isempty(x0)
		x0 = 1 : numel(dm);
	else
		x0 = x0(:)';
	end
	CFlag = false;
	if nargin < 4 || isempty(C)
		C = 'b';
	elseif ~ischar(C)
		Cn = C;
		C = 'b';
		CFlag = true;
	end
	if nargin < 5 || isempty(av)
		av = 0.2;
	end
	x = [x0, fliplr(x0)];
	h = patch(x, y, C,'EdgeColor','none', 'FaceAlpha', av);
	if CFlag
		set(h, 'FaceColor', Cn);
	end
% 	tmp = get(h, 'FaceColor');
% 	set(h, 'EdgeColor', tmp);
% 	patch
end %end of function