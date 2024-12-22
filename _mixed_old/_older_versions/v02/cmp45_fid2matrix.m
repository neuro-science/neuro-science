function A = cmp45_fid2matrix(X, Z)
% % % 14/03/2015	written by wp
% % % Z - [x1 y1 z1 x2 y2 z2 x3 y3 z3]

	%% 1. prepare
	% % % initialize
	str = {'x1', 'y1', 'z1', 'x2', 'y2', 'z2', 'x3', 'y3', 'z3'};
	for k = 1 : 9
		cmd = [str{k}, '= Z(', num2str(k), ');'];
		eval(cmd);
	end

	% % % origin
	x0 = (x2 + x3)/2;
	y0 = (y2 + y3)/2;
	z0 = (z2 + z3)/2;
	
	% % %
	
end %end of function
