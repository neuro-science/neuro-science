function plt24_merge4BrainPlots(fein, faus, cutFlag, flipColorFlag)

	% % % av5 @25/05/2018
	D1 = 131;
	D2 = 70;
	D3 = 160;
	D4 = 160;
	D5 = 40;

% 	% % % sy5 @30/05/2018
% 	D1 = 60;
% 	D2 = 10;
% 	D3 = 80;
% 	D4 = 70;
% 	D5 = 70;

	% % % check input
	if nargin < 4
		flipColorFlag = 0;
	end
	if nargin < 3
		cutFlag = true;
	end
	if numel(fein) ~= 4
		error('Now I support only four figures');
	end
	
	% % % get data
	v = cell(2, 2);
	sz = zeros(3, 4);
	for k = 1 : 4
		v{k} = imread(fein{k});
		sz(:, k) = size(v{k});
	end
	SZ = max(sz, [], 2);
	for k = 1 : 4
		if any(SZ - sz(:, k))
			v{k}(end:SZ(1), end:SZ(2), :) = 255;
		end
	end
	
	% % % crop if needed before 21/12/2017 sy05
	if cutFlag
		for k = 1 : 4
			v{k}([1:D1, end-D2:end], :, :) = [];
			v{k}(:, [1:D3, end-D4:end], :) = [];
		end
		v{1}(:, 1:D5, :) = [];
		v{4}(:, 1:D5, :) = [];
		v{2}(:, end-D5+1 : end, :) = [];
		v{3}(:, end-D5+1 : end, :) = [];
	end
	
% % % 	% % % crop if needed
% % % 	if cutFlag
% % % 		for k = 1 : 4
% % % 			v{k}([1:135 435:end], :, :) = [];
% % % 			v{k}(:, [1:175, 595:end], :) = [];
% % % 		end
% % % 		v{1}(:, 1:35, :) = [];
% % % 		v{4}(:, 1:35, :) = [];
% % % 		v{2}(:, end-34 : end, :) = [];
% % % 		v{3}(:, end-34 : end, :) = [];
% % % 	end
	% % % merge data
	vv = [v{1, 1}, v{1, 2}; v{2, 1}, v{2, 2}];
	
	% % % flip background color black->white
	if flipColorFlag
		sz = size(vv);
		vv = reshape(vv, [prod(sz(1:2)), sz(3)]);
		if flipColorFlag > 0
			X = any(vv, 2);
			vv(~X, :) = 255;
		else
			X = any(255 - vv, 2);
			vv(~X, :) = 0;
		end			
		vv = reshape(vv, sz);
	end
	
	% % % write to file
	imwrite(vv, faus, 'png');
	
% 	% % % online test
% 	clf;image(vv);axis equal;size(vv)
end %end of function