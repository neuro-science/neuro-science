function f605_merge_plots(dpath, cutFlag, flipColorFlag, theSuffix, Ds)
	% % % re-write 04/02/2019 to be more convenient
	
	%% 1. prepare 
	if nargin < 1
		fprintf('We need at least a folder name!\n');
	end
	
	if nargin < 2 || isempty(cutFlag)
		cutFlag = true;
	end
		
	if nargin < 3 || isempty(flipColorFlag)
		flipColorFlag = false;
	end
		
	if nargin < 4	|| isempty(theSuffix)
		sh = {'_l_lh', '_r_rh'};
		ss = {'_inner', '_outer'}; 
	else
		sh = theSuffix(1:2);
		ss = theSuffix(3:4);
	end
	
	% % % av5 @04/02/2019@pro80, default
	D1 = 105;
	D2 = 85;
	D3 = 280;
	D4 = 280;
	D5 = 0;
	if nargin < 5	|| isempty(Ds)
		try tmp = char(java.net.InetAddress.getLocalHost.getHostName);
			if strcmp(tmp(1:6), 'iMac80')	%newer @04/03/2019
				D1 = 241;
				D2 = 200;
				D3 = 660;
				D4 = 660;
				D5 = 0;		
	% % % 				D1 = 181;
	% % % 				D2 = 150;
	% % % 				D3 = 500;
	% % % 				D4 = 500;
	% % % 				D5 = 0;		
			end
		catch
		end
	else
		for ik = 1 : 5
			cmd = ['D', num2str(ik), ' = Ds(', num2str(ik), ');'];
			eval(cmd);
		end
	end
	
	%% 2. read
	% % % list of files	
	dr = dir([dpath, '*', sh{1}, ss{1}, '.png']);
	for ic = 1 : length(dr)
		% % % prefix		
		preFix = [dpath, dr(ic).name(1:end-(numel(sh{1}) + numel(ss{1}) + 4))];
		% % % 	initialize	
		v = cell(2, 2);
		sz = zeros(3, 2, 2);
		% % % 	read	data
		for ih = 1 : 2
			for is = 1 : 2
				v{is, ih} = imread([preFix, sh{ih}, ss{is}, '.png']);
				sz(:, is, ih) = size(v{is, ih});
			end
		end
		% % % 	fill blank	
		SZ = max(reshape(sz, [3 4]), [], 2);
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
		try
			preFix = strrep(preFix, '/sm_', '/../bbb/');
		catch
			preFix = strrep(preFix, '/sm', '/../bbb/');
		end
		imwrite(vv, [preFix, '.png'], 'png');
	end

end %end of function