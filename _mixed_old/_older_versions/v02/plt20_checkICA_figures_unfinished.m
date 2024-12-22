function v = plt20_checkICA (myDir, myHeader, v)
% % % written 18/11/16 by wp: this function is to semi-automatic ICA component rejection

	%% 0. preparison
	% % % indication of starting time
	tmp = clock;
	fprintf('\n===========\nThe ICA semi-manual check start at %04d-%02d-%02d %02d:%02d:%02d.\n', round(tmp));
	fprintf('The location of ICA figures are: %s.\n', myDir);
	v.fname = [myDir, '/', myHeader, num2str(tmp(2) * 1000000 + tmp(3) * 10000 + ...
		tmp(4) * 100 + tmp(5), '%8.8d'), '_ica.mat'];
	
	if nargin < 3 || isempty(v)
		fprintf('We start from the beginning: %s.\n', v.dir);
		v.dir = myDir;
		v.hdr = myHeader;
		vFlag = true;
	elseif strcmp(v.dir, myDir) && strcmp(v.hdr, myHeader)
		fprintf('We continue with previous work.');
		vFlag = false;
	else
		fprintf('The previous work was not consistent with this, check again!');
		return;
	end
	
	%% 1. check the folders
	if vFlag
		tmp1 = dir(myDir);
		tmp2 = struct2cell(tmp1);
		f1 = false(size(tmp2, 2), 1);
		for k = 1 : size(tmp2, 2)
			if length(tmp2{1, k}) >= length(myHeader)
				f1(k) = strcmp(myHeader, tmp2{1, k}(1:length(myHeader)));
			end
		end
		v.dn = tmp2(1, f1);
		v.nd = length(v.dn);
		v.nf = zeros(v.nd, 1);
		v.fn = cell(v.nd, 1);
		v.ff = cell(v.nd, 1);
		v.fl = cell(v.nd, 1);
		clear tmp1 tmp2 f1;
	end	
	
	%% 2. check the files
	if vFlag
		for id = 1 : v.nd
			tmp1 = dir([myDir, '/', v.dn{id}, '/*.png']);
			tmp2 = struct2cell(tmp1);
			v.nf(id) = size(tmp2, 2);
			v.fn{id} = strrep(tmp2(1, :), myHeader, [myDir, '/', v.dn{id}, '/', myHeader]);
			v.ff{id} = false(v.nf(id), 1);
			v.fl{id} = char(zeros(v.nf(id), 1, 'uint8'));
		end
	end	
	
	%% 3. work on the data
	if vFlag
		v.cID = 1;
		v.cFD = 1;
	end
	
	fid = figure(1);
	xs = 80; ys = 10;
	set(fid, 'WindowStyle', 'normal', 'visible', 'on', ...
		'position', [xs ys xs +  1000 ys + 700]);
% 	for d = v.cID : v.nd
	d = 1;
		for f = 1 : 10
% 		for f = v.cFD : v.nf{id}
			str = ['open -a Preview ', v.fn{d}{f}];
			unix(str);
			
		end
% 	end
	
end %end of function
