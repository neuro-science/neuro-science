function v = plt_checkICA_matfile (ein, aus, para)
% % % updated 27/07/18 by wp: modify the program
% % % updated 02/12/16 by wp: modify the program to include MS data
% % % written 20/11/16 by wp: this function is to semi-automatic ICA component rejection

	%% 1. preparison
	% % % some paras
	if nargin < 3 || isempty(para)
		icaWaveDispRange = [-1200 3600];
		skip_id = [];
		reductionFlag = 1;
		reductionZ = 1;
	else
		icaWaveDispRange = para.icaWaveDispRange;
		skip_id = para.skip_id;
		try
			reductionFlag = para.reductionFlag;
			reductionZ = para.reductionZ;
		catch
			reductionFlag = 1;
			reductionZ = 1;
		end
	end
	% % % PSYCHTOOLBOX
	KbName('UnifyKeyNames');
	ListenChar(2);
	% % % indication of starting time
	tmp = clock;
	fprintf('\n===========\nThe ICA semi-manual check start at %04d-%02d-%02d %02d:%02d:%02d.\n', round(tmp));
	% % % get the ica data file structure, but not read in
	try
		fn = matfile(ein);
	catch
		fprintf('error in loading data!');
		return;
	end
	% % % check the file existence
	if exist(aus, 'file')
		fprintf('The previous checking exists, trying to load ...\n');
		try
			tmp = load(aus);
			v = tmp.v;
			v_old = v;
			vFlag = false;
			fprintf(' ...now loaded!\n');
		catch
			fprintf(' ...loading failed!\n');
			vFlag = true;
		end
	else
		vFlag = true;
	end
	% % % construct the file if not there
	if vFlag
		tmp = size(fn, 'ica');
		v.nFs = max(tmp);
		if v.nFs ~= prod(tmp)
			fprintf('I wish the ica data to be 1-D!\n');
			return;
		end
		clear tmp;
		v.oc = cell(v.nFs, 1);
		v.label = {'VEOG', 'HEOG', 'EKG', 'LINE', 'MUSCLE', 'JUMP', 'OTHER'};
		v.cf = 1;
		v.cc = 1;
		v_old = [];
	end
	v.ein = ein;
	v.aus = aus;

	%% 2. plot and check
	for is = v.cf : v.nFs
		% % % whether move on
		while KbCheck;	end;
		cFlag = 0;
		fprintf('We are going to work on data %02d of %02d, press <space> to continue and <q> to quit.\n', is, v.nFs);
		while ~cFlag
			[tmp1, tmp2, tmp3] = KbCheck;
			if tmp1
				tmp4 = find(tmp3);
				myKey = KbName(tmp4(1));
				fprintf([myKey, 'pressed!\n']);
				switch myKey
					case 'q'
						v.cf = is;
						v.cc = 1;
						save(aus, 'v', 'v_old');
						ListenChar;
						return;
					case 'space'
						cFlag = 1;
					otherwise
						fprintf('\n');
				end
			end
		end
		% % % load ica data
		tmp = fn.ica(is, :);
		if numel(tmp) > 1
			fprintf('The data size was not as expected, please check!');
			return;
		elseif reductionFlag
			ic = tmp{1};
			x0 = 0.5 : 0.1 : 100;
			[tmp1, tmp2] = min(abs(bsxfun(@minus, ic.xfs', x0)));
			x = ic.xfs(tmp2);
			tmp3 = tmp2(2 : end);
			tmp4 = [tmp2(1 : end - 1); tmp3];
			tmp5 = round(mean(tmp4));
			y = ic.yfs(tmp2, :);
			for k = 2 : length(x) - 1
				y(k, :) = mean(ic.yfs(tmp5(k - 1) : tmp5(k), :));
			end
			ic.xfs = x;
			ic.yfs = y;
			clear tmp1 tmp2 tmp3 tmp5 tmp4 x y x0;
			fprintf('ICA data loaded and power data reduced!\n');
		else
			ic = tmp{1};
			fprintf('ICA data loaded and no reduction was made!\n');
		end

		% % % initialize the output data if not yet
		oc = v.oc{is};
		if ~isfield(oc, 'flag')
			oc.n = ic.nICs;
			if isfield(ic, 'fname')
				oc.fname = ic.fname;
			elseif isfield(ic, 'iPath')
				oc.fname = ic.iPath;
			else
				oc.fname = ['File#', num2str(is)];
% 				error('unknown file ID!');
			end
			oc.infoRatio = ic.infoRatio;
			oc.flag = zeros(oc.n, 1, 'single'); %default space
		else
			fprintf('The data was checked partially before, trying to resume...\n');
			if isfield(ic, 'fname') && strcmp(oc.fname, ic.fname)
				fprintf('consistency check succeed!\n');
			elseif isfield(ic, 'iPath') && strcmp(oc.fname, ic.iPath)
				fprintf('consistency check succeed!\n');
			else
				fprintf('consistency check failed!\n');
				fprintf('previous input file name %s, \n', oc.fname);
				fprintf('loaded input file name %s, \n', ic.fname);
				return;
			end
		end
		% % % plot the component templates
		myIn = {'fColor', [1 1 1], 'cvRangeX', icaWaveDispRange, ...
			'cvTmpPtsNum', ic.nxts, 'cvNum', size(ic.yts, 2), 'tpRangeXY', 0.5};
		fp = plt_icaFigParaSet(myIn{:});
		fp.tpTmpPtsNum = round(fp.tpTmpPtsNum / reductionZ);
		ff = plt_icaFigureTemplate(fp);		% plot template
		set(ff.fid, 'visible', 'on');
		ttt = clock;
		fprintf('ICA template prepared @%02d:%02d:%02.0f!\n\n', ttt(4:6));

		% % % plot for each ic and check
		if is == v.cf
			ii = v.cc;
		else
			ii = 1;
		end

		while ii <= oc.n
			tic;
			if ismember(oc.flag(ii), skip_id)
				if oc.flag(ii)
					fprintf('According to pre-defined criteria, we would skip IC #%02d as<%s>.\n', ii, v.label{oc.flag(ii)});
				else
					fprintf('According to pre-defined criteria, we would skip IC #%02d as<CLEAN>.\n', ii);
				end
				ii = ii + 1;
				continue;
			end
			% % % plot the ICA component
			if isfield(para, 'loc2d')
				plt_plotICA_matfile (ff, fp, ic, ii, para.loc2d);
			else
				plt_plotICA_matfile (ff, fp, ic, ii);
			end			
			% % % do check
			while KbCheck(-1);
			end
			cFlag = 0;
			tFlag = 1;
			while ~cFlag
				[tmp1, tmp2, tmp3] = KbCheck;
				if tmp1
					tmp4 = find(tmp3);
					myKey = KbName(tmp4(1));
					switch myKey
						case 'q'
							v.cf = is;
							v.cc = ii;
							v.oc{is} = oc;
							save(aus, 'v', 'v_old');
							ListenChar;
							return;
						case 'b'
							ii = ii - 2;
							cFlag = 1;
							tFlag = 0;
							if ii + 1 > 0
								fprintf('We would go back to IC #%02d!\n', ii + 1);
							else
								ii = 0;
								fprintf('We could not go back to stone age, no go to IC #%02d!\n', ii + 1);
							end
						case 'v'
							oc.flag(ii) = 1;
							cFlag = 1;
						case 'h'
							oc.flag(ii) = 2;
							cFlag = 1;
						case 'k'
							oc.flag(ii) = 3;
							cFlag = 1;
						case 'l'
							oc.flag(ii) = 4;
							cFlag = 1;
						case 'm'
							oc.flag(ii) = 5;
							cFlag = 1;
						case 'j'
							oc.flag(ii) = 6;
							cFlag = 1;
						case 'n'
							oc.flag(ii) = 7;
							cFlag = 1;
						case 'space'
							oc.flag(ii) = 0;
							cFlag = 1;
						otherwise
							fprintf('\n');
					end
				end
			end
			ttt = clock;
			if tFlag && oc.flag(ii)
				fprintf('IC #%02d of %d was labeled <%s> on @%02d:%02d:%02.0f!\n', ii, ic.nICs, v.label{oc.flag(ii)}, ttt(4:6));
			elseif tFlag
				fprintf('IC #%02d of %d was labeled <CLEAN> on @%02d:%02d:%02.0f!\n', ii, ic.nICs, ttt(4:6));
			end
			ii = ii + 1;
		end % end of all components
		v.oc{is} = oc;
		close all;
		save(aus, 'v', 'v_old');
	end % end of all data files

	%% 3. save the data
	v.cf = is;
	v.cc = ii - 1;
	save(aus, 'v', 'v_old');
	ListenChar;
end %end of function
