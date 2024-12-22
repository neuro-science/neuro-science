function a = fcv06_fs2mat(subDir, subName, srcNums)

% % % 28/11/14	wp write for sources

	%% 1. prepare file names and paths etc.
	nSrcs = length(srcNums);
	
	%% 2. read the src locations
	for ic = 1 : nSrcs
		% % % read source data
		fName = [subDir, subName, '/bem/', subName, '-ico-', num2str(srcNums(ic)), '-src.fif'];
		tmp1 = ft_read_headshape(fName, 'format', 'mne_source');
		clear fName;
		% % % get it
		if strcmp(tmp1.unit, 'm')
			tmp2 = tmp1.pnt * 1000;
		else
			error('wrong format of source space!');
		end
		cmd = ['a.pnt', num2str(srcNums(ic)), ' = tmp2;', ...
		'a.tri', num2str(srcNums(ic)), ' = int32(tmp1.tri);', ...
		'a.idx', num2str(srcNums(ic)), ' = logical(tmp1.orig.inuse);'];
		eval(cmd);
	end
	
	%% 3. read white surface
	a.vcWhite = tmp1.orig.pnt * 1000;
	a.fcWhite = int32(tmp1.orig.tri);
	a.nnWhite = tmp1.orig.nn;

end

