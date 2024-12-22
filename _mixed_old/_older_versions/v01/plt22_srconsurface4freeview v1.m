function plt22_srconsurface4freeview(fpath, fName_suffix, data, ana, para)
% % % written 26/01/17 by wp: generate activation file for freeview
% % % need <plt07_gaussian_dist>

	%% 1. check inputs and prepare
	if nargin < 4 || nargin > 5
		fprintf('Input Error!\nWe need 4-5 inputs: FilePath, FileSuffix, SourceData, Anatomy and distPara[2]!\n');
		return;
	end
	
	if nargin == 4
		para = 2;
	end
	
	if abs(max(size(data) - [300 1])) > 0.1
		fprintf('Input Error!\nI support only 300x1 voxels!\n');
		return;
	end
		
	lname = [fpath, 'a_lh_', fName_suffix, '.act'];
	rname = [fpath, 'a_rh_', fName_suffix, '.act'];

	
	%% 2. put 300x1 data into 324x1 source locations
	idx = find(ana.idx20);
	if abs(length(idx) - 324) > 0.5
		error('Ana data size mismatch!');
	end
	d0 = zeros(324, 1);
	id02 = true(324, 1);
	id02(ana.null2_idx) = false;
	d0(id02) = data;
	ld0 = d0(1 : 162);
	rd0 = d0(163 : end);
	lsrc = ana.vcFL(idx(1 : 162), :);
	rsrc = ana.vcFR(idx(163 : end) - ana.leftLastVertex, :);
	
	%% 3. compute the data
	ld = plt07_gaussian_dist (ana.vcFL, lsrc, para, ld0);
	rd = plt07_gaussian_dist (ana.vcFR, rsrc, para, rd0);
	
	%% 4. write to file
	write_curv(lname, ld, size(ana.fcFL, 1));
	write_curv(rname, rd, size(ana.fcFR, 1));
	
	%% 5. clean up
	clear lname rname ld ld0 rd rd0 para d0 id02 idx lsrc rsrc;
end

function [curv] = write_curv(fname, curv, fnum)
% writes a curvature vector into a binary file [Copy From FreeSurfer]
%				fname - name of file to write to
%				curv  - vector of curvatures
%				fnum  - # of faces in surface.
%
% Original Author: Bruce Fischl
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:13 $
%    $Revision: 1.3 $
%
% Copyright ?? 2011 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu

	% open it as a big-endian file
	fid = fopen(fname, 'wb', 'b') ;
	vnum = length(curv) ;
	NEW_VERSION_MAGIC_NUMBER = 16777215;
	fwrite3(fid, NEW_VERSION_MAGIC_NUMBER ) ;
	fwrite(fid, vnum,'int32') ;
	fwrite(fid, fnum,'int32') ;
	fwrite(fid, 1, 'int32');
	fwrite(fid, curv, 'float') ;
	fclose(fid) ;
end
