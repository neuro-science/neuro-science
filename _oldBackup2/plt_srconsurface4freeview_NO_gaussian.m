function d = plt_srconsurface4freeview_NO_gaussian(fname, data, source_idx, vc, nFC)
% % % updated 12/02/24 by wp: default to freesurfer than sy05
% % % updated 12/06/18 by wp: remove redundant paras
% % % updated 15/08/17 by wp: rename without numbers
% % % updated 07/02/17 by wp: we use more generic parameters
% % % written 26/01/17 by wp: generate activation file for freeview
% % % need <plt07_gaussian_dist>

	%% 1. check inputs and prepare
	% % % input number check
	if nargin < 5 || nargin > 6
		fprintf('Input Error!\nWe need 5-6 inputs: Filename, SourceData, sourceIndex, vertices, numOfFaces and distPara[2]!\n');
		return;
	end
	
	% % % default path and name
	if ismac
		addpath('/Applications/freesurfer/matlab/');
	elseif isunix
		addpath('/home/pwang/8TL/freesurfer/matlab/');
	end
	
	if isempty(fname)
		fname = '/Users/wang/Documents/FreeSurfer/fsaverage/surf/_tmp.act';
	elseif isnumeric(fname)
		fname = ['/Users/wang/Documents/FreeSurfer/fsaverage/surf/_tmp', num2str(fname, '%02d'), '.act'];
	end
	
	%% 2. check data
	if islogical(source_idx)
		idx = find(source_idx);
	else
		idx = source_idx;
	end
	if abs(length(idx) - length(data)) > 0.5
		error('Ana data size mismatch!');
	end
% 	src = vc(idx, :);
	
	%% 3. compute the data
	d = zeros(size(vc, 1), 1) + nan;
	d(idx) = data;
	
	%% 4. write to file
	write_curv(fname, d, nFC);
	
	%% 5. clean up
	clear d0 idx src;
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
