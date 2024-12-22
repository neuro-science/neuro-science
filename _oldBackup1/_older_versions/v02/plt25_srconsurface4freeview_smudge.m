function d = plt24_srconsurface4freeview_smudge(fpath, fName, data, sourceIndices, fc, depth, thresh)
% % % updated 07/02/17 by wp: we use more generic parameters
% % % written 26/01/17 by wp: generate activation file for freeview
% % % need <plt07_gaussian_dist>

	%% 1. check inputs and prepare
	if nargin < 5 || nargin > 7
		fprintf('Input Error!\nWe need 4-5 inputs: FilePath, Filename, SourceData, sourceIndices, faces and diepth[6], thresh[0]!\n');
		return;
	end
	
	if nargin < 7
		thresh = 0;
	end
	if nargin < 6
		depth = 6;
	end
			
	fname = [fpath, fName, '.act'];
	
	%% 3. compute the data
	nFC = length(fc);
	d = zeros(size(sourceIndices));
	d(sourceIndices) = data;
	d = smudge(d, fc, depth, thresh);
	
	%% 4. write to file
	write_curv(fname, d, nFC);
	
	%% 5. clean up
	clear depth nFC;
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
