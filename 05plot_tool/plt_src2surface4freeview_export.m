function plt_src2surface4freeview_export(varargin)
% % % written 22/08/24 by wp: generate activation file for freeview with weighted average

	%% 1. input check and prepare
	
	% Initialize default values
	fname = '/Users/wang/Documents/FreeSurfer/fsaverage/surf/_tmp.act';
	faces = [];
	src = [];
	theMatrix = [];
	depth = 21;
	alpha = 1;
	method = 'Barycentric';

	% Parse input arguments
	for i = 1 : 2 : length(varargin)
		switch lower(varargin{i})
			case 'src'
				src = varargin{i+1};
			case 'thematrix'
				theMatrix = varargin{i+1};
			case 'faces'
				faces = varargin{i+1};
			case 'node_ids'
				node_ids = varargin{i+1};
				if islogical(node_ids)
					node_ids = find(node_ids);
				end
			case 'fname'
				fname = varargin{i+1};
			case 'method'
				 method = varargin{i+1};
			case 'depth'
				 depth = varargin{i+1};
			case 'alpha'
				 alpha = varargin{i+1};
			otherwise
				error('Unknown parameter name: %s', varargin{i});
		end
	end

	% % % default path and name
	if ismac
		addpath('/Applications/freesurfer/matlab/');
		fpath = '/Users/wang/Documents/FreeSurfer/fsaverage/';
	elseif isunix
		addpath('/home/pwang/8TL/freesurfer/matlab/');
		fpath = '/home/pwang/0FS/fsaverage/';
	else
		error('We support only mac and linux');
	end

	if isempty(src)
		error('I need input data!');
	end

	if isnumeric(fname)
		fname = [fpath, 'surf/_tmp', num2str(fname, '%02d'), '.act'];
	else
		fname = [fpath, 'surf/', fname, '.act'];
	end
	
	%% 2. do computation
	% % % get M if not there	
	if ~isempty(theMatrix)

		[num_faces, num_nodes] = size(theMatrix);
		
	else
		if isempty(faces) || isempty(node_ids)
			error('I need input triagulation and source ids, while no matrix data!');
		end
		
		num_nodes = numel(node_ids);
		num_faces = size(faces, 1);
		
		theMatrix = plt_src2surface4freeview_matrix(faces, node_ids, method, depth, alpha);
	end
	
	% % % check data size	
	if abs(numel(src) - num_nodes) > 0.5 
		if abs(numel(src) - 2 * num_nodes) > 0.5
			error('input data mismatch in size!');
		else
			fprintf('spliting data ...\n');
			dFlag = true;
			src_left = src(1 : num_nodes);
			src_right = src(num_nodes + 1 : 2 * num_nodes);
		end
	else
		dFlag = false;
		fprintf('single data.\n');
	end
	
	% % % here it is
	if strcmpi(method(1), 'b')	%	Barycentric
		sumM = sum(theMatrix, 2);

		if dFlag
			act_left = theMatrix * src_left ./ sumM;
			act_left(~sumM) = 0;

			write_curv(strrep(fname, '.act', '_l.act'), act_left, num_faces);

			act_right = theMatrix * src_right ./ sumM;
			act_right(~sumM) = 0;

			write_curv(strrep(fname, '.act', '_r.act'), act_right, num_faces);
		else
			act = theMatrix * src ./ sumM;
			act(~sumM) = 0;

			write_curv(fname, act, num_faces);
		end
	elseif strcmpi(method(1), 'b')	%	Barycentric
		act = theMatrix * src;
		write_curv(fname, act, num_faces);
	else
		error('Unsupported interpolation method.');
	end
		
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
