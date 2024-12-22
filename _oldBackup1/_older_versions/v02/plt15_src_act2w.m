function plt15_src_act2w(fname, src_data, src_loc, vc, varargin)
% % % written 23/12/15 by wp: for better visualization in freesurfer.

	%% 1. prepare
	% % % default parameters
	p.intMethod = 'smudge';
	p.color_map_type = 'jet';
	p.intplPara = 20;
   p.subThresh = 0.2;
	p.maxValue = max(src_data(:));
	% % % input parameter overwrite defaults
	nVarArgIn = length(varargin);
	kCounter = 1; %arg in counter
	nUnknownPara = 0;	%number of unkown input, display only
	while kCounter < nVarArgIn
		switch lower(varargin{kCounter})
			case 'ip'
				p.intplPara = varargin{kCounter + 1};
				kCounter = kCounter + 1;
			case 'th'
				p.subThresh = varargin{kCounter + 1};
				kCounter = kCounter + 1;
			case 'mx'
				p.maxValue = varargin{kCounter + 1};
				kCounter = kCounter + 1;
			case 'cm'
				p.color_map_type = varargin{kCounter + 1};
				kCounter = kCounter + 1;
			case 'im'
				p.intMethod = varargin{kCounter + 1};
				kCounter = kCounter + 1;
			case 'fc'
				p.fc = varargin{kCounter + 1};
				kCounter = kCounter + 1;
			otherwise
				nUnknownPara = nUnknownPara + 1;
				fprintf('\nUnknown parameter #%d: " %s " ...\n', varargin{kCounter});
		end
		kCounter = kCounter + 1;
	end % readin end
	
	%% 2. compute the data
	% % % interpolate raw data to all vc
	if nargin < 3 || (length(src_data) == length(vc))
		vout = src_data;
	else
		switch p.intMethod
			case 'gaussian'
				vout = plt07_gaussian_dist(vc, src_loc, p.intplPara, src_data);
			case 'nearest'
				vout = interp_ungridded(src_loc, vc, 'data', src_data, 'projmethod', 'nearest');
			case 'sphere_avg'
				vout = interp_ungridded(src_loc, vc, 'data', src_data, 'projmethod', 'sphere_avg', 'sphereradius', p.intplPara);
			case 'sphere_weighteddistance'
				vout = interp_ungridded(src_loc, vc, 'data', src_data, 'projmethod', 'sphere_weighteddistance', 'sphereradius', p.intplPara);
			case 'smudge'
				vout = interp_ungridded(src_loc, vc, 'data', src_data, 'projmethod', 'smudge', 'triout', p.fc);
			otherwise
				error('unknown projection method!\n');
		end
	end
	% % % scale the values to fit the table
	
	
	%% 3. write the w file
	% % % construct the color table
	err = write_wfile(fname, vout);
end





% % % ================= below is from freesurfer ============== % % %
% write_wfile.m
%
% Original Author: Doug Greve
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
%

function err = write_wfile(fname, w, v)
	err = 1;

	if(nargin ~= 2 & nargin ~= 3)
	  fprintf('USAGE: err = write_wfile(fname, w, <v>) \n');
	  return;
	end

	vnum = length(w) ;

	% Handle when v not given or is empty %
	if(exist('v') ~= 1) v = []; end
	if(isempty(v)) v = [0:vnum-1]; end

	% open it as a big-endian file
	fid = fopen(fname, 'wb', 'b') ;
	if(fid == -1)
	  fprintf('ERROR: could not open %s\n',fname);
	  return;
	end

	fwrite(fid, 0, 'int16') ;
	fwrite3(fid, vnum) ;
	for i=1:vnum
	  fwrite3(fid, v(i)) ;          % vertex number (0-based)
	  fwrite(fid,  w(i), 'float') ; % vertex value
	end

	fclose(fid) ;
	err = 0;
	return;
end



