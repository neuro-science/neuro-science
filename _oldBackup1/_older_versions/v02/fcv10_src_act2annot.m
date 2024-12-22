function fcv10_src_act2annot(fname, vc, fc, src_loc, src_data, varargin)
% function p = fcv10_src_act2annot(vc, tri, cv, src, ax, varargin)
% % % written 21/12/15 by wp: for better visualization in freesurfer.

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
			otherwise
				nUnknownPara = nUnknownPara + 1;
				fprintf('\nUnknown parameter #%d: " %s " ...\n', varargin{kCounter});
		end
		kCounter = kCounter + 1;
	end % readin end
	
	%% 2. compute the data
	% % % interpolate raw data to all vc
	if length(src_data) == length(vc)
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
				vout = interp_ungridded(src_loc, vc, 'data', src_data, 'projmethod', 'smudge', 'triout', 'fc');
			otherwise
				error('unknown projection method!\n');
		end
	end
	% % % scale the values to fit the table
	d = vout / p.maxValue;
	clear vout;
	sid = abs(d) < p.subThresh;
	d(d > 1) = 1;
	d(d < -1) = -1;
	d = (d + 1) * 100 + 1;
	d = round(d);
	d(sid) = 0;
	
	
	%% 3. write the file
	% % % construct the color table
	N = 100; %color scale
	ct.numEntries = 2 * N + 1;
	ct.orig_tab = p.color_map_type;
	ct.struct_names = cell(ct.numEntries, 1);
	% % % color map type to 1-letter id
	color_map_id = lower(p.color_map_type(1));
	% % % get the color table values
	cMap = plt03_scale2color(-N : N, 0, [], color_map_id);
	% % % finish the table
	tmp = round(cMap * 255);
	clear cMap;
	tmp(:, 4) = 0; %alpha?
	tmp(:, 5) = tmp(:, 1) + tmp(:, 2) * 2^8 + tmp(:, 3) * 2^16;
	ct.table = tmp;
	clear tmp;
	% % % write it
	write_annotation_modified(fname, d, ct);
end



