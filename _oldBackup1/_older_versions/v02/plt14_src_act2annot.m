function plt14_src_act2annot(fname, vc, fc, src_loc, src_data, varargin)
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



%%%%%%%%%%%%%%%%%% Below is copywrite info from original function %%%%%%%%%%%%%%%%%%%%%%
%
% Contact ythomas@csail.mit.edu or msabuncu@csail.mit.edu for bugs or questions 
%
%=========================================================================
%
%  Copyright (c) 2008 Thomas Yeo and Mert Sabuncu
%  All rights reserved.
%
%=========================================================================
function write_annotation_modified(filename, label, ct)

	% write_annotation(filename, vertices, label, ct)
	%
	% Only writes version 2...
	%
	% label is the vector of annotation
	% ct is a struct
	% ct.numEntries = number of Entries
	% ct.orig_tab = name of original ct
	% ct.struct_names = list of structure names (e.g. central sulcus and so on)
	% ct.table = n x 5 matrix. 1st column is r, 2nd column is g, 3rd column
	% is b, 4th column is flag, 5th column is resultant integer values
	% calculated from r + g*2^8 + b*2^16 + flag*2^24. flag expected to be all 0

	%% 1. first write vertices and label
	% % % open file
	fp = fopen(filename, 'w', 'b');

	% % % test the file/size
	nVs = length(label); %number of vertices
	count = fwrite(fp, int32(nVs), 'int');
	if(count~=1)
		error('write_annotation: Writing #vertices/labels not successful!!');
	end

	% % % prepare stuff to write
	temp = zeros(nVs * 2, 1);
	temp(1:2:end) = 0 : nVs - 1;
	temp(2:2:end) = label;
	temp = int32(temp);

	% % % do write 
	count = fwrite(fp, int32(temp), 'int');
	if(count~=length(temp))
		error('write_annotation: Writing labels/vertices not successful!!');
	end

	%% 2. ct and more
	% % % Write that ct exists
	count = fwrite(fp, int32(1), 'int');
	if(count~=1)
		error('write_annotation: Unable to write flag that ct exists!!');
	end

	% % % write version number
	count = fwrite(fp, int32(-2), 'int');
	if(count~=1)
		 error('write_annotation: Unable to write version number!!');
	end

	% % % write number of entries
	count = fwrite(fp, int32(ct.numEntries), 'int');
	if(count~=1)
		 error('write_annotation: Unable to write number of entries in ct!!');
	end

	% % % write original table
	orig_tab = [ct.orig_tab char(0)];
	count = fwrite(fp, int32(length(orig_tab)), 'int');
	if(count~=1)
		 error('write_annotation: Unable to write length of ct source!!');
	end

	count = fwrite(fp, orig_tab, 'char');
	if(count~=length(orig_tab))
		 error('write_annotation: Unable to write orig_tab!!');
	end

	% % % write number of entries
	count = fwrite(fp, int32(ct.numEntries), 'int');
	if(count~=1)
		 error('write_annotation: Unable to write number of entries in ct!!');
	end

	% % % write ct
	for i = 1:ct.numEntries
		 count = fwrite(fp, int32(i-1), 'int');
		 if(count~=1)
			  error('write_annotation: Unable to write structure number!!');
		 end

		 structure_name = [ct.struct_names{i} char(0)];
		 count = fwrite(fp, int32(length(structure_name)), 'int');
		 if(count~=1)
			  error('write_annotation: Unable to write length of structure name!!');
		 end
		 count = fwrite(fp, structure_name, 'char');
		 if(count~=length(structure_name))
			  error('write_annotation: Unable to write structure name!!');
		 end

		 count = fwrite(fp, int32(ct.table(i, 1)), 'int');
		 if(count~=1)
			 error('write_annotation: Unable to write red color'); 
		 end

		 count = fwrite(fp, int32(ct.table(i, 2)), 'int');
		 if(count~=1)
			 error('write_annotation: Unable to write blue color'); 
		 end

		 count = fwrite(fp, int32(ct.table(i, 3)), 'int');
		 if(count~=1)
			 error('write_annotation: Unable to write green color'); 
		 end

		 count = fwrite(fp, int32(ct.table(i, 4)), 'int');
		 if(count~=1)
			 error('write_annotation: Unable to write padded color'); 
		 end 

	end
	fclose(fp);
end




