function [lbl, clt, snm] = fcv01_fs_annt2mat(subjectDir, ann_Name)
% % % 04/07/14 wp re-write for more general use
% % % 19/06/14	wp write for source and labels

	%% 1. prepare data
	fname = [subjectDir, '/label/lh.', ann_Name, '.annot'];
	[v1, l1, c1] = read_annotation(fname);
	if isempty(l1 == c1.table(1, 5))
		l1(~l1) = c1.table(1, 5);
	end
	clear fname;
	fname = [subjectDir, '/label/rh.', ann_Name, '.annot'];
	[v2, l2, c2] = read_annotation(fname);
	if isempty(l2 == c2.table(1, 5))
		l2(~l2) = c2.table(1, 5);
	end
	clear fname;
	
	%% 2. conversion
	nLeft = size(c1.table, 1);
	lb_left = zeros(nLeft, 1);
	for k = 1 : nLeft
		lb_left(l1 == c1.table(k, 5)) = k;
		c1.table(k, 4) = k;
	end
	nRight = size(c2.table, 1);
	lb_right = zeros(nRight, 1);
	for k = 1 : nRight
		lb_right(l2 == c2.table(k, 5)) = k + nLeft;
		c2.table(k, 4) = k + nLeft;
	end
	lbl = [lb_left; lb_right];
	snm = [strcat('lh_', c1.struct_names); strcat('rh_', c2.struct_names)];
	clt = [c1.table; c2.table];

end