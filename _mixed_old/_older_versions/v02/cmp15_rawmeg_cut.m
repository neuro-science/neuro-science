function [data, idx, idt] = cmp15_rawmeg_cut(data, trg, theRange)
% % % 25/07/14	written by wp: 

	% % % data size
	[npts, nchs] = size(data);
	if npts < nchs
		fprintf(fid, '\nWarning: I assume your input data is "chans x points".\n');
		fprintf(fid, 'I will change it to "points x chans"...\n');
		data = data';
		fprintf(fid, 'Done.\n\n');
	end
	% % % trg size
	trg = trg(:);
	
	% % % idx generate
	id1 = theRange(1) : theRange(2);
	id2 = bsxfun(@plus, trg, id1);
	idx = unique(id2(:));
	[tmp, idt] = ismember(id2', idx);
	
	% % % data reduction	
	data = data(idx, :);
end %end of function

