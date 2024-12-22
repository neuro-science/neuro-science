function roi = fcv07_lb2roi(src_label, roi_label)
% % % 01/12/14 wp write for roi from ids by manual selection

	roi = [];
	for k = 1 : numel(roi_label);
		roi = [roi; find(src_label == roi_label(k))];
	end
end