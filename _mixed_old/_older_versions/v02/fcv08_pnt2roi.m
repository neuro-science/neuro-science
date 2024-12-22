function [r_id, lb_name] = fcv08_pnt2roi(src_idx, roi_label, roi_name)
% % % 01/12/14 wp write for roi from ids by manual selection

	if nargin < 2
		error('Both source indices and roi labels needed!');
	else
		r_id = roi_label(src_idx);
		if nargout > 1 && nargin >= 3
			lb_name = repmat({'unknown'}, size(r_id));
			lb_name(logical(r_id)) = roi_name(r_id(logical(r_id)));
		end
	end
end