function p = cmp24_value2p (H0, d)
% % % written 03/09/2014 by wp

	H = sort([H0(:); d], 'descend');
	p = find(H==d)./length(H);
	
end
