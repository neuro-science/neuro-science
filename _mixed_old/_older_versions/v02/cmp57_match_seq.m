function [i1, i2, di1, di2] = cmp57_match_seq(d1, d2, vf)
% % % 13/04/2016	written by wp
	d1 = d1(:);
	d2 = d2(:);
	n1 = length(d1);
	n2 = length(d2);
	dd = abs(bsxfun(@minus, d1, d2')); %[n1, n2]
	M = max(dd(:)) + 1;
	
	if nargin < 3
		vf = 0;
	end
	
	if n1 < n2
		i1 = 1 : n1;
		for ii =1 : n1
			[tmp1, tmp2] = min(dd(ii, :));
			i2(ii) = tmp2;
			dd(:, tmp2) = M;
		end
	else
		i2 = 1 : n2;
		for ii =1 : n2
			[tmp1, tmp2] = min(dd(:, ii));
			i1(ii) = tmp2;
			dd(tmp2, :) = M;
		end
	end	
	
	if vf
		if n1 < n2
			[tmp3, tmp4] = sort(d2);
			[la, lb] = ismember(tmp4, i2);
			lb = lb(logical(lb));
			tmp4(la) = [];
			tmp5 = mean(d2(i2));
			tmp6 = mean(d1);
			ct = 0;
			while (tmp5 < tmp6) && (vf < 2) && vf
				ct = ct + 1;
				tmp7 = i2(lb(ct));
				i2(lb(ct)) = tmp4(end - ct + 1);
				if (mean(d2(i2))) < tmp5
					vf = 0;
					i2(lb(ct)) = tmp7;
				else
					tmp5 = mean(d2(i2));
				end
			end
			if ct > 0
				vf = 0;
				ct = 0;
			end
			while (tmp5 > tmp6) && (vf > -2) && vf
				ct = ct + 1;
				tmp7 = i2(lb(end - ct + 1));
				i2(lb(end - ct + 1)) = tmp4(ct);
				if (mean(d2(i2))) > tmp5
					vf = 0;
					i2(lb(end - ct + 1)) = tmp7;
				else
					tmp5 = mean(d2(i2));
				end
			end
		else
			[tmp3, tmp4] = sort(d1);
			[la, lb] = ismember(tmp4, i1);
			tmp4(la) = [];
			tmp5 = mean(d1(i1));
			tmp6 = mean(d2);
			ct = 0;
			while (tmp5 < tmp6) && (vf > -2) && vf
				ct = ct + 1;
				tmp7 = i1(lb(ct));
				i1(lb(ct)) = tmp4(end - ct + 1);
				if (mean(d1(i1))) < tmp5
					vf = 0;
					i1(lb(ct)) = tmp7;
				else
					tmp5 = mean(d1(i1));
				end
			end
			if ct > 0
				vf = 0;
				ct = 0;
			end
			while (tmp5 > tmp6) && (vf < 2) && vf
				ct = ct + 1;
				tmp7 = i1(lb(end - ct + 1));
				i1(lb(end - ct + 1)) = tmp4(ct);
				if (mean(d1(i1))) > tmp5
					vf = 0;
					i1(lb(end - ct + 1)) = tmp7;
				else
					tmp5 = mean(d1(i1));
				end
			end
		end
	end
	
	if nargout > 2
		di1 = d1(i1);
		di2 = d2(i2);
	end
	clear tmp1 tmp2 tmp4 tmp 3 tmp5 tmp6 tmp7 ct la lb vf dd M n1 n2 d1 d2;
	
end % end of function
