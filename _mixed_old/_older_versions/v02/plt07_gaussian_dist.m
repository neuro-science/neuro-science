function vout = plt07_gaussian_dist (vertices, sources, dThresh, data)
% % % updated 08/08/2014 by wp, for more general purpose
% % % a function to interpolate data for surface graph
% % % checked 9/5/14 wp
% % % It was adopted from Guido's program, the original was attached below for respect
% % % Compare to the orginal form, it need more memory but a little faster
% % % an option to output matrix is also provided 

	% % % check input	
	if nargin > 4 || nargin < 2
		error('input shall be 2 (vertices, grids) or 3(+ distance threshold) or 4 (+data)!');
	elseif nargin == 2
		dThresh = 2;
	end
	
	% % % compute the transform matrix	
	alpha2 = (dThresh / sqrt(log(2))) .^ 2;
	dist2 = sum(bsxfun(@minus, permute(vertices, [1 3 2]), permute(sources, [3 1 2])).^2, 3);
	w = exp(-dist2 / alpha2);
	
	% % % decide output
	if nargin == 3
		vout = w;
	else
		vout = bsxfun(@rdivide, w * data, sum(w, 2));
	end

end %end of function

    

% % % The following is from guido
% function vout = spatfiltergauss (vin, grid, d, grid2)
% % makes a spatialfilter with a Gaussian function
% % usage: vout=spatfiltergauss(vin,grid,d)
% % or:   vout=spatfiltergauss(vin,grid,d,grid2)
% % 
% % if only one grid is provide the smoothed field is calcuated on the 
% % the same grid as for original field
% %
% % input: vin   NxK matrix, for N voxels and K fields, each column is a field (e.g. power) 
% %              defined on grid provided as the second argument
% %        grid  Nx3 matrix, grid locations of original data, 
% %        d     width of Gaussian, at distance d the Gaussian drops to 50%
% %                of its maximum
% %        grid2 (optional) Mx3 matrix of grid locations where the smoothed field is calculated
% %               if not provided, it is set to grid.
% %
% % output: vout NxK (or MxK, if grid2 is provided) matrix of smoothed
% %               fields. 
% %         
% 
% if nargin==3;
%     
% [ng ns]=size(vin);
% alpha=d/sqrt(-log(.5));
% vout=vin;
% for i=1:ng
%     r0=grid(i,:);
%     rd=grid-repmat(r0,ng,1);
%     dis=sqrt(sum(rd.^2,2));
%     w=exp(-dis.^2/alpha^2);
%     for j=1:ns
%     vout(i,j)=sum(w.*vin(:,j))/sum(w);
%     end
% end
% 
% elseif nargin==4
%     
%     [ng ns]=size(vin);
%     alpha=d/sqrt(-log(.5));
% 
%     [ng2 ndum]=size(grid2);
%     vout=zeros(ng2,ns);
% 
% 
%     for i=1:ng2
%         r0=grid2(i,:);
%         rd=grid-repmat(r0,ng,1);
%         dis=sqrt(sum(rd.^2,2));
%         w=exp(-dis.^2/alpha^2);
%         for j=1:ns
%         vout(i,j)=sum(w.*vin(:,j))/sum(w);
%         end
%     end
%     
% end
% 
% 
% 
% return
    