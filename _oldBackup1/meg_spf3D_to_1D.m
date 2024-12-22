function A1 = meg_spf3D_to_1D(A3, C)
% % % 13/06/18	updated for more filters, good for multiple sessions
% % % 04/07/14	updated for only filter ori
% % % 24/06/14 updated by wp for complex inputs
% % % 20/06/14 modified a little by wp for parallel computation
% % % original code by guido as [F1,p,dipori] = getdipdir(cs, A), listed
% % % below to show respect

	% % % C - N x N			power at certain time points/condition etc. N = nchans
	% % % A3 - N x M x 3		spatial filter M = nvxs

	%% 1. check inputs
	if iscell(A3)
		fprintf('The inputs are cell arrays, multiple filters may exist!\n');
		sFlag = false;
	else
		fprintf('The inputs are numeric arrays, only one filter exists!\n');
		sFlag = true;
	end
	
	%% 2. prepare data
	if sFlag
		C = real(C);
		[nchs, nvxs, ndms] = size(A3); %get size
		if ndms ~= 3
			error('spatial filter shall be nChs x nVXs x 3 in size!');
		else
			A3 = permute(A3, [1, 3, 2]);	%Nx3xM
		end
		A1 = zeros(nvxs, nchs);
	else
		n = numel(A3);
		for ic = 1 : n
			C{ic} = real(C{ic});
			[nchs, nvxs, ndms] = size(A3{ic}); %get size
			if ndms ~= 3
				error('spatial filter shall be nChs x nVXs x 3 in size!');
			else
				A3{ic} = permute(A3{ic}, [1, 3, 2]);	%Nx3xM
			end
			A1{ic} = zeros(nvxs, nchs);
		end
	end	
	%% 2. do computation
	if sFlag
		parfor k = 1 : nvxs;
			c1 = A3(:, :, k)' * C * A3(:, :, k);
			[u, s, v] = svd(c1);
			A1(k, :) = A3(:, :, k) * u(:, 1);
		end
	else
		for k = 1 : nvxs
			cc = zeros(3, 3, n);
			for ic = 1 : n
				c(:, :, ic) = A3{ic}(:, :, k)' * C{ic} * A3{ic}(:, :, k);
			end
			c1 = mean(c, 3);
			[u, s, v] = svd(c1);
			for ic = 1 : n
				A1{ic}(k, :) = A3{ic}(:, :, k) * u(:, 1);
			end
		end
	end
end
    
%% original code by Guido Nolte
% % % function [F1,p,dipori]=getdipdir(cs,A);
% % % 
% % % cs=real(cs);
% % % [nchan ng ndum]=size(A);
% % % p=zeros(ng,1);
% % % F1=zeros(nchan,ng);
% % % dipori=zeros(ndum,ng);
% % % for i=1:ng;
% % %     Aloc=squeeze(A(:,i,:));
% % %     csloc=Aloc'*cs*Aloc;
% % %     [u s v]=svd(csloc);
% % %     F1(:,i)=Aloc*u(:,1);
% % %     p(i)=s(1,1);
% % %     dipori(:,i)=u(:,1);
% % % end
% % % 
% % % return;
% % %     