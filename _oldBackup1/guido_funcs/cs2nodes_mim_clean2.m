
function [res]=cs2nodes_mim_clean2(cs,F1,F2);
% maximizes absolute value of coherence for all pairs 
% of voxels with respect to both dipole orientations for any 
% pair. 
% usage [res1,res2]=cs2nodes_mim_complex(cs,F1,F2);
% input: 
% cs:  NxN matrix for N channels, complex cross-spectrum, 
% F1: NxMx3 tensor, spatial filters for M voxels, specifies all first sources of the pairs  
% F2: NxPx3 trensor, same as F1 for P voxels, specifies all second sources
%       of the pairs. (To calculate all combinations in a grid, choose F2=F1 
% 
% output: 
% res: MxP matrix of complex coherences
%

regu=.000001;

[nchan ng ndipdir]=size(F1);
ng2=ng;
% if nargin>2
    [nchan2 ng2 ndipdir2]=size(F2);
% end


    res=zeros(ng,ng2);
 % res2=zeros(ng,ng2);

csvoxrealinv=zeros(ndipdir,ndipdir,ng);
csvoxrealinvsqrt=zeros(ndipdir,ndipdir,ng);
csvoxvec=zeros(ndipdir,nchan,ng);
c1=zeros(ng,nchan);
for i=1:ng
    Floc=squeeze(F1(:,i,:));
    csloc=Floc'*real(cs)*Floc;
    csvoxrealinv(:,:,i)=inv(csloc+regu*eye(ndipdir)*mean(diag(csloc)));
    csvoxrealinvsqrt(:,:,i)=sqrt(csvoxrealinv(:,:,i));
    csvoxvec(:,:,i)=Floc'*cs;
end
if nargin>2
csvoxrealinv2=zeros(ndipdir2,ndipdir2,ng2);
csvoxrealinvsqrt2=zeros(ndipdir2,ndipdir2,ng2);
csvoxvec2=zeros(ndipdir2,nchan,ng2);
    for i=1:ng2
        Floc2=squeeze(F2(:,i,:));
%         if i==1;
%             xxx3=size(Floc2)
%             xxx4=size(cs)
%         end
csloc2=Floc2'*real(cs)*Floc2;
    csvoxrealinv2(:,:,i)=inv(csloc2+regu*eye(ndipdir2)*mean(diag(csloc2)));
        %csvoxrealinv2(:,:,i)=inv(Floc2'*real(cs)*Floc2);
        csvoxrealinvsqrt2(:,:,i)=sqrtm(csvoxrealinv2(:,:,i));
        csvoxvec2(:,:,i)=Floc2'*cs;
    end
end


for i=1:ng
%     if round(i/100)*100==i;disp(i);end
    csvoxvecloc=csvoxvec(:,:,i);
%     caainv=csvoxrealinv(:,:,i);
    caainvsqrt=sqrtm(csvoxrealinv(:,:,i));
    
    for j=1:ng2;
%         if nargin==2;
%             Floc2=squeeze(F1(:,j,:));
%         else
            %ng2=ng2
            Floc2=squeeze(F2(:,j,:));
%         end
        cab_comp=csvoxvecloc*Floc2;
%         nphi=101;
        %resloc2=zeros(nphi,1);
          cbbinv=csvoxrealinv2(:,:,j);
          nite=5;
          nite1=5;
          restest=zeros(nite1,1);
          for ite=1:nite1;
              phi=(ite-1)/nite1*pi/2;
           restest(ite)=cohmax(phi,cab_comp,cbbinv,caainvsqrt);
          end
          [resloc iresmax]=max(restest);
          phi=(iresmax-1)/nite1*pi/4;
          dphi=.000001;
           akont=1;
          %restest=restest
              phi=mod(phi+pi/2,pi)-pi/2;
          % phi=pi/4;
        for ite=1:nite;
            %akont=10;
        
            resloc=cohmax(phi,cab_comp,cbbinv,caainvsqrt);
            reslocu=cohmax(phi+dphi,cab_comp,cbbinv,caainvsqrt);
            reslocd=cohmax(phi-dphi,cab_comp,cbbinv,caainvsqrt);
            fprime=(reslocu-reslocd)/(2*dphi);
            fprimeprime=(reslocu+reslocd-2*resloc)/(dphi^2);
            deltaphi=-fprime/(fprimeprime-akont);
            phin=phi+deltaphi;
             phin=mod(phin+pi/2,pi)-pi/2;
            reslocn=cohmax(phin,cab_comp,cbbinv,caainvsqrt);
            %[deltaphi,fprime,fprimeprime]
            if reslocn>resloc;
                akont=akont/2;
                phi=phin;
                resloc=reslocn;
            else
                akont=akont*2;
            end
            %cbbinvsqrt=sqrtm(csvoxrealinv(:,:,j));
           % [ite,phi,resloc,log(akont)/log(2)]
        end 
        
            res(i,j)=abs(resloc)*exp(sqrt(-1)*phi);
%         res2(i,j)=phi;
       
    end
    
      
 end


return;

function resloc=cohmax(phi,cab_comp,cbbinv,caainvsqrt)
  cab=real(exp(sqrt(-1)*phi)*cab_comp);
   X=cab*cbbinv*cab';
            %resloc1(iphi)=trace(caainv*X);
            Y=caainvsqrt*X*caainvsqrt;
            [u,s,v]=svd(Y);
            resloc=sqrt(s(1,1));
            
return
