function [F1,p,dipori]=getdipdir(cs,A);

cs=real(cs);
[nchan ng ndum]=size(A);
p=zeros(ng,1);
F1=zeros(nchan,ng);
dipori=zeros(ndum,ng);
for i=1:ng;
    Aloc=squeeze(A(:,i,:));
    csloc=Aloc'*cs*Aloc;
    [u s v]=svd(csloc);
    F1(:,i)=Aloc*u(:,1);
    p(i)=s(1,1);
    dipori(:,i)=u(:,1);
end

return;
    