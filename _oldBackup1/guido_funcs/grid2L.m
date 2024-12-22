function L=grid2L(grid,fp);

vtest=forward_general(zeros(1,6),fp);
nchan=length(vtest);
[ns ndum]=size(grid);
L=zeros(nchan,ns,3);

u=eye(3);
for i=1:3;
    uloc=u(i,:);
    dips=[grid,repmat(uloc,ns,1)];
    L(:,:,i)=forward_general(dips,fp);
end

return


