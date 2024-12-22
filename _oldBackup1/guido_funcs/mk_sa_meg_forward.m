function sa_out=mk_sa_meg_forward(sa,fn);


p1=15;
p2=10;

sens = ft_read_sens(fn);
   
sa_out=sa;

n=length(sens.label);
inds=[];

%save sens sens
for i=1:n;
    x=sens.label{i};
    if x(1)=='M';
        inds=[inds,i];
    end
end



%inds=inds
%length(inds)
%length(sa.inds)
% if norm(sa.inds-inds)>.01
%     error('you use a different MEG system')
% end

nchan=length(inds);
sens1=[sens.coilpos(inds,:),sens.coilori(inds,:)];
refs=[sens.coilpos(inds+nchan+min(inds)-1,:),-sens.coilori(inds+nchan+min(inds)-1,:)];

[vc,center,radius,coeffs]=pointsonsurface(sa.vc_indi.vc(:,1:3),sa.vc_indi.vc(:,1:3),p1);

sa_out.fp_indi=meg_ini(vc,center',p2,sens1,refs); 

sa_out.locs_3D_indi=sens.coilpos(inds,:);
%sa_out.ori_3D_indi=sens.coilori(inds,:);
sa_out.sens_indi=sens;
sa_out.coils_indi=sens1;
sa_out.inds=inds;

if nchan==-271;
    sa_out.locs_2D=sa_out.locs_2D_271;
    sa_out.locs_2D_sparse=sa_out.locs_2D_sparse_271;
elseif nchan==-273;
    sa_out.locs_2D=sa_out.locs_2D_273;
    sa_out.locs_2D_sparse=sa_out.locs_2D_sparse_273;
else
    para1.rot=90;
     locs_2D=mk_sensors_plane(sa_out.locs_3D_indi,para1);
     para2.rot=90;
     para2.nin=50;
     locs_2D_sparse=mk_sensors_plane(sa_out.locs_3D_indi,para2);
     sa_out.locs_2D=locs_2D;
     sa_out.locs_2D_sparse=locs_2D_sparse;
end
sa_out=rmfield(sa_out,'locs_2D_271');
sa_out=rmfield(sa_out,'locs_2D_sparse_271');
sa_out=rmfield(sa_out,'locs_2D_273');
sa_out=rmfield(sa_out,'locs_2D_sparse_273');

%sa_out.ctf_sens=sens;    
    
return;



