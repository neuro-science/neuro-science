function sa_out=mk_sa_meg_mri(sa,fn);
sa_out=sa;
mri = ft_read_mri(fn);

cfg           = [];
cfg.coordsys  = 'ctf';
cfg.output    = 'brain';
segmentedmri  = ft_volumesegment(cfg, mri);

cfg = [];
cfg.method='singleshell';
vol = ft_prepare_headmodel(cfg, segmentedmri);

%save vol vol segmentedmri 

xx=vol.bnd.pnt;
for i=1:3;xx(:,i)=xx(:,i)-mean(xx(:,i));end
xrad=mean(sqrt(sum(xx.^2,2)));
if xrad>20;
    warning('units mm?')
    vol.bnd.pnt=vol.bnd.pnt/10;
elseif xrad<.2
     warning('units m?')
    vol.bnd.pnt=vol.bnd.pnt*100;
end



% % % ======== here I made a change ========= % % %
% % %	template = ft_read_mri('T1.nii');	%original
% % % template = ft_read_mri('single_subj_T1.nii');	%worked on local computer
template = ft_read_mri('~/8TL/fieldtrip/external/spm8/templates/T1.nii');	%worked on server
% % % ======== here I made a change ========= % % %

template.coordsys = 'spm';

mygrid=sa.grid_medium;

[xgrid,ygrid,zgrid]=grid2cube(mygrid);

%mygrid=[-mygrid(:,1),-mygrid(:,2),mygrid(:,3)];
[xgrid,ygrid,zgrid]=grid2cube(mygrid);

cfg = [];
cfg.grid.warpmni   = 'yes';


mytemplate_grid.pos=mygrid;
mytemplate_grid.xgrid=xgrid;
mytemplate_grid.ygrid=ygrid;
mytemplate_grid.zgrid=zgrid;
mytemplate_grid.dim=[length(xgrid) length(ygrid) length(zgrid)];


mytemplate_grid.inside=1:length(mygrid);
mytemplate_grid.outside=[];
cfg.grid.template=mytemplate_grid;
cfg.grid.nonlinear = 'yes'; % use non-linear normalization
cfg.mri            = mri;
grid               = ft_prepare_sourcemodel(cfg);

grid1=grid.pos;
grid2=mygrid;

r1=mean(grid1);
r2=mean(grid2);
ng=length(grid1);
grid1m=grid1-repmat(r1,ng,1);
grid2m=grid2-repmat(r2,ng,1);

A=inv(grid1m'*grid1m)*(grid1m'*grid2m)
Ainv=inv(A);
r=r1-r2*Ainv;

ng=length(sa.grid_coarse);
sa_out.grid_coarse_indi=sa_out.grid_coarse*Ainv+repmat(r,ng,1);

ng=length(sa.grid_medium);
sa_out.grid_medium_indi=sa_out.grid_medium*Ainv+repmat(r,ng,1);

ng=length(sa.grid_fine);
sa_out.grid_fine_indi=sa_out.grid_fine*Ainv+repmat(r,ng,1);

ng=length(sa.grid_xcoarse);
sa_out.grid_xcoarse_indi=sa_out.grid_xcoarse*Ainv+repmat(r,ng,1);

ng=length(sa.grid_cortexhippo);
sa_out.grid_cortexhippo_indi=sa_out.grid_cortexhippo*Ainv+repmat(r,ng,1);

ng=length(sa.grid_cortex3000);
sa_out.grid_cortex3000_indi=sa_out.grid_cortex3000*Ainv+repmat(r,ng,1);

sa_out.vc_indi.vc=vol.bnd.pnt;
sa_out.vc_indi.tri=vol.bnd.tri;

sa_out.trafo.u_indi2template=A;
sa_out.trafo.r_indi2template=-r1*A+r2;
sa_out.trafo.readme='xtemplate=xindi*u+r, row-vectors';

return;



