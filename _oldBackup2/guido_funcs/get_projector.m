function [P,pats]=get_projector(data,fs,artifact_freqs)
% calculated a projector to remove this weird artifact
%
% usuage [P,pats]=get_projector(data,fs,artifact_freqs)
% Input: 
% data: NxM data matrix for N time points and M channels
% fs: sampling rate in Hz
% artifact_freqs: row vector containing all frequencies (in Hz) used to estimate artifact
%                  pattern. E.g. artifact_freqq=[55:95] uses all
%                  frequencies between 55 and 95 Hz. 
%
% Output: 
% P: the projector, an MxM matrix for M channels. Then data*P is the clean
%     data
% pats: Mx2 matrix, two spatial patterns of of the artifact, which were used to construct P.
%



[n,nchan]=size(data);
ff=artifact_freqs+1;
ffmax=max(ff);

data=detrend(data,'constant');
segleng=fs;segshift=segleng/2;epleng=n;
[cs coh]=data2cs_event(data,segleng,segshift,epleng,ffmax);

csm=mean(cs(:,:,artifact_freqs),3);

[u,s,v]=svd(csm);
pats=[real(u(:,1)),imag(u(:,1))];

pro=pats*inv(pats'*pats)*pats';
P=eye(nchan)-pro;

return;
