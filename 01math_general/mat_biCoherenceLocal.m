% % % The code was adopted from Guido's function data2bs_univar, which was
% % % quoted below to show respect. The modification is to increase efficency
% % % for the format/structure of my own data - wp@2018-12-20

function [cs, pwr] = mat_biCoherenceLocal (data, nf)
	%% 1. prepare the data
	% % % power4power
	NP = 3;
	% % % get the size and reshape
	[nch, npt, ntr] = size(data);
	data = reshape(permute(data, [2 3 1]), [npt, ntr*nch]);
	% % % output size
	nf0 = floor(npt / 4);
	if nargin < 2 || isempty(nf)
		nf = nf0;
	elseif nf > nf0
		fprintf('I don''t think it a good idea: you ask too more frequencies!');
		return;
	end
	
	%% 2. real computation
	% % % 	detrend and fft
	dfft = fft(bsxfun(@times, detrend(data), hanning(npt)));
	dfft = dfft(1 : 2*nf - 1, :);
	% % % below check spectrum 
	ctmp = myHK(conj(dfft), nf, ntr, nch);
	% % % get cs	
	tmp = bsxfun(@times, permute(dfft(1 : nf, :), [1 3 2]), permute(dfft(1 : nf, :), [3 1 2]));
	cs = permute((mean(reshape(tmp, [nf, nf, ntr, nch]) .* ctmp, 3)), [1 2 4 3]);
	clear tmp ctmp;
	% % % get power for f1/f2	
	pwr0 = (abs(dfft).^ NP);
	p1 = mean(reshape(pwr0(1:nf, :), [nf, ntr, nch]), 2) .^ (1/NP);
	% % % get power for f1+f2-1
	ptmp = myHK(pwr0, nf, ntr, nch);
	p2 = permute(mean(ptmp, 3) .^ (1/NP), [1 2 4 3]);
	% % % merge
	pwr = bsxfun(@times, p1, permute(p1, [2 1 3])) .* p2;
end

% % % function to get (f1 + f2 - 1)
function out = myHK(data, nf, ntr, nch)
	% % % realize with indices - faster!!
	ij = bsxfun(@plus, (1 : nf).', (0 : nf - 1));
	out = reshape(data(ij, :), [nf, nf, ntr, nch]);
	clear tmp tmp1 data;	
% % % 	% % % realize with loop - slower
% % % 	for k1 = nf : -1 : 1
% % % 		for k2 = nf : -1 : 1
% % % 			out(k1, k2, :) = data(k1 + k2 - 1, :);
% % % 		end
% % % 	end
% % % 	out = reshape(out, [nf nf ntr nch]);
end

% % % ============================================ Guido's function ============================================
% % % function [cs,csnr,nave]=data2bs_univar(data,segleng,segshift,epleng,maxfreqbins,para)
% % % % calculates bispectrum  from data in general for event-related measurement
% % % % as univariate measures, i.e. always within each sensor. 
% % % %
% % % % usage: [cs,csnr,nave]=data2bs_univar(data,segleng,segshift,epleng,maxfreqbins,para);
% % % %
% % % % input: 
% % % % data: ndat times nchan matrix each colum is the time-series in one
% % % %             channel;
% % % % segleng: length of each segment in bins, e.g. segleng=1000;  
% % % % segshift: numer of bins by which neighboring segments are shifted;
% % % %           e.g. segshift=segleng/2 makes overlapping segments
% % % % epleng: leng of each epoch
% % % % maxfreqbins: maximum frequency in bins, starting at zeros Hertz. 
% % % %             The frequency resolution, df, is given by the physical length of a
% % % %             segment, say T. Then df=1/T. E.g. if T=2 seconds, the maxfreqbins=101
% % % %             means that the maximum physical frequency is 50 Hertz.
% % % % para: structure which is eventually used later
% % % %
% % % % output: 
% % % % cs: nchan  by nf by nf tensor for nf frequencies (i.e. nf=maxfreqbins)   
% % % %  cs(i,f1,f2)=<x(f1)_i*x(f2)_i*conj(x(f1+f2-1)_i)>
% % % %  where  x is the Fourier-transform of the data of each segment
% % % %
% % % % csn:  corresponding normalization factor defined by 
% % % %       csn(i,f1,f2)=N_i(f1) N_i(f2) N_i(f1+f2-1);
% % % %       where N_p(f) is defined as (<abs(x(f)_p)^3>)^(1/3) 
% % % %       Bicoherence can be calculated as cs./csn 
% % % %
% % % % nave: number of averages
% % % 
% % % [ndat,nchan]=size(data);
% % % nf=maxfreqbins;
% % % 
% % % 
% % % mywindow=repmat(hanning(segleng),1,nchan);
% % % if nargin>5
% % %     if isfield(para,'fave')
% % %      % fave=para.fave;
% % %      end;
% % %     if isfield(para,'mywindow');
% % %        mywindow=repmat(para.mywindow,1,nchan);
% % %     end
% % % end
% % % 
% % % 
% % % 
% % % 
% % % 
% % % nep=floor(ndat/epleng);
% % % 
% % % nseg=floor((epleng-segleng)/segshift)+1; %total number of segments
% % % 
% % % 
% % %  cs=zeros(nchan,nf,nf);
% % %  csnr=zeros(nchan,nf,nf);
% % %  csn=zeros(nchan,2*nf-1);
% % %  %csloc=zeros(nf,nf);
% % % 
% % % %figure;plot(mywindow);
% % % nave=0;
% % % for j=1:nep;
% % %     %disp(j)
% % %     dataep=data((j-1)*epleng+1:j*epleng,:);
% % %     for i=1:nseg; %average over all segments;
% % %       dataloc=dataep((i-1)*segshift+1:(i-1)*segshift+segleng,:);
% % %       datalocfft=fft(detrend(dataloc).*mywindow);
% % %       datalocfft=datalocfft(1:2*nf-1,:);
% % %       cslocn=((abs(datalocfft)).^3)';
% % %       for ichan=1:nchan;
% % %           xx=hankel(conj(datalocfft(1:2*nf-1,ichan)));
% % %           csloc(ichan,:,:)=(datalocfft(1:nf,ichan)*transpose(datalocfft(1:nf,ichan))).*xx(1:nf,1:nf);
% % %           cs(ichan,:,:)=cs(ichan,:,:)+csloc(ichan,:,:);
% % %       end
% % %       %xxx=squeeze(csloc(1,1:3,1:3))
% % %       
% % %       nave=nave+1;
% % %       csn=csn+cslocn;
% % %     end;
% % % end
% % % 
% % % cs=cs/nave;
% % % csn=csn/nave;
% % % csn=power(csn,1/3);
% % % for i=1:nchan;for f1=1:nf;for f2=1:nf;
% % %     csnr(i,f1,f2)=(csn(i,f1)*csn(i,f2)*csn(i,f1+f2-1));
% % % end;end;end
% % % 
% % % 
% % %   
% % %     
% % % 
% % % return;
    
