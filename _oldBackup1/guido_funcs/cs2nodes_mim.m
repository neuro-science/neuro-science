
function [res1,res2,res3,res4]=cs2nodes_mim(cs,F1,F2);



[nchan ng ndipdir]=size(F1);
ng2=ng;
if nargin>2
    [nchan2 ng2 ndipdir2]=size(F2)
end

if nargin==2;
    res1=zeros(ng,1);
    res2=zeros(ng,1);
    res3=zeros(ng,1);
    res4=zeros(ng,1);
else
    res1=zeros(ng,ng2);
    res2=zeros(ng,ng2);
end


csvoxrealinv=zeros(ndipdir,ndipdir,ng);
csvoxrealinvsqrt=zeros(ndipdir,ndipdir,ng);
csvoxvec=zeros(ndipdir,nchan,ng);
c1=zeros(ng,nchan);
for i=1:ng
    Floc=squeeze(F1(:,i,:));
    csvoxrealinv(:,:,i)=inv(Floc'*real(cs)*Floc);
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
        csvoxrealinv2(:,:,i)=inv(Floc2'*real(cs)*Floc2);
        csvoxrealinvsqrt2(:,:,i)=sqrtm(csvoxrealinv2(:,:,i));
        csvoxvec2(:,:,i)=Floc2'*cs;
    end
end


for i=1:ng
    if round(i/1000)*1000==i;disp(i);end
    csvoxvecloc=csvoxvec(:,:,i);
    caainv=csvoxrealinv(:,:,i);
    caainvsqrt=sqrtm(csvoxrealinv(:,:,i));
    
    for j=1:ng2;
        if nargin==2;
            Floc2=squeeze(F1(:,j,:));
        else
            %ng2=ng2
            Floc2=squeeze(F2(:,j,:));
        end
        cab=imag(csvoxvecloc*Floc2);
        cbbinv=csvoxrealinv2(:,:,j);
        %cbbinvsqrt=sqrtm(csvoxrealinv(:,:,j));
        X=cab*cbbinv*cab';
        resloc1(j)=trace(caainv*X);
        Y=caainvsqrt*X*caainvsqrt;
        [u,s,v]=svd(Y);
        resloc2(j)=s(1,1);
        if nargin>2
            res1(i,j)=resloc1(j);
            res2(i,j)=resloc2(j);
        end
    end
     if nargin==2;
       res1(i)=mean(resloc1);
       res2(i)=mean(resloc2);
       res3(i)=max(resloc1);
       res4(i)=max(resloc2);
     end
      
 end


return;
	  
	  
	  
%其实呢，国际贸易本身没有错，错在以压榨本国国民（直接压榨企业员工和通过人为汇率压榨全国人民）换取市场的办法。有了GDP和税收，毁了整个国家的产业结构。早先我老家当地很多企业生产的机械装备是大量出口到国外甚至包括欧洲一些机械工业不错的国家的。但是自从血汗工厂兴起后，以对内压榨员工，对外以低价格，对采购人员行贿为利器的私营企业迅速挤压这些国营大厂在国内的生存空间（毕竟主要市场还在国内），然后朱镕基主导下国企大量破产：厂领导贷款组建私人企业，廉价买下国有资产，然后倒手卖给南方私人老板，转产无技术含量的低端产品，宝贵的技术工人下岗去扫马路或者当小贩。这个模式在北方很多工业城市很常见。如今的中国制造已经成为劣质品的代名词。不是我们不能造好东西，而是造好东西的企业被灭了。以前的国企工人收入尚可，可以支持大量的市场。现在的私企工人收入过低基本无法支持市场。	  
% 	好吧，很多人在说和我原意无关的东西。



						  