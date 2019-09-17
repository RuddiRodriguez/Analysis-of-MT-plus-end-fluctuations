function [muc ,mu_std1c ,mur,mu_std1r,vecc,vecr]=main


clear;

% xp=[57 100 367 53 368 193 395 133];
% yp=[283 315 291 337 279 268 217 214];
% index=[0 2 3 5 6 7 8 12];

%control48h
 xp=[328 102 103 293 55 298 353];
 yp=[254 277 267 287 207 213 236];
 index=[8 6 5 3 0 2 9];
 
    majora=[266.490 255.732 296.482 298.194 339.218 280.903 260.300]./2;
   minora=[194.443 174.343 184.204 191.310 161.417 188.536 199.512]./2;
  
   tiempo = [169 162 30 195 88 129 123 ]
   
   [muc,mu_std1c,vecc]= analysis_temp (xp,yp,index,majora,minora,tiempo,0)
%  
%gpr48h
  
  xp=[354 92 398  395  102 370 345 115 132 137 142];
  yp=[249 275 235  309  292 267 258 285 204 312 226];
  index=[13 12 10  5  3 2 0 14 4 11 6];
  
  majora=[252.066 294.307 298.66 264.842 271.902 305.143 290.985 197.055 243.929 274.790 286.266]./2;
  minora=[176.004 198.700 188.893 170.793 190.610 187.087 180.000 154.297 164.828 194.885 152.286]./2;
  
  tiempo = [185 334 342 140 80 326 162 130 55 75 62];
  
  [mur,mu_std1r,vecr]= analysis_temp (xp,yp,index,majora,minora,tiempo,1);
  
 muc;
 mur;
end

%   


%  xp=[359];
%    yp=[115];
%    index=[3];

%density
%gpr48
% index= [14 13 ]
% areaa = 
% areap
% index=[1 2]
%   xp=[420 59];
%   yp=[333 93];
    function [mu,mu_std1,vector]= analysis_temp (xp,yp,index,majora,minora,tiempo,control)
lengthd=0;
lengthdif=0;
histcell=cell(1,11);
histcellc=cell(1,11);
number=0;
gx=majora(1)/2:majora(1)/2:2*majora(1);

bb=zeros(length(gx),length(index));nb=zeros(length(gx),length(index));sb=zeros(length(gx),length(index));
for l=1:length(index)
di=0;
dif=0;
gxx(:,l)=majora(l)/2:majora(l)/2:2*majora(l);
asymParam=0;
asymFlag=0;
c=0;
distdf=0;
cx=0;
j=index(l);
dif=0;
di=0;
cxdif=0;
cxd=0;
if control==0
path=sprintf('/Users/ruddirodriguez/Documents/GPR_new_analysis/control/%ifull.xml',j);
else
path=sprintf('/Users/ruddirodriguez/Documents/GPR_new_analysis/gpr48/%ifull.xml',j);
end
[tracks, md] = importTrackMateTracks(path);
for i=1:length(tracks)
 Xb=tracks{i, 1}(1:end,2);
 Yb=tracks{i, 1}(1:end,3);
 if length (Xb)>3
 position= [Xb Yb];
 c=c+1;
 
%    distdf(i)=sqrt((mean(Xb)-xp(l)).^2+((mean(Yb)-yp(l)).^2));
  cx(i)=length(Xb);
  cy=length(Yb);
  
 [asymParam(c),asymFlag(c)] = asymDeterm2D3D(position,0.1);
 if asymFlag(c)==1
      %figure (1); plot (Xb,Yb,'r');
     di=di+1;
      lengthd(l,di)=length(Xb);
%      distdi=
   %distdf(di)=sqrt((mean(Xb)-xp(l)).^2+((mean(Yb)-yp(l)).^2));
  cxd(di)=length(Xb);
  cy=length(Yb);
 else
      %figure (1); plot (Xb,Yb,'b');
     dif=dif+1;
     lengthdif(l,dif)=length(Xb);
         distdf(dif)=sqrt((mean(Xb)-xp(l)).^2+((mean(Yb)-yp(l)).^2));
    cxdif(dif)=length(Xb);
  cy=length(Yb);
 end
 hold on
 end
end
gx=40:40:4*40;
% bin_data_my( distdf,cxdif,[0 300] );
[bb(:,l),nb(:,l),sb(:,l),yt]=bin_data_myy(distdf,cxdif,gx);%,gxx(:,l));
% edges = [4 4:0.1:20 20];
% hh=histogram(yt,edges)
% eq1=@(theta) (majora(l)*cos(theta) + minora(l)*sin(theta));
% eq2=@(theta) eq1*eq1;
% 
% theta0=atan((2.*minora(l))./(majora(l)));
% areap(l)=integral(eq2,0,theta0);
% minora(l)=80;
% majora(l)=61;
 theta =sym ('theta');
 a =sym ('a');
 b =sym ('b');
 f= a*cos(theta)+b*sin(theta);
 theta0=atan((2.*minora(l))./(majora(l)));
 inter=int(f^2,theta,0,theta0);

a=majora(l);
b=minora(l);
areap(l)=0.5.*eval(inter);

% %areap(l)= ((majora(l)*minora(l)*(majora(l)^2+(3*minora(l)^2)))/(majora(l)^2+4*minora(l)^2))+(((majora(l)^2+minora(l)^2)/2)*atan(2*minora(l)/majora(l)));
% 
areato(l)=majora(l).*minora(l).*pi;
areamidpost(l)=(areato(l)./2)-areap(l);

vector1 =[areap(l) areamidpost(l) areamidpost(l) areap(l)]';

% density(:,l)= nb(:,l)./(vector1.*0.130.*(tiempo(l).*0.1));
% gxxx(:,l)=gxx(:,l);

[nn,chist]=hist(yt,4:0.1:20);
histcell{1,l}=nn;
histcellc{1,l}=chist;
%  [ bb(:,l),h,sb(:,l),vectoranova] = bin_version2( distdf,cx,40,320,40 )
 %figure ;errorbar (gx,bb,sb);
 number(l)=di+dif;
end
lengthd(lengthd==0)=NaN;
lengthdif(lengthdif==0)=NaN;
lengthd=lengthd.*0.1;
lengthdif=lengthdif.*0.1;
  % figure(3456) ;errorbar ((gx).*0.130,nanmean(bb,2).*0.1,(nanstd(bb,[],2).*0.1)./11);
      %%
%    stdm= nanstd(density,[],2)./11;
%    meant= nanmean(density,2);
%    figure ; errorbar (meant,stdm)
   if control==0
    [mu,mu_std1,vector]= get_results (histcell,chist,7);
   else
       [mu,mu_std1,vector]= get_results (histcell,chist,11);
   end
    figure(67) ; errorbar (gx.*0.130,mu,mu_std1);
    end
    
