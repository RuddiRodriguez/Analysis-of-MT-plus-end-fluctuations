
close all;
for j=2:2
pathopen=fullfile('Z:\users\Rodri029\All plus end accumulation events\Plus end\Un-FRAP\data 7\analysis\roi1\tracks',...
    sprintf('tracksc_%i.xml',j));
% path=sprintf('%s\tracksc_%i.xml',pathopen,j);
[tracks, md] = importTrackMateTracks(pathopen);
Xb=tracks{1, 1}(1:end-5,2);
Yb=tracks{1, 1}(1:end-5,3);
[ xrr,yrr,gof,fitresult ] = coordinate_rotation1( Xb,Yb,length(Xb) );
yrrt(1:length(yrr),j )=yrr;
[ psdx,freq,relaxationtime] = spectr_calcu( yrr );
end

