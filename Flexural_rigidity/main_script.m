
close all;
dbt = rdir(['Z:\users\Rodri029\All plus end accumulation events\EB3 20 nM + 100 nM Fchitax3 2spf_10min_200ms_50laser_2*', '\**\tracksca*.xml']);
for j=1:length(dbt)
pathopen=dbt(j);%fullfile('Z:\users\Rodri029\All plus end accumulation events\Plus end\Un-FRAP\data 4\analysis\roi1\tracks',...
    %sprintf('tracksca_%i.xml',j));
% path=sprintf('%s\tracksc_%i.xml',pathopen,j);
[~,namef,~] = fileparts(pathopen.name)
if strncmp (namef,'tracksc',7) &&  strncmp (namef,'tracksca_1',10) 
[tracks, md] = importTrackMateTracks(pathopen.name);
Xb=(tracks{1, 1}(60:end,2))./16;
Yb=(tracks{1, 1}(60:end,3))./16;
[ xrr,yrr,gof,fitresult ] = coordinate_rotation1( Xb,Yb,length(Xb),1 );
vectorr= [xrr yrr];
totalxy= sqrt((xrr.^2)+(yrr.^2));
yrrn= yrr-mean(yrr);
varyrr = var (yrrn);
 yrrt(1:length(yrrn),j )=yrrn;
%  [ psdx,freq,relaxationtime] = spectr_calcu( yrr );
end
end

