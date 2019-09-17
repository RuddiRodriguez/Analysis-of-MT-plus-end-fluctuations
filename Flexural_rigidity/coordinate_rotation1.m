function [ xrr,yrr,gof,fitresult ] = coordinate_rotation1( x,y, ldata,inid )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

  x=abs(x-x(1));
  y=abs(y-y(1));
[fitresult, gof] = createFit_linear_tracks(x,y,ldata,inid);
angle1=atand(fitresult.p1)
 %angle1 = 89;
xr=x*cosd(angle1)+y*sind(angle1);
yr=-x*sind(angle1)+y*cosd(angle1);
% angle=1:1:360;



%Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 'y vs. x', 'untitled fit 1', 'Location', 'NorthEast' );
% % Label axes
% xlabel( 'x' );
% ylabel( 'y' );
% grid on


% % for i=1:length(angle)
% % xr=x*cosd(angle(i))-y*sind(angle(i));
% % yr=x*sind(angle(i))+y*cos(angle(i));
% %
% %
% % [fitresult, gof] = createFit_linear_tracks(xr, yr);
% %
% % if abs(fitresult.p1)<0.01
% %     figure ; plot (xr-xr(1),yr-yr(1),x-x(1),y-y(1))
    xrr=xr;
    yrr=yr;%-nanmean(yr);
    
    figure ; plot (x,y,xrr,yrr);
% %     timevector=xrr(1):0.05:xrr(end);
% %     [yGrid] = interp1 ( xrr, yrr ,timevector);
% %     break;
% % else
% % end
% dis=abs((y-((fitresult.p1).*x+fitresult.p2))./(sqrt(((fitresult.p1).^2)+1)));
% % %  figure(1900) ; plot (angle(i),fitresult.p1,'o');hold on
% % end

end