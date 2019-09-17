function [ output_args ] = nlinearfitting( xdata1,ydata1 )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
a_y = []; fo_y = [];
xdata = log(xdata1);
    ydata = log(ydata1);
    start_point = [1e-13 1e-6]%rand(1, 2); %Generates 2 points to start minimization
    model = @logfun1;
    opts1 = optimset ('MaxFunEvals',1e5);
    opts2 = optimset('maxiter',1e5);
    options = optimset(opts1,opts2);
    estimates1 = fminsearch(model, start_point,options);
    %The plotting assumes an offset of 1 so that each plot in the loop will
    %not overlap the previous plots. The term 'col(mod((i+1)/2,8)))' 
    %indicates the color used that is guaranteed to be different from each
    %of the 7 plots above or below.
    figure ;
    plot(log10(exp(xdata)), log10(exp(ydata)))
    hold on
    [sse1, FittedCurve1] = model(estimates1);
    plot(log10(exp(xdata)), log10(exp(FittedCurve1)));
    %The next four lines help generate the legend of the graph (my guess is
    %there is an easier way to do this, but this way works)
%     str = ['Power ', num2str(pwr((i+1)/2),2)];
%     str2 = ['fitted curve a =',num2str(Alpha1, 4),' fo =', num2str(rolloff1, 4)];
%     leg(i,1:length(str)) = str;
%     leg((i+1),1:length(str2)) = str2;
    %Finally, the next two lines concatonates the alpha and rolloff values
    %into two vectors, a_x and fo_x, respectively.
    a_y = [a_y Alpha1];
    fo_y = [fo_y rolloff1];

% logfun accepts curve parameters as inputs, and outputs sse,
% the sum of squares error for A + .5*log(B^2+xdata.^2) - ydata, 
% and the FittedCurve. 
    function [sse1, FittedCurve1] = logfun1(params1)
        Alpha1 = params1(1);
        rolloff1 = params1(2);
        FittedCurve1 = log(Alpha1) - log(rolloff1.^2+exp(xdata.*2)); %1.0? or 0.5
        ErrorVector1 = FittedCurve1 - ydata;
        sse1 = sum(ErrorVector1 .^ 2);
    end
% legend(leg);
xlabel('Log Frequency (Hz)');
ylabel('PSD (V^2/Hz)');
title('Log Power Spectral Density (Y)');
end






