function [tr, tr_std1, A, pdf] = fit_spectrum_nonli(ti, ni)


opts = optimset('Jacobian', 'off', ...
    'MaxFunEvals', 1e4, ...
    'MaxIter', 1e4, ...
    'Display', 'off', ...
    'TolX', 1e-8, ...
    'Tolfun', 1e-8);

mu0 = nansum(ni.*ti)/nansum(ni);

[p,resnorm,R,~,~,~,J] = lsqnonlin(@costFct, [1e8 1], [], [], opts, ti, ni);
A = p(1);
tr = p(2);

J = full(J);
mu_std = sqrt( resnorm/(numel(ni)-2) * inv(J'*J) ); %#ok<MINV>

% [ypred,delta] = nlpredci(@costFct,ti,p,R,'Jacobian',J)
ci = nlparci(p,R,'jacobian',J)
tr_std1=abs(tr-(ci(1,2)));
pdf = 4.*A.*tr./(1+((6.28.*ti.*tr).^2));

figure( 'Name', 'untitled fit 1' );
h = plot( ti, ni,'o', ti, pdf );hold on
% bar (ti,ni)

legend( h, 'n vs. xout', 'untitled fit 1', 'Location', 'NorthEast' );
% Label axes
xlabel( 'xout' );
ylabel( 'n' );
grid on

% DOF = numel(x)-numel(beta0); % Degrees of Freedom 
% % This function is called after lsqcurvefit 
% mse = resnorm/(DOF); % Mean square error
% COVBw = (Jw'*Jw)./mse; 
% 

pdf(isnan(ni)) = NaN;


function v = costFct(p, ti, ni)
tr = p(1);
A = p(2);

pdf = 4.*A.*tr./(1+((6.28.*ti.*tr).^2));
v = pdf-ni;
v(isnan(v)) = [];