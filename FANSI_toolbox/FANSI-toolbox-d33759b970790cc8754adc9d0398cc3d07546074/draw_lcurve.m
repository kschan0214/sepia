% Draw a L-curve graph and calculate the curvature
function [ Kappa ] = draw_lcurve( Lambda, regularization, consistency, fig )
% input
% Lambda: vector with the regularization parameters
% regularization: cost associated to the regularization term
% consistency: cost associated to the data fidelity term
% fig: figure number for output
%
% output - Kappa: point of maximum curvature, standard optimization criteria
%
% Based on the code by Bilgic Berkin at http://martinos.org/~berkin/software.html
% Las modified by Carlos Milovic 2017.12.27


figure(fig), subplot(1,2,1), plot((consistency), (regularization), 'marker', '*')

% cubic spline differentiation to find Kappa (largest curvature) 

eta = log(regularization.^2);
rho = log(consistency.^2);

M = [0 3 0 0;0 0 2 0;0 0 0 1;0 0 0 0];

pp = spline(Lambda, eta);
ppd = pp;

ppd.coefs = ppd.coefs*M;
eta_del = ppval(ppd, Lambda);
ppd.coefs = ppd.coefs*M;
eta_del2 = ppval(ppd, Lambda);


pp = spline(Lambda, rho);
ppd = pp;

ppd.coefs = ppd.coefs*M;
rho_del = ppval(ppd, Lambda);
ppd.coefs = ppd.coefs*M;
rho_del2 = ppval(ppd, Lambda);


Kappa = 2 * (rho_del2 .* eta_del - eta_del2 .* rho_del) ./ (rho_del.^2 + eta_del.^2).^1.5;

index_opt = find(Kappa == max(Kappa));
disp(['Optimal lambda, consistency, regularization: ', num2str([Lambda(index_opt), consistency(index_opt), regularization(index_opt)])])

figure(fig), subplot(1,2,2), semilogx(Lambda, Kappa, 'marker', '*')
end

