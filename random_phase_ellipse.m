%% Initial phase space for geant4 simulations
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%close all; 
clear; 
%clc;
%format long;
x0  = 0.0003443;     % m
xp0 = 0.000029023;     % rad 
alpha = 0.59018;
R1 = normrnd(0,x0,[1 1000]);
R2 = normrnd(0,xp0,[1 1000]);


for i=1:max(size(R1))
    x(i) = R1(i);
    xp(i) = R2(i) - alpha*x(i);
end


figure(20)
subplot(1,2,1)
h1 = histfit(x,50);
legend('x')
subplot(1,2,2)
h2= histfit(xp,50);
legend('xp')


[~,ind0_1] = min(abs(x-x0));
[~,ind0_2] = min(abs(xp-xp0));

% Find corresponding coordinates
xmax_x = x0;     
xmax_y = xp(ind0_1);
%
xpmax_y = xp0;
xpmax_x = x(ind0_2); 
%


% %%%%%%%% Calculate the phase space %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Method 1
%emitt0 = sqrt( (xmax_x*xpmax_y)^2 - (xmax_x*xpmax_y)*(xmax_y*xpmax_x) );

% %%%%%%%%%%%%%%%%%%%%%%
% Method 2
% rms values
xrms   = rms(x-mean(x));
xprms  = rms(xp-mean(xp));
% second terms
xxp = sum(   (x-mean(x)).*(xp-mean(xp)) ) / 1000;
emitt0 = sqrt(  xrms^2*xprms^2-xxp^2  );
% %%%%%%%%%%%%%%%%%%%%%%

% Twiss parameters
beta0 = (xmax_x^2) / emitt0;
gamma0 = (xpmax_y^2)/ emitt0;
alpha0 = -(xpmax_x)*sqrt(gamma0/emitt0);
% Display
disp(['Geometric Emittance (m-rad): ' num2str(emitt0)]);
disp(['Beta: ' num2str(beta0)]);
disp(['Alpha: ' num2str(alpha0)]);
disp(['Gamma: ' num2str(gamma0)]);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(30)
h20 = plot(x,xp,'o');
legend(['\epsilon: ' num2str(emitt0) ', \beta: ' num2str(beta0) ', \alpha: ' num2str(alpha0)]);
title('Initial emittance generated for Geant4 simulations.')
%%
saveas(h20,['initial_emitt_' num2str(x0*1000) '_' num2str(xp0*1000) '.fig']);
saveas(h20,['initial_emitt_' num2str(x0*1000) '_' num2str(xp0*1000) '.eps'],'epsc');