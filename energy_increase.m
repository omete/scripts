% Energy increase in the plasma
close all;
clear;
clc;

E0 = 10; %GeV
Ee = 0.000511; %GeV
g = 0.495; % GeV/m  

% Energy high resolution
for i=1:500
    x(i) = i;
    E(i) = g*x(i) + E0;
end

%Energy low resolution
%dx_ = [10 10 10 50 50 50 50 50 50 50 50 70];
dx_ = [10 10 10 50 50 50 50 50 50 50 50 50 20];
x_(1) = 0;
for i=1:length(dx_)
    x_(i+1) = x_(i)+ dx_(i); 
end
E_ = g*x_ + E0;

% Values for Set2
%emitt = [0.1073 0.1205    0.2395    0.4217    1.4388    5.8866    9.8954   15.3130   26.4340]; %um
%emitt_growth = [0 0.1245 0.1435 0.3549 3.8121 2.8409 6.3908 7.7950 7.9480 17.0080];
%%emitt_growth = emitt_growth./emitt(1:length(emitt_growth));   % as (%)
%err = [0 0.0338 0.0261 0.0667 1.8471 1.0370 3.7198 1.8668 0.2837];
% Values for Set3
%emitt        = [0.1 1.6325    4.2366    8.3892   42.5210  112.4600  167.8500  267.5100  399.8100];
% Values for Set4
emitt = [1.4681e-12 9.6878e-10 3.2906e-9 6.5801e-9 8.5537e-8 1.5521e-7 3.6629e-7 6.1420e-7 1.1027e-6 1.5815e-6 2.4494e-6 3.2541e-6 5.1314e-6 6.3244e-6];
err = [7.2943e-13 6.8647e-10 1.5287e-9 1.6830e-9 7.9035e-9 8.1337e-8 1.1518e-7 1.0877e-7 1.1105e-7 4.0877e-7 4.5346e-7 5.9929e-7 8.9951e-7 7.0960e-7];
beta = [10.652 15.069 16.7716 23.2933 84.8684 71.5894 74.2311 67.3084 78.445 78.7795 93.558 84.5688 74.8135];

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%% Implementation of the theoretical curve %%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r0 = 2.8179403267*10^-15; % Classical electron radius, m
Q  = 1;         % np/n
Rb = 2.5e-5;    % Plasma radius (ion column radius), m
Ra = 1e-10;     % Atomic radius of Lithium
Z  = 3;         % Atomic number of Lithium

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PAC07 Dr. Kirby %%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1
% S term
S = Q*( log(Rb/Ra) + (1.78*Z*(Z+1)*log(287/sqrt(Z)))/(Q^2)  );
% Energy term
%e_term(1) = sqrt(E(1)/Ee);
e_term(1) = 0;
for i=1:length(E)-1
    e_term(i+1) = sqrt(E(i+1)/Ee) - sqrt(E(i)/Ee);
end
demitt_theory = sqrt(2)*r0*S*e_term*(180)/(E(i)/Ee); %un-normalise
% un-normalise
%demitt_theory = demitt_theory.*(pi./(E./Ee));

% Emittance curve
%emitt_theory(1) = 0;
emitt_theory(1) = emitt(1);
for i=1:length(E)-1
    emitt_theory(i+1) = sqrt(emitt_theory(i)^2 + demitt_theory(i)^2);  %geometric
end

end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% O. Mete  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Constants
ngas = 6*10^20;   % m^-3
c = 3*10^8; % m/s
theta_max = 2*pi;   % degrees
r0 = 2.8179403267*10^-15; % Classical electron radius, m
hbar = 6.58211928*10^-16; %eV.s

% Average along the plasma channel
for i=1:length(x_)-1
    Av(i) = c*ngas*beta(i);
end

% Set theta_min
%for i=1:length(x_)-1
%theta_min(i) = hbar/((E_(i+1)/Ee)*me*c*152e-12);
%end
for i=1:length(x_)-1
theta_min(i) = Z^(1/3)*(1/(1.4*hbar*c))*(1/(E_(i+1)/Ee));
end

% Conversion to degrees
%theta_min = (180/pi)*theta_min;

% Angular part
for i=1:length(x_)-1
%Ang(i) = (pi/(2*theta_max))*(3*theta_min(i)*atan(theta_max/theta_min(i)) + theta_max*(log(theta_min(i)^2+theta_max^2)-2));
Ang(i) = log((theta_min(i)^2 + theta_max^2) / theta_min(i)^2) - (theta_max^2/(theta_min(i)^2 + theta_max^2)); 
end

  
% Emittance
for i=1:length(x_)-1
demitt_theory2(i) = (2*pi*(Z*r0)^2*Av(i)*Ang(i))/(E_(i+1)/Ee);
demitt_theory2(i) = demitt_theory2(i)*pi/(E_(i+1)/Ee); % *pi/(E_(i+1)/Ee) for un-normalise (pi/gamma);
end


emitt_theory2(1)=emitt(1);
for i=1:length(E_)-1
    emitt_theory2(i+1) = sqrt(emitt_theory2(i)^2 + demitt_theory2(i)^2);  %geometric
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Emittance with a factor
for i=1:length(x_)-1
%demitt_theory2(i) = ((2*Z*r0)^2/(g*x_(i+1)*1e9+gamma0))*(tau/2)*Av(i)*Ang(i); 
demitt_theory2n(i) = (2*pi*(Z*r0)^2*Av(i)*Ang(i))/(E_(i+1)/Ee);
demitt_theory2n(i) = demitt_theory2n(i)*600*pi/(E_(i+1)/Ee); % *pi/(E_(i+1)/Ee) for un-normalise (pi/gamma);
end

emitt_theory2n(1)=emitt(1);
for i=1:length(E_)-1
    emitt_theory2n(i+1) = sqrt(emitt_theory2n(i)^2 + demitt_theory2n(i)^2);  %geometric
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%grey = [0.4 0.3 0.6];
figure(1)
h1 = errorbar(x_,emitt*1e6, err*1e6,'ob');
hold on;
%errorbar(x_,emitt_growth*1e6,err*1e6,'mo');
%h2 = plot(x,emitt_theory*1e6,'--k','linewidth',2);
h3 = plot(x_,emitt_theory2*1e6,'--og','linewidth',1);
h4 = plot(x_,emitt_theory2n*1e6,'--or','linewidth',1);
errorbar(x_,emitt*1e6, err*1e6,'ob');
plot(500, (0.02*1e-9)*1e6,'ok','linewidth',2);
%set(gca,'yscale','log','YColor',grey)
set(gca,'yscale','log')
hold off;
grid on;
xlabel(gca,'s (m)','fontsize',14);
ylabel(gca,'\epsilon_{geometric} (\mu m)','fontsize',14);
xlim([0 max(x_)]);
%legend([h1 h2 h3], 'GEANT4 results','PAC07, N. Kirby et al.', 'Section 2, Eq. 7');
legend([h1 h3 h4], 'GEANT4 results', 'Theory', 'Theory with normalised \Delta\epsilon');
%ylim([1e-25 1e5]);





