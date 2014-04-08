% Energy increase in the plasma

E0 = 10; %GeV
Ee = 0.000511; %GeV
g = 0.495; % GeV/m  

% Energy high resolution
for i=1:500
    x(i) = i;
    E(i) = g*x(i) + E0;
end

%Energy low resolution
dx_ = [10 10 10 50 50 50 50 50 50 50 50 70];
x_(1) = 0;
for i=1:length(dx_)
    x_(i+1) = x_(i)+ dx_(i); 
end
E_ = g*x_ + E0;

% Values for Set2
emitt = [0.1073 0.1205    0.2395    0.4217    1.4388    5.8866    9.8954   15.3130   26.4340]; %um
%emitt_growth = [0 0.1245 0.1435 0.3549 3.8121 2.8409 6.3908 7.7950 7.9480 17.0080];
%%emitt_growth = emitt_growth./emitt(1:length(emitt_growth));   % as (%)
%err = [0 0.0338 0.0261 0.0667 1.8471 1.0370 3.7198 1.8668 0.2837];

% Values for Set3
%emitt        = [0.1 1.6325    4.2366    8.3892   42.5210  112.4600  167.8500  267.5100  399.8100];
%err = [];

% Emittance growth
emitt_growth(1) = 0;
for i=1:length(emitt)-1
    emitt_growth(i+1) = emitt(i+1)-emitt(i);
end


figure(1)
[haxes, hline1, hline2] = plotyy(x_(1:length(emitt)),emitt,x,E,'semilogy','plot');
hold(haxes(2),'on');
h1 = plot(haxes(2),x_,E_,'go','linewidth',1);
hold(haxes(1),'on');
h2 = plot(haxes(1),x_(1),emitt(1),'ob');
h3 = plot(haxes(1),x_(1:length(emitt)),emitt_growth(1:length(emitt)),'mo');
legend([h2 h3], '\epsilon','\Delta \epsilon');
% %%%%%%%%%%%%%
xlabel(haxes(2),'s (m)','fontsize',14);
ylabel(haxes(1),'\epsilon (geometric) (\mu m)','fontsize',14);
ylabel(haxes(2),'E (GeV) ','fontsize',14);
set(hline1,'linestyle','o','linewidth',2);
set(hline2,'linestyle','-','linewidth',1);
grid on;
xlim([0 max(x_)]);
ylim([0 1000]);