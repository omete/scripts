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
%emitt = [0.1025 0.2268 0.3702 0.7252 4.5373 7.3782 13.7690 21.5640 29.5120 46.5200];
%emitt_growth = [0 0.1245 0.1435 0.3549 3.8121 2.8409 6.3908 7.7950 7.9480 17.0080];
%%emitt_growth = emitt_growth./emitt(1:length(emitt_growth));   % as (%)
%err = [0 0.0338 0.0261 0.0667 1.8471 1.0370 3.7198 1.8668 0.2837];

% Values for Set3
emitt        = [0.1 3.9535 8.5656 14.0460 66.6870 150.6700 238.16];
emitt_growth = [0 3.8535 4.6121 5.4804 52.6410 83.9830 87.49 ];
err          = [0 0.7245 5.0274 3.1753 3.4663 49.9001];

figure(1)
[haxes, hline1, hline2] = plotyy(x_(1:length(emitt)),emitt,x,E,'semilogy','plot');
hold(haxes(2),'on');
h1 = plot(haxes(2),x_,E_,'go','linewidth',2);
hold(haxes(1),'on');
h2 = errorbar(haxes(1),x_(1:length(err)),emitt(1:length(err)),err,'o');
% Growth %%%%%%
h3 = errorbar(haxes(1),x_(1:length(err)),emitt_growth(1:length(err)),err,'mo');
plot(haxes(1),x_(1:length(emitt)),emitt_growth(1:length(emitt)),'mo');
legend([h2 h3], '\epsilon','\Delta \epsilon');
% %%%%%%%%%%%%%
xlabel(haxes(2),'s (m)','fontsize',14);
ylabel(haxes(1),'\epsilon (geometric) (nm)','fontsize',14);
ylabel(haxes(2),'E (GeV) ','fontsize',14);
set(hline1,'linestyle','o','linewidth',2);
set(hline2,'linestyle',':','linewidth',2);
grid on;
xlim([0 max(x_)]);
ylim([0 10000]);