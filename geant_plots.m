% Variables:
% Step#    X(mm)    Y(mm)    Z(mm) KinE(MeV)  dE(MeV) StepLeng TrackLeng  NextVolume ProcName

fileID = fopen('run5_Li_primary.txt');
formatSpec = '%f %f %f %f %f %f %f %f %s %s';
Data = textscan(fileID, formatSpec);

%VarName1 = Data{1};
%VarName2 = Data{2};
%VarName3 = Data{3};
%VarName4 = Data{4};
%VarName7=Data{7};
%Shape1 = Data{9};

% Assign some parameters
part_num = length(Data{1});

%Step = VarName1(1:part_num);
%StepLength = VarName7(1:part_num)*1e-3; % m
%Xmm  = VarName2(1:part_num)*1e-3;       % m
%Ymm  = VarName3(1:part_num)*1e-3;       % m 
%Zmm  = VarName4(1:part_num)*1e-3;       % m 

Step       = Data{1};
Xmm        = Data{2}*1e-3; %m
Ymm        = Data{3}*1e-3; %m
Zmm        = Data{4}*1e-3; %m
StepLength = Data{7}*1e-3; %m
Shape1     = Data{9};

sfig = 1; %  to save figures
%% Real and phase space

indx = max(size(StepLength));
% give the longitudinal coordinate of the selected slice -- the last z point of the simulation.
slice_ini = find(Step == 0 & (abs(Xmm) <= 0.1) & (abs(Ymm) <= 0.1));    % indices of primaries at initial position
slice_f   = find(ismember(Shape1,'OutOfWorld') & (abs(Xmm) <= 0.1) & (abs(Ymm) <= 0.1));    % indices of primaries at final position
%slice_f   = find(Zmm == 5e+02 & (abs(Xmm) <= 0.1) & (abs(Ymm) <= 0.1));    % indices of primaries at final position
% Calculate the final angles of the primary particles at the initial and
% the final StepLength
%num_prim = max(size(slice_ini));  %number of initial particles
num_ini = max(size(slice_ini));  %number of initial particles
num_final = max(size(slice_f));  %number of final particles 

% Final angle of primaries at initial position
for i=1:num_ini
    angle_ini(i) = (Xmm(slice_ini(i)+1)-Xmm(slice_ini(i))) / StepLength(slice_ini(i)+1);
    x_ini(i) = Xmm(slice_ini(i)+1);
    y_ini(i) = Ymm(slice_ini(i)+1);
end
disp('ok')
% Final angle of primaries at final position
for i=1:num_final
    angle_f(i) = (Xmm(slice_f(i))-Xmm(slice_f(i)-1)) / StepLength(slice_f(i)); % go one step back than "envelope" to get the value at the exit of "Shape1"
    x_f(i) = Xmm(slice_f(i));
    y_f(i) = Ymm(slice_f(i));
end


figure(1);
plot3(Zmm*1e3,Xmm*1e3,Ymm*1e3,'-o')
xlabel('z (mm)')
ylabel('x (mm)')
zlabel('y (mm)')
grid on;
pbaspect([10 1 1])
if (sfig == 1)
    saveas(gca,'xyz.eps','epsc')
end

figure(2)
h2 = plot(x_f*1e3,y_f*1e3,'ob','linewidth',2)
hold on;
h1 = plot(x_ini*1e3,y_ini*1e3,'or','linewidth',2)
xlabel('x position (mm)')
ylabel('y position (mm)')
legend([h1 h2],'Initial','Final')
xlim([-100 100])
ylim([-100 100])
if (sfig == 1)
    saveas(gca,'xyif.eps','epsc')
end

figure(3)
plot(x_ini*1e3,angle_ini,'ro')
hold on;
plot(x_f*1e3,angle_f,'bo')
plot(x_ini*1e3,angle_ini,'ro')
hold off;
xlabel('x  (mm)')
ylabel('xp  (rad)')
legend('Initial','Final','Initial')
grid on;
xlim([-110 110])
%ylim([-1 1])
if (sfig == 1)
    saveas(gca,'emitt_if.eps','epsc')
end


%% Projections and fits
% Projections
% Fetch data
numbins1 = 40;
numbins2 = 40;
figure(5)
h1 = histfit(x_ini,numbins1);
h1_x = get(h1(2),'XData');
h1_y = get(h1(2),'YData');
h2 = histfit(x_f,numbins1);
h2_x = get(h2(2),'XData');
h2_y = get(h2(2),'YData');
h3 = histfit(angle_ini,numbins2);
h3_x = get(h3(2),'XData');
h3_y = get(h3(2),'YData');
h4 = histfit(angle_f,numbins2);
h4_x = get(h4(2),'XData');
h4_y = get(h4(2),'YData');

% Model: Gaussian distribution 
F = @(x,xdata)x(1)*exp(-((xdata-x(2))/x(3)).^2);
% Initial values
% Perform the fittings  
x0 =[max(h1_y); 0; (h1_x(100)-h1_x(1))/4];
[x1,resnorm1,~,exitflag1,output1] = lsqcurvefit(F,x0,h1_x,h1_y);
x0 =[max(h2_y); 0; (h2_x(100)-h2_x(1))/4];
[x2,resnorm2,~,exitflag2,output2] = lsqcurvefit(F,x0,h2_x,h2_y);
x0 =[max(h3_y); 0; (h3_x(100)-h3_x(1))/4];
[x3,resnorm3,~,exitflag3,output3] = lsqcurvefit(F,x0,h3_x,h3_y);
x0 =[max(h4_y); 0; (h4_x(100)-h4_x(1))/4];
[x4,resnorm4,~,exitflag4,output4] = lsqcurvefit(F,x0,h4_x,h4_y);

% Retrieve fit parameters
xsize1 = x1(3)/sqrt(2);   % 1sigma = 68.27% x_ini
xsize2 = x2(3)/sqrt(2);   % x_f
xsize3 = x3(3)/sqrt(2);   % angle_ini
xsize4 = x4(3)/sqrt(2);   % angle_f
posx1 = x1(2);
posx2 = x2(2);
posx3 = x3(2);
posx4 = x4(2);
%gofx = gof1.rsquare;


figure(6)
subplot(2,2,1)
hist(x_ini,numbins1)
hold on;
plot(h1_x,h1_y,'b')
plot(h1_x,F(x1,h1_x),'r')
hold off;
title('Real Space')
xlabel('x_{ini} (m)')
%
subplot(2,2,2)
hist(x_f,numbins1)
hold on;
plot(h2_x,h2_y,'b')
plot(h2_x,F(x2,h2_x),'r')
hold off;
xlabel('x_f (m)')
%
subplot(2,2,3)
hist(angle_ini,numbins2)
hold on;
plot(h3_x,h3_y,'b')
plot(h3_x,F(x3,h3_x),'r')
hold off;
title('Phase Space')
xlabel('xp_{ini} (rad)')
%
subplot(2,2,4)
hist(angle_f,numbins2)
hold on;
plot(h4_x,h4_y,'b')
plot(h4_x,F(x4,h4_x),'r')
hold off;
xlabel('xp_f (rad)')
if (sfig == 1)
    saveas(gca,'histos.eps','epsc')
end

%% Calculate the emittance and the Twiss parameters
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use measured, x_ini, x_f, angle_ini, angle_f
% Retrieve the coordinates in the phase space

[~,ind_1] = min(abs(x_ini-xsize1));
[~,ind_2] = min(abs(x_f-xsize2));
[~,ind_3] = min(abs(angle_ini-xsize3));
[~,ind_4] = min(abs(angle_f-xsize4));

% Confidence 
sig = 1;

xmax_x = sig*[xsize1 xsize2]; %x1, measured
xmax_y(1) = sig*angle_ini(ind_1); %x2 initial, retrieved 
xmax_y(2) = sig*angle_f(ind_2);   %x2 final, retrieved

xpmax_y = sig*[xsize3 xsize4];%xp2, measured
xpmax_x(1) = sig*x_ini(ind_3);    %xp1 initial, retrieved   
xpmax_x(2) = sig*x_f(ind_4);      %xp1 final, retrieved
                                          %xp2, measured

%emittance_i = (xmax_x(1)*xpmax_y(1))*(xmax_y(1)*xpmax_x(1))+(xmax_x(1)*xpmax_y(1))^2;
%emittance_f = (xmax_x(2)*xpmax_y(2))*(xmax_y(2)*xpmax_x(2))+(xmax_x(2)*xpmax_y(2))^2;

xmax_x = abs(xmax_x);
xmax_y = abs(xmax_y);
xpmax_x = abs(xpmax_x);
xpmax_y = abs(xpmax_y);

emittance_i = sqrt( (xmax_x(1)*xpmax_y(1))^2 + (xmax_x(1)*xpmax_y(1))*(xmax_y(1)*xpmax_x(1)) );
emitt_in = (10000/0.511)*emittance_i;
emittance_f = sqrt( (xmax_x(2)*xpmax_y(2))^2 + (xmax_x(2)*xpmax_y(2))*(xmax_y(2)*xpmax_x(2))  );
emitt_fn = (10000/0.511)*emittance_f;

% Twiss parameters

alpha_ini = (xmax_x(1)*xmax_y(1))/emittance_i;
alpha_f   = (xmax_x(2)*xmax_y(2))/emittance_f;

beta_ini = (xmax_x(1)^2) / emittance_i;
beta_f   = (xmax_x(2)^2) / emittance_f;

gamma_ini = (xpmax_y(1)^2)/ emittance_i;
gamma_f   = (xpmax_y(2)^2)/ emittance_f;

demitt = emitt_fn-emitt_in;


disp(['Geometric Emittance (m-rad): ' num2str(emittance_i) ', ' num2str(emittance_f)])
disp(['Normalised Emittance (m-rad): ' num2str(emitt_in) ', ' num2str(emitt_fn)])
disp(['Emittance growth (m-rad): ' num2str(demitt)])
disp(['Beta: ' num2str(beta_ini) ', ' num2str(beta_f)])
disp(['Alpha: ' num2str(alpha_ini) ', ' num2str(alpha_f)])
disp(['Gamma: ' num2str(gamma_ini) ', ' num2str(gamma_f)])



%% Normalise the phase space
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_ini_norm = x_ini/sqrt(beta_ini);
x_f_norm = x_f/sqrt(beta_f);
%x_ini_norm = x_ini;
%x_f_norm = x_f;

angle_ini_norm = (alpha_ini/sqrt(beta_ini))*x_ini + (sqrt(beta_ini)*angle_ini);
angle_f_norm = (alpha_f/sqrt(beta_f))*x_f + (sqrt(beta_f)*angle_f);

figure(7)
plot(x_f_norm,angle_f_norm,'bo')
hold on;
plot(x_ini_norm,angle_ini_norm,'ro')
hold off;
%xlabel('x  (mm)')
%ylabel('xp  (mrad)')
legend('Final','Initial')
grid on;
xlim([-20 20]*1e-3)
ylim([-20 20]*1e-3)
if (sfig == 1)
    saveas(gca,'norm,_phasespace.eps','epsc')
end














