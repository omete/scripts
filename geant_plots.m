%% Variables:
% Step#    X(mm)    Y(mm)    Z(mm) KinE(MeV)  dE(MeV) StepLeng TrackLeng  NextVolume ProcName
% Command line import 
close all;
clear;
clc;

run_num = 8;
fileID = fopen(['run_distnew3_step0' num2str(run_num) '_primary.txt']);
formatSpec = '%f %f %f %f %f %f %f %f %s %s';
Data = textscan(fileID, formatSpec);
fclose(fileID);

Step       = Data{1};
Xmm        = Data{2}*1e-3; %m
Ymm        = Data{3}*1e-3; %m
Zmm        = Data{4}*1e-3; %m
KinE       = Data{5};
dE         = Data{6};
StepLength = Data{7}*1e-3; %m
Shape1     = Data{9};

sfig = 0; %  to save figures
%% Use after interactive data import
% Assign some parameters
part_num = length(Shape1);

Step = VarName1(1:part_num);
StepLength = VarName7(1:part_num)*1e-3; % m
Xmm  = VarName2(1:part_num)*1e-3;       % m
Ymm  = VarName3(1:part_num)*1e-3;       % m 
Zmm  = VarName4(1:part_num)*1e-3;       % m 

sfig = 1; %  to save figures
%% Real and phase space

indx = length(StepLength);
% give the longitudinal coordinate of the selected slice -- the last z point of the simulation.
%slice_ini = find(Step == 0 & (abs(Xmm) <= 0.1) & (abs(Ymm) <= 0.1));    % indices of primaries at initial position
slice_ini = find((Step == 0) & (Zmm == 0) & (abs(Xmm.^2 + Ymm.^2) <= 0.1));    % indices of primaries at initial position
%slice_f   = find(ismember(Shape1,'OutOfWorld') & (abs(Xmm) <= 0.1) & (abs(Ymm) <= 0.1));    % indices of primaries at final position
slice_f   = find(ismember(Shape1,'OutOfWorld') & (abs(Xmm.^2 + Ymm.^2) <= 0.1));    % indices of primaries at final position
%slice_f   = find(Zmm == 5e+02 & (abs(Xmm) <= 0.1) & (abs(Ymm) <= 0.1));    % indices of primaries at final position
% Calculate the final angles of the primary particles at the initial and
% the final StepLength
%num_prim = max(size(slice_ini));  %number of initial particles
num_ini = max(size(slice_ini));  %number of initial particles
num_final = max(size(slice_f));  %number of final particles 

% Calculate mean initial angles from first snum sample step
% Final angle of primaries at initial position
snum = 50; 
for i=1:num_ini
    for j=1:snum %collect samples to calculate mean angle for each particle
        sangle_ini(i,j) = (Xmm(slice_ini(i)+1+j-1)-Xmm(slice_ini(i)+j-1)) / StepLength(slice_ini(i)+1+j-1);
        sx_ini(i,j)     = Xmm(slice_ini(i)+1+j-1);
        sy_ini(i,j)     = Ymm(slice_ini(i)+1+j-1);
    end
    angle_ini(i) = mean(sangle_ini(i,:));
    x_ini(i)     = mean(sx_ini(i,:));
    y_ini(i)     = mean(sy_ini(i,:));
end
disp('x ok')

% Final angle of primaries at final position
for i=1:num_final
    for j=1:snum
        sangle_f(i,j) = (Xmm(slice_f(i)-j+1)-Xmm(slice_f(i)-1-j+1)) / StepLength(slice_f(i)-j+1); % go one step back than "envelope" to get the value at the exit of "Shape1"
        sx_f(i,j) = Xmm(slice_f(i)-j+1);
        sy_f(i,j) = Ymm(slice_f(i)-j+1);
    end
    angle_f(i) = mean(sangle_f(i,:));
    x_f(i)     = mean(sx_f(i,:));
    y_f(i)     = mean(sy_f(i,:));
end
disp('xp ok')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Check against Inf numbers --- BUG to be corrected.%%%%%%%%%%%%%%%
%angle_ini = angle_ini(find(~isinf(angle_ini)));
%x_ini = x_ini(find(~isinf(angle_ini)));
%y_ini = y_ini(find(~isinf(angle_ini)));

%angle_f = angle_f(find(~isinf(angle_f)));
%x_f = x_f(find(~isinf(angle_f)));
%y_f = y_f(find(~isinf(angle_f)));
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1);
plot3(Zmm*1e3,Xmm*1e3,Ymm*1e3,'-o')
xlabel('z (mm)')
ylabel('x (mm)')
zlabel('y (mm)')
grid on;
pbaspect([10 1 1])
%ylim([-30 30]);
%zlim([-30 30]);
%ylim([-100 100]);
if (sfig == 1)
    saveas(gca,['xyz' num2str(run_num) '.eps'],'epsc')
end

figure(2)
h2 = plot(x_f*1e3,y_f*1e3,'ob','linewidth',2)
hold on;
h1 = plot(x_ini*1e3,y_ini*1e3,'or','linewidth',2)
xlabel('x position (mm)')
ylabel('y position (mm)')
legend([h1 h2],'Initial','Final')
if (sfig == 1)
    saveas(gca,['xyif' num2str(run_num) '.eps'],'epsc')
end

figure(3)
plot(x_f*1e3,angle_f*1e3,'bo')
hold on;
plot(x_ini*1e3,angle_ini*1e3,'ro')
hold off;
xlabel('x  (mm)')
ylabel('xp  (mrad)')
legend('Final','Initial')
grid on;
%xlim([-4 4])
%ylim([-.4 .4])
if (sfig == 1)
    saveas(gca,['emitt_if ' num2str(run_num) '.eps'],'epsc')
end


%% Projections and fits
% Projections
% Fetch data
numbins1 = 30;
numbins2 = 30;
figure(5);
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
plot(h1_x,h1_y,'ob')
plot(h1_x,F(x1,h1_x),'r')
hold off;
title('Real Space')
xlabel('x_{ini} (m)')
%xlim([-6 6]*1e-3);
%
subplot(2,2,2)
hist(x_f,numbins1)
hold on;
plot(h2_x,h2_y,'ob')
plot(h2_x,F(x2,h2_x),'r')
hold off;
xlabel('x_f (m)')
%xlim([-6 6]*1e-3);
%
subplot(2,2,3)
hist(angle_ini,numbins2)
hold on;
plot(h3_x,h3_y,'ob')
plot(h3_x,F(x3,h3_x),'r')
hold off;
title('Phase Space')
xlabel('xp_{ini} (rad)')
%xlim([-2 2]*1e-3);
%
subplot(2,2,4)
hist(angle_f,numbins2)
hold on;
plot(h4_x,h4_y,'ob')
plot(h4_x,F(x4,h4_x),'r')
hold off;
xlabel('xp_f (rad)')
%xlim([-2 2]*1e-3);
if (sfig == 1)
    saveas(gca,['histos' num2str(run_num) '.eps'],'epsc')
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

xmax_x = sig*[xsize1 xsize2];     %x1, measured
xmax_y(1) = sig*angle_ini(ind_1); %x2 initial, retrieved 
xmax_y(2) = sig*angle_f(ind_2);   %x2 final, retrieved

xpmax_y = sig*[xsize3 xsize4];    %xp2, measured
xpmax_x(1) = sig*x_ini(ind_3);    %xp1 initial, retrieved   
xpmax_x(2) = sig*x_f(ind_4);      %xp1 final, retrieved
                                         

%emittance_i = (xmax_x(1)*xpmax_y(1))*(xmax_y(1)*xpmax_x(1))+(xmax_x(1)*xpmax_y(1))^2;
%emittance_f = (xmax_x(2)*xpmax_y(2))*(xmax_y(2)*xpmax_x(2))+(xmax_x(2)*xpmax_y(2))^2;

%xmax_x = abs(xmax_x);
%xmax_y = abs(xmax_y);
%xpmax_x = abs(xpmax_x);
%xpmax_y = abs(xpmax_y);

emittance_i = sqrt( (xmax_x(1)*xpmax_y(1))^2 + (xmax_x(1)*xpmax_y(1))*(xmax_y(1)*xpmax_x(1)) );
%emitt_in = (10000/0.511)*emittance_i;
emittance_f = sqrt( (xmax_x(2)*xpmax_y(2))^2 + (xmax_x(2)*xpmax_y(2))*(xmax_y(2)*xpmax_x(2))  );
%emitt_fn = (10000/0.511)*emittance_f;

% Twiss parameters
beta_ini = (xmax_x(1)^2) / emittance_i;
beta_f   = (xmax_x(2)^2) / emittance_f;

gamma_ini = (xpmax_y(1)^2)/ emittance_i;
gamma_f   = (xpmax_y(2)^2)/ emittance_f;

alpha_ini = xpmax_x(1)*sqrt(gamma_ini/emittance_i);
alpha_f   = xpmax_x(2)*sqrt(gamma_f/emittance_f);

demitt = emittance_f-emittance_i;

disp(['Geometric Emittance (m-rad): ' num2str(emittance_i) ', ' num2str(emittance_f)])
%disp(['Normalised Emittance (m-rad): ' num2str(emitt_in) ', ' num2str(emitt_fn)])
disp(['Emittance growth (m-rad): ' num2str(demitt)]);
disp(['Beta: ' num2str(beta_ini) ', ' num2str(beta_f)]);
disp(['Alpha: ' num2str(alpha_ini) ', ' num2str(alpha_f)])
disp(['Gamma: ' num2str(gamma_ini) ', ' num2str(gamma_f)])
disp(['Initialise next run with x_ini = ' num2str(xsize2*1000) 'mm, xp_ini = ' num2str(xsize4*1000) ' mrad.']);


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
%xlim([-20 20]*1e-5)
%ylim([-20 20]*1e-5)
if (sfig == 1)
    saveas(gca,['norm_phasespace' num2str(run_num) '.eps'],'epsc')
end



%% WHEN SORTED GAUSSIAN ENERGY DISTRIBUTION -- TO BE DONE IN SUMMER 2014
% Energy distribution

figure(8)
subplot(1,2,1)
he1 = histfit(KinE(slice_ini),10);
he1x = get(he1(2),'XData');
he1y = get(he1(2),'YData');
xlabel('Initial')
subplot(1,2,2)
he2 = histfit(KinE(slice_f),10);
he2x = get(he1(2),'XData');
he2y = get(he1(2),'YData');
xlabel('Final')

% Model: Gaussian distribution 
F = @(x,xdata)x(1)*exp(-((xdata-x(2))/x(3)).^2);
% Initial values
% Perform the fittings  
x0 =[max(he1y); 0; (he1x(100)-he1x(1))/4];
[xe1,resnorm1,~,exitflag1,output1] = lsqcurvefit(F,x0,he1x,he1y);
x0 =[max(he2y); 0; (he2x(100)-he2x(1))/4];
[xe2,resnorm2,~,exitflag2,output2] = lsqcurvefit(F,x0,he2x,he2y);


figure(8)
subplot(1,2,1)
h1 = histfit(KinE(slice_ini),20);
hold on;
plot(he1x,he1y,'ob')
plot(he1x,F(xe1,he1x),'r')
hold off;
xlabel('Initial')
subplot(1,2,2)
h2 = histfit(KinE(slice_f),20);

hold on;
plot(he2x,he2y,'ob')
plot(he2x,F(xe2,he2x),'r')
hold off;
xlabel('Final')
