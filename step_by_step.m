%%
%num1=slice_ini;
num1 = slice_f;
figure(1);
plot3(Zmm(num1)*1e3,Xmm(num1)*1e3,Ymm(num1)*1e3,'-ob')
hold on;
%plot3(Zmm(num2)*1e3,Xmm(num2)*1e3,Ymm(num2)*1e3,'-mo')
grid on;
xlabel('z (mm)')
ylabel('x (mm)')
zlabel('y (mm)')
grid on;
pbaspect([10 1 1])

%%
figure(10)
plot(Xmm(slice_f),Ymm(slice_f),'.b')


%%
xmax_xt = 0.001;
xmax_yt = 0.0;
xpmax_xt = 0.0;
xpmax_yt = 0.001;
emittance_test = sqrt( (xmax_xt*xpmax_yt)^2 + (xmax_xt*xpmax_yt)*(xmax_yt*xpmax_xt) );
%%
