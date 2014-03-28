% Random phase space
x0  = 5;
xp0 = 2*pi;
alpha = -10;
R1 = normrnd(0,x0,[1 1000]);
R2 = normrnd(0,xp0,[1 1000]);

figure(1)
subplot(1,2,1)
histfit(R1,40)
legend('x')
subplot(1,2,2)
histfit(R2,40)
legend('xp')


for i=1:max(size(R1))
    x(i) = R1(i);
    xp(i) = R2(i) - alpha*x(i);
end

figure(2)
plot(x,xp,'o')