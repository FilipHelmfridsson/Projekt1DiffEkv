clear all;
close all;

N = 200;
tslut = 2*pi;
u = [1;0;0;1];

for i = [1:(N-1)]
    h = tslut / N;
    u = [u, u(:,end) + h * planetrhs(u(:,end))];
end

theta = 0:0.0001:2*pi;
plot(cos(theta),sin(theta));
hold on;
plot(u(1,:),u(2,:));
axis equal;
hold off;