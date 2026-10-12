% U2b: Framåt Euler, en period
T  = 2*pi;
N  = 200;               % antal tidssteg (välj själv)
dt = T/N;

u = zeros(4, N+1);
u(:,1) = [1; 0; 0; 1];

for n = 1:N
    u(:,n+1) = u(:,n) + dt*planetrhs(u(:,n));
end

th = linspace(0, 2*pi, 400);
plot(cos(th), sin(th), 'k--', 'LineWidth', 1.2); hold on
plot(u(1,:), u(2,:), 'r-', 'LineWidth', 1.2)
axis equal; grid on
legend('Exakt cirkel', 'Framåt Euler')
xlabel('x'); ylabel('y')
title(sprintf('Framåt Euler, N = %d steg', N))