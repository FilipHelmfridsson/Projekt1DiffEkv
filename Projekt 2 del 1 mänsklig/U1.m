clear all; close all;


u = @(x) sin(exp(x));
du = @(x) exp(x).*cos(exp(x));    % exakt derivata
x  = 0.75;
k = [1:8];
h  = 2.^(-k);

% Formel (2) centrerad
Derivata2 = (u(x+h)-u(x-h))./(2*h);

% Formel (6) framåt
Derivata6 = (-u(x+2*h)+4*u(x+h)-3*u(x))./(2*h);

% Formel (7) bakåt
Derivata7 = (3*u(x)-4*u(x-h)+u(x-2*h))./(2*h);

% De olika felens noggranhetsordning
E2 = (abs(Derivata2-du(x)))
E6 = (abs(Derivata6-du(x)))
E7 = (abs(Derivata7-du(x)))

text(k(end), E2(end), '  (2)')
text(k(end), E6(end), '  (6)')
text(k(end), E7(end), '  (7)')

p2 = polyfit(log(h), log(E2), 1)
p6 = polyfit(log(h), log(E6), 1)
p7 = polyfit(log(h), log(E7), 1)

ref = h.^2;


loglog(h, E2, h, E6, h, E7, h, ref)
legend('Formel (2)', 'Formel (6)','Formel (7)', 'Referens h^2', 'Location','best')
axis square;
grid on;
axis equal;
% hold on;
% loglog(h, E6)
% hold on;
% loglog(h, E7)
% hold on;
% loglog(h, ref)


