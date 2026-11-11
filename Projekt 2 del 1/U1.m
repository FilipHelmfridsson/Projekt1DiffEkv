% U1b: Fel i differensformlerna (2), (6) och (7) för u = sin(e^x) i x = 0.75
clear; close all;

u  = @(x) sin(exp(x));
du = @(x) exp(x).*cos(exp(x));    % exakt derivata
x  = 0.75;
h  = 2.^-(1:8);

central = (u(x+h) - u(x-h)) ./ (2*h);                  % (2)
framat  = (-u(x+2*h) + 4*u(x+h) - 3*u(x)) ./ (2*h);    % (6)
bakat   = (3*u(x) - 4*u(x-h) + u(x-2*h)) ./ (2*h);     % (7)

loglog(h, abs(central - du(x)), 'o-', ...
       h, abs(framat  - du(x)), 's-', ...
       h, abs(bakat   - du(x)), 'd-', ...
       h, h.^2, 'k--')
legend('Central (2)', 'Framåt (6)', 'Bakåt (7)', 'Referens h^2', 'Location', 'northwest')
xlabel('h'), ylabel('|fel|'), grid on
