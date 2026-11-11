function [x, w] = heat_euler(N, M, T, D, ffcn, uinitfcn, TLfcn, TRfcn)
% Löser u_t = D u_xx + f(x,t) på [0,1] x [0,T] med Dirichletrandvillkor.
% Centrala differenser i rummet, Euler framåt i tiden.
% Returnerar x (N+1 punkter) och lösningen w vid t = T.

dx = 1/N;
dt = T/M;
x  = (0:N)' * dx;
w  = uinitfcn(x);
j  = 2:N;                       % inre punkter

for n = 0:M-1
    t = n*dt;
    w(j) = w(j) + dt*( D*(w(j+1) - 2*w(j) + w(j-1))/dx^2 + ffcn(x(j), t) );
    w(1)   = TLfcn(t + dt);
    w(end) = TRfcn(t + dt);
end
end
