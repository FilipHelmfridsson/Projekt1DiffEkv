clear all;
close all;

% Vi väljer STÖRMER-Verlet

N = 100;

tslut = 2*pi;

u = [1,0,0,1];

resultatVektor = u;
for i = [1:N]

    h = tslut/N;
    du = planetrhs(u);

    u = [u(1) + h*(du(1)) + (1/2)*(h^2)*du(3), u(2) + h*(du(2)) + (1/2)*(h^2)*du(4), u(3), u(4)];

    gamla_ax = du(3);
    gamla_ay = du(4);

    du = planetrhs(u);

    u = [u(1), u(2), u(3) + (1/2)*h*(gamla_ax + du(3)),  u(4) + (1/2)*h*(gamla_ay + du(4))];

    resultatVektor = [resultatVektor;u];
end

plot(resultatVektor(:,1),resultatVektor(:,2));
axis equal;

%% Resonemang:

% Vi valde STÖRMER-Verlet för att den har andra ordningens noggranhet, 
% men har så gott som samma beräkningskrav som euler.

% 