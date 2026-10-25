clear all;
close all;

% (2) Newton
function J = planetJacobian(u)
    x = u(1);
    y = u(2);

    r2 = x^2 + y^2;
    r5 = r2^(5/2);

    J31 = (3*x^2 - r2) / r5;
    J32 = 3*x*y / r5;
    J41 = J32;
    J42 = (3*y^2 - r2) / r5;

    J = [0,   0,   1, 0;
         0,   0,   0, 1;
         J31, J32, 0, 0;
         J41, J42, 0, 0];
end

function u_ny = trapetsmetoden(u, h)
    u = u(:);
    tol = 1e-6;
    z = u;

    while 1

        % Residual för trapetsmetoden
        G = z - u - h/2 * (planetrhs(u) + planetrhs(z));

        % Jacobian för högerledet vid aktuellt z
        J = planetJacobian(z);

        % Jacobian för residualen
        DG = eye(4) - h/2 * J;

        % Newton-uppdatering
        delta = DG \ G;
        z_ny = z - delta;

        % Konvergenskontroll
        if norm(z_ny - z, inf) < tol
            u_ny = z_ny;
            return
        end

        z = z_ny;
    end
end

N = 100;
u = [1;0;0;1]; % MÅSTE GISSA DENNA
tslut = 2*pi;
h = tslut/N;

resultatVektor = [];
for i = [0:h:tslut]
    u = trapetsmetoden(u, h);
    resultatVektor = [resultatVektor, u];
end

plot(resultatVektor(1,:), resultatVektor(2,:));

axis equal;

tolerans = 1e-6;

% while 1
%     g_varde = g(aR, h);
%     dg_varde = eye(4) - h/2 * J;

%     u_ny = u - (g_varde / dg_varde);
%     fel = abs(u_ny - u);
    
%     disp("Iteration:")
%     disp(i)
%     disp("Nytt aR värde:")
%     disp(aR_ny)
%     disp("Felet är:")
%     disp(fel)
%     if i >= 3
%         konvergensOrdning = log(fel / felGammal) / log(felGammal / felGammalGammal);
%         disp("Konvergensordningen är:")
%         disp(konvergensOrdning)
%     else
%         konvergensOrdning = 0;
%         disp("Det krävs ett par iterationer för att beräkna konvergensordning")
%     end

%     if max(fel) < tolerans
%         disp('hittat aR!')
%         disp(aR_ny)
%         break
%     end

%     aR = aR_ny;
%     felGammalGammal=felGammal;
%     felGammal = fel;
%     i = i+1;
% end