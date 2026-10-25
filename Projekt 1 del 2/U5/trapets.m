function resultat = trapets(u,N,tslut)

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

    h = tslut/N;

    resultatVektor = u';
    for i = 1:N
        u = trapetsmetoden(u, h);
        resultatVektor = [resultatVektor; u'];
    end
    resultat = resultatVektor;

end