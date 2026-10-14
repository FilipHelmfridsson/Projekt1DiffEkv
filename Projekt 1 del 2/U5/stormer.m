function resultat = stormer(u,N,tslut)

    resultatVektor = u';
    for i = [1:(N-1)]

        h = tslut/N;
        du = planetrhs(u);

        u = [u(1) + h*(du(1)) + (1/2)*(h^2)*du(3), u(2) + h*(du(2)) + (1/2)*(h^2)*du(4), u(3), u(4)];

        gamla_ax = du(3);
        gamla_ay = du(4);

        du = planetrhs(u);

        u = [u(1), u(2), u(3) + (1/2)*h*(gamla_ax + du(3)),  u(4) + (1/2)*h*(gamla_ay + du(4))];

        resultatVektor = [resultatVektor;u];
    end
    resultat = resultatVektor;

end