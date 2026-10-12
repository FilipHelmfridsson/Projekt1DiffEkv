function f = planetrhs(u)
% PLANETRHS  Högerled f(u) för Kepler-problemet, u = [x; y; vx; vy]
    x  = u(1);
    y  = u(2);
    vx = u(3);
    vy = u(4);

    r3 = (x^2 + y^2)^(3/2);

    f = [ vx;
          vy;
         -x/r3;
         -y/r3 ];
end