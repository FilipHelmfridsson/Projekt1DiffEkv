

function dudt = planetrhs(u)
    x = u(1);
    y = u(2);
    xhastighet = u(3);
    yhastighet = u(4);
    
    r = sqrt(x^2 + y^2);
    
    xacceleration = -x/(r^3);
    yacceleration = -y/(r^3);

    dudt = [xhastighet; yhastighet; xacceleration; yacceleration];
end