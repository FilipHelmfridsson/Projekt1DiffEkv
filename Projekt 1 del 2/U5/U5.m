clear all;
close all;
tslut = 2*pi/4;
N = 100;
u = [1;0;0;1];
exakt_slut = [cos(tslut),sin(tslut),-sin(tslut),cos(tslut)]

u_euler = euler(u,N,tslut);
u_trapets = trapets(u,N,tslut);
u_stormer = stormer(u,N,tslut);

fel_euler = abs(u_euler(end,:) - exakt_slut);
fel_trapets = abs(u_trapets(end,:) - exakt_slut);
fel_stormer = abs(u_stormer(end,:) - exakt_slut);

felvektor_euler = [];
felvektor_stormer = [];
felvektor_trapets = [];

for i = [0:1]
    N = N*(2^i);

    u_euler = euler(u,N,tslut);
    u_trapets = trapets(u,N,tslut);
    u_stormer = stormer(u,N,tslut);

    fel_euler = abs(u_euler(end,:) - exakt_slut);
    fel_trapets = abs(u_trapets(end,:) - exakt_slut);
    fel_stormer = abs(u_stormer(end,:) - exakt_slut);

    felvektor_euler = [felvektor_euler;fel_euler];
    felvektor_trapets = [felvektor_trapets;fel_trapets];
    felvektor_stormer = [felvektor_stormer;fel_stormer];
end

konvergensordning_euler = (log10(felvektor_euler(1,:) ./ felvektor_euler(2,:))) ./ log10(2)

konvergensordning_stormer = (log10(felvektor_stormer(1,:) ./ felvektor_stormer(2,:))./log10(2))

konvergensordning_trapets = (log10(felvektor_trapets(1,:) ./ felvektor_trapets(2,:))./log10(2))

%konvergensordning_euler = log(abs(felvektor_euler(1,:) - felvektor_euler(2,:))) ./ 