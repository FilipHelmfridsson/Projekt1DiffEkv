clear all;
close all;
tslut = 2*pi/4;
N = 100;
u = [1;0;0;1];
exakt_slut = [cos(tslut),sin(tslut),-sin(tslut),cos(tslut)]

felvektor_euler = [];
felvektor_stormer = [];
felvektor_trapets = [];

for i = [0:1]
    Ni = N * 2^i;

    u_euler = euler(u,Ni,tslut);
    u_trapets = trapets(u,Ni,tslut);
    u_stormer = stormer(u,Ni,tslut);

    fel_euler = norm(u_euler(end,:) - exakt_slut);
    fel_trapets = norm(u_trapets(end,:) - exakt_slut);
    fel_stormer = norm(u_stormer(end,:) - exakt_slut);

    felvektor_euler = [felvektor_euler;fel_euler];
    felvektor_trapets = [felvektor_trapets;fel_trapets];
    felvektor_stormer = [felvektor_stormer;fel_stormer];
end

konvergensordning_euler = log10(felvektor_euler(1) / felvektor_euler(2)) / log10(2)

konvergensordning_trapets = log10(felvektor_trapets(1) / felvektor_trapets(2)) / log10(2)

konvergensordning_stormer = log10(felvektor_stormer(1) / felvektor_stormer(2)) / log10(2)

loglog([tslut/(2*N);tslut/N],felvektor_euler,'r')
hold on;
loglog([tslut/(2*N);tslut/N],felvektor_trapets,'g')
hold on;
loglog([tslut/(2*N);tslut/N],felvektor_stormer,'b')