clear all;
close all;
T_perioder = 40 %Antal omloppsbanor
tslut = T_perioder * 2*pi;
h = 0.01;
N = tslut / h;
u = [1;0;0;1];
exakt_slut = [cos(tslut),sin(tslut),-sin(tslut),cos(tslut)]

for i = [1]
    Ni = N * 2^i;

    u_euler = euler(u,Ni,tslut);
    u_trapets = trapets(u,Ni,tslut);
    u_stormer = stormer(u,Ni,tslut);
end

E0 = - 1/2; %Beräknat utifrån begynnelsevärden


E_euler = 1/2*(u_euler(:,3).^2 + u_euler(:,4).^2) - 1./sqrt(u_euler(:,1).^2 + u_euler(:,2).^2);

E_trapets = 1/2*(u_trapets(:,3).^2 + u_trapets(:,4).^2) - 1./sqrt(u_trapets(:,1).^2 + u_trapets(:,2).^2);

E_stormer = 1/2*(u_stormer(:,3).^2 + u_stormer(:,4).^2) - 1./sqrt(u_stormer(:,1).^2 + u_stormer(:,2).^2);


deltaE_euler = E_euler - E0;
deltaE_trapets = E_trapets - E0;
deltaE_stormer = E_stormer - E0;

% Tidsvektor
h2 = tslut/Ni;
t = [0:h2:tslut];


plot(t,deltaE_euler,'r')
figure;
plot(t,deltaE_trapets,'b')
figure;
plot(t,deltaE_stormer,'g')

% Energifel
deltaE_euler = E_euler(end) - E0
deltaE_trapets = E_trapets(end) - E0
deltaE_stormer = E_stormer(end) - E0

% b)

T_perioder = 2500; %Antal omloppsbanor
tslut = T_perioder * 2*pi;
N = tslut / h;

u_stormer2 = stormer(u,Ni,tslut);
E_stormer2 = 1/2*((u_stormer(end,3)^2) + u_stormer(end,4)^2) - 1/sqrt((u_stormer(end,1)^2) + (u_stormer(end,2)^2));
deltaE_stormer2 = E_stormer2 - E0

%"Felet" i energin växer inte alls med fler perioder! Den hålls konstant,  
%vilket matchar dess beskrivning med att "kretsa" kring energin.

% c) 
% Störmer-Verlet rekommenderas för miljontals perioder för att den begränsar energifelet.
% För väldigt få eller enskilda perioder fungerar trapetsmetoden bra, med marginellt större beräkningskostnad
% jämtemot Euler framåt, då den har högre konvergensordning och bevarar banans form, då dess R(z) = 1 för många steglängder
% Euler har få användsningsområden mer än för utbildningssyfte då den är dyr (pga kravet av mycket små tidssteg),
% låg precision samt stabilitet och är usel på att bevara energi. Men, väldigt få rader kod, nice!