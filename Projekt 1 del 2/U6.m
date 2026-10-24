% framåt euler R(z) = 1+z
% Funktionen är absolut stabil då abs(R(z))<1 efterssom yi --> 0
% 1+z < 1
% z < 0
% z = (lambda * h) < 0, så lambda måste vara negativt för euler frammåt.


% bakåt euler R(z)= 1/(1-z)

% 1/(1-z) < 1
% 1-z > 1



% trapetsmetoden R(z) = (2+z)/(2-z)



omega_i = 1i;
h = 0.1;
z = omega_i*h;

% Framåt Euler
R_euler = 1 + z;
belopp_euler = abs(R_euler)

R_trapets = (2+z)/(2-z);
belopp_trapets = abs(R_trapets)

%För euler spiralar cirkeln utåt, vilket stämmer med våra resultat, då R är något större än 1.

%För trapetsmetoden håller sig cirkeln i princip rund, vilket också stämmer med våra resultat, då R = 1!

%Hade R varit mindre än 1 hade banan spiralat in mot centrum.