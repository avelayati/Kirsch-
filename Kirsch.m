%%
%%% Arian Velayati, PhD
%%%% This script is used to find normal and shear stresses around a circular cavity in a homogeneous linear elastic solid . The complete Kirsch solution assumes independent action of multiple factors,
%%%% namely far-field isotropic stress, deviatoric stress, wellbore pressure and pore pressure.

clc; clear; close;

%% Inputs

% teta = ; % Angle between the direction of SHmax and the point in which stress is considered
% r = ; % The distance r is measured from the center of the wellbore. at the wellbore wall r=a
a = 1; % At the wellbore wall r=a
SH = 22; Sh = 13; Sv = 25;  % Effective stresses (MPa) "SigH-aPp"
Pw = 10; % Wellbore pressure
Pp = 10; % Pore pressure 
v = 0.25; %Poisson ratio

%% Calculations: The Kirsch solution for a vertical wellbore with radius a within a linear elastic and isotropic solid is:

r = 1:0.1:5;
teta = 0:5:360;


for i = 1:length(r)
    for j = 1:length(teta)
Srr(i,j) = (Pw-Pp)*(power(a,2)/power(r(i),2)) + 0.5*(SH + Sh)*(1-(power(a,2)/power(r(i),2))) + 0.5*(SH-Sh)*(1-4*a^2/r(i)^2 + 3*a^4/r(i)^4)*cosd(2*teta(j)) ; % Radial effective stress
Stt(i,j) = -(Pw-Pp)*(power(a,2)/power(r(i),2)) + 0.5*(SH + Sh)*(1+(power(a,2)/power(r(i),2))) - 0.5*(SH-Sh)*(1+3*a^4/r(i)^4)*cosd(2*teta(j)) ; % Tangential (hoop) Stress
Srt(i,j) = 0.5*(SH-Sh)*(1 + 2*a^2/r(i)^2 - 3*a^4/r(i)^4)*sind(2*teta(j)); % Shear stress in a plane perpendicular to rin tangential direction teta
Szz(i,j) = Sv - 2*v*(SH-Sh)*(a^2/r(i)^2)*cosd(2*teta(j)) ; % Vertical effective stress in direction z
    end
end

%% plot

% At the wellbore wall r = a 

figure(1)
plot(teta, Stt(1,:),'r','LineWidth',3)
hold on
plot(teta, Srr(1,:),'b','LineWidth',3)
plot(teta, Srt(1,:),'g','LineWidth',3)
plot(teta, Szz(1,:),'y','LineWidth',3)
xlabel('teta (deg)')
ylabel('Stress (MPa)')
legend('Stt','Srr','Srt','Szz')


figure(2) %  deg
plot(r, Stt(:,10),'r','LineWidth',3)
hold on
plot(r, Srr(:,10), 'b','LineWidth',3)
plot(r, Srt(:,10), 'g','LineWidth',3)
plot(r, Szz(:,10), 'y','LineWidth',3)

xlabel('r (m)')
ylabel('Stress (MPa)')
legend('Stt','Srr','Srt','Szz')
