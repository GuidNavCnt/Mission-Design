%Name: Mohammad A. Edaibat ----- 2/19/2021
%Johns Hopkins University - APL
%Space Mission Design and Navigation
%This function takes a position and velocity vector and central body mass
%parameter and calculates the six classical orbital elements
%Inputs
%---------------r: position vector
%---------------v: velocity vector
%---------------muo: central body mass of the sun
%Outputs
%---------------a: semimajor axis in Km
%---------------e: eccentricity
%---------------i: inclination
%---------------w: argument of periapsis
%---------------Omega: Longitude of ascending node
%---------------Theta: True Anomaly
%verify this function using the following values
%r=[227939282.200749 -11219880.0592502 2764663.06791779]
%v=[-307970.911257186 1894120.02681853 22144.8242010879]
%muo=132712440041.94
%r=[227939282.200749 -11219880.0592502 2764663.06791779];v=[-3.56447813955076 21.9226854955848 0.25630583566074];muo=132712440041.94;
function [a,e,i,w,Omega,Theta] = problem1_sixOrbitalElements(r,v,muo)
I_hat = [1 0 0];
K_hat = [0 0 1];
Energy = (((dot(v,v))/2)-(muo/norm(r)));
a = -muo/(2*Energy) %Km, this is the semimajor axis
h = cross(r,v); %Km^2/sec, this is the angular momentum
N_hat = cross(K_hat,h)/norm(cross(K_hat,h));
e_vector = (cross(v,h)/muo)-(r/norm(r)); %this is the eccentricity vector
e = norm(e_vector) %this is the eccentricity value
%e = sqrt(1-((dot(h,h))./(muo*a))) Another equation for eccentricity
i = acos(dot(h,K_hat)/norm(h)) %rad, this is the inclination angle
Omega = acos(dot(N_hat,I_hat)); %rad, this is the ascending node Omega
if N_hat(2)<0
    Omega=2*pi-Omega
end
w = acos(dot(N_hat,e_vector)/e) %rad, this is the argument of periapsis angle
if e_vector(3)<0
    w=2*pi-w
end
Theta=acos(dot(e_vector,r)/(e*norm(r)));
if dot(r,v)<0
    Theta=2*pi-Theta
end