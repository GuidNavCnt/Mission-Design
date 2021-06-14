%Name: Mohammad A. Edaibat ----- 2/19/2021
%Johns Hopkins University - APL
%Space Mission Design and Navigation
%This function calculates position and velocity vectors by taking the six
%classical orbital elements and mass parameter.
%Inputs
%---------------a: semimajor axis in Km
%---------------e: eccentricity
%---------------i: inclination
%---------------w: argument of periapsis
%---------------Omega: Longitude of ascending node
%---------------Theta: True Anomaly
%---------------muo: central body mass of the sun
%Outputs
%---------------rp: position in perifocal frame
%---------------vp: velocity in perifocal frame
%---------------r: position vector
%---------------v: velocity vector
%verify this function using the following values
%r=[227939282.200749 -11219880.0592502 2764663.06791779] Km
%v=[-3.56447813955076 21.9226854955848 0.25630583566074] Km/sec
%Use the following inputs to verify the above position and velocity values
%a=198200000
%e=0.2559
%i=0.0188 Radians
%Omega=5.5341 Radians
%w=3.0985 Radians
%Theta=3.8847 Radians
%muo=132712440041.94
%a=198200882.566171;e=0.2559;i=0.0188;w=3.0985;Omega=5.534;Theta=3.884;muo=132712440041.94;
function [r,v]=problem2_PositionVelocity(a,e,i,w,Omega,Theta,muo)
p = a*(1-(e^2)); %the semi-parameter
rp = p/(1+e*cos(Theta))*[cos(Theta);sin(Theta);0]; %km, position vector in perifocal frame
vp = sqrt(muo/p)*[-sin(Theta);e+cos(Theta);0]; %km/sec velcoity vector in perifocal frame
C1_i=[1,0,0;0,cos(i),sin(i);0,-sin(i),cos(i)] %used for frame rotation
C3_Omega=[cos(Omega),sin(Omega),0;-sin(Omega),cos(Omega),0;0,0,1] %used for frame rotation
C3_w=[cos(w),sin(w),0;-sin(w),cos(w),0;0,0,1] %used for frame rotation
C_rot=(C3_w*C1_i*C3_Omega)' %matrix rotated from perifocal frame to geocentric frame
format long g
r=C_rot*rp; %km, position vector in geocentric frame (after rotating from perifocal frame)
v=C_rot*vp; %km/sec, velocity vector in geocentric frame (after rotating from perifocal frame)
r=r'
v=v'