%MAE 162A Final Project Spring 2021
%Duncan Di Mauro
%UID: 805163177
%Due: 6/11/2021 10:00 AM

clear all; close all; clc;

%% Establishing Constants

r1 = 0.805; %meters
r2 = r1*(5/4);
r3 = r1*(1/2);
r4 = r1*(5/4);
r5 = r1*(5/8);
r6 = r1*(1/2);
r26 = r1*(1/4);
r15 = r1*(5/8);

%Note that t2 means theta2. Same with other thetas.
t2_min_degrees = 36.88;
t2_min = t2_min_degrees*(pi/180);
t2_max_degrees = 96.8;
t2_max = t2_max_degrees*(pi/180);

%Making the t2 input array, goes from the minimum to maximum in radians
t2 = linspace(t2_min, t2_max, 1000);


%% Position Analysis - find t3, t4, t5, and t6 for all t2

[t3,t4] = pos_L1(r1, r2, r3, r4, t2, 1);
[t5,t6] = pos_L2(r2, r26, r6, r5, r15, r1, t2, t3, t4, 1);

figure(2);
plot(t2, t6)
title('Rotation of Tabletop in Z-Direction');
%axis([0.3,2,-10^-13,10^-13])
xlabel('Input Angle, Theta 2 (radians)');
ylabel('Theta 6 (radians)');


%% Link 6 Center of Mass Displacement Calculations

%Finding center of mass locations relative to the tail of r6
r6_cm_x = r6*(1/2); 
r6_cm_y = 0.14915; %in meters, found from SolidWorks
r6_cm_mag = sqrt(r6_cm_x^2 + r6_cm_y^2);

t_r6_cm = atan(r6_cm_y/r6_cm_x); % = radians counterclockwise from r6

%Finding center of mass locations relative to the origin of the system
xg6_1 = r2*cos(t2) + r26*cos(t3) + r6_cm_mag*cos(t6 + t_r6_cm);
yg6_1 = r2*sin(t2) + r26*sin(t3) + r6_cm_mag*sin(t6 + t_r6_cm);

%Finding the ranges of xg6_1 and yg6_1
range_xg6_1 = max(xg6_1) - min(xg6_1);
range_yg6_1 = max(yg6_1) - min(yg6_1);

% Making the displacement vectors have a minimum value of 0 
xg6_1 = xg6_1 - min(xg6_1);
yg6_1 = yg6_1 - min(yg6_1);

figure(3); 
plot(t2,xg6_1*1000);
title('Tabletop Center of Mass Horizontal Displacement');
xlabel('Input Angle, Theta 2 (radians)');
ylabel('Displacement (mm)');

figure(4);
plot(t2,yg6_1*1000);
title('Tabletop Center of Mass Vertical Displacement');
xlabel('Input Angle, Theta 2 (radians)');
ylabel('Displacement (mm)');


%% Motor Driving Link 2 Setup, Position Analysis, and Kinematic Analysis

w = 1*(pi/180); %degrees/sec

%Creating time array for 3 full cycles
%seconds needed = 360*3 = 180
t = linspace(0, 1080, 1081);

%Modeling theta 2 based on the motor (t2m means theta 2 motor)
t2m = (t2_max - t2_min)*0.5*sin(w*t) + (t2_max + t2_min)*0.5;

%Finding 1st and 2nd derivatives of t2m (w2 and a2)
w2 = (1/2)*w*(t2_max - t2_min)*cos(w*t);
a2 = -(1/2)*w^2*(t2_max - t2_min)*sin(w*t);

%Finding angles and first/second order coefficients when motor is applied
[t3m, t4m] = pos_L1(r1,r2,r3,r4,t2m,2);
[t5m, t6m] = pos_L2(r2,r26,r6,r5,r15,r1,t2m,t3m,t4m,2);

[t3mp,t4mp] = FirstOrder_L1(r2,r3,r4,t2m,t3m,t4m);
[t5mp,t6mp] = FirstOrder_L2(r2,r26,r6,r5,r15,t2m,t3m,t4m,t5m,t6m,t3mp,t4mp);

[t3mpp,t4mpp] = SecondOrder_L1(r2,r3,r4,t2m,t3m,t4m,t3mp,t4mp);
[t5mpp,t6mpp] = SecondOrder_L2(r2,r26,r6,r5,r15,t2m,t3m,t4m,t5m,t6m,t3mp,t4mp,t5mp,t6mp,t3mpp,t4mpp);

%% Center of mass location, velocity, and acceleration of each link

%Finding center of mass for each link
xg2 = 0.5*r2*cos(t2m);
yg2 = 0.5*r2*sin(t2m);

xg3 = r2*cos(t2m) + 0.5*r3*cos(t3m);
yg3 = r2*sin(t2m) + 0.5*r3*sin(t3m);

xg4 = r2*cos(t2m) + r3*cos(t3m) + 0.5*r4*cos(t4m);
yg4 = r2*sin(t2m) + r3*sin(t3m) + 0.5*r4*sin(t4m);

xg5 = xg3 + r6*cos(t6m) + 0.5*r5*cos(t5m);
yg5 = yg3 + r6*sin(t6m) + 0.5*r5*sin(t5m);

xg6 = r2*cos(t2m) + r26*cos(t3m) + r6_cm_mag*cos(t_r6_cm + t6m);
yg6 = r2*sin(t2m) + r26*sin(t3m) + r6_cm_mag*sin(t_r6_cm + t6m);

% Finding velocity for each center of mass
xg2p = -0.5*r2*sin(t2m);
yg2p = 0.5*r2*cos(t2m);

xg3p = -r2*sin(t2m) - 0.5*r3*sin(t3m).*t3mp;
yg3p = r2*cos(t2m) + 0.5*r3*cos(t3m).*t3mp;

xg4p = -r2*sin(t2m) - r3*cos(t3m).*t3mp - 0.5*r4*sin(t4m).*t4mp;
yg4p = r2*cos(t2m) + r3*cos(t3m).*t3mp + 0.5*r4*cos(t4m).*t4mp;

xg5p = xg3p - r6*sin(t6m).*t6mp - 0.5*r5*sin(t5m).*t5mp;
yg5p = yg3p + r6*cos(t6m).*t6mp + 0.5*r5*cos(t5m).*t5mp;

xg6p = -r2*sin(t2m) - r26*sin(t3m).*t3mp - r6_cm_mag*sin(t_r6_cm + t6m).*t6mp;
yg6p = r2*cos(t2m) + r26*cos(t3m).*t3mp + r6_cm_mag*cos(t_r6_cm + t6m).*t6mp;

% Finding accelerations for each center of mass
xg2pp = -0.5*r2*cos(t2m);
yg2pp = -0.5*r2*sin(t2m);

xg3pp = -r2*cos(t2m) - 0.5*r3*cos(t3m).*t3mp.^2 - 0.5*r3*sin(t3m).*t3mpp;
yg3pp = -r2*sin(t2m) - 0.5*r3*sin(t3m).*t3mp.^2 + 0.5*r3*cos(t3m).*t3mpp;

xg4pp = -r2*cos(t2m) - r3*cos(t3m).*t3mp.^2 - r3*sin(t3m).*t3mpp - 0.5*r4*cos(t4m).*t4mp.^2 - 0.5*r4*sin(t4m).*t4mpp;
yg4pp = -r2*sin(t2m) - r3*sin(t3m).*t3mp.^2 + r3*cos(t3m).*t3mpp - 0.5*r4*sin(t4m).*t4mp.^2 + 0.5*r4*cos(t4m).*t4mpp;

xg5pp = xg3pp - r6*cos(t6m).*t6mp.^2 - r6*sin(t6m).*t6mpp - 0.5*r5*cos(t5m).*t5mp.^2 - 0.5*r5*sin(t5m).*t5mpp;
yg5pp = yg3pp - r6*sin(t6m).*t6mp.^2 + r6*cos(t6m).*t6mpp - 0.5*r5*sin(t5m).*t5mp.^2 + 0.5*r5*cos(t5m).*t5mpp;

xg6pp = -r2*cos(t2m) - r26*cos(t3m).*t3mp.^2 - r26*sin(t3m).*t3mpp - r6_cm_mag*cos(t_r6_cm + t6m).*t6mp.^2 - r6_cm_mag*sin(t_r6_cm + t6m).*t6mpp;
yg6pp = -r2*sin(t2m) - r26*sin(t3m).*t3mp.^2 + r26*cos(t3m).*t3mpp - r6_cm_mag*sin(t_r6_cm + t6m).*t6mp.^2  + r6_cm_mag*cos(t_r6_cm + t6m).*t6mpp;


%Calculating mass of each link
rho = 2698.9;  % Density of aluminum in kg/m^3
m1 = r1 * rho * 0.03^2;    
m2 = r2 * rho * 0.03^2;   
m3 = r3 * rho * 0.03^2;  
m4 = r4 * rho * 0.03^2;   
m5 = r5 * rho * 0.03^2;   
m6 = r6 * rho * 0.03^2;    

I1 = m1/12*(0.03^2 + r1^2); 
Ig2 = m2/12*(0.03^2 + r2^2);
Ig3 = m3/12*(0.03^2 + r3^2);
Ig4 = m4/12*(0.03^2 + r4^2);
Ig5 = m5/12*(0.03^2 + r5^2);
Ig6 = m6/12*(0.03^2 + r6^2);

A = m2*(xg2p.^2+yg2p.^2) + Ig2 + m3*(xg3p.^2+yg3p.^2) + Ig3*t3mp.^2 + ...
    m4*(xg4p.^2+yg4p.^2) + Ig4*t4mp.^2 + m5*(xg5p.^2+yg5p.^2) + ...
    Ig5*t5mp.^2 + m6*(xg6p.^2 + yg6p.^2) + Ig6*t6mp.^2;

B = m2*(xg2p.*xg2pp + yg2p.*yg2pp) + m3*(xg3p.*xg3pp + yg3p.*yg3pp) + ...
    Ig3*t3mp.*t3mpp + m4*(xg4p.*xg4pp + yg4p.*yg4pp) + Ig4*t4mp.*t4mpp + ...
    m5*(xg5p.*xg5pp + yg5p.*yg5pp) + Ig5*t5mp.*t5mpp +...
    m6*(xg6p.*xg6pp + yg6p.*yg6pp) + Ig6*t6mp.*t6mpp;

% Kinetic Energy
KE = A.*w2.*a2 + B.*(w2.^3);

% Potential Energy
Ugr = 9.81*w2.*(m2*yg2p + m3*yg3p +m4*yg4p + m5*yg5p + m6*yg6p);

Torque = (KE + Ugr)./w2;
figure(5);

plot(t,Torque);
title('Input Torque of Motor vs. Time');
axis([0,1080,-3,6])
xlabel('Time (s)');
ylabel('Torque (Nm)');