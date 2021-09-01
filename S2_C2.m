clear;
clc;

%% Time Setting

sec2day = 1/(3600 * 24);                    % 1 day = 3600 * 24 secs
t_ref = 2451545.0;                          % J2000 reference time point
t_initial = 365.25 * 21;                    % departure time(2021.01.01. 00:00)
t_step = 1;                                 % time step(sec)
t_final = 405890;                           % Before the Second-kick(sec)
t = 1:t_step:t_final;                       % Whole time interval

%% Orbit Preprocessing

d = 384400;                                 % Radius of the Moon orbit(km)
mu_E = 3.98595624*10^5;                     % Standard gravitaional parameter of the Earth(km^3/sec^2)
mu_M = 4.90278576*10^3;                     % Standard gravitaional parameter of the Moon(km^3/sec^2)
w = sqrt((mu_E+mu_M)/d^3);                  % Angular velocity of the Moon(rad/sec)

theta = rad2deg((t-1) * w);                 % Transformation angle(deg)

Observer = load("Orbit.mat");
Observer = Observer.Orbit;
Observer = Observer(:, t);                  % Position of the observer at the target time interval

% Coordinate transformation (rotational -> fixed)
for i = 1:length(t)
    Observer(:, i) = RotMx(4, theta(i)) * Observer(:, i);
end

%% Orbit Transformation

% Axis generate
% X : Initial direction of the Moon
% Z : cross product of the X and the initial velocity vector
% Y : cross product of the Z and X
[Moon0_p, Moon0_v] = planetEphemeris([t_ref t_initial], 'Earth', 'Moon', '405');
X = ToUnit(Moon0_p');
Z = ToUnit(cross(X, Moon0_v'));
Y = cross(Z, X);

% Get Planets position (from DE405)
VenusEph = planetEphemeris(t'*sec2day+t_initial+t_ref, 'Earth', 'Venus', '405');
MarsEph = planetEphemeris(t'*sec2day+t_initial+t_ref, 'Earth', 'Mars', '405');

% 2D projection (to the orbital plane of the Moon, XY plane)
Venus.P = [VenusEph*X VenusEph*Y]';
Mars.P = [MarsEph*X MarsEph*Y]';
Moon.P = d * [cosd(theta);sind(theta)];
Earth.P = zeros(2, length(theta));

%% Additional information

clear VenusEph MarsEph Moon0_p Moon0_v

% Observer to planets direction vector(true value)
Venus.E = (Venus.P-Observer)./vecnorm((Venus.P-Observer), 2, 1);
Mars.E = (Mars.P-Observer)./vecnorm((Mars.P-Observer), 2, 1);
Moon.E = (Moon.P-Observer)./vecnorm((Moon.P-Observer), 2, 1);
Earth.E = (Earth.P-Observer)./vecnorm(Earth.P-Observer, 2, 1);

% Angle between the direction vector of the planets and X axis
Venus.Ang = acosd(Venus.E(1, :)) .* (2*(Venus.E(2, :)>0)-1);
Mars.Ang = acosd(Mars.E(1, :)) .* (2*(Mars.E(2, :)>0)-1);
Moon.Ang = acosd(Moon.E(1, :)) .* (2*(Moon.E(2, :)>0)-1);
Earth.Ang = acosd(Earth.E(1, :)) .* (2*(Earth.E(2, :)>0)-1);

%% Case #2

% Setting parameters
Case2.sigma = [1;1;1;1] * 10/3600;
Case2.Observer = Observer;
Case2.t = t;

% Making empty matrix for the results of the simulation
Case2.P_hat_LS = zeros(2, length(t));
Case2.P_hat_WLS = zeros(2, length(t));
Case2.P_hat_2P = zeros(2, length(t));
Case2.Error_LS = zeros(2, length(t));
Case2.Error_WLS = zeros(2, length(t));
Case2.Error_2P = zeros(2, length(t));
Case2.RMS_LS = zeros(3, 2);
Case2.RMS_WLS = zeros(3, 2);
Case2.RMS_2P = zeros(3, 2);
Case2.Ang = abs(abs(Earth.Ang-Moon.Ang) - (360 * (abs(Earth.Ang-Moon.Ang)>180)));

% Generation noise
Noise = randn(4, length(t)) .* Case2.sigma;
e1 = [cosd(Moon.Ang+Noise(1, :));sind(Moon.Ang+Noise(1, :))];
e2 = [cosd(Earth.Ang+Noise(2, :));sind(Earth.Ang+Noise(2, :))];
e3 = [cosd(Venus.Ang+Noise(3, :));sind(Venus.Ang+Noise(3, :))];
e4 = [cosd(Mars.Ang+Noise(4, :));sind(Mars.Ang+Noise(4, :))];

% Simulation loop
for i = 1:length(t)
    Planets = [Moon.P(:, i) Earth.P(:, i) Venus.P(:, i) Mars.P(:, i)];
    E = [e1(:, i) e2(:, i) e3(:, i) e4(:, i)];
    Case2.P_hat_LS(:, i) = GetPosition2D_LS(Planets, E);
    if i==1
        Case2.P_hat_WLS(:, i) = GetPosition2D_LS(Planets, E);
    else
        Case2.P_hat_WLS(:, i) = GetPosition2D_WLS(Planets, E, P_prev);
    end
    Case2.P_hat_2P(:, i) = GetPosition2D_LS(Planets(:, 1:2), E(:, 1:2));
    
    Case2.Error_LS(:, i) = Case2.P_hat_LS(:, i)-Case2.Observer(:, i);
    Case2.Error_WLS(:, i) = Case2.P_hat_WLS(:, i)-Case2.Observer(:, i);
    Case2.Error_2P(:, i) = Case2.P_hat_2P(:, i)-Case2.Observer(:, i);
    P_prev = Case2.P_hat_WLS(:, i);
end

Case2.RMS_LS(1:2, 1) = rms(Case2.Error_LS, 2);
Case2.RMS_LS(3, 1) = norm(Case2.RMS_LS(1:2, 1));
Case2.RMS_WLS(1:2, 1) = rms(Case2.Error_WLS, 2);
Case2.RMS_WLS(3, 1) = norm(Case2.RMS_WLS(1:2, 1));
Case2.RMS_2P(1:2, 1) = rms(Case2.Error_2P, 2);
Case2.RMS_2P(3, 1) = norm(Case2.RMS_2P(1:2, 1));

Case2.CriticalTime1 = find(Case2.Ang>170);
Case2.CriticalTime2 = find(Case2.Ang<10);

Case2.RMS_LS(1:2, 2) = sqrt((length(t)*rms(Case2.Error_LS, 2).^2 - length([Case2.CriticalTime1 Case2.CriticalTime2])*rms(Case2.Error_LS(:, [Case2.CriticalTime1 Case2.CriticalTime2]), 2).^2)/(length(t)-length([Case2.CriticalTime1 Case2.CriticalTime2])));
Case2.RMS_LS(3, 2) = norm(Case2.RMS_LS(1:2, 2));
Case2.RMS_WLS(1:2, 2) = sqrt((length(t)*rms(Case2.Error_WLS, 2).^2 - length([Case2.CriticalTime1 Case2.CriticalTime2])*rms(Case2.Error_WLS(:, [Case2.CriticalTime1 Case2.CriticalTime2]), 2).^2)/(length(t)-length([Case2.CriticalTime1 Case2.CriticalTime2])));
Case2.RMS_WLS(3, 2) = norm(Case2.RMS_WLS(1:2, 2));
Case2.RMS_2P(1:2, 2) = sqrt((length(t)*rms(Case2.Error_2P, 2).^2 - length([Case2.CriticalTime1 Case2.CriticalTime2])*rms(Case2.Error_2P(:, [Case2.CriticalTime1 Case2.CriticalTime2]), 2).^2)/(length(t)-length([Case2.CriticalTime1 Case2.CriticalTime2])));
Case2.RMS_2P(3, 2) = norm(Case2.RMS_2P(1:2, 2));

%% Save the Results

save('S2_C2.mat', 'Case2');