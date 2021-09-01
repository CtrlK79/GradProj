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
%% Additional Information

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

%% Case #1

% Setting parameters
Distance = vecnorm(Observer, 2, 1);
Case1.T = [length(t)-length(find(Distance>(d/5)));length(t)-length(find(Distance>(2*d/5)));
    length(t)-length(find(Distance>(3*d/5)));length(t)-length(find(Distance>(4*d/5)))];

Case1.sigma = [1;1;1;1] * 10/3600;
Case1.N = 50000;
Case1.Observer = Observer;

% Making empty matrix for the results of the simulation
Case1.P_hat_LS = zeros(2, Case1.N, length(Case1.T));
Case1.P_hat_WLS = zeros(2, Case1.N, length(Case1.T));
Case1.P_hat_2P = zeros(2, Case1.N, length(Case1.T));
Case1.Error_LS = zeros(2, Case1.N, length(Case1.T));
Case1.Error_WLS = zeros(2, Case1.N, length(Case1.T));
Case1.Error_2P = zeros(2, Case1.N, length(Case1.T));
Case1.RMS_LS = zeros(3, length(Case1.T));
Case1.RMS_WLS = zeros(3, length(Case1.T));
Case1.RMS_2P = zeros(3, length(Case1.T));
Case1.Ang = zeros(1, length(Case1.T));

% Simulation loop
for i = 1:length(Case1.T)
    Noise = randn(4, Case1.N) .* Case1.sigma;
    Planets = [Venus.P(:, Case1.T(i)) Mars.P(:, Case1.T(i)) Moon.P(:, Case1.T(i)) Earth.P(:, Case1.T(i))];
    e1 = [cosd(Venus.Ang(Case1.T(i))+Noise(1, :));sind(Venus.Ang(Case1.T(i))+Noise(1, :))];
    e2 = [cosd(Mars.Ang(Case1.T(i))+Noise(2, :));sind(Mars.Ang(Case1.T(i))+Noise(2, :))];
    e3 = [cosd(Moon.Ang(Case1.T(i))+Noise(3, :));sind(Moon.Ang(Case1.T(i))+Noise(3, :))];
    e4 = [cosd(Earth.Ang(Case1.T(i))+Noise(4, :));sind(Earth.Ang(Case1.T(i))+Noise(4, :))];
    for j = 1:Case1.N
        E = [e1(:, j) e2(:, j) e3(:, j) e4(:, j)];
        Case1.P_hat_LS(:, j, i) = GetPosition2D_LS(Planets, E);
        Case1.P_hat_WLS(:, j, i) = GetPosition2D_WLS(Planets, E, Case1.Observer(:, Case1.T(i)));
        Case1.P_hat_2P(:, j, i) = GetPosition2D_LS(Planets(:, 3:4), E(:, 3:4));
    end
    Case1.Error_LS(:, :, i) = Case1.P_hat_LS(:, :, i)-Case1.Observer(:, Case1.T(i));
    Case1.Error_WLS(:, :, i) = Case1.P_hat_WLS(:, :, i)-Case1.Observer(:, Case1.T(i));
    Case1.Error_2P(:, :, i) = Case1.P_hat_2P(:, :, i)-Case1.Observer(:, Case1.T(i));
    Case1.RMS_LS(1:2, i) = rms(Case1.Error_LS(:, :, i), 2);
    Case1.RMS_LS(3, i) = norm(Case1.RMS_LS(1:2, i));
    Case1.RMS_WLS(1:2, i) = rms(Case1.Error_WLS(:, :, i), 2);
    Case1.RMS_WLS(3, i) = norm(Case1.RMS_WLS(1:2, i));
    Case1.RMS_2P(1:2, i) = rms(Case1.Error_2P(:, :, i), 2);
    Case1.RMS_2P(3, i) = norm(Case1.RMS_2P(1:2, i));
    Case1.Ang(i) = acosd([1 0] * (Earth.E(:, Case1.T(i)).*Moon.E(:, Case1.T(i))));
end
%% Save the Results

save('S2_C1.mat', 'Case1');