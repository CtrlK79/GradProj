clear;
clc;

%% Load Data

load('S2_C1.mat');
load('S2_C2.mat');
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

%% Case #1 Visualization

x = (0:0.01:1) * 10^8;
for n = 1:length(Case1.T)
    figure();
    p1 = plot(0, 0, 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 15);
    hold on;
    p2 = plot(Moon.P(1, Case1.T(n)), Moon.P(2, Case1.T(n)), 'co', 'MarkerFaceColor', 'c', 'MarkerSize', 15);
    p3 = plot(Observer(1, Case1.T(n)), Observer(2, Case1.T(n)), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 15);
    
    y_Venus = (x-Observer(1, Case1.T(n)))*(Venus.P(2, Case1.T(n))-Observer(2, Case1.T(n)))/(Venus.P(1, Case1.T(n))-Observer(1, Case1.T(n))) + Observer(2, Case1.T(n));
    y_Mars = (x-Observer(1, Case1.T(n)))*(Mars.P(2, Case1.T(n))-Observer(2, Case1.T(n)))/(Mars.P(1, Case1.T(n))-Observer(1, Case1.T(n))) + Observer(2, Case1.T(n));
    
    y_Moon = (x-Observer(1, Case1.T(n)))*(Moon.P(2, Case1.T(n))-Observer(2, Case1.T(n)))/(Moon.P(1, Case1.T(n))-Observer(1, Case1.T(n))) + Observer(2, Case1.T(n));
    y_Earth = (x-Observer(1, Case1.T(n)))*(Earth.P(2, Case1.T(n))-Observer(2, Case1.T(n)))/(Earth.P(1, Case1.T(n))-Observer(1, Case1.T(n))) + Observer(2, Case1.T(n));
    
    
    
    p4 = plot(-x+Observer(1, Case1.T(n)), -y_Earth + Observer(2, Case1.T(n)), '--g', 'LineWidth', 2);
    p5 = plot(x+Observer(1, Case1.T(n)), y_Moon + Observer(1, Case1.T(n))*(Moon.P(2, Case1.T(n))-Observer(2, Case1.T(n)))/(Moon.P(1, Case1.T(n))-Observer(1, Case1.T(n))), '--c', 'LineWidth', 2);
    p6 = plot(-x + Observer(1, Case1.T(n)), 2 * Observer(2, Case1.T(n)) - (y_Venus + Observer(1, Case1.T(n))*(Venus.P(2, Case1.T(n))-Observer(2, Case1.T(n)))/(Venus.P(1, Case1.T(n))-Observer(1, Case1.T(n)))), '--m', 'LineWidth', 2);
    p7 = plot(-x + Observer(1, Case1.T(n)), 2 * Observer(2, Case1.T(n)) - (y_Mars + Observer(1, Case1.T(n))*(Mars.P(2, Case1.T(n))-Observer(2, Case1.T(n)))/(Mars.P(1, Case1.T(n))-Observer(1, Case1.T(n)))), '--r', 'LineWidth', 2);
    %title(['$Estimated\ Position\ Error\ ($', num2str(n*0.2),'$R)$'], 'Interpreter', 'latex');
    legend([p1 p2 p3 p6 p7], {'$Earth$', '$Moon$', '$Observer$', '$Direction\ of\ the\ Venus$', '$Direction\ of\ the\ Mars$'}, 'Location', 'southeast', 'Interpreter', 'latex');
    xlabel('$X-direction(km)$', 'Interpreter', 'latex');
    ylabel('$Y-direction(km)$', 'Interpreter', 'latex');
    axis(1.2 * 384400*[-0.5 1 -0.5 1]);
    hold off;
end

%% Case #1 Results

x = (0:0.01:1) * 10^8;
for n = 1:length(Case1.T)
    figure();
    p1 = plot(Case1.Error_LS(1, :, n), Case1.Error_LS(2, :, n), '.');
    hold on;
    p2 = plot(Case1.Error_WLS(1, :, n), Case1.Error_WLS(2, :, n), '.');
    p3 = plot(Case1.Error_2P(1, :, n), Case1.Error_2P(2, :, n), '.');
    
    y_Venus = -x*(Venus.P(2, Case1.T(n))-Observer(2, Case1.T(n)))/(Venus.P(1, Case1.T(n))-Observer(1, Case1.T(n)));
    y_Mars = -x*(Mars.P(2, Case1.T(n))-Observer(2, Case1.T(n)))/(Mars.P(1, Case1.T(n))-Observer(1, Case1.T(n)));
    
    y_Moon = x*(Moon.P(2, Case1.T(n))-Observer(2, Case1.T(n)))/(Moon.P(1, Case1.T(n))-Observer(1, Case1.T(n)));
    y_Earth = -x*(Earth.P(2, Case1.T(n))-Observer(2, Case1.T(n)))/(Earth.P(1, Case1.T(n))-Observer(1, Case1.T(n)));
    
    p4 = plot(-x, y_Earth, '--g', 'LineWidth', 2);
    p5 = plot(-x, y_Moon, '--c', 'LineWidth', 2);
    p6 = plot(x, y_Venus, '--m', 'LineWidth', 2);
    p7 = plot(-x, y_Mars, '--r', 'LineWidth', 2);
    title(['$Estimated\ Position\ Error\ ($', num2str(n*0.2),'$R)$'], 'Interpreter', 'latex');
    legend([p1 p2 p3 p4 p5 p6 p7], {'$LS$', '$WLS$', '$2-Planet$', '$Direction of the Earth$', '$Direction of the Moon$', '$Direction of the Venus$', '$Direction of the Mars$'}, 'Location', 'southeast', 'Interpreter', 'latex');
    xlabel('$X-direction(km)$', 'Interpreter', 'latex');
    ylabel('$Y-direction(km)$', 'Interpreter', 'latex');
    axis(2.5e4*[-1 1 -1 1]);
    hold off;
end

for n = 1:length(Case1.T)
    figure();
    p1 = plot(Case1.Error_WLS(1, :, n), Case1.Error_WLS(2, :, n), '.');
    hold on;
    p2 = plot(Case1.Error_2P(1, :, n), Case1.Error_2P(2, :, n), '.');
    
    y_Venus = -x*(Venus.P(2, Case1.T(n))-Observer(2, Case1.T(n)))/(Venus.P(1, Case1.T(n))-Observer(1, Case1.T(n)));
    y_Mars = -x*(Mars.P(2, Case1.T(n))-Observer(2, Case1.T(n)))/(Mars.P(1, Case1.T(n))-Observer(1, Case1.T(n)));
    
    y_Moon = x*(Moon.P(2, Case1.T(n))-Observer(2, Case1.T(n)))/(Moon.P(1, Case1.T(n))-Observer(1, Case1.T(n)));
    y_Earth = -x*(Earth.P(2, Case1.T(n))-Observer(2, Case1.T(n)))/(Earth.P(1, Case1.T(n))-Observer(1, Case1.T(n)));
    
    p3 = plot(-x, y_Earth, '--g', 'LineWidth', 2);
    p4 = plot(x, y_Moon, '--c', 'LineWidth', 2);
    p5 = plot(-x, y_Venus, '--m', 'LineWidth', 2);
    p6 = plot(-x, y_Mars, '--r', 'LineWidth', 2);
    title(['$Estimated\ Position\ Error\ ($', num2str(n*0.2),'$R)$'], 'Interpreter', 'latex');
    legend([p1 p2 p3 p4 p5 p6], {'$WLS$', '$2-Planet$', '$Direction of the Earth$', '$Direction of the Moon$', '$Direction of the Venus$', '$Direction of the Mars$'}, 'Location', 'southeast', 'Interpreter', 'latex');
    xlabel('$X-direction(km)$', 'Interpreter', 'latex');
    ylabel('$Y-direction(km)$', 'Interpreter', 'latex');
    axis(100*[-1 1 -1 1]);
    hold off;
end

%% Case #2 Visualization - 1

d = 384400;

figure();
p1 = plot(Case2.t, Case2.Ang, 'LineWidth', 2);
hold on;
p2 = plot(Case2.t, 10 * ones(1, length(Case2.t)), 'r', 'LineWidth', 2);
p3 = plot(Case2.t, 170 * ones(1, length(Case2.t)), 'r', 'LineWidth', 2);
xlabel('$time(sec)$', 'Interpreter', 'latex');
ylabel('$Earth-Observer-Moon\ Angle(^\circ)$', 'Interpreter', 'latex');
axis([1 length(t) 0 180]);
%title('$time(sec)\ vs\ Earth-Observer-Moon\ angle(^\circ)\ plot$', 'Interpreter', 'latex');
legend([p1 p2], {'$E-O-M\ Angle$', '$10^\circ,\ 170^\circ\ line$'}, 'Interpreter', 'latex');
hold off;

figure();
p1 = plot(Case2.P_hat_LS(1, :), Case2.P_hat_LS(2, :), 'b');
hold on;
xlabel('$X-direction(km)$', 'Interpreter', 'latex');
ylabel('$Y-direction(km)$', 'Interpreter', 'latex');
p2 = plot(Case2.Observer(1, :), Case2.Observer(2, :), 'r', 'Linewidth', 2);
title('$True\ Orbit\ vs\ LS\ Estimated\ Position$', 'Interpreter', 'latex');
axis(1.1*[-d d -d d]);
legend([p2 p1 p3], {'$True$', '$LS$'}, 'Interpreter', 'latex');
hold off;

figure();
p1 = plot(Case2.P_hat_WLS(1, :), Case2.P_hat_WLS(2, :), 'b');
hold on;
xlabel('$X-direction(km)$', 'Interpreter', 'latex');
ylabel('$Y-direction(km)$', 'Interpreter', 'latex');
p2 = plot(Case2.Observer(1, :), Case2.Observer(2, :), 'r', 'Linewidth', 2);
title('$True\ Orbit\ vs\ WLS\ Estimated\ Position$', 'Interpreter', 'latex');
axis(1.1*[-d d -d d]);
legend([p2 p1 p3], {'$True$', '$WLS$'}, 'Interpreter', 'latex');
hold off;

figure();
p1 = plot(Case2.P_hat_2P(1, :), Case2.P_hat_2P(2, :), 'b');
hold on;
xlabel('$X-direction(km)$', 'Interpreter', 'latex');
ylabel('$Y-direction(km)$', 'Interpreter', 'latex');
p2 = plot(Case2.Observer(1, :), Case2.Observer(2, :), 'r', 'Linewidth', 2);
title('$True\ Orbit\ vs\ 2-Planet\ Estimated\ Position$', 'Interpreter', 'latex');
axis(1.1*[-d d -d d]);
legend([p2 p1 p3], {'$True$', '$2-Planet$'}, 'Interpreter', 'latex');
hold off;

%% Case #2 Visualization - 2

h = figure();

axis_max = 1.05 * d;
t_1000 = 1:1000:length(t);

for i = 1:length(t_1000)
   hold off
   p0 = plot(0, 0, 'bo', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'MarkerSize', 10);
   
   xlim(axis_max * [-1, 1]);
   ylim(axis_max * [-1, 1]);
   zlim(axis_max * [-1, 1]);
   xlabel('$X-direction(km)$', 'Interpreter', 'latex');
   ylabel('$Y-direction(km)$', 'Interpreter', 'latex');
   hold on;
   p1 = plot(Moon.P(1,t_1000(i)), Moon.P(2,t_1000(i)), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 10);
   p2 = plot(Case2.Observer(1,t_1000(i)), Case2.Observer(2,t_1000(i)), 'o', 'Color', [0 1 0], 'MarkerFaceColor', [0 1 0]);
   
   p3 = plot(Moon.P(1,1:1000:t_1000(i)), Moon.P(2,1:1000:t_1000(i)), '-r', 'linewidth', 2);
   p4 = plot(Case2.Observer(1,1:1000:t_1000(i)), Case2.Observer(2,1:1000:t_1000(i)), '-', 'Color', [0 1 0], 'linewidth', 2);
   
   x = 0:0.05:1 * axis_max * 2;
   x_Venus = -x+Observer(1, t_1000(i));
   y_Venus = (x_Venus-Observer(1, t_1000(i)))*(Venus.P(2, t_1000(i))-Observer(2, t_1000(i)))/(Venus.P(1, t_1000(i))-Observer(1, t_1000(i))) + Observer(2, t_1000(i));
   x_Mars = -x+Observer(1, t_1000(i));
   y_Mars = (x_Mars-Observer(1, t_1000(i)))*(Mars.P(2, t_1000(i))-Observer(2, t_1000(i)))/(Mars.P(1, t_1000(i))-Observer(1, t_1000(i))) + Observer(2, t_1000(i));
   
   p5 = plot(x_Venus, y_Venus, '--m', 'LineWidth', 2);
   p6 = plot(x_Mars, y_Mars, '--r', 'LineWidth', 2);
   
   p3.Annotation.LegendInformation.IconDisplayStyle = 'off';
   p4.Annotation.LegendInformation.IconDisplayStyle = 'off';
   
   legend([p0 p1 p2 p5 p6], {'$Earth$', '$Moon$', '$Observer$', '$Direction\ of\ the\ Venus$', '$Direction\ of\ the\ Mars$'}, 'location', 'southwest', 'Interpreter', 'latex');
   
   frame = getframe;
   im = frame2im(frame);
   [A,map] = rgb2ind(im,256); 

   if (i == 1)       
       imwrite(A,map,'Orbit.gif','gif','LoopCount',Inf,'DelayTime',1);
   else       
       imwrite(A,map,'Orbit.gif','gif','WriteMode','append','DelayTime', 0);
   end
end