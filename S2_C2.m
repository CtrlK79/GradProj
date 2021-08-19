clear;
clc;

%% Time Setting
sec2day = 1/(3600 * 24);                    % 1 일 = 3600 * 24 초 

t_ref = 2451545.0;                          % J2000 기준점
t_initial = 365.25 * 21;                    % 출발 시간
t_step = 1;                                 % 측정치 생성 간격(초)
t_final = 405890;                           % Second-kick 이전까지 시간(초)
t = 1:t_step:t_final;                       % 전체 시간 구간(초)

%% Orbit Preprocessing
d = 384400;
mu_E = 3.98595624*10^5;                     % km^3/sec^2
mu_M = 4.90278576*10^3;                     % km^3/sec^2
w = sqrt((mu_E+mu_M)/d^3);                  % rad/sec

theta = rad2deg((t-1) * w);                 % 회전좌표계 변환각(deg)

Observer = load("Orbit.mat");
Observer = Observer.Orbit;
Observer = Observer(:, t);                  % 측정치 생성 시간에서 위치      

for i = 1:length(t)
    Observer(:, i) = RotMx(4, theta(i)) * Observer(:, i);
end

%% Orbit Transformation
% axis generate
[Moon0_p, Moon0_v] = planetEphemeris([t_ref t_initial], 'Earth', 'Moon', '405');
X = ToUnit(Moon0_p');
Z = ToUnit(cross(X, Moon0_v'));
Y = cross(Z, X);

% get Planets position
VenusEph = planetEphemeris(t'*sec2day+t_initial+t_ref, 'Earth', 'Venus', '405');
MarsEph = planetEphemeris(t'*sec2day+t_initial+t_ref, 'Earth', 'Mars', '405');

% 2D projection
Venus.P = [VenusEph*X VenusEph*Y]';
Mars.P = [MarsEph*X MarsEph*Y]';
Moon.P = d * [cosd(theta);sind(theta)];
Earth.P = zeros(2, length(theta));
%% Data Preprocessing
clear VenusEph MarsEph Moon0_p Moon0_v

% Results : oberser->planets 방향벡터(true 값)
Venus.E = (Venus.P-Observer)./vecnorm((Venus.P-Observer), 2, 1);
Mars.E = (Mars.P-Observer)./vecnorm((Mars.P-Observer), 2, 1);
Moon.E = (Moon.P-Observer)./vecnorm((Moon.P-Observer), 2, 1);
Earth.E = (Earth.P-Observer)./vecnorm(Earth.P-Observer, 2, 1);

Venus.Ang = acosd(Venus.E(1, :)) .* (2*(Venus.E(2, :)>0)-1);
Mars.Ang = acosd(Mars.E(1, :)) .* (2*(Mars.E(2, :)>0)-1);
Moon.Ang = acosd(Moon.E(1, :)) .* (2*(Moon.E(2, :)>0)-1);
Earth.Ang = acosd(Earth.E(1, :)) .* (2*(Earth.E(2, :)>0)-1);

%% Case #2
Case2.sigma = [1;1;1;1] * 10/3600;
Case2.Observer = Observer;
Case2.t = t;

Case2.P_hat_LS = zeros(2, length(t));
Case2.P_hat_WLS = zeros(2, length(t));
Case2.P_hat_2P = zeros(2, length(t));
Case2.Error_LS = zeros(2, length(t));
Case2.Error_WLS = zeros(2, length(t));
Case2.Error_2P = zeros(2, length(t));
Case2.RMS_LS = zeros(3, 1);
Case2.RMS_WLS = zeros(3, 1);
Case2.RMS_2P = zeros(3, 1);
Case2.Ang = abs(abs(Earth.Ang-Moon.Ang) - (360 * (abs(Earth.Ang-Moon.Ang)>180)));

Noise = randn(4, length(t)) .* Case2.sigma;
e1 = [cosd(Moon.Ang+Noise(1, :));sind(Moon.Ang+Noise(1, :))];
e2 = [cosd(Earth.Ang+Noise(2, :));sind(Earth.Ang+Noise(2, :))];
e3 = [cosd(Venus.Ang+Noise(3, :));sind(Venus.Ang+Noise(3, :))];
e4 = [cosd(Mars.Ang+Noise(4, :));sind(Mars.Ang+Noise(4, :))];
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
Case2.RMS_2P(1:2, 1) = rms(Case2.Error_LS, 2);
Case2.RMS_2P(3, 1) = norm(Case2.RMS_LS(1:2, 1));

%% Case #2 Results
d = 384400;
figure();
p1 = plot(Case2.t, Case2.Ang);
hold on;
p2 = plot(Case2.t, 10 * ones(1, length(Case2.t)), 'r');
p3 = plot(Case2.t, 170 * ones(1, length(Case2.t)), 'r');
xlabel('time(sec)');
ylabel(['Earth-Observer-Moon Angle(', '$^\circ$', ')'], 'Interpreter', 'latex');
axis([1 length(t)*1.4 0 180]);
title('time(sec) vs Earth-Observer-Moon angle(^\circ) plot');
legend([p1 p2], {'E-O-M Angle', '10^\circ, 170^\circ line'});
hold off;

figure();
p1 = plot(Case2.P_hat_LS(1, :), Case2.P_hat_LS(2, :), 'y');
hold on;

p2 = plot(Case2.Observer(1, :), Case2.Observer(2, :), 'r', 'Linewidth', 2);
title('True Orbit vs LS Estimated Position', 'Color', 'w');
set(gca, 'color', [0 0 0]);
axis off
set(gcf, 'color', [0 0 0]);
axis(1.1*[-0.3*d d -0.3*d d]);
legend([p2 p1 p3], {'True', 'LS'}, 'TextColor', 'w');
hold off;

figure();
p1 = plot(Case2.P_hat_WLS(1, :), Case2.P_hat_WLS(2, :), 'y');
hold on;

p2 = plot(Case2.Observer(1, :), Case2.Observer(2, :), 'r', 'Linewidth', 2);
title('True Orbit vs WLS Estimated Position', 'Color', 'w');
set(gca, 'color', [0 0 0]);
axis off
set(gcf, 'color', [0 0 0]);
axis(1.1*[-0.3*d d -0.3*d d]);
legend([p2 p1 p3], {'True', 'WLS'}, 'TextColor', 'w');
hold off;

figure();

p1 = plot(Case2.P_hat_2P(1, :), Case2.P_hat_2P(2, :), 'y');
hold on;

p2 = plot(Case2.Observer(1, :), Case2.Observer(2, :), 'r', 'Linewidth', 2);
title('True Orbit vs 2-Planet Estimated Position', 'Color', 'w');
set(gca, 'color', [0 0 0]);
axis off
set(gcf, 'color', [0 0 0]);
axis(1.1*[-0.3*d d -0.3*d d]);
legend([p2 p1 p3], {'True', '2-Planet'}, 'TextColor', 'w');
hold off;

%% Case #2 Visualization
h = figure();

axis_max = 1.05 * d;
t_1000 = 1:1000:length(t);

for i = 1:length(t_1000)
   hold off
   p0 = plot(0, 0, 'bo', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
   
   xlim(axis_max * [-1, 1]);
   ylim(axis_max * [-1, 1]);
   zlim(axis_max * [-1, 1]);
   hold on;
   p1 = plot(Moon.P(1,t_1000(i)), Moon.P(2,t_1000(i)), 'wo');
   p2 = plot(Case2.Observer(1,t_1000(i)), Case2.Observer(2,t_1000(i)), 'ro');

   legend([p0 p1 p2], {'Earth', 'Moon', 'Observer'}, 'TextColor', 'w', 'location', 'northwest', 'Color', 'k');
   
   p3 = plot(Moon.P(1,1:1000:t_1000(i)), Moon.P(2,1:1000:t_1000(i)), '-w', 'linewidth', 2);
   p4 = plot(Case2.Observer(1,1:1000:t_1000(i)), Case2.Observer(2,1:1000:t_1000(i)), '-r', 'linewidth', 2);
   
   set(gca, 'color', [0 0 0]);
   axis off
   set(gcf, 'color', [0 0 0]);
   
   p3.Annotation.LegendInformation.IconDisplayStyle = 'off';
   p4.Annotation.LegendInformation.IconDisplayStyle = 'off';
   
   frame = getframe;
   im = frame2im(frame);
   [A,map] = rgb2ind(im,256); 

   if (i == 1)       
       imwrite(A,map,'Orbit.gif','gif','LoopCount',Inf,'DelayTime',1);
   else       
       imwrite(A,map,'Orbit.gif','gif','WriteMode','append','DelayTime', 0);
   end
end
%% Save the Results
save('S2_C2.mat', 'Case2');