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

%% Case #1
Distance = vecnorm(Observer, 2, 1);
T = [length(t)-length(find(Distance>(d/5)));length(t)-length(find(Distance>(2*d/5)));
    length(t)-length(find(Distance>(3*d/5)));length(t)-length(find(Distance>(4*d/5)))];

Case1.sigma = [1;1;1;1] * 10/3600;
Case1.N = 50000;
Case1.Observer = Observer;
Case1.t = t;
Case1.T = T;

Case1.P_hat_LS = zeros(2, Case1.N, length(T));
Case1.P_hat_WLS = zeros(2, Case1.N, length(T));
Case1.P_hat_2P = zeros(2, Case1.N, length(T));
Case1.Error_LS = zeros(2, Case1.N, length(T));
Case1.Error_WLS = zeros(2, Case1.N, length(T));
Case1.Error_2P = zeros(2, Case1.N, length(T));
Case1.RMS_LS = zeros(3, length(T));
Case1.RMS_WLS = zeros(3, length(T));
Case1.RMS_2P = zeros(3, length(T));
Case1.Ang = zeros(1, length(T));

for i = 1:length(T)
    Noise = randn(4, Case1.N) .* Case1.sigma;
    Planets = [Venus.P(:, T(i)) Mars.P(:, T(i)) Moon.P(:, T(i)) Earth.P(:, T(i))];
    e1 = [cosd(Venus.Ang(T(i))+Noise(1, :));sind(Venus.Ang(T(i))+Noise(1, :))];
    e2 = [cosd(Mars.Ang(T(i))+Noise(2, :));sind(Mars.Ang(T(i))+Noise(2, :))];
    e3 = [cosd(Moon.Ang(T(i))+Noise(3, :));sind(Moon.Ang(T(i))+Noise(3, :))];
    e4 = [cosd(Earth.Ang(T(i))+Noise(4, :));sind(Earth.Ang(T(i))+Noise(4, :))];
    for j = 1:Case1.N
        E = [e1(:, j) e2(:, j) e3(:, j) e4(:, j)];
        Case1.P_hat_LS(:, j, i) = GetPosition2D_LS(Planets, E);
        Case1.P_hat_WLS(:, j, i) = GetPosition2D_WLS(Planets, E, Case1.Observer(:, T(i)));
        Case1.P_hat_2P(:, j, i) = GetPosition2D_LS(Planets(:, 3:4), E(:, 3:4));
    end
    Case1.Error_LS(:, :, i) = Case1.P_hat_LS(:, :, i)-Case1.Observer(:, T(i));
    Case1.Error_WLS(:, :, i) = Case1.P_hat_WLS(:, :, i)-Case1.Observer(:, T(i));
    Case1.Error_2P(:, :, i) = Case1.P_hat_2P(:, :, i)-Case1.Observer(:, T(i));
    Case1.RMS_LS(1:2, i) = rms(Case1.Error_LS(:, :, i), 2);
    Case1.RMS_LS(3, i) = norm(Case1.RMS_LS(1:2, i));
    Case1.RMS_WLS(1:2, i) = rms(Case1.Error_WLS(:, :, i), 2);
    Case1.RMS_WLS(3, i) = norm(Case1.RMS_WLS(1:2, i));
    Case1.RMS_2P(1:2, i) = rms(Case1.Error_2P(:, :, i), 2);
    Case1.RMS_2P(3, i) = norm(Case1.RMS_2P(1:2, i));
    Case1.Ang(i) = acosd([1 0] * (Earth.E(:, T(i)).*Moon.E(:, T(i))));
end

%% Case #1 Results
for n = 1:length(Case1.T)
    figure();
    plot(Case1.Error_LS(1, :, n), Case1.Error_LS(2, :, n), '.');
    hold on;
    plot(Case1.Error_WLS(1, :, n), Case1.Error_WLS(2, :, n), '.');
    plot(Case1.Error_2P(1, :, n), Case1.Error_2P(2, :, n), '.');
    title(['Estimated Position Error (', num2str(n*0.2),'R)']);
    legend('LS', 'WLS', '2-Planet', 'Location', 'southeast');
    xlabel('X-direction(km)');
    ylabel('Y-direction(km)');
    axis(2.5e4*[-1 1 -1 1]);
    hold off;
end

for n = 1:length(Case1.T)
    figure();
    plot(Case1.Error_WLS(1, :, n), Case1.Error_WLS(2, :, n), '.');
    hold on;
    plot(Case1.Error_2P(1, :, n), Case1.Error_2P(2, :, n), '.');
    title(['Estimated Position Error (', num2str(n*0.2),'R)']);
    legend('WLS', '2-Planet', 'Location', 'southeast');
    xlabel('X-direction(km)');
    ylabel('Y-direction(km)');
    axis(100*[-1 1 -1 1]);
    hold off;
end

%% Save the Results
save('S2_C1.mat', 'Case1');