clear;
clc;

%% Case #2

Case2.sigma = [1;1];
Case2.Observer = [0;0];
Case2.theta = [0.01:0.01:10 170:0.01:180];
Case2.N = 5000;

Case2.P_hat = zeros(2, Case2.N, length(Case2.theta));
Case2.Error = zeros(2, Case2.N, length(Case2.theta));
Case2.Counts = zeros(1, length(Case2.theta));
Case2.Freq = zeros(1, length(Case2.theta));

for j = 1:length(Case2.theta)
    Planets = [1 cosd(Case2.theta(j)); 0 sind(Case2.theta(j))];
    Noise = randn(2, Case2.N) .* Case2.sigma;
    e1 = [cosd(Noise(1, :));sind(Noise(1, :))];
    e2 = [cosd(Case2.theta(j)+Noise(2, :));sind(Case2.theta(j)+Noise(2, :))];
    for k = 1:Case2.N
        EMatrix = [e1(:, k) e2(:, k)];
        Case2.P_hat(:, k, j) = GetPosition2D_LS(Planets, EMatrix);
    end
    Case2.Error(:, :, j) = squeeze(Case2.P_hat(:, :, j))-Case2.Observer;
    Case2.Counts(j) = sum(vecnorm(squeeze(Case2.Error(:, :, j)))>1);
end
Case2.Freq = Case2.Counts/Case2.N;

%% Case #2 Results
n = sum(Case2.theta<90);
subplot(1, 2, 1);
plot(Case2.theta(1:n), Case2.Freq(1:n));
title('$\ \theta(^\circ)\ vs\ Frequency\ of\ the\ Critical\ Error\ plot$', 'Interpreter', 'latex');
xlabel('$\theta(^\circ)$', 'Interpreter', 'latex');
ylabel('Frequency of the Critical Error', 'Interpreter', 'latex');
axis([0 Case2.theta(n) 0 1]);
hold off;

subplot(1, 2, 2);
plot(Case2.theta(n+1:length(Case2.theta)), Case2.Freq(n+1:length(Case2.theta)));
title('$\ \theta(^\circ)\ vs\ Frequency\ of\ the\ Critical\ Error\ plot$', 'Interpreter', 'latex');
xlabel('$\theta(^\circ)$', 'Interpreter', 'latex');
ylabel('Frequency of the Critical Error', 'Interpreter', 'latex');
axis([Case2.theta(n+1) Case2.theta(length(Case2.theta)) 0 1]);
hold off;

Freq_small = Case2.Freq.*(Case2.theta<90);
Freq_large = Case2.Freq.*(Case2.theta>90);

Case2.CritS = 0.01 * (1+[sum(Freq_small>0.1) sum(Freq_small>0.01) sum(Freq_small>0.001) sum(Freq_small>0.0001)]);
Case2.CritL = 180 - 0.01 * ([sum(Freq_large>0.1) sum(Freq_large>0.01) sum(Freq_large>0.001) sum(Freq_small>0.0001)]);
%% Case #2 Visualization

angle_n = 180;
R_n = 1;

figure();

p1 = plot(squeeze(Case2.P_hat(1, :, R_n, angle_n)), squeeze(Case2.P_hat(2, :, R_n, angle_n )), '.');
hold on;
p2 = plot(Case2.Observer(1), Case2.Observer(2), 'ro', 'MarkerSize', 10, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');
p3 = plot(1, 0, 'go', 'MarkerSize', 10, 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g');
p4 = plot(Case2.R(R_n)*cosd(Case2.theta(angle_n)), Case2.R(R_n)*sind(Case2.theta(angle_n)), 'go', 'MarkerSize', 10, 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g');
x = Case2.R(R_n)*(-1.1:0.01:1.1);
y1 = (x-Case2.R(R_n)) * tand(3);
y2 = (x-Case2.R(R_n)) * tand(-3);
y3 = (x-Case2.R(R_n)*cosd(Case2.theta(angle_n))) * tand(Case2.theta(angle_n)+3)+Case2.R(R_n)*sind(Case2.theta(angle_n));
y4 = (x-Case2.R(R_n)*cosd(Case2.theta(angle_n))) * tand(Case2.theta(angle_n)-3)+Case2.R(R_n)*sind(Case2.theta(angle_n));

p5 = plot(x, y1, '-g');
p6 = plot(x, y2, '-g');
p7 = plot(x, y3, '-g');
p8 = plot(x, y4, '-g');

axis(Case2.R(R_n)*[-5 5 -1 1]);
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
hold off;
legend([p1 p2 p3], {'Estimated Position','Observer', 'Planets'});
title(['R = ', num2str(squeeze(Case2.R(R_n))), '$,\ \theta^{(\circ)} = $', num2str(squeeze(Case2.theta(angle_n)))], 'Interpreter', 'latex');

%% Save the Results
save('S1_C2.mat', 'Case2');