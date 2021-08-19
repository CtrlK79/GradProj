clear;
clc;

%% Case #1

Case1.sigma = [1;1];
Case1.Observer = [0;0];
Case1.theta = 1:180;
Case1.N = 5000;
Case1.R = [1 10 100 1000];

Case1.P_hat = zeros(2, Case1.N, length(Case1.R), length(Case1.theta));
Case1.Error = zeros(2, Case1.N, length(Case1.R), length(Case1.theta));
Case1.RMS = zeros(3, length(Case1.R), length(Case1.theta));

for i = 1:length(Case1.R)
    for j = 1:length(Case1.theta)
        Planets = [1 Case1.R(i)*cosd(Case1.theta(j)); 0 Case1.R(i)*sind(Case1.theta(j))];
        Noise = randn(2, Case1.N) .* Case1.sigma;
        e1 = [cosd(Noise(1, :));sind(Noise(1, :))];
        e2 = [cosd(Case1.theta(j)+Noise(2, :));sind(Case1.theta(j)+Noise(2, :))];
        for k = 1:Case1.N
            EMatrix = [e1(:, k) e2(:, k)];
            Case1.P_hat(:, k, i, j) = GetPosition2D_LS(Planets, EMatrix);
        end
        Case1.Error(:, :, i, j) = squeeze(Case1.P_hat(:, :, i, j))-Case1.Observer;
        Case1.RMS(1:2, i, j) = rms(squeeze(Case1.Error(:, :, i, j)), 2);
        Case1.RMS(3, i, j) = norm(squeeze(Case1.RMS(1:2, i, j)));
    end
end


%% Case #1 Results
for n = 1:length(Case1.R)
    figure();
    plot(Case1.theta, squeeze(Case1.RMS(1, n, :)));
    hold on;
    plot(Case1.theta, squeeze(Case1.RMS(2, n, :)));
    plot(Case1.theta, squeeze(Case1.RMS(3, n, :)));
    legend('RMSE_X', 'RMSE_Y', 'RMSE_{tot}');
    title(['R = ',num2str(Case1.R(n)), ',$\ \theta(^\circ)\ vs\ RMSE\ plot$'], 'Interpreter', 'latex');
    xlabel('$\theta(^\circ)$', 'Interpreter', 'latex');
    ylabel('RMSE', 'Interpreter', 'latex');
    axis([0 180 0 Case1.R(n)]);
    hold off;
end
%% Case #1 Visualization

angle_n = 150;
R_n = 4;

figure();

p1 = plot(squeeze(Case1.P_hat(1, :, R_n, angle_n)), squeeze(Case1.P_hat(2, :, R_n, angle_n )), '.');
hold on;
p2 = plot(Case1.Observer(1), Case1.Observer(2), 'ro', 'MarkerSize', 10, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');
p3 = plot(1, 0, 'go', 'MarkerSize', 10, 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g');
p4 = plot(Case1.R(R_n)*cosd(Case1.theta(angle_n)), Case1.R(R_n)*sind(Case1.theta(angle_n)), 'go', 'MarkerSize', 10, 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g');
x = Case1.R(R_n)*(-1.1:0.01:1.1);
y1 = (x-1) * tand(3);
y2 = (x-1) * tand(-3);
y3 = (x-Case1.R(R_n)*cosd(Case1.theta(angle_n))) * tand(Case1.theta(angle_n)+3)+Case1.R(R_n)*sind(Case1.theta(angle_n));
y4 = (x-Case1.R(R_n)*cosd(Case1.theta(angle_n))) * tand(Case1.theta(angle_n)-3)+Case1.R(R_n)*sind(Case1.theta(angle_n));

p5 = plot(x, y1, '-g');
p6 = plot(x, y2, '-g');
p7 = plot(x, y3, '-g');
p8 = plot(x, y4, '-g');

axis(1.5 * Case1.R(R_n)*[-1 1 -1 1]);
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
hold off;
legend([p1 p2 p3], {'Estimated Position','Observer', 'Planets'});
title(['R = ', num2str(squeeze(Case1.R(R_n))), '$,\ \theta^{(\circ)} = $', num2str(squeeze(Case1.theta(angle_n)))], 'Interpreter', 'latex');

%% Save the Results
save('S1_C1.mat', 'Case1');