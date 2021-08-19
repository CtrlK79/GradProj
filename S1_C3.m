clear;
clc;

%% Case #3

Case3.sigma = [1;1;1];
Case3.Observer = [0;0];
Case3.theta = 45:225;
Case3.N = 5000;
Case3.R = [1 10 100 1000];

Case3.P_hat_LS = zeros(2, Case3.N, length(Case3.R), length(Case3.theta));
Case3.P_hat_WLS = zeros(2, Case3.N, length(Case3.R), length(Case3.theta));
Case3.Error_LS = zeros(2, Case3.N, length(Case3.R), length(Case3.theta));
Case3.Error_WLS = zeros(2, Case3.N, length(Case3.R), length(Case3.theta));
Case3.RMS_LS = zeros(3, length(Case3.R), length(Case3.theta));
Case3.RMS_WLS = zeros(3, length(Case3.R), length(Case3.theta));

for i = 1:length(Case3.R)
    for j = 1:length(Case3.theta)
        Planets = [1 0 Case3.R(i)*cosd(Case3.theta(j)); 0 1 Case3.R(i)*sind(Case3.theta(j))];
        Noise = randn(3, Case3.N) .* Case3.sigma;
        e1 = [cosd(Noise(1, :));sind(Noise(1, :))];
        e2 = [cosd(Noise(2, :)+90);sind(Noise(2, :)+90)];
        e3 = [cosd(Case3.theta(j)+Noise(3, :));sind(Case3.theta(j)+Noise(3, :))];
        for k = 1:Case3.N
            EMatrix = [e1(:, k) e2(:, k) e3(:, k)];
            Case3.P_hat_LS(:, k, i, j) = GetPosition2D_LS(Planets, EMatrix);
            Case3.P_hat_WLS(:, k, i, j) = GetPosition2D_WLS(Planets, EMatrix, Case3.Observer);
        end
        Case3.Error_LS(:, :, i, j) = squeeze(Case3.P_hat_LS(:, :, i, j))-Case3.Observer;
        Case3.Error_WLS(:, :, i, j) = squeeze(Case3.P_hat_WLS(:, :, i, j))-Case3.Observer;
        Case3.RMS_LS(:, i, j) = [rms(squeeze(Case3.Error_LS(:, :, i, j)), 2);norm(rms(squeeze(Case3.Error_LS(:, :, i, j)), 2))];
        Case3.RMS_WLS(:, i, j) = [rms(squeeze(Case3.Error_WLS(:, :, i, j)), 2);norm(rms(squeeze(Case3.Error_WLS(:, :, i, j)), 2))];
    end
end


%% Case #3 Results

for n = 1:length(Case3.R)
        figure();
        plot(Case3.theta, squeeze(Case3.RMS_LS(3, n, :)));
        hold on;
        plot(Case3.theta, squeeze(Case3.RMS_WLS(3, n, :)));
        plot(Case3.theta, ones(1, length(Case3.theta)).*squeeze(Case1.RMS(3, 1, 90)));
        legend('RMSE_{tot, LS}', 'RMSE_{tot, WLS}', 'RMSE_{tot, 2-Planet}');
        title(['LS/WLS/2-Planet Comparision, R = ',num2str(Case3.R(n)), ',$\ \theta(^\circ)\ vs\ RMSE\ plot$'], 'Interpreter', 'latex');
        xlabel('$\theta(^\circ)$', 'Interpreter', 'latex');
        ylabel('RMSE', 'Interpreter', 'latex');
        xlim([45 225]);
        ylim([0 1.5 * max(Case3.RMS_LS(3, n, :))]);
        hold off;
end

for n = 1:length(Case3.R)
        figure();
        plot(Case3.theta, squeeze(Case3.RMS_WLS(3, n, :)));
        hold on;
        plot(Case3.theta, ones(1, length(Case3.theta)).*squeeze(Case1.RMS(3, 1, 90)));
        legend('RMSE_{tot, WLS}', 'RMSE_{tot, 2-Planet}');
        title(['WLS/2-Planet Comparision, R = ',num2str(Case3.R(n)), ',$\ \theta(^\circ)\ vs\ RMSE\ plot$'], 'Interpreter', 'latex');
        xlabel('$\theta(^\circ)$', 'Interpreter', 'latex');
        ylabel('RMSE', 'Interpreter', 'latex');
        xlim([45 225]);
        ylim([0 1.5 * max(Case3.RMS_WLS(3, n, :))]);
        hold off;
end
%% Case #3 Visualization

angle_n = 91;
R_n = 2;

figure();

p1 = plot(squeeze(Case3.P_hat_LS(1, :, R_n, angle_n)), squeeze(Case3.P_hat_LS(2, :, R_n, angle_n )), 'b.');
hold on;
p2 = plot(squeeze(Case3.P_hat_WLS(1, :, R_n, angle_n)), squeeze(Case3.P_hat_WLS(2, :, R_n, angle_n )), 'm.');
p3 = plot(squeeze(Case1.P_hat(1, :, 1, angle_n)), squeeze(Case1.P_hat(2, :, 1, angle_n )), 'k.');
p4 = plot(Case3.Observer(1), Case3.Observer(2), 'ro', 'MarkerSize', 10, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');
p5 = plot(1, 0, 'go', 'MarkerSize', 10, 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g');
p6 = plot(0, 1, 'go', 'MarkerSize', 10, 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g');
p7 = plot(Case3.R(R_n)*cosd(Case3.theta(angle_n)), Case3.R(R_n)*sind(Case3.theta(angle_n)), 'go', 'MarkerSize', 10, 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g');
x = Case3.R(R_n)*(-1.1:0.01:1.1);
y1 = (x-1) * tand(3);
y2 = (x-1) * tand(-3);
y3 = x * tand(90+3)+1;
y4 = x * tand(90-3)+1;
y5 = (x-Case3.R(R_n)*cosd(Case3.theta(angle_n))) * tand(Case3.theta(angle_n)+3) + Case3.R(R_n) * sind(Case3.theta(angle_n));
y6 = (x-Case3.R(R_n)*cosd(Case3.theta(angle_n))) * tand(Case3.theta(angle_n)-3) + Case3.R(R_n) * sind(Case3.theta(angle_n));

p8 = plot(x, y1, '-g');
p9 = plot(x, y2, '-g');
p10 = plot(x, y3, '-g');
p11 = plot(x, y4, '-g');
p12 = plot(x, y5, '-c');
p13 = plot(x, y6, '-c');


axis([-0.3 0.3 -0.3 0.3]);
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
hold off;
%legend([p1 p2 p3 p4 p5 p12], {'Estimated Position(LS)', 'Estimated Position(WLS)', 'Estimated Position(2-Planets)', 'Observer', 'Planets', '3-sigma line of the 3rd Planet'});
title(['LS, R = ', num2str(squeeze(Case3.R(R_n))), '$,\ \theta^{(\circ)} = $', num2str(squeeze(Case3.theta(angle_n)))], 'Interpreter', 'latex');

%% Save the Results
save('S1_C3.mat', 'Case3');