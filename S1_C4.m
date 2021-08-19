clear;
clc;

%% Case #4

Case4.sigma = [1;1;1];
Case4.Observer = [0;0];
Case4.theta = 0:10:90;
Case4.R = 100:10:1000;
Case4.N = 5000;

Case4.P_hat_LS = zeros(2, Case4.N, length(Case4.R), length(Case4.theta));
Case4.P_hat_WLS = zeros(2, Case4.N, length(Case4.R), length(Case4.theta));
Case4.Error_LS = zeros(2, Case4.N, length(Case4.R), length(Case4.theta));
Case4.Error_WLS = zeros(2, Case4.N, length(Case4.R), length(Case4.theta));
Case4.Counts_LS = zeros(length(Case4.R), length(Case4.theta));
Case4.Counts_WLS = zeros(length(Case4.R), length(Case4.theta));
Case4.Freq_LS = zeros(length(Case4.R), length(Case4.theta));
Case4.Freq_WLS = zeros(length(Case4.R), length(Case4.theta));

for i = 1:length(Case4.R)
    for j = 1:length(Case4.theta)
        Planets = [1 -1 Case4.R(i) * cosd(Case4.theta(j)); 0 0 Case4.R(i) * sind(Case4.theta(j))];
        Noise = randn(3, Case4.N) .* Case4.sigma;
        e1 = [cosd(Noise(1, :));sind(Noise(1, :))];
        e2 = -[cosd(Noise(2, :));sind(Noise(2, :))];
        e3 = [cosd(Case4.theta(j)+Noise(2, :));sind(Case4.theta(j)+Noise(2, :))];
        for k = 1:Case4.N
            EMatrix = [e1(:, k) e2(:, k) e3(:, k)];
            Case4.P_hat_LS(:, k, i, j) = GetPosition2D_LS(Planets, EMatrix);
            Case4.P_hat_WLS(:, k, i, j) = GetPosition2D_WLS(Planets, EMatrix, Case4.Observer);
        end
        Case4.Error_LS(:, :, i, j) = squeeze(Case4.P_hat_LS(:, :, i, j))-Case4.Observer;
        Case4.Error_WLS(:, :, i, j) = squeeze(Case4.P_hat_WLS(:, :, i, j))-Case4.Observer;
        Case4.Counts_LS(i, j) = sum(vecnorm(squeeze(Case4.Error_LS(:, :, i, j)))>1);
        Case4.Counts_WLS(i, j) = sum(vecnorm(squeeze(Case4.Error_WLS(:, :, i, j)))>1);
    end
end
Case4.Freq_LS = Case4.Counts_LS/Case4.N;
Case4.Freq_WLS = Case4.Counts_WLS/Case4.N;

%% Case #4 Results
figure();
plot(Case4.R, ones(1, length(Case4.R)) * Case2.Freq(length(Case2.Freq)), 'LineWidth', 2);
hold on;
for n = 1:length(Case4.theta)
    plot(Case4.R, Case4.Freq_LS(:, n));
end
title({'LS/2-Planet Comparision';'R vs Frequency of the Critical Error'}, 'Interpreter', 'latex');
xlabel('R', 'Interpreter', 'latex');
ylabel('Frequency of the Critical Error', 'Interpreter', 'latex');
axis([100 1000 0 1]);
legend({'2-Planet', '$\theta=0$', '$\theta=10$', '$\theta=20$', '$\theta=30$', '$\theta=40$', '$\theta=50$', '$\theta=60$', '$\theta=70$', '$\theta=80$', '$\theta=90$'}, 'Location', 'southeast', 'NumColumns', 2,'Interpreter', 'latex');
hold off;


figure();
plot(Case4.R, ones(1, length(Case4.R)) * Case2.Freq(length(Case2.Freq)), 'LineWidth', 2);
hold on;
for n = 1:length(Case4.theta)
    plot(Case4.R, Case4.Freq_WLS(:, n));
end
title({'WLS/2-Planet Comparision';'R vs Frequency of the Critical Error'}, 'Interpreter', 'latex');
xlabel('R', 'Interpreter', 'latex');
ylabel('Frequency of the Critical Error', 'Interpreter', 'latex');
axis([100 1000 0 1]);
legend({'2-Planet', '$\theta=0$', '$\theta=10$', '$\theta=20$', '$\theta=30$', '$\theta=40$', '$\theta=50$', '$\theta=60$', '$\theta=70$', '$\theta=80$', '$\theta=90$'}, 'Location', 'southeast', 'NumColumns', 2,'Interpreter', 'latex');
hold off;
%% Case #4 Visualization

figure();
plot(squeeze(Case2.P_hat(1, :, length(Case2.theta))), squeeze(Case2.P_hat(2, :, length(Case2.theta))), '.');
title('Estimated Position, 2-Planet Method');
axis([-40 40 -2 2]);

figure();
plot(squeeze(Case4.P_hat_LS(1, :, 51, 10)), squeeze(Case4.P_hat_LS(2, :, 51, 10)), '.');
title('Estimated Position, LS Method');
axis([-40 40 -2 2]);

figure();
plot(squeeze(Case4.P_hat_WLS(1, :, 51, 10)), squeeze(Case4.P_hat_WLS(2, :, 51, 10)), '.');
title('Estimated Position, WLS Method');
axis([-40 40 -2 2]);
%% Save the Results
save('S1_C4.mat', 'Case4');