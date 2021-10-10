clear;
clc;

%% Load Data

load('S1_C1_1.mat');
load('S1_C1.mat');
load('S1_C2.mat');
load('S1_C3.mat');
load('S1_C4.mat');

%% Case #1 Results

for n = 1:length(Case1.R)
    figure();
    plot(Case1.theta, squeeze(Case1.RMS(1, n, :)), 'LineWidth', 2);
    hold on;
    plot(Case1.theta, squeeze(Case1.RMS(2, n, :)), 'LineWidth', 2);
    plot(Case1.theta, squeeze(Case1.RMS(3, n, :)), 'LineWidth', 2);
    legend('$RMSE_X$', '$RMSE_Y$', '$RMSE_{tot}$', 'Interpreter', 'latex');
    %title(['$R = $', num2str(Case1.R(n)), '$,\ \theta(^\circ)\ vs\ RMSE\ plot$'], 'Interpreter', 'latex');
    xlabel('$\theta(^\circ)$', 'Interpreter', 'latex');
    ylabel('$RMSE$', 'Interpreter', 'latex');
    axis([0 180 0 Case1.R(n)]);
    hold off;
end

%% Case #1 Visualization

% Visuals change as setting these variables
angle_n = 180;
R_n = 1;

% Figure generation
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

axis(1.2 * Case1.R(R_n)*[-10 10 -1 1]);
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
hold off;
legend([p1 p2 p3 p5], {'$Estimated Position$','$Observer$', '$Planets$', '$3-\sigma\ lines$'}, 'Interpreter', 'latex');
title(['$R = $', num2str(squeeze(Case1.R(R_n))), '$,\ \theta^{(\circ)} = $', num2str(squeeze(Case1.theta(angle_n)))], 'Interpreter', 'latex');

%% Case #1_1 Results

figure();
plot(Case1_1.R, squeeze(Case1_1.RMS(3, :)), 'LineWidth', 1.5);
%title('$\theta = 90^\circ,\ R\ vs\ RMSE\ plot$', 'Interpreter', 'latex');
xlabel('$R$', 'Interpreter', 'latex');
ylabel('$RMSE$', 'Interpreter', 'latex');
hold off;

%% Case #2 Results

n = sum(Case2.theta<90);
subplot(1, 2, 1);
plot(Case2.theta(1:n), Case2.Freq(1:n), 'LineWidth', 1.5);
%title('$\theta(^\circ)\ vs\ Frequency\ of\ the\ Critical\ Error\ plot$', 'Interpreter', 'latex');
xlabel('$\theta(^\circ)$', 'Interpreter', 'latex');
ylabel('$Frequency\ of\ the\ Critical\ Error$', 'Interpreter', 'latex');
axis([0 Case2.theta(n) 0 1]);
hold off;

subplot(1, 2, 2);
plot(Case2.theta(n+1:length(Case2.theta)), Case2.Freq(n+1:length(Case2.theta)), 'LineWidth', 1.5);
%title('$\theta(^\circ)\ vs\ Frequency\ of\ the\ Critical\ Error\ plot$', 'Interpreter', 'latex');
xlabel('$\theta(^\circ)$', 'Interpreter', 'latex');
ylabel('$Frequency\ of\ the\ Critical\ Error$', 'Interpreter', 'latex');
axis([Case2.theta(n+1) Case2.theta(length(Case2.theta)) 0 1]);
hold off;

%% Case #3 Results

for n = 1:length(Case3.R)
        figure();
        plot(Case3.theta, squeeze(Case3.RMS_LS(3, n, :)), 'LineWidth', 1.5);
        hold on;
        plot(Case3.theta, squeeze(Case3.RMS_WLS(3, n, :)), 'LineWidth', 1.5);
        plot(Case3.theta, ones(1, length(Case3.theta)).*squeeze(Case1.RMS(3, 1, 90)), 'LineWidth', 1.5);
        legend('$RMSE_{tot, LS}$', '$RMSE_{tot, WLS}$', '$RMSE_{tot, 2-Planet}$', 'Interpreter', 'latex');
        %title(['$LS/WLS/2-Planet\ Comparision,\ R = $', num2str(Case3.R(n)), ',$\ \theta(^\circ)\ vs\ RMSE\ plot$'], 'Interpreter', 'latex');
        xlabel('$\theta(^\circ)$', 'Interpreter', 'latex');
        ylabel('$RMSE$', 'Interpreter', 'latex');
        xlim([45 225]);
        ylim([0 1.5 * max(Case3.RMS_LS(3, n, :))]);
        hold off;
end

%for n = 1:length(Case3.R)
%        figure();
%        plot(Case3.theta, squeeze(Case3.RMS_WLS(3, n, :)));
%        hold on;
%        plot(Case3.theta, ones(1, length(Case3.theta)).*squeeze(Case1.RMS(3, 1, 90)));
%        legend('$RMSE_{tot, WLS}$', '$RMSE_{tot, 2-Planet}$', 'Interpreter', 'latex');
%        title(['$WLS/2-Planet\ Comparision,\ R = $', num2str(Case3.R(n)), ',$\ \theta(^\circ)\ vs\ RMSE\ plot$'], 'Interpreter', 'latex');
%        xlabel('$\theta(^\circ)$', 'Interpreter', 'latex');
%        ylabel('$RMSE$', 'Interpreter', 'latex');
%        xlim([45 225]);
%        ylim([0 1.5 * max(Case3.RMS_WLS(3, n, :))]);
%        hold off;
%end

%subplot(2, 1, 1);
%n=1;
%plot(Case3.theta, squeeze(Case3.RMS_LS(3, n, :)));
%hold on;
%plot(Case3.theta, squeeze(Case3.RMS_WLS(3, n, :)));
%plot(Case3.theta, ones(1, length(Case3.theta)).*squeeze(Case1.RMS(3, 1, 90)));
%legend('$RMSE_{tot, LS}$', '$RMSE_{tot, WLS}$', '$RMSE_{tot, 2-Planet}$', 'Interpreter', 'latex');
%title(['$LS/WLS/2-Planet\ Comparision,\ R = $', num2str(Case3.R(n)), ',$\ \theta(^\circ)\ vs\ RMSE\ plot$'], 'Interpreter', 'latex');
%xlabel('$\theta(^\circ)$', 'Interpreter', 'latex');
%ylabel('$RMSE$', 'Interpreter', 'latex');
%xlim([45 225]);
%ylim([0 1.5 * max(Case3.RMS_LS(3, n, :))]);

%subplot(2, 1, 2);
%n=2;
%plot(Case3.theta, squeeze(Case3.RMS_LS(3, n, :)));
%hold on;
%plot(Case3.theta, squeeze(Case3.RMS_WLS(3, n, :)));
%plot(Case3.theta, ones(1, length(Case3.theta)).*squeeze(Case1.RMS(3, 1, 90)));
%legend('$RMSE_{tot, LS}$', '$RMSE_{tot, WLS}$', '$RMSE_{tot, 2-Planet}$', 'Interpreter', 'latex');
%title(['$LS/WLS/2-Planet\ Comparision,\ R = $', num2str(Case3.R(n)), '$,\ \theta(^\circ)\ vs\ RMSE\ plot$'], 'Interpreter', 'latex');
%xlabel('$\theta(^\circ)$', 'Interpreter', 'latex');
%ylabel('$RMSE$', 'Interpreter', 'latex');
%xlim([45 225]);
%ylim([0 1.5 * max(Case3.RMS_LS(3, n, :))]);

%% Case #3 Visualization

% Visuals change as setting these variables
angle_n = 91;
R_n = 2;

% Figure generation
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
y5 = (x-Case3.R(R_n)*cosd(Case3.theta(angle_n))) * tand(Case3.theta(angle_n)) + Case3.R(R_n) * sind(Case3.theta(angle_n));

p8 = plot(x, y1, '-g');
p9 = plot(x, y2, '-g');
p10 = plot(x, y3, '-g');
p11 = plot(x, y4, '-g');

l = find(x>0);
l = l(1);
p12 = plot(x(1:(l-1)), y5(1:(l-1)), '--r', 'LineWidth', 2);


axis([-0.3 0.3 -0.3 0.3]);
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
hold off;
legend([p1 p2 p3 p4 p5 p12], {'$Estimated\ Position(LS)$', '$Estimated\ Position(WLS)$', '$Estimated\ Position(2-Planets)$', '$Observer$', '$Planets$', '$Direction\ of\ the\ 3rd\ Planet$'}, 'Interpreter', 'latex', 'Location', 'southeast');
title(['$R = $', num2str(squeeze(Case3.R(R_n))), '$,\ \theta^{(\circ)} = $', num2str(squeeze(Case3.theta(angle_n)))], 'Interpreter', 'latex');

%% Case #4 Results

subplot(1, 2, 1);
plot(Case4.R, ones(1, length(Case4.R)) * Case2.Freq(length(Case2.Freq)), '--r', 'LineWidth', 3);
hold on;
for n = 1:length(Case4.theta)
    plot(Case4.R, Case4.Freq_LS(:, n));
end
%title({'$LS/2-Planet\ Comparision$';'$R\ vs\ Frequency\ of\ the\ Critical\ Error$'}, 'Interpreter', 'latex');
xlabel('$R$', 'Interpreter', 'latex');
ylabel('$Frequency\ of\ the\ Critical\ Error$', 'Interpreter', 'latex');
axis([100 1000 0 1]);
legend({'$2-Planet$', '$\theta=0$', '$\theta=10$', '$\theta=20$', '$\theta=30$', '$\theta=40$', '$\theta=50$', '$\theta=60$', '$\theta=70$', '$\theta=80$', '$\theta=90$'}, 'Location', 'southeast', 'NumColumns', 3,'Interpreter', 'latex');
hold off;


subplot(1, 2, 2);
plot(Case4.R, ones(1, length(Case4.R)) * Case2.Freq(length(Case2.Freq)),'--r',  'LineWidth', 3);
hold on;
for n = 1:length(Case4.theta)
    plot(Case4.R, Case4.Freq_WLS(:, n));
end
%title({'$WLS/2-Planet\ Comparision$';'$R\ vs\ Frequency\ of\ the\ Critical\ Error$'}, 'Interpreter', 'latex');
xlabel('$R$', 'Interpreter', 'latex');
ylabel('$Frequency\ of\ the\ Critical\ Error$', 'Interpreter', 'latex');
axis([100 1000 0 1]);
legend({'$2-Planet$', '$\theta=0$', '$\theta=10$', '$\theta=20$', '$\theta=30$', '$\theta=40$', '$\theta=50$', '$\theta=60$', '$\theta=70$', '$\theta=80$', '$\theta=90$'}, 'Location', 'southeast', 'NumColumns', 3,'Interpreter', 'latex');
hold off;

%% Case #4 Visualization

figure();
plot(squeeze(Case2.P_hat(1, :, length(Case2.theta))), squeeze(Case2.P_hat(2, :, length(Case2.theta))), '.');
title('$Estimated\ Position,\ 2-Planet$', 'Interpreter', 'latex');
xlabel('$X$', 'Interpreter', 'latex');
ylabel('$Y$', 'Interpreter', 'latex');
axis([-100 100 -2 2]);

figure();
plot(squeeze(Case4.P_hat_LS(1, :, 41, 10)), squeeze(Case4.P_hat_LS(2, :, 41, 10)), '.');
title('$Estimated\ Position,\ LS$', 'Interpreter', 'latex');
xlabel('$X$', 'Interpreter', 'latex');
ylabel('$Y$', 'Interpreter', 'latex');
axis([-100 100 -2 2]);

figure();
plot(squeeze(Case4.P_hat_WLS(1, :, 41, 10)), squeeze(Case4.P_hat_WLS(2, :, 41, 10)), '.');
title('$Estimated\ Position,\ WLS$', 'Interpreter', 'latex');
xlabel('$X$', 'Interpreter', 'latex');
ylabel('$Y$', 'Interpreter', 'latex');
axis([-100 100 -2 2]);