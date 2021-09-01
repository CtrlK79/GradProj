clear;
clc;

%% Case #2

% Setting parameters
Case2.sigma = [1;1];
Case2.Observer = [0;0];
Case2.theta = [0.01:0.01:10 170:0.01:180];
Case2.N = 5000;

% Making empty matrix for the results of the simulation
Case2.P_hat = zeros(2, Case2.N, length(Case2.theta));
Case2.Error = zeros(2, Case2.N, length(Case2.theta));
Case2.Counts = zeros(1, length(Case2.theta));
Case2.Freq = zeros(1, length(Case2.theta));

% Simulation loop
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

Freq_small = Case2.Freq.*(Case2.theta<90);
Freq_large = Case2.Freq.*(Case2.theta>90);

Case2.CritS = 0.01 * (1+[sum(Freq_small>0.1) sum(Freq_small>0.01) sum(Freq_small>0.001) sum(Freq_small>0.0001)]);
Case2.CritL = 180 - 0.01 * ([sum(Freq_large>0.1) sum(Freq_large>0.01) sum(Freq_large>0.001) sum(Freq_small>0.0001)]);

%% Save the Results
save('S1_C2.mat', 'Case2');