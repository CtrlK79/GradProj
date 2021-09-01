clear;
clc;

%% Case #1

% Setting parameters
Case1.sigma = [1;1];
Case1.Observer = [0;0];
Case1.theta = 1:180;
Case1.N = 5000;
Case1.R = [1 10 100 1000];

% Making empty matrix for the results of the simulation
Case1.P_hat = zeros(2, Case1.N, length(Case1.R), length(Case1.theta));
Case1.Error = zeros(2, Case1.N, length(Case1.R), length(Case1.theta));
Case1.RMS = zeros(3, length(Case1.R), length(Case1.theta));

% Simulation loop
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

%% Save the Results
save('S1_C1.mat', 'Case1');