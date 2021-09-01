clear;
clc;

%% Case #1

% Setting parameters
Case1_1.sigma = [1;1];
Case1_1.Observer = [0;0];
Case1_1.theta = 90;
Case1_1.N = 5000;
Case1_1.R = 1:1000;

% Making empty matrix for the results of the simulation
Case1_1.P_hat = zeros(2, Case1_1.N, length(Case1_1.R));
Case1_1.Error = zeros(2, Case1_1.N, length(Case1_1.R));
Case1_1.RMS = zeros(3, length(Case1_1.R));

% Simulation loop
for i = 1:length(Case1_1.R)
        Planets = [1 Case1_1.R(i)*cosd(Case1_1.theta); 0 Case1_1.R(i)*sind(Case1_1.theta)];
        Noise = randn(2, Case1_1.N) .* Case1_1.sigma;
        e1 = [cosd(Noise(1, :));sind(Noise(1, :))];
        e2 = [cosd(Case1_1.theta+Noise(2, :));sind(Case1_1.theta+Noise(2, :))];
        for k = 1:Case1_1.N
            EMatrix = [e1(:, k) e2(:, k)];
            Case1_1.P_hat(:, k, i) = GetPosition2D_LS(Planets, EMatrix);
        end
        Case1_1.Error(:, :, i) = squeeze(Case1_1.P_hat(:, :, i))-Case1_1.Observer;
        Case1_1.RMS(1:2, i) = rms(squeeze(Case1_1.Error(:, :, i)), 2);
        Case1_1.RMS(3, i) = norm(squeeze(Case1_1.RMS(1:2, i)));
end

%% Save the Results
save('S1_C1_1.mat', 'Case1_1');