clear;
clc;

%% Case #3

% Setting parameters
Case3.sigma = [1;1;1];
Case3.Observer = [0;0];
Case3.theta = 45:225;
Case3.N = 5000;
Case3.R = [1 10 100 1000];

% Making empty matrix for the results of the simulation
Case3.P_hat_LS = zeros(2, Case3.N, length(Case3.R), length(Case3.theta));
Case3.P_hat_WLS = zeros(2, Case3.N, length(Case3.R), length(Case3.theta));
Case3.Error_LS = zeros(2, Case3.N, length(Case3.R), length(Case3.theta));
Case3.Error_WLS = zeros(2, Case3.N, length(Case3.R), length(Case3.theta));
Case3.RMS_LS = zeros(3, length(Case3.R), length(Case3.theta));
Case3.RMS_WLS = zeros(3, length(Case3.R), length(Case3.theta));

% Simulation loop
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

%% Save the Results
save('S1_C3.mat', 'Case3');