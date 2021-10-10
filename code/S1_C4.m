clear;
clc;

%% Case #4

% Setting parameters
Case4.sigma = [1;1;1];
Case4.Observer = [0;0];
Case4.theta = 0:10:90;
Case4.R = 100:10:1000;
Case4.N = 5000;

% Making empty matrix for the results of the simulation
Case4.P_hat_LS = zeros(2, Case4.N, length(Case4.R), length(Case4.theta));
Case4.P_hat_WLS = zeros(2, Case4.N, length(Case4.R), length(Case4.theta));
Case4.Error_LS = zeros(2, Case4.N, length(Case4.R), length(Case4.theta));
Case4.Error_WLS = zeros(2, Case4.N, length(Case4.R), length(Case4.theta));
Case4.Counts_LS = zeros(length(Case4.R), length(Case4.theta));
Case4.Counts_WLS = zeros(length(Case4.R), length(Case4.theta));
Case4.Freq_LS = zeros(length(Case4.R), length(Case4.theta));
Case4.Freq_WLS = zeros(length(Case4.R), length(Case4.theta));

% Simulation loop
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
        Case4.RMS_LS(:, i, j) = [rms(squeeze(Case4.Error_LS(:, :, i, j)), 2);norm(rms(squeeze(Case4.Error_LS(:, :, i, j)), 2))];
        Case4.RMS_WLS(:, i, j) = [rms(squeeze(Case4.Error_WLS(:, :, i, j)), 2);norm(rms(squeeze(Case4.Error_WLS(:, :, i, j)), 2))];
        
    end
end

Case4.Freq_LS = Case4.Counts_LS/Case4.N;
Case4.Freq_WLS = Case4.Counts_WLS/Case4.N;

%% Save the Results

save('S1_C4.mat', 'Case4');