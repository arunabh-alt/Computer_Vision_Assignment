x_true = csvread('x.csv');
y_true = csvread('y.csv');
na_noisy = csvread('na.csv');
nb_noisy = csvread('nb.csv');

delta_t = 0.5; % Time interval
F_baseline = [1 delta_t 0 0; 0 1 0 0; 0 0 1 delta_t; 0 0 0 1]; % State transition matrix
H_baseline = [1 0 0 0; 0 0 1 0]; % Observation matrix
Q_baseline = [0.16 0 0 0; 0 0.36 0 0; 0 0 0.16 0; 0 0 0 0.36]; % Process noise covariance
%Q_baseline = [0.03 0 0 0; 0 0.20 0 0; 0 0 0.03 0; 0 0 0 0.20];
R_baseline = [0.25 0; 0 0.25]; % Measurement noise covariance

% Initialize state and covariance matrices
x_est_baseline = [na_noisy(1); 0; nb_noisy(1); 0]; % Initial state estimate
P_baseline = eye(4); % Initial error covariance matrix

% Initialize variables for storing estimated trajectory
x_est_traj_baseline = zeros(size(x_true));
y_est_traj_baseline = zeros(size(y_true));

% Kalman filter loop
for i = 1:length(x_true)
    % Prediction step
    x_pred = F_baseline * x_est_baseline;
    P_pred = F_baseline * P_baseline * F_baseline' + Q_baseline;
    
    % Update step
    K = P_pred * H_baseline' * inv(H_baseline * P_pred * H_baseline' + R_baseline);
    z = [na_noisy(i); nb_noisy(i)]; % Observation
    x_est = x_pred + K * (z - H_baseline * x_pred);
    P_baseline = (eye(4) - K * H_baseline) * P_pred;
    
    % Store estimated trajectory
    x_est_traj_baseline(i) = x_est(1);
    y_est_traj_baseline(i) = x_est(3);
end
% Calculate RMSE
rmse_noisy_baseline = sqrt(mean((x_true - na_noisy).^2 + (y_true - nb_noisy).^2));
rmse_estimated_baseline = sqrt(mean((x_true - x_est_traj_baseline).^2 + (y_true - y_est_traj_baseline).^2));

% Finetuned Kalman Filter
% Define parameters
delta_t = 0.4; % Time interval
F = [1 delta_t 0 0; 0 1 0 0; 0 0 1 delta_t; 0 0 0 1]; % State transition matrix
H = [1 0 0 0; 0 0 1 0]; % Observation matrix
% Define parameter search ranges
Q_range = linspace(0.01, 0.019, 2); % Coarse search range for Q
R_range = linspace(0.1,0.2 , 100);    % Coarse search range for R

best_rmse = Inf; % Initialize best RMSE
best_Q = [];
best_R = [];
for q = Q_range
    for r = R_range
        % Define parameters
        Q = [q 0 0 0; 0 q 0 0; 0 0 q 0; 0 0 0 q]; % Process noise covariance
        R = [r 0; 0 r]; % Measurement noise covariance

        % Initialize state and covariance matrices
        x_est = [na_noisy(1); 0; nb_noisy(1); 0]; % Initial state estimate
        P = eye(4); % Initial error covariance matrix

        % Initialize variables for storing estimated trajectory
        x_est_traj = zeros(size(x_true));
        y_est_traj = zeros(size(y_true));

        % Kalman filter loop
        for i = 1:length(x_true)
            % Prediction step
            x_pred = F * x_est;
            P_pred = F * P * F' + Q;

            % Update step
            K = P_pred * H' * inv(H * P_pred * H' + R);
            z = [na_noisy(i); nb_noisy(i)]; % Observation
            x_est = x_pred + K * (z - H * x_pred);
            P = (eye(4) - K * H) * P_pred;

            % Store estimated trajectory
            x_est_traj(i) = x_est(1);
            y_est_traj(i) = x_est(3);
        end

        % Calculate RMSE
        rmse_estimated = sqrt(mean((x_true - x_est_traj).^2 + (y_true - y_est_traj).^2));

        % Update best RMSE and parameters if current RMSE is lower
        if rmse_estimated < best_rmse
            best_rmse = rmse_estimated;
            best_Q = Q;
            best_R = R;
        end
    end
end

fprintf('Best Finetuned RMSE: %.4f\n', best_rmse);
fprintf('Best Finetuned Q matrix:\n');
disp(best_Q);
fprintf('Best Finetuned R matrix:\n');
disp(best_R);
rmse_noisy = sqrt(mean((x_true - na_noisy).^2 + (y_true - nb_noisy).^2));
rmse_estimated = sqrt(mean((x_true - x_est_traj).^2 + (y_true - y_est_traj).^2));


fprintf('Ground Truth RMSE between true and noisy measurements: %.4f\n', rmse_noisy_baseline);
fprintf('Ground Truth RMSE between true and estimated coordinates: %.4f\n', rmse_estimated_baseline);
fprintf('Finetuned RMSE between true and noisy measurements: %.4f\n', rmse_noisy);
fprintf('Finetuned RMSE between true and estimated coordinates: %.4f\n', rmse_estimated);

figure;
subplot(2,1,1);
plot(x_true, y_true, 'b', 'LineWidth', 2); hold on;
plot(na_noisy, nb_noisy, 'r--', 'LineWidth', 1.5);
plot(x_est_traj_baseline, y_est_traj_baseline, 'g', 'LineWidth', 1.5);
xlabel('x');
ylabel('y');
legend('True trajectory', 'Noisy measurements', 'Estimated trajectory');
title('Ground Truth Trajectories of True, Noisy, and Estimated Coordinates');

subplot(2,1,2);
plot(x_true, y_true, 'b', 'LineWidth', 2); hold on;
plot(na_noisy, nb_noisy, 'r--', 'LineWidth', 1.5);
plot(x_est_traj, y_est_traj, 'g', 'LineWidth', 1.5);
xlabel('x');
ylabel('y');
legend('True trajectory', 'Noisy measurements', 'Estimated trajectory');
title('Finetuned Trajectories of True, Noisy, and Estimated Coordinates');


figure;
b = bar([rmse_estimated_baseline, rmse_estimated], 'BarWidth', 0.2, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none');
hold on;
% Add a red line connecting the tops of the bars
x = b(1).XEndPoints + b(1).XOffset; % X positions of the bars
y = b(1).YEndPoints; % Y positions of the tops of the bars
plot(x, y, 'r', 'LineWidth', 2, 'LineStyle', '--'); % Plot red dashed line
xticks(1:2);
xticklabels({'Baseline', 'Fine-tuned'});
ylabel('RMSE', 'FontSize', 12);
title('Comparison of RMSE: Baseline vs. Fine-tuned Kalman Filter', 'FontSize', 14);
set(gca, 'FontSize', 12); % Set font size for axis ticks
legend('Estimated RMSE', 'Location', 'northwest');
box off; % Remove box around the plot