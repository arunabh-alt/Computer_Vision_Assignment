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
errors_baseline = zeros(size(x_true));
errors_baseline_noisy=zeros(size(x_true));
% Define gate threshold
threshold = 1.345977052032410e+06;
% Kalman filter loop
for i = 1:length(x_true)
    % Prediction step
    x_pred = F_baseline * x_est_baseline;
    P_pred = F_baseline * P_baseline * F_baseline' + Q_baseline;
    
    % Calculate predicted measurement
    z_pred = H_baseline * x_pred;
    
    % Calculate innovation covariance
    S = H_baseline * P_pred * H_baseline' + R_baseline;
    z = [na_noisy(i); nb_noisy(i)];
    % Define Gate Logic
    gate = (z-z_pred)'*inv(S) *(z-z_pred);
    % Check if measurement is within validation gate
    if gate <= threshold
        disp('Measurement is validated By gate');
        % Update step
        K = P_pred * H_baseline' * inv(S);
        x_est = x_pred + K * (z - H_baseline * x_pred);
        P_baseline = (eye(4) - K * H_baseline) * P_pred;
        
        % Store estimated trajectory
        x_est_traj_baseline(i) = x_est(1);
        y_est_traj_baseline(i) = x_est(3);
        %Calculate The Noisy RMSE for Every Coordinate  
        errors_baseline_noisy(i) = sqrt(mean((x_true(i) - na_noisy(i))^2 + (y_true(i) - nb_noisy(i))^2));
        % Calculate RMSE for Every coordinate
        errors_baseline(i) = sqrt(mean((x_true(i) - x_est_traj_baseline(i))^2 + (y_true(i) - y_est_traj_baseline(i))^2));
    else
        % If measurement is outside validation gate, do not update state estimate
        disp('Measurement is outside validation gate');
        x_est = x_pred;
        % Store estimated trajectory
        x_est_traj_baseline(i) = x_est(1);
        y_est_traj_baseline(i) = x_est(3);
        %Calculate The Noisy RMSE for Every Coordinate  
        errors_baseline_noisy(i) = sqrt(mean((x_true(i) - na_noisy(i))^2 + (y_true(i) - nb_noisy(i))^2));
        % Calculate RMSE for Every coordinate without updating the state estimate
        errors_baseline(i) = sqrt(mean((x_true(i) - x_est_traj_baseline(i))^2 + (y_true(i) - y_est_traj_baseline(i))^2));
    end
end
% Calculate RMSE

rmse_estimated_baseline_mean = mean(errors_baseline);
rmse_estimated_baseline_std = std(errors_baseline);

rmse_estimated_baseline_mean_noisy = mean(errors_baseline_noisy);
rmse_estimated_baseline_std_noisy = std(errors_baseline_noisy);

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
errors_finetune_estimated = zeros(size(x_true));
errors_finetune_noisy=zeros(size(x_true));
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
            % Calculate RMSE For Noisy Coordinate
            errors_finetune_noisy(i) = sqrt(mean((x_true(i) - na_noisy(i))^2 + (y_true(i) - nb_noisy(i))^2));
            % Calculate RMSE For Estimated Coordinate
            errors_finetune_estimated(i) = sqrt(mean((x_true(i) - x_est_traj(i))^2 + (y_true(i) - y_est_traj(i))^2));
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

rmse_estimated_finetuned_mean = mean(errors_finetune_estimated);
rmse_estimated_finetuned_std = std(errors_finetune_estimated);

rmse_estimated_finetuned_mean_noisy = mean(errors_finetune_noisy);
rmse_estimated_finetuned_std_noisy = std(errors_finetune_noisy);

std_error_x = std(x_true - x_est_traj);
std_error_y = std(y_true - y_est_traj);
std_error_x_baseline = std(x_true - x_est_traj_baseline );
std_error_y_baseline = std(y_true- y_est_traj_baseline);
fprintf('Baseline Filter RMSE in between true and Noisy coordinates Mean : %.4f, Standard Deviation: %.4f\n', rmse_estimated_baseline_mean_noisy, rmse_estimated_baseline_std_noisy);
fprintf('Finetuned Filter RMSE in between true and Noisy coordinates Mean : %.4f, Standard Deviation: %.4f\n', rmse_estimated_finetuned_mean_noisy,rmse_estimated_finetuned_std_noisy);
fprintf('Baseline Filter RMSE between true and estimated coordinates Mean : %.4f, Standard Deviation: %.4f\n', rmse_estimated_baseline_mean, rmse_estimated_baseline_std);
fprintf('Finetuned Filter RMSE between true and estimated coordinates Mean : %.4f, Standard Deviation: %.4f\n', rmse_estimated_finetuned_mean,rmse_estimated_finetuned_std);
fprintf('Baseline Kalman Filter Standard Deviation Error in X direction Error: %.4f, Y Direction Error : %.4f\n', std_error_x_baseline, std_error_y_baseline);
fprintf('Fine Tuned Kalman Filter Standard Deviation Error in X direction Error: %.4f, Y Direction Error : %.4f\n', std_error_x, std_error_y);


% Plotting
figure;

% Subplot 1 - Baseline Trajectory
subplot(2, 2, 1);
plot(x_true, y_true, 'b', 'LineWidth', 2); hold on;
plot(na_noisy, nb_noisy, 'r--', 'LineWidth', 1.5);
plot(x_est_traj_baseline, y_est_traj_baseline, 'g', 'LineWidth', 1.5);
xlabel('x');
ylabel('y');
legend('True trajectory', 'Noisy measurements', 'Estimated trajectory');
title('Baseline Trajectories');

% Subplot 2 - Fine-tuned Trajectory
subplot(2, 2, 2);
plot(x_true, y_true, 'b', 'LineWidth', 2); hold on;
plot(na_noisy, nb_noisy, 'r--', 'LineWidth', 1.5);
plot(x_est_traj, y_est_traj, 'g', 'LineWidth', 1.5);
xlabel('x');
ylabel('y');
legend('True trajectory', 'Noisy measurements', 'Estimated trajectory');
title('Fine-tuned Trajectories');

% Subplot 3 - Mean RMSE Comparison Bar Graph
subplot(2, 2, 3);
bar([rmse_estimated_baseline_mean, rmse_estimated_finetuned_mean], 'BarWidth', 0.4, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none');
hold on;
errorbar(1:2, [rmse_estimated_baseline_mean, rmse_estimated_finetuned_mean], [rmse_estimated_baseline_std, rmse_estimated_finetuned_std], '.', 'Color', 'r', 'LineWidth', 2);
ylabel('RMSE');
xticks(1:2);
xticklabels({'Baseline', 'Fine-tuned'});
title('Mean RMSE Comparison');
legend('Mean RMSE', 'Standard Deviation');

% Subplot 4 - X and Y Direction Error Bar Graph
subplot(2, 2, 4);
bar([std_error_x_baseline, std_error_y_baseline; std_error_x, std_error_y], 'BarWidth', 0.4, 'FaceColor', [0.5 0.5 0.5]);
hold on;
ylabel('Standard Deviation Error');
xticks(1:2);
xticklabels({'Baseline', 'Fine-tuned'});
legend({'X Direction', 'Y Direction'});
title('X and Y Direction Error Comparison');


sgtitle('Kalman Filter Analysis');
