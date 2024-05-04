% Define the paths to the folders containing predicted and ground-truth segmented images
predictedFolder = 'C:\Users\Computing\Downloads\OneDrive_1_04-05-2024\Segmented_Image';
groundTruthFolder = 'C:\Users\Computing\Downloads\OneDrive_1_04-05-2024\Mask';

% Get a list of all files in the predicted and ground-truth folders
predictedFiles = dir(fullfile(predictedFolder, '*.png'));
groundTruthFiles = dir(fullfile(groundTruthFolder, '*.png'));

% Check if the number of files in both folders is the same
if numel(predictedFiles) ~= numel(groundTruthFiles)
    error('Number of files in predicted and ground truth folders do not match.');
end

% Initialize arrays to store Dice scores and image numbers
diceScores = zeros(numel(predictedFiles), 1);
imageNumbers = 1:numel(predictedFiles);

% Initialize arrays to store Dice scores and corresponding image names
diceAndNames = cell(numel(predictedFiles), 2);

% Loop through each pair of images
for i = 1:numel(predictedFiles)
    % Read the predicted and ground-truth images
    predictedImg = imread(fullfile(predictedFolder, predictedFiles(i).name));
    groundTruthImg = imread(fullfile(groundTruthFolder, groundTruthFiles(i).name));
    
    % Convert predicted mask to a single-channel binary mask
    if size(predictedImg, 3) == 3
        predictedMask = rgb2gray(predictedImg) > 0; % Convert to grayscale and create binary mask
    else
        predictedMask = predictedImg > 0;
    end
    
    % Convert ground truth mask to a binary mask
    groundTruthMask = groundTruthImg > 0;
    
    % Check if the dimensions of the masks match
    if ~isequal(size(predictedMask), size(groundTruthMask))
        fprintf('Dimensions of predicted mask: %s\n', mat2str(size(predictedMask)));
        fprintf('Dimensions of ground truth mask: %s\n', mat2str(size(groundTruthMask)));
        error('Dimensions of predicted and ground truth masks do not match.');
    end
    
    % Calculate the intersection and union of the masks
    intersection = sum(predictedMask(:) & groundTruthMask(:));
    union = sum(predictedMask(:)) + sum(groundTruthMask(:));
    
    % Calculate the Dice coefficient
    diceCoefficient = 2 * intersection / union;
    
    % Store the Dice coefficient for this image pair
    diceScores(i) = diceCoefficient;
    
    % Store the Dice coefficient and image name for this image
    diceAndNames{i, 1} = diceCoefficient;
    diceAndNames{i, 2} = predictedFiles(i).name;
end

% Sort Dice scores and get indices for best and worst 5 images
[sortedDiceScores, sortedIndices] = sort(diceScores, 'descend');
bestIndices = sortedIndices(1:5);
worstIndices = sortedIndices(end-4:end);

% Display the best and worst 5 Dice scores with their image numbers
fprintf('Best 5 Dice Scores:\n');
for i = 1:numel(bestIndices)
    fprintf('Image %d: %.4f\n', bestIndices(i), sortedDiceScores(i));
end
fprintf('\nWorst 5 Dice Scores:\n');
for i = 1:numel(worstIndices)
    fprintf('Image %d: %.4f\n', worstIndices(i), sortedDiceScores(end-i+1));
end

% Display the best and worst 5 segmented images
% Display the best and worst 5 segmented images
figure;

% Display best segmented images with corresponding ground truth masks
for i = 1:5
    % Display best segmented images
    subplot(4, 5, i);
    imshow(fullfile(predictedFolder, diceAndNames{bestIndices(i), 2}));
    title(['Best Dice Score: ', num2str(sortedDiceScores(i))]);
    
    % Load and display corresponding ground truth mask
    subplot(4, 5, i+5);
    groundTruthName = strrep(diceAndNames{bestIndices(i), 2}, '_segmented.png', '_GT.png');
    imshow(fullfile(groundTruthFolder, groundTruthName));
    title('Ground Truth Mask');
end

% Display worst segmented images with corresponding ground truth masks
for i = 1:5
    % Display worst segmented images
    subplot(4, 5, i+10);
    imshow(fullfile(predictedFolder, diceAndNames{worstIndices(i), 2}));
    title(['Worst Dice Score: ', num2str(sortedDiceScores(end-i+1))]);
    
    % Load and display corresponding ground truth mask
    subplot(4, 5, i+15);
    groundTruthName = strrep(diceAndNames{worstIndices(i), 2}, '_segmented.png', '_GT.png');
    imshow(fullfile(groundTruthFolder, groundTruthName));
    title('Ground Truth Mask');
end
sgtitle('Top Five Best and Worst Dice Score Images with Ground Truth Masks');
figure;
bar(imageNumbers, diceScores);
xlabel('Image Number');
ylabel('Dice Score');
title('Dice Scores for Segmented Images');
xticks(1:numel(predictedFiles));
xticklabels(imageNumbers);
grid on;

% Calculate mean and standard deviation of total Dice scores
meanDiceScore = mean(diceScores);
stdDiceScore = std(diceScores);

% Print mean and standard deviation
fprintf('Mean Dice Score: %.4f\n', meanDiceScore);
fprintf('Standard Deviation of Dice Scores: %.4f\n', stdDiceScore);
