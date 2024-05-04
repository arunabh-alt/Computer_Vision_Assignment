rgb_folder = 'C:\Users\Computing\Downloads\Final\Final\Color Image';
mask_folder = 'C:\Users\Computing\Downloads\Final\Final\Mask';

% Get the list of files in the RGB and mask folders
rgb_files = dir(fullfile(rgb_folder, '*.png')); % Assuming images are in PNG format
mask_files = dir(fullfile(mask_folder, '*.png')); % Assuming masks are in PNG format

% Check if the number of RGB images matches the number of mask images
if numel(rgb_files) ~= numel(mask_files)
    error('Number of RGB images does not match number of mask images');
end

% Define color names for the types of balls
colorNames = {'Yellow', 'White', 'Orange'};

% Initialize arrays to store features for each color ball
yellowSolidity = [];
yellowEccentricity = [];
yellowNonCompactness = [];
yellowCircularity = [];

whiteSolidity = [];
whiteEccentricity = [];
whiteNonCompactness = [];
whiteCircularity = [];

orangeSolidity = [];
orangeEccentricity = [];
orangeNonCompactness = [];
orangeCircularity = [];

% Iterate through each pair of RGB image and mask
for idx = 1:numel(rgb_files)
    % Read RGB image and corresponding mask
    rgb_image = imread(fullfile(rgb_folder, rgb_files(idx).name));
    mask_image = imread(fullfile(mask_folder, mask_files(idx).name));
    binary_image = imbinarize(rgb2gray(rgb_image));

    % Find connected components in the ground truth mask
    cc = bwconncomp(mask_image);
    stats = regionprops(cc, 'BoundingBox');
    
    % Extract ball patches using bounding boxes
    ball_patches = cell(1, cc.NumObjects);
    for j = 1:cc.NumObjects
        bb = stats(j).BoundingBox;
        ball_patches{j} = imcrop(binary_image, bb);
    end
    
    % Sort the patches based on area
    areas = zeros(1, cc.NumObjects);
    for k = 1:cc.NumObjects
        bb = stats(k).BoundingBox;
        areas(k) = bb(3) * bb(4); % width * height
    end
    [~, sorted_indices] = sort(areas);
    sorted_patches = ball_patches(sorted_indices);
    
    % Store the patches: smallest is yellow, medium is white, large is orange
    num_patches = numel(sorted_patches);
    
    for j = 1:num_patches
        if j <= ceil(num_patches / 3)
            color = 'Yellow';
            idx_color_1 = numel(yellowSolidity) + 1;
            idx_color_2= numel(yellowCircularity)+1;
            idx_color_3= numel(yellowNonCompactness)+1;
            idx_color_4= numel(yellowEccentricity)+1;
        elseif j <= ceil(2 * num_patches / 3)
            color = 'White';
            idx_color_5 = numel(whiteSolidity) + 1;
            idx_color_6= numel(whiteCircularity)+1;
            idx_color_7= numel(whiteNonCompactness)+1;
            idx_color_8= numel(whiteEccentricity)+1;
        else
            color = 'Orange';
            idx_color_9 = numel(orangeSolidity) + 1;
            idx_color_10= numel(orangeCircularity)+1;
            idx_color_11= numel(orangeNonCompactness)+1;
            idx_color_12= numel(orangeEccentricity)+1;
        end
        
        % Calculate additional shape features
        stats_region = regionprops('table', sorted_patches{j}, ...
            'Solidity', 'Eccentricity', 'Area', 'Perimeter','Circularity');
        area=stats_region.Area;
       
        % circularity = (stats_region.Perimeter.^2) ./ (4 * pi * stats_region.Area);
        circularity = stats_region.Circularity;
        non_compactness = (1/ circularity)+((2*6*pi)/area);
        
        % Store the features based on color
        switch color
            case 'Yellow'
                yellowSolidity{idx_color_1} = mean(stats_region.Solidity);
                yellowCircularity{idx_color_2} = mean(circularity);
                yellowNonCompactness{idx_color_3} = mean(non_compactness);
                yellowEccentricity{idx_color_4} = mean(stats_region.Eccentricity);
            case 'White'
                whiteSolidity{idx_color_5} = mean(stats_region.Solidity);
                whiteCircularity{idx_color_6} = mean(circularity);
                whiteNonCompactness{idx_color_7} = mean(non_compactness);
                whiteEccentricity{idx_color_8} = mean(stats_region.Eccentricity);
            case 'Orange'
                orangeSolidity{idx_color_9} = mean(stats_region.Solidity);
                orangeCircularity{idx_color_10} = mean(circularity);
                orangeNonCompactness{idx_color_11} = mean(non_compactness);
                orangeEccentricity{idx_color_12} = mean(stats_region.Eccentricity);
        end
    end
end
% Convert cell arrays to matrices for yellow, white, and orange color balls
yellowSolidity_all = cat(1, yellowSolidity{:}); 
yellowNonCompactness_all = cat(1, yellowNonCompactness{:}); 
yellowCircularity_all = cat(1, yellowCircularity{:}); 
yellowEccentricity_all = cat(1, yellowEccentricity{:}); 

whiteSolidity_all = cat(1, whiteSolidity{:}); 
whiteNonCompactness_all = cat(1, whiteNonCompactness{:}); 
whiteCircularity_all = cat(1, whiteCircularity{:});
whiteEccentricity_all = cat(1, whiteEccentricity{:})';

orangeSolidity_all = cat(1, orangeSolidity{:}); 
orangeNonCompactness_all = cat(1, orangeNonCompactness{:}); 
orangeCircularity_all = cat(1, orangeCircularity{:}); 
orangeEccentricity_all = cat(1, orangeEccentricity{:}); 
num_bins=10;

%Plot histograms for shape features of yellow color ball
figure;
subplot(2, 2, 1);
histogram(yellowSolidity_all,num_bins, 'FaceColor', [1 1 0]);
title('Yellow Ball - Solidity');
xlabel('Solidity');
ylabel('Frequency');

subplot(2, 2, 2);
histogram(yellowNonCompactness_all,num_bins, 'FaceColor', [1 1 0]);
title('Yellow Ball - Non-Compactness');
xlabel('Non-Compactness');
ylabel('Frequency');

subplot(2, 2, 3);
histogram(yellowCircularity_all,num_bins, 'FaceColor', [1 1 0]);
title('Yellow Ball - Circularity');
xlabel('Circularity');
ylabel('Frequency');

subplot(2, 2, 4);
histogram(yellowEccentricity_all,num_bins, 'FaceColor', [1 1 0]);
title('Yellow Ball - Eccentricity');
xlabel('Eccentricity');
ylabel('Frequency');
sgtitle('Shape Features Distribution For Yellow Ball');

% Plot histograms for shape features of white color ball
figure;
subplot(2, 2, 1);
histogram(whiteSolidity_all,num_bins, 'FaceColor', [0/255, 21/255, 71/255]);
title('White Ball - Solidity');
xlabel('Solidity');
ylabel('Frequency');

subplot(2, 2, 2);
histogram(whiteNonCompactness_all,num_bins, 'FaceColor', [0/255, 21/255, 71/255]);
title('White Ball - Non-Compactness');
xlabel('Non-Compactness');
ylabel('Frequency');

subplot(2, 2, 3);
histogram(whiteCircularity_all,num_bins, 'FaceColor', [0/255, 21/255, 71/255]);
title('White Ball - Circularity');
xlabel('Circularity');
ylabel('Frequency');

subplot(2, 2, 4);
histogram(whiteEccentricity_all,num_bins, 'FaceColor', [0/255, 21/255, 71/255]);
title('White Ball - Eccentricity');
ylabel('Frequency');
xlabel('Eccentricity');
sgtitle('Shape Features Distribution For White Ball');

% Plot histograms for shape features of orange color ball
figure;
subplot(2, 2, 1);
histogram(orangeSolidity_all, num_bins,'FaceColor', [1 0 0]);
title('Orange Ball - Solidity');
xlabel('Solidity');
ylabel('Frequency');

subplot(2, 2, 2);
histogram(orangeNonCompactness_all,num_bins, 'FaceColor', [1 0 0]);
title('Orange Ball - Non-Compactness');
xlabel('Non-Compactness');
ylabel('Frequency');

subplot(2, 2, 3);
histogram(orangeCircularity_all,num_bins, 'FaceColor', [1 0 0]);
title('Orange Ball - Circularity');
xlabel('Circularity');
ylabel('Frequency');

subplot(2, 2, 4);
histogram(orangeEccentricity_all,num_bins, 'FaceColor', [1 0 0]);
title('Orange Ball - Eccentricity');
xlabel('Eccentricity');
ylabel('Frequency');
sgtitle('Shape Features Distribution For Orange Ball');

figure;

% Set number of bins
%num_bins = 20;

% Plot histograms for Solidity Features
subplot(2,2,1);
hold on;
num_bins1 = getBinNumber(yellowSolidity_all,whiteSolidity_all,orangeSolidity_all);
histogram(yellowSolidity_all, num_bins1, 'FaceColor', [1 1 0], 'EdgeColor', 'none');
histogram(whiteSolidity_all, num_bins1, 'FaceColor', [0/255, 21/255, 71/255], 'EdgeColor', 'none');
histogram(orangeSolidity_all, num_bins1, 'FaceColor', [1 0 0], 'EdgeColor', 'none');
xlabel('Solidity');
ylabel('Frequency');
title('Solidity Features Distribution For All Balls');
legend('Yellow Ball', 'White Ball', 'Orange Ball');
hold off;

% Plot histograms for Non-Compactness Features
subplot(2,2,2);
hold on;
num_bins2=getBinNumber(yellowNonCompactness_all,whiteNonCompactness_all,orangeNonCompactness_all);
histogram(yellowNonCompactness_all, num_bins2, 'FaceColor', [1 1 0], 'EdgeColor', 'none');
histogram(whiteNonCompactness_all + 0.02, num_bins2, 'FaceColor', [0/255, 21/255, 71/255], 'EdgeColor', 'none');
histogram(orangeNonCompactness_all + 0.04, num_bins2, 'FaceColor', [1 0 0], 'EdgeColor', 'none');
xlabel('Non-Compactness');
ylabel('Frequency');
title('Non-Compactness Features Distribution For All Balls');
legend('Yellow Ball', 'White Ball', 'Orange Ball');
hold off;

% Plot histograms for Circularity Features
subplot(2,2,3);
hold on;
num_bins3 = getBinNumber(yellowCircularity_all,whiteCircularity_all,orangeCircularity_all);
histogram(yellowCircularity_all, num_bins3, 'FaceColor', [1 1 0], 'EdgeColor', 'none');
histogram(whiteCircularity_all + 0.02, num_bins3, 'FaceColor', [0/255, 21/255, 71/255], 'EdgeColor', 'none');
histogram(orangeCircularity_all + 0.04, num_bins3, 'FaceColor', [1 0 0], 'EdgeColor', 'none');
xlabel('Circularity');
ylabel('Frequency');
title('Circularity Features Distribution For All Balls');
legend('Yellow Ball', 'White Ball', 'Orange Ball');
hold off;

% Plot histograms for Eccentricity Features
subplot(2,2,4);
hold on;
num_bins4=getBinNumber(yellowEccentricity_all,whiteEccentricity_all,orangeEccentricity_all);
histogram(yellowEccentricity_all, num_bins4, 'FaceColor', [1 1 0], 'EdgeColor', 'none');
histogram(whiteEccentricity_all + 0.02, num_bins4, 'FaceColor', [0/255, 21/255, 71/255], 'EdgeColor', 'none');
histogram(orangeEccentricity_all + 0.04, num_bins4, 'FaceColor', [1 0 0], 'EdgeColor', 'none');
xlabel('Eccentricity');
ylabel('Frequency');
title('Eccentricity Features Distribution For All Balls');
legend('Yellow Ball', 'White Ball', 'Orange Ball');
hold off;
sgtitle('All Shape Features for All Ball Patches');
function binEdges = getBinNumber(val_1, val_2, val_3)
    % Combine all values into one array
    allValues = [val_1(:); val_2(:); val_3(:)];
    
    % Find the minimum and maximum values
    min_val = min(allValues);
    max_val = max(allValues);
    
    
    Number = 10; % Assuming each bin has a range of 10 units
    
    % Calculate the bin edges
    binEdges = linspace(min_val, max_val, Number+1);
end
