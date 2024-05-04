% Define the folders containing RGB images and mask images
image_folder = 'C:\Users\Computing\Downloads\Final\Final\Color Image';
mask_folder = 'C:\Users\Computing\Downloads\Final\Final\Mask';

% Get list of image and mask files
image_files = dir(fullfile(image_folder, '*.png'));
mask_files = dir(fullfile(mask_folder, '*.png'));

% Initialize cell arrays to store extracted ball patches, contrast, correlation, and ASM values
red_yellow_balls_contrast = cell(1, 63);
red_yellow_balls_correlation = cell(1, 63);
red_yellow_balls_asm = cell(1, 63);
red_white_balls_contrast = cell(1, 63);
red_white_balls_correlation = cell(1, 63);
red_white_balls_asm = cell(1, 63);
red_orange_balls_contrast = cell(1, 63);
red_orange_balls_correlation = cell(1, 63);
red_orange_balls_asm = cell(1, 63);

green_yellow_balls_contrast = cell(1, 63);
green_yellow_balls_correlation = cell(1, 63);
green_yellow_balls_asm = cell(1, 63);
green_white_balls_contrast = cell(1, 63);
green_white_balls_correlation = cell(1, 63);
green_white_balls_asm = cell(1, 63);
green_orange_balls_contrast = cell(1, 63);
green_orange_balls_correlation = cell(1, 63);
green_orange_balls_asm = cell(1, 63);

blue_yellow_balls_contrast = cell(1, 63);
blue_yellow_balls_correlation = cell(1, 63);
blue_yellow_balls_asm = cell(1, 63);
blue_white_balls_contrast = cell(1, 63);
blue_white_balls_correlation = cell(1, 63);
blue_white_balls_asm = cell(1, 63);
blue_orange_balls_contrast = cell(1, 63);
blue_orange_balls_correlation = cell(1, 63);
blue_orange_balls_asm = cell(1, 63);

red_yellow_balls_contrast_range = cell(1, 63);
red_white_balls_contrast_range = cell(1, 63);
red_orange_balls_contrast_range = cell(1, 63);

green_yellow_balls_correlation_range = cell(1, 63);
green_white_balls_correlation_range = cell(1, 63);
green_orange_balls_correlation_range = cell(1, 63);

blue_yellow_balls_asm_range = cell(1, 63);
blue_white_balls_asm_range = cell(1, 63);
blue_orange_balls_asm_range = cell(1, 63);

for i = 1:length(image_files)
    % Read image and mask
    image = imread(fullfile(image_folder, image_files(i).name));
    mask = imread(fullfile(mask_folder, mask_files(i).name));
    binary_image = imbinarize(rgb2gray(image));
    
    % Find connected components in the ground truth mask
    cc = bwconncomp(mask);
    stats = regionprops(cc, 'BoundingBox','Area');
    
    % Extract ball patches using bounding boxes from the binary image
    ball_patches = cell(1, cc.NumObjects);
    for j = 1:cc.NumObjects
        bb = stats(j).BoundingBox;
        ball_patches{j} = imcrop(binary_image, bb);
    end

    % Sort the patches based on area
    areas = zeros(1, cc.NumObjects);
    for j = 1:cc.NumObjects
        bb = stats(j).BoundingBox;
        areas(j) = bb(3) * bb(4); % width * height
    end
    [~, sorted_indices] = sort(areas);
    sorted_patches = ball_patches(sorted_indices);
    
    % Convert sorted patches back to RGB using original image
    num_patches = numel(sorted_patches);
    for j = 1:num_patches
        % Crop corresponding region from original RGB image
        stats = regionprops(cc, 'BoundingBox');
        bb = stats(sorted_indices(j)).BoundingBox;
        ball_rgb_patch = imcrop(image, bb);
        
        % Calculate GLCM and texture features only for yellow balls
        if j <= ceil(num_patches / 3)
            % Extract red channel for GLCM calculation
            for channel = 1:3
                channel_patch = ball_rgb_patch(:, :, channel);
                offsets = [0 1; -1 1; -1 0; -1 -1];

                glcms = graycomatrix(channel_patch, 'Offset', offsets, 'Symmetric', true, 'NumLevels', 256);
            
                % Normalize the GLCM (assuming each slice represents a different orientation)
                total_sum = sum(glcms(:));
                scaling_factor =550;
                scaled_glcm = glcms*scaling_factor;
                normalized_glcm = scaled_glcm / total_sum;
                normalized_glcm=uint8(normalized_glcm);
                contrast_values = zeros(4, 1);
                correlation_values = zeros(4, 1);
                asm_values = zeros(4, 1);
                for k = 1:4
                    stats = graycoprops(normalized_glcm(:,:,k), {'Contrast', 'Correlation', 'Energy'});
                    contrast_values(k) = mean(stats.Contrast);
                    correlation_values(k) = mean(stats.Correlation);
                    asm_values(k) = mean(stats.Energy);
                end
                if channel == 1
                    % Red channel
                    red_yellow_balls_contrast{i}{j} = mean(contrast_values);
                    red_yellow_balls_correlation{i}{j} = mean(correlation_values);
                    red_yellow_balls_asm{i}{j} = mean(asm_values);
                    red_yellow_balls_contrast_range{i}{j} = max(contrast_values)-min(contrast_values);
                    
                elseif channel == 2
                    % Green channel
                    green_yellow_balls_contrast{i}{j} = mean(contrast_values);
                    green_yellow_balls_correlation{i}{j} = mean(correlation_values);
                    green_yellow_balls_asm{i}{j} = mean(asm_values);
                    
                    green_yellow_balls_correlation_range{i}{j} = max(correlation_values)-min(correlation_values);
                    
                else
                    % Blue channel
                    blue_yellow_balls_contrast{i}{j} = mean(contrast_values);
                    blue_yellow_balls_correlation{i}{j} = mean(correlation_values);
                    blue_yellow_balls_asm{i}{j} = mean(asm_values);
                    
                    blue_yellow_balls_asm_range{i}{j} = max(asm_values)-min(asm_values);
                end    
            end
        elseif j <= ceil(2 * num_patches / 3)
                % Extract red channel for GLCM calculation
                for channel = 1:3
                    channel_patch = ball_rgb_patch(:, :, channel);
                    offsets = [0 1; -1 1; -1 0; -1 -1];

                    glcms = graycomatrix(channel_patch, 'Offset', offsets, 'Symmetric', true, 'NumLevels', 256);
                
                    % Normalize the GLCM (assuming each slice represents a different orientation)
                    total_sum = sum(glcms(:));
                    scaling_factor = 550;
                    scaled_glcm = glcms*scaling_factor;
                    normalized_glcm = scaled_glcm / total_sum;
                    normalized_glcm=uint8(normalized_glcm);
                    contrast_values = zeros(4, 1);
                    correlation_values = zeros(4, 1);
                    asm_values = zeros(4, 1);
                    for k = 1:4
                        stats = graycoprops(normalized_glcm(:,:,k), {'Contrast', 'Correlation', 'Energy'});
                        contrast_values(k) = mean(stats.Contrast);
                        correlation_values(k) = mean(stats.Correlation);
                        asm_values(k) = mean(stats.Energy);
                    end
                    if channel == 1
                        % Red channel
                        red_white_balls_contrast{i}{j} = mean(contrast_values);
                        red_white_balls_correlation{i}{j} = mean(correlation_values);
                        red_white_balls_asm{i}{j} = mean(asm_values);
                        red_white_balls_contrast_range{i}{j} = max(contrast_values)-min(contrast_values);
                       
                    elseif channel == 2
                        % Green channel
                        green_white_balls_contrast{i}{j} = mean(contrast_values);
                        green_white_balls_correlation{i}{j} = mean(correlation_values);
                        green_white_balls_asm{i}{j} = mean(asm_values);
                       
                        green_white_balls_correlation_range{i}{j} = max(correlation_values)- min(correlation_values);
                        
                    else
                        % Blue channel
                        blue_white_balls_contrast{i}{j} = mean(contrast_values);
                        blue_white_balls_correlation{i}{j} = mean(correlation_values);
                        blue_white_balls_asm{i}{j} = mean(asm_values);
                        
                        blue_white_balls_asm_range{i}{j} = max(asm_values)- min(asm_values);
                    end  
                end
        else 
            for channel = 1:3
                    channel_patch = ball_rgb_patch(:, :, channel);
                    offsets = [0 1; -1 1; -1 0; -1 -1];

                    glcms = graycomatrix(channel_patch, 'Offset', offsets, 'Symmetric', true, 'NumLevels', 256);
                
                    % Normalize the GLCM (assuming each slice represents a different orientation)
                    total_sum = sum(glcms(:));
                    scaling_factor = 550;
                    scaled_glcm = glcms*scaling_factor;
                    normalized_glcm = scaled_glcm / total_sum;
                    normalized_glcm=uint8(normalized_glcm);
                    contrast_values = zeros(4, 1);
                    correlation_values = zeros(4, 1);
                    asm_values = zeros(4, 1);
                    for k = 1:4
                        stats = graycoprops(normalized_glcm(:,:,k), {'Contrast', 'Correlation', 'Energy'});
                        contrast_values(k) = mean(stats.Contrast);
                        correlation_values(k) = mean(stats.Correlation);
                        asm_values(k) = mean(stats.Energy);
                    end
                    if channel == 1
                        % Red channel
                        red_orange_balls_contrast{i}{j} = mean(contrast_values);
                        red_orange_balls_correlation{i}{j} = mean(correlation_values);
                        red_orange_balls_asm{i}{j} = mean(asm_values);
                        red_orange_balls_contrast_range{i}{j} = max(contrast_values)- min(contrast_values);
                       
                    elseif channel == 2
                        % Green channel
                        green_orange_balls_contrast{i}{j} = mean(contrast_values);
                        green_orange_balls_correlation{i}{j} = mean(correlation_values);
                        green_orange_balls_asm{i}{j} = mean(asm_values);
                       
                        green_orange_balls_correlation_range{i}{j} = max(correlation_values)- min(correlation_values);
                        
                    else
                        % Blue channel
                        blue_orange_balls_contrast{i}{j} = mean(contrast_values);
                        blue_orange_balls_correlation{i}{j} = mean(correlation_values);
                        blue_orange_balls_asm{i}{j} = mean(asm_values);
                       
                        blue_orange_balls_asm_range{i}{j} = max(asm_values)- min(asm_values);
                    end    
              end
        end
    end
end    

% Define the number of bins for the histograms
num_bins = 10;


% Flatten the nested cell arrays
red_contrast_all = cell2mat(cellfun(@(x) cell2mat(x), red_yellow_balls_contrast, 'UniformOutput', false));
red_white_contrast_all = cell2mat(cellfun(@(x) cell2mat(x), red_white_balls_contrast, 'UniformOutput', false));
red_orange_contrast_all = cell2mat(cellfun(@(x) cell2mat(x), red_orange_balls_contrast, 'UniformOutput', false));
red_contrast_all_range = cell2mat(cellfun(@(x) cell2mat(x), red_yellow_balls_contrast_range, 'UniformOutput', false));
red_white_contrast_all_range = cell2mat(cellfun(@(x) cell2mat(x), red_white_balls_contrast_range, 'UniformOutput', false));
red_orange_contrast_all_range = cell2mat(cellfun(@(x) cell2mat(x), red_orange_balls_contrast_range, 'UniformOutput', false));



green_correlation_all = cell2mat(cellfun(@(x) cell2mat(x), green_yellow_balls_correlation, 'UniformOutput', false));
green_white_correlation_all = cell2mat(cellfun(@(x) cell2mat(x), green_white_balls_correlation, 'UniformOutput', false));
green_orange_correlation_all = cell2mat(cellfun(@(x) cell2mat(x), green_orange_balls_correlation, 'UniformOutput', false));
green_correlation_all_range = cell2mat(cellfun(@(x) cell2mat(x), green_yellow_balls_correlation_range, 'UniformOutput', false));
green_white_correlation_all_range = cell2mat(cellfun(@(x) cell2mat(x), green_white_balls_correlation_range, 'UniformOutput', false));
green_orange_correlation_all_range = cell2mat(cellfun(@(x) cell2mat(x), green_orange_balls_correlation_range, 'UniformOutput', false));

blue_asm_all = cell2mat(cellfun(@(x) cell2mat(x), blue_yellow_balls_asm, 'UniformOutput', false));
blue_white_asm_all = cell2mat(cellfun(@(x) cell2mat(x), blue_white_balls_asm, 'UniformOutput', false));
blue_orange_asm_all = cell2mat(cellfun(@(x) cell2mat(x), blue_orange_balls_asm, 'UniformOutput', false));
blue_asm_all_range = cell2mat(cellfun(@(x) cell2mat(x), blue_yellow_balls_asm_range, 'UniformOutput', false));
blue_white_asm_all_range = cell2mat(cellfun(@(x) cell2mat(x), blue_white_balls_asm_range, 'UniformOutput', false));
blue_orange_asm_all_range = cell2mat(cellfun(@(x) cell2mat(x), blue_orange_balls_asm_range, 'UniformOutput', false));



% Figure 4: Red Channel
figure;
subplot(2,3,1);
histogram(red_contrast_all,num_bins,'FaceColor', [1 1 0]);
title('Histogram of red channel\_yellow ball\_contrast');
xlabel('Contrast');
ylabel('Frequency');
subplot(2,3,2);
histogram(red_white_contrast_all,num_bins,'FaceColor', [0/255, 21/255, 71/255]);
title('Histogram of red\_white\_contrast');
xlabel('Contrast');
ylabel('Frequency');
subplot(2,3,3);
histogram(red_orange_contrast_all,num_bins,'FaceColor', [1 0 0]);
title('Histogram of red\_orange\_contrast');
xlabel('Contrast');
ylabel('Frequency');
subplot(2,3,4);
histogram(red_contrast_all_range,num_bins,'FaceColor', [1 1 0]);
title('Histogram of red\_yellow\_contrast\_range');
xlabel('Contrast');
ylabel('Frequency');
subplot(2,3,5);
histogram(red_white_contrast_all_range,num_bins,'FaceColor', [0/255, 21/255, 71/255]);
title('Histogram of red\_white\_contrast\_range');
xlabel('Contrast');
ylabel('Frequency');
subplot(2,3,6);
histogram(red_orange_contrast_all_range,num_bins,'FaceColor', [1 0 0]);
title('Histogram of red\_orange\_contrast\_range');
xlabel('Contrast');
ylabel('Frequency');
sgtitle('Red Channel Ball Contrast');

% Figure 5: Green Channel
figure;
subplot(2,3,1);
histogram(green_correlation_all,num_bins,'FaceColor', [1 1 0]);
title('Histogram of green\_yellow\_correlation');
xlabel('Correlation');
ylabel('Frequency');
subplot(2,3,2);
histogram(green_white_correlation_all,num_bins,'FaceColor', [0/255, 21/255, 71/255]);
title('Histogram of green\_white\_correlation');
xlabel('Correlation');
ylabel('Frequency');
subplot(2,3,3);
histogram(green_orange_correlation_all,num_bins,'FaceColor', [1 0 0]);
title('Histogram of green\_orange\_correlation');
xlabel('Correlation');
ylabel('Frequency');
subplot(2,3,4);
histogram(green_correlation_all_range,num_bins,'FaceColor', [1 1 0]);
title('Histogram of green\_yellow\_correlation\_range');
xlabel('Correlation');
ylabel('Frequency');
subplot(2,3,5);
histogram(green_white_correlation_all_range,num_bins,'FaceColor', [0/255, 21/255, 71/255]);
title('Histogram of green\_white\_correlation\_range');
xlabel('Correlation');
ylabel('Frequency');
subplot(2,3,6);
histogram(green_orange_correlation_all_range,num_bins,'FaceColor', [1 0 0]);
title('Histogram of green\_orange\_correlation\_range');
xlabel('Correlation');
ylabel('Frequency');
sgtitle('Green Channel Ball Correlation');

% Figure 6: Blue Channel
figure;
subplot(2,3,1);
histogram(blue_asm_all,num_bins,'FaceColor', [1 1 0]);
title('Histogram of blue\_yellow\_asm');
xlabel('ASM');
ylabel('Frequency');
subplot(2,3,2);
histogram(blue_white_asm_all,num_bins,'FaceColor', [0/255, 21/255, 71/255]);
title('Histogram of blue\_white\_asm');
xlabel('ASM');
ylabel('Frequency');
subplot(2,3,3);
histogram(blue_orange_asm_all,num_bins,'FaceColor', [1 0 0]);
title('Histogram of blue\_orange\_asm');
xlabel('ASM');
ylabel('Frequency');
subplot(2,3,4);
histogram(blue_asm_all_range,num_bins,'FaceColor', [1 1 0]);
title('Histogram of blue\_yellow\_asm\_range');
xlabel('ASM');
ylabel('Frequency');
subplot(2,3,5);
histogram(blue_white_asm_all_range,num_bins,'FaceColor', [0/255, 21/255, 71/255]);
title('Histogram of blue\_white\_asm\_range');
xlabel('ASM');
ylabel('Frequency');
subplot(2,3,6);
histogram(blue_orange_asm_all_range,num_bins,'FaceColor', [1 0 0]);
title('Histogram of blue\_orange\_asm\_range');
xlabel('ASM');
ylabel('Frequency');
sgtitle('Blue Channel Ball ASM');


figure;

% Subplot 1: Histograms of Red Channel Contrast Mean
subplot(2, 3, 1);
hold on;
num_bins1 = getBinNumber(red_contrast_all, red_white_contrast_all, red_orange_contrast_all);
histogram(red_contrast_all, num_bins1, 'FaceColor', [1 1 0], 'EdgeColor', 'none');
histogram(red_white_contrast_all+0.2, num_bins1, 'FaceColor', [0/255, 21/255, 71/255], 'EdgeColor', 'none');
histogram(red_orange_contrast_all+0.4, num_bins1, 'FaceColor', [1 0 0], 'EdgeColor', 'none');
hold off;
legend('Yellow Ball', 'White Ball', 'Orange Ball');
title('Histograms of Red Channel Contrast Mean');
xlabel('Contrast');
ylabel('Frequency');
% Subplot 2: Histograms of Green Channel Correlation
subplot(2, 3, 2);
hold on;
num_bins2 = getBinNumber(green_correlation_all, green_white_correlation_all, green_orange_correlation_all);
histogram(green_correlation_all, num_bins2, 'FaceColor', [1 1 0], 'EdgeColor', 'none');
histogram(green_white_correlation_all, num_bins2, 'FaceColor', [0/255, 21/255, 71/255], 'EdgeColor', 'none');
histogram(green_orange_correlation_all, num_bins2, 'FaceColor', [1 0 0], 'EdgeColor', 'none');
hold off;
legend('Yellow Ball', 'White Ball', 'Orange Ball');
title('Histograms of Green Channel Correlation Mean');
xlabel('Correlation');
ylabel('Frequency');
% Subplot 3: Histograms of Blue Channel ASM
subplot(2, 3, 3);
hold on;
num_bins3 = getBinNumber(blue_asm_all, blue_white_asm_all, blue_orange_asm_all);
histogram(blue_asm_all, num_bins3, 'FaceColor', [1 1 0], 'EdgeColor', 'none');
histogram(blue_white_asm_all, num_bins3, 'FaceColor', [0/255, 21/255, 71/255], 'EdgeColor', 'none');
histogram(blue_orange_asm_all, num_bins3, 'FaceColor', [1 0 0], 'EdgeColor', 'none');
hold off;
legend('Yellow Ball', 'White Ball', 'Orange Ball');
title('Histograms of Blue Channel ASM Mean');
xlabel('ASM');
ylabel('Frequency');
subplot(2, 3, 4);
hold on;
num_bins4 = getBinNumber(red_contrast_all_range, red_white_contrast_all_range, red_orange_contrast_all_range);
histogram(red_contrast_all_range, num_bins4, 'FaceColor', [1 1 0], 'EdgeColor', 'none');
histogram(red_white_contrast_all_range+0.2, num_bins4, 'FaceColor', [0/255, 21/255, 71/255], 'EdgeColor', 'none');
histogram(red_orange_contrast_all_range+0.4, num_bins4, 'FaceColor', [1 0 0], 'EdgeColor', 'none');
hold off;
legend('Yellow Ball', 'White Ball', 'Orange Ball');
title('Histograms of Red Channel Contrast Range');
xlabel('Contrast');
ylabel('Frequency');
subplot(2, 3, 5);
hold on;
num_bins5 = getBinNumber(green_correlation_all_range, green_white_correlation_all_range, green_orange_correlation_all_range);
histogram(green_correlation_all_range, num_bins5, 'FaceColor', [1 1 0], 'EdgeColor', 'none');
histogram(green_white_correlation_all_range, num_bins5, 'FaceColor', [0/255, 21/255, 71/255], 'EdgeColor', 'none');
histogram(green_orange_correlation_all_range, num_bins5, 'FaceColor', [1 0 0], 'EdgeColor', 'none');
hold off;
legend('Yellow Ball', 'White Ball', 'Orange Ball');
title('Histograms of Green Channel Correlation Range');
xlabel('Correlation');
ylabel('Frequency');
subplot(2, 3, 6);
hold on;
num_bins6 = getBinNumber(blue_asm_all_range, blue_white_asm_all_range, blue_orange_asm_all_range);
histogram(blue_asm_all_range, num_bins6, 'FaceColor', [1 1 0], 'EdgeColor', 'none');
histogram(blue_white_asm_all_range, num_bins6, 'FaceColor', [0/255, 21/255, 71/255], 'EdgeColor', 'none');
histogram(blue_orange_asm_all_range, num_bins6, 'FaceColor', [1 0 0], 'EdgeColor', 'none');
hold off;
legend('Yellow Ball', 'White Ball', 'Orange Ball');
title('Histograms of Blue Channel ASM Range');
xlabel('ASM');
ylabel('Frequency');
sgtitle('Histograms of Color Channel Features');

disp('Red Channel Contrast:');
disp('Mean Values:');
disp(['Yellow Ball: ', num2str(mean(red_contrast_all))]);
disp(['White Ball: ', num2str(mean(red_white_contrast_all))]);
disp(['Orange Ball: ', num2str(mean(red_orange_contrast_all))]);
disp('Range Values:');
disp(['Yellow Ball: ', num2str(mean(red_contrast_all_range))]);
disp(['White Ball: ', num2str(mean(red_white_contrast_all_range))]);
disp(['Orange Ball: ', num2str(mean(red_orange_contrast_all_range))]);

disp('Green Channel Correlation:');
disp('Mean Values:');
disp(['Yellow Ball: ', num2str(nanmean(green_correlation_all))]);
disp(['White Ball: ', num2str(nanmean(green_white_correlation_all))]);
disp(['Orange Ball: ', num2str(nanmean(green_orange_correlation_all))]);
disp('Range Values:');
disp(['Yellow Ball: ', num2str(mean(green_correlation_all_range))]);
disp(['White Ball: ', num2str(mean(green_white_correlation_all_range))]);
disp(['Orange Ball: ', num2str(mean(green_orange_correlation_all_range))]);

disp('Blue Channel ASM:');
disp('Mean Values:');
disp(['Yellow Ball: ', num2str(mean(blue_asm_all))]);
disp(['White Ball: ', num2str(mean(blue_white_asm_all))]);
disp(['Orange Ball: ', num2str(mean(blue_orange_asm_all))]);
disp('Range Values:');
disp(['Yellow Ball: ', num2str(mean(blue_asm_all_range))]);
disp(['White Ball: ', num2str(mean(blue_white_asm_all_range))]);
disp(['Orange Ball: ', num2str(mean(blue_orange_asm_all_range))]);
% 
function binNumber = getBinNumber(val_1, val_2, val_3)
    % Calculate the average value
    allValues = [val_1(:); val_2(:); val_3(:)];
    
    % Find the minimum and maximum values
    min_val = min(allValues);
    max_val = max(allValues);
    
    % Determine the number of bins
    Number = 15; % Assuming each bin has a range of 10 units
    binNumber = linspace(min_val, max_val, Number+1);
    % Display the determined bin number
    %disp(['The bin number for the input values is: ', num2str(binNumber)]);
end
