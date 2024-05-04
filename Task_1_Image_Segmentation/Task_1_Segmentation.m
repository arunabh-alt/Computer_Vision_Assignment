% Specify input and output folders
inputFolder = 'C:\Users\Computing\Downloads\OneDrive_1_04-05-2024\Color Image';
outputFolder = 'C:\Users\Computing\Downloads\OneDrive_1_04-05-2024\Segmented_Image';
% Check if the output folder exists, if not, create it
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
    disp('Output folder created.');
end
% Get a list of all PNG files in the input folder
fileList = dir(fullfile(inputFolder, '*.png'));
edgeThreshold = 0.2;
% Iterate over each image file
for i = 1:numel(fileList)
    % Read the image
    filePath = fullfile(inputFolder, fileList(i).name);
    RGB = imread(filePath);
    
    % Mask 1 - Find circles
    [centers, radii, ~] = imfindcircles(RGB, [10 75], 'ObjectPolarity', 'bright', 'Sensitivity', 0.90);
    BW1 = false(size(RGB, 1), size(RGB, 2));
    [Xgrid, Ygrid] = meshgrid(1:size(BW1, 2), 1:size(BW1, 1));
    BW1 = hypot(Xgrid - centers(1, 1), Ygrid - centers(1, 2)) <= radii(1);

    % Mask 2 - Convert RGB image to chosen color space
    I = rgb2lab(RGB);
    % Define thresholds for LAB channels
    channel1Min = 46.026;
    channel1Max = 78.974;
    channel2Min = 19.241;
    channel2Max = 50.203;
    channel3Min = 18.205;
    channel3Max = 64.115;
    
    % Create mask based on LAB thresholds
    sliderBW2 = (I(:,:,1) >= channel1Min) & (I(:,:,1) <= channel1Max) & ...
        (I(:,:,2) >= channel2Min) & (I(:,:,2) <= channel2Max) & ...
        (I(:,:,3) >= channel3Min) & (I(:,:,3) <= channel3Max);
    BW2 = sliderBW2;
    BW2 = edge(BW2, 'prewitt', edgeThreshold);
    se_ellipse = strel('line', 2, 60);
    BW2 = imdilate(sliderBW2, se_ellipse);
    se1= strel('square',10);
    BW2 = imopen(BW2,se1);
   
    se3 =strel('diamond',8);
    BW2 = imdilate(BW2, se3);
  
    % Mask 3 - Convert RGB image to YCbCr color space
    I = rgb2ycbcr(RGB);
    
     
    channel1Min = 120.000;
    channel1Max = 255.000;
    channel2Min = 0.000;
    channel2Max = 90.000;
    channel3Min = 0.000;
    channel3Max = 255.000;


    sliderBW3 = (I(:,:,1) >= channel1Min) & (I(:,:,1) <= channel1Max) & ...
        (I(:,:,2) >= channel2Min) & (I(:,:,2) <= channel2Max) & ...
        (I(:,:,3) >= channel3Min) & (I(:,:,3) <= channel3Max);
    BW3 = sliderBW3;
    BW3 = imfill(BW3, 'holes');
    BW3 = imfill(BW3, 'holes');
    se2 = strel('disk',2);
    BW3 = imdilate(sliderBW3, se2);
   
    % Resize masks to match the size of BW1
    if ~isequal(size(BW1), size(BW2))
        BW2 = imresize(BW2, size(BW1));
    end
    
    if ~isequal(size(BW1), size(BW3))
        BW3 = imresize(BW3, size(BW1));
    end
    
    % Combine masks
    combinedMask = BW1 | BW2 | BW3;
    se4 = strel('disk',1);
    combinedMask = imopen(combinedMask, se4);
   
    % Save the combined mask as segmented image
    outputFilePath = fullfile(outputFolder, [fileList(i).name(1:end-4), '_segmented.png']); % Adjust file extension if needed
    imwrite(combinedMask, outputFilePath);
end
