# Computer_Vision_Assignment

This assignment is consist of primarily 3 tasks.
## Task 1 
Task Statment - Imaage Segmentation and Compare the segmented image with Ground Truth. Determine the Dice score Mean and Standard Deviation.
Task Result - 
              
              Best 5 Dice Scores

              Image 55: 0.9095
              Image 54: 0.8801
              Image 27: 0.8798
              Image 26: 0.8792
              Image 46: 0.8786
    
              Worst 5 Dice Scores:
              Image 16: 0.5569
              Image 15: 0.5840
              Image 18: 0.5995
              Image 14: 0.6187
              Image 17: 0.6191
              
              Mean Dice Score: 0.8049
              Standard Deviation of Dice Scores: 0.0861
## Task 2
Task Statment - From the RGB Images and Ground Truth Determine the Shape Features and Texture Features for the Ball Patches. 
Task Result -  Shape Feature Result (Part A)

                Mean Yellow Solidity: 0.98
                Mean Yellow Non-Compactness: 1.15
                Mean Yellow Circularity: 0.97
                Mean Yellow Eccentricity: 0.46
                Mean White Solidity: 0.88
                Mean White Non-Compactness: 1.60
                Mean White Circularity: 0.58
                Mean White Eccentricity: 0.61
                Mean Orange Solidity: 0.79
                Mean Orange Non-Compactness: 1.61
                Mean Orange Circularity: 0.50
                Mean Orange Eccentricity: 0.71
Task Result - Textural Features Result (Part B)

              Red Channel Contrast:
                Mean Values:
                Yellow Ball: 5.8597
                White Ball: 0.75889
                Orange Ball: 0.11612
                Range Values:
                Yellow Ball: 9.1622
                White Ball: 1.946
                Orange Ball: 0.20066
              Green Channel Correlation:
                Mean Values:
                Yellow Ball: 0.97125
                White Ball: 0.88523
                Orange Ball: 0.96211
                Range Values:
                Yellow Ball: 0.10129
                White Ball: 0.25243
                Orange Ball: 0.10097
              Blue Channel ASM:
                Mean Values:
                Yellow Ball: 0.20043
                White Ball: 0.29113
                Orange Ball: 0.23043
                Range Values:
                Yellow Ball: 0.25866
                White Ball: 0.38031
                Orange Ball: 0.33471
## Task 3 
Task Statement - Implement Kalman Filter for Given Ground Truth Coordinates and Noisy coordinates of a Moving Ball. For given matrix calculate the RMSE and standard deviation and compare with finetuned Kalman filter. 
Task Result - 
- Best Finetuned RMSE: 4.9627
- Best Finetuned Q matrix:
  
      0.0100         0         0         0
           0    0.0100         0         0
           0         0    0.0100         0
           0         0         0    0.0100

- Best Finetuned R matrix:

      0.1000         0
             0    0.1000

- Baseline Filter RMSE between true and estimated coordinates Mean : 168.4126, Standard Deviation: 83.0663
- Finetuned Filter RMSE between true and estimated coordinates Mean : 4.6792, Standard Deviation: 2.1713
- Baseline Kalman Filter Standard Deviation Error in X direction Error: 82.7309, Y Direction Error : 11.0045
- Fine Tuned Kalman Filter Standard Deviation Error in X direction Error: 3.6485, Y Direction Error : 3.0581
