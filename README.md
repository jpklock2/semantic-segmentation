# UPDATE 21/07/2021

Updating project to latest version, please wait...



# semantic-segmentation
Semantic Segmentation Project 

Instructions to run:

## Step 1

* Unzip the train image in Images/Train (twice!)
* Unzip the train image class file in Images/Train/Classes

## Step 2

Make sure your C++ compiler is set up in MATLAB
If not, run the following and configure:

`mex -setup C++`

Then go to Classifier/SVM and run the `make.m` file to compile `libsvm` with the following:

`make`

## Step 3

This step is divided between the Semantic Segmentation and Position Estimation pipelines.

### Step 3.1

Semantic Segmentation Pipeline:

In the main folder, edit the `main.m` file flags for your use case.
All flags are explained whithin the code.
Then run the project with:

`main`

### Step 3.2

Position Estimation Pipeline:

In the main folder, edit the `semanticSegmentation.m` file flags for your use case.
All flags are explained whithin the code.
Then run the project with:

`positionEstimation`


Enjoy!
