# semantic-segmentation
Semantic Segmentation Project 

This project contains the implementation for the paper:

B. R. A. Jaimes, J. P. K. Ferreira and C. L. Castro, *"Unsupervised Semantic Segmentation of Aerial Images With Application to UAV Localization,"* in **IEEE Geoscience and Remote Sensing Letters**, doi: 10.1109/LGRS.2021.3113878.

Instructions to run:

## Step 1

* Download the georeferenced map and place it in CnnComparison\Images\Train:
https://drive.google.com/file/d/1nYzS5sosIneXw0YELjQkVZUGy9UoYq7I/view?usp=sharing

* Download the map semantic segmentation and place it in PositionEstimation:
https://drive.google.com/file/d/1o2VmEuGHtU__y6IbbzXFngTRbVUd2o2a/view?usp=sharing

* Download the pre-trained networks and place in CnnComparison\TrainedNetworks:
https://drive.google.com/drive/folders/1FWG16mBCjY4j4f_P3YLe1gf3MASMDvbr?usp=sharing

## Step 2

Make sure your C++ compiler is set up in MATLAB
If not, run the following and configure:

`mex -setup C++`

Then go to PositionEstimation\SVM and run the `make.m` file to compile `libsvm` with the following:

`make`

## Step 3

This step is divided between each experiment from the paper.

Firt of all, in the main folder, run the initial function to add folders to path:

`initCode.m`

### Step 3.1

Semantic Segmentation:

In order to obtain the semantic segmentation of a single image, run the `SemanticSegmentation\semanticSegmentation.m` as follows:

`semanticSegmentation(imagePath, 1, 3);`

If you want to change the semantic segmentation parameters, all flags are explained whithin the code.

### Step 3.2

Position Estimation Experiment:

To replicate the paper experiments, for the Brazilian datast simply run:

`PositionEstimation\BR_experiment.m`

And for the Sweeden dataset, run:

`PositionEstimation\SW_experiment.m`

### Step 3.3

To run the CNN comparison experiment, run:

`CnnComparison\CNN_Experiment.m`

Enjoy!


(Early Access) Cite as:

@article{Jaimes2021,
  author={Jaimes, Brayan Rene Acevedo and Ferreira, João Pedro Klock and Castro, Cristiano Leite},
  journal={IEEE Geoscience and Remote Sensing Letters}, 
  title={Unsupervised Semantic Segmentation of Aerial Images With Application to UAV Localization}, 
  year={2021},
  volume={},
  number={},
  pages={1-5},
  doi={10.1109/LGRS.2021.3113878}}
