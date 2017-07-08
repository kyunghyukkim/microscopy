# Cell Lineage Tracking and Time Delay Correlation

Copyright (c) 2017 Kiri Choi

## Introduction
Python code to track cell lineage and apply time delay correlation using the output of [SuperSegger](https://github.com/wiggins-lab/SuperSegger/), a MATLAB-based high-throughput image cell segmentation software suite.

The input file is an output of cell segmentation process done by SuperSegger software called `uclist.mat`. A short description of what the code does is provided in each file. 

## Dependencies
The script requires [scipy stack](https://www.scipy.org/), [seaborn](https://seaborn.pydata.org/), and [h5py](http://www.h5py.org/)