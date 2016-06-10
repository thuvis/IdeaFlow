Introduction
============

This directory includes sources used in the following paper:

    Tracking Idea Flows between Social Groups
    Yangxin Zhong, Shixia Liu, Xiting Wang, Jiannan Xiao, Yangqiu Song
    In Proceedings of AAAI, 2016
    
This code has been tested under 64-bit Windows environment using MATLAB 8.4.0.150421 (R2014b).

In this directory, we include part of the implementation of our method: tensor-based algorithm 
to cluster words into ideas and aggregate the word-level relationships into idea flows (with the 
third kind of tensor representation in our paper). This part of implementation can be run after 
calculating the lead-lag relationship of each word pair between two social groups (using the other 
part of our implementation).


Dependency Package
==================

Before running the code, you have to install a MATLAB Tensor Toolbox, which was released by Tamara 
G. Kolda, et al. In our experiment, we use version 2.6 (released Feb. 6, 2015). The website to 
download the toolbox is as below:

    http://www.sandia.gov/~tgkolda/TensorToolbox/index-2.6.html
    
    
File Description
================

ideaFlow.m  : The implementation of tensor-based clustering algorithm. Please refer to the comments 
              in this source file for more details about its usage (including input and output).
            
main.m      : This file includes a simple example to illustrate how to use our implementation. 
              Detailed guides for the example are also included in the comments in this file.
              
test.txt    : This file is the index matrix for the word-level lead-lag relationships between two 
              social groups. It is used as the input file for main.m and we will descripe more about
              it in the next section.
              
test.jpg    : This file is used to show the ground truth of data in test.txt. It tells about which words 
              in each group belong to the same idea as well as when and how the ideas flow to each other. 
              After running main.m, you can compare the result with this file.


Test Data Description
=====================

In the example of main.m, we need to input an index matrix for the sparse tensor of word-level 
lead-lag relationships. The index matrix is in the file "test.txt" in this directory and we would
like to illustrate more about its format and meaning here.

The first line of the index matrix gives the information about the size of the sparse tensor: 
"7 10 12 5" means this is a 4-order tensor with a size of 7 x 10 x 12 x 5, where 7 is the number 
of words in the first social group, 10 is the number of words in the second social group, 12 is 
the number of timepoints, and 5 means there are five possible integer values for lead-lag time 
(from -2 to +2, which means taoMax == 2).

The rest part of the matrix gives the information about the values in the sparse tensor: every
line means a 1 at the corresponding position of the tensor. For instance, "1 3 4 2" in the 13th 
line indicates there is a 1 at the index of (1, 3, 4, 2) in the 4-order tensor, which means 
the 1st word of the first social group correlated with the 3rd word of the second social group 
at the 4th timepoints with a lead-lag time -1. Here lead-lag time is calculated by adding a 
-(taoMax+1) on the last index. In this case, since taoMax == 2, the last index 2 means lead-lag 
time 2 - (2+1) == -1. If a certain position in the tensor is not included in the index matrix, 
it indicates a 0 at the corresponding place.

Since the lead-lag relationship of word pair is calculated by an augmented algorithm combining 
Bayesian Conditional Cointegration and Dynamic Time Warping, the sparse tensor actually indicates 
the cointegration relationshiip between two words along the time axis. For instance, since the 
lead-lag time in the case "1 6 7 5" is 2, it means the 7th timepoint in the time series of the 
1st word in the first social group cointegrate with the 7+2 == 9th timepoint in the time series 
of the 6th word in the second social group. And in the case of "1 3 4 2", since the lead-lag time  
is -1, which is a NEGATIVE number, it means the 4th timepoint in the time series of the 3rd word 
in the SECOND social group cointegrate with the 7+|-1| == 8th timepoint in the time series of the 
1st word in the FIRST social group. In a word, the sign of the lead-lag time indicates which social
group leads the flow and the absolute value means the propogation delay time in this flow.

The ground truth of test.txt is in test.jpg, including the which idea a word belong to as well as 
when the flows happen and the idea-level lead-lag times.


Additional Information
======================

If you find this tool useful, please cite it as

    Tracking Idea Flows between Social Groups
    Yangxin Zhong, Shixia Liu, Xiting Wang, Jiannan Xiao, and Yangqiu Song
    In Proceedings of AAAI, 2016

The bibtex format is

    @inproceedings{zhong2016,
        title = {Tracking Idea Flows between Social Groups},
        author = {Yangxin Zhong and Shixia Liu and Xiting Wang and Jiannan Xiao and Yangqiu Song},
        booktitle = {AAAI},
        year = {2016}
    }

For any question, please contact Yangxin Zhong <zhongyx4869@163.com>.