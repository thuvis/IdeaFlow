Introduction
============

This directory includes sources used in the following paper:

    Tracking Idea Flows between Social Groups
    Yangxin Zhong, Shixia Liu, Xiting Wang, Jiannan Xiao, Yangqiu Song
    In Proceedings of AAAI, 2016
    
This code has been tested under 64-bit Windows environment using Eclipse 4.3.0 with JRE6.

In this directory, we include part of the implementation of our method: an augmented algorithm 
combining Bayesian Conditional Cointegration and Dynamic Time Warping to calculate the lead-
lag relationship of word pair between two social group, using their time series of the word 
frequency. This part of implementation is supposed to run first in our method, that is to say, 
the output of this part will be the input of the other part of our implementation (the MATLAB code).

File Description
================

GaussianIntegral.java        : The code in this file is used to calculate integral for Gaussian function 
                               limited in a certain interval. This part is necessary when we calculate
                               the Bayesian Conditional Cointegration between two time series.

BayesianCointegration.java   : This file implement Bayesian Conditional Cointegration, Dynamic Time Warping 
                               and the augmented algorithm combining them presented in our paper. This file
                               contain a few important function, so we would like to introduce more about them.
                               You can also refer to the comments in this file for more details.
                               
        public double[] calculateBayesianCointegration(double x[],double y[]):
                The function to calculate whether series x and y cointegrated.
                
        public ResultOfDTW DTW_MeanNormalized(double x[], double y[], int normalizeFlag, int trendFlag):
                Align series x and y with DTW.
                
        public double[] cointegration_DTW_MeanNormalized(double x[], double y[]),
        public double[] cointegration_DTW_MeanNormalized(double x[], double y[], int DTW_win, double gamma_thr):
                Augmented BCC, combined with DTW.
                
        public void outputIndexMatrix(String fileName, double[][][] leadLagTime, int taoMax) throws IOException:
                Output an index matrix according to the leadLagTime tensor between words of two social group.
                This output file can be used as the input file for the tensor-based clustering algorithm (MATLAB code).
                 
        public static void main(String[] args) throws IOException:
                This function gives examples about how to use the functions above.


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