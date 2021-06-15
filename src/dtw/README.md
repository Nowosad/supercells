# DTW_cpp
A small Dynamic Time Warping (DTW) single header library for C++

[DTW.hpp](https://github.com/cjekel/DTW_cpp/blob/master/include/DTW.hpp) computes the DTW distance between two c++ vectors ```a``` and ```b```! 

[![Build Status](https://travis-ci.com/cjekel/DTW_cpp.svg?branch=master)](https://travis-ci.com/cjekel/DTW_cpp) [![Coverage Status](https://coveralls.io/repos/github/cjekel/DTW_cpp/badge.svg?branch=master)](https://coveralls.io/github/cjekel/DTW_cpp?branch=master)

# Features

- Supports N-Dimensional data
- ```a``` and ```b```can have different number of data points
- Compute the distance using any ```p```-norm
- ```DTW::dtw_distance_only(a, b, p);``` function which returns only the DTW distance
- ```DTW::DTW MyDtw (a, b, p);``` class contains the pairwise distance vector, DTW distance vector, DTW distance, and a function to calculate the DTW alignment path

# What is Dynamic Time Warping ? 

Dynamic Time Warping (DTW) is an algorithm to measure the similarity between two temporal curves. The wiki page on [DTW](https://en.wikipedia.org/wiki/Dynamic_time_warping) is a great place to learn more.

![Image of two different curves](https://raw.githubusercontent.com/cjekel/similarity_measures/master/images/TwoCurves.png)

Consider the above Numerical and Experimental curves in 2D space. DTW can be used to measure the similarity between the two curves. A DTW distance of zero would mean that the warped curves match exactly.

The order of data points matters. Each curve is a sequence of data points, with a known beginning and ending.

# Examples

Check out the two [examples](https://github.com/cjekel/DTW_cpp/tree/master/examples).

# Tests

Run [run_tests.sh](https://github.com/cjekel/DTW_cpp/blob/master/run_tests.sh) in a linux environment. 
- travisci tests using Ubuntu Xenial and g++ version 5.4.0
- also tested on openSUSE Leap 15.1 and g++ version 7.4.0

# Requirements

- C++11 standard or later