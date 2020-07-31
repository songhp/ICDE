# ICDE
This is a MATLAB packege for Cosparse Signal Recovery via Iterative Cosupport Detection-Estimation (ICDE) Algorithm


This packeage is written and maintained by Heping Song
Email: songhp@usj.edu.cn

## Requirements and Dependencies
- MATLAB (we test with MATLAB R2018b on Windows 7)
- [CVX](https://github.com/cvxr/CVX)

## Installation

	% Start MATLAB
	% First, set up CVX package
	cd cvx
	cvx_setup
	cd ..
	% select SeDuMi solver for optimization
	cvx_solver sedumi


## Demo

To test ICDE:

    >> demo_ICDEL1
    >> demo_ICDEL2

