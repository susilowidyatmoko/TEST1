% MAIN_DEMO Demonstration script for Simplified MITC4+ Shell Element MATLAB Module
%
% Paper Reference:
% "The simplified MITC4+ shell element and its performance in linear and nonlinear analysis"
% Hyung-Gyu Choi, Phill-Seung Lee (Computers & Structures 290, 2024, 107177)

clear; clc; close all;

fprintf('========================================================================\n');
fprintf('  SIMPLIFIED MITC4+ SHELL ELEMENT MATLAB IMPLEMENTATION DEMO\n');
fprintf('  Computers & Structures (2024)\n');
fprintf('========================================================================\n\n');

% Add paths
addpath('src');
addpath('utils');
addpath('benchmarks');

%% Step 1: Basic Numerical Tests (Section 4)
test_basic_tests();

%% Step 2: Linear Benchmark Problems (Section 5)
run_linear_benchmarks();

%% Step 3: Geometric Nonlinear Benchmark Problems (Section 5)
run_nonlinear_benchmarks();

fprintf('\n========================================================================\n');
fprintf('  DEMO COMPLETE: All linear, nonlinear, and benchmark tests executed!\n');
fprintf('========================================================================\n');
