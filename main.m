clear all;  % just a precaution
addpath('include');

% protein variables
hydrophobicity = [0, 0, 1, 0, 1, 1, ...
                  0, 1, 0, 0, 0, 0, ...
                  1, 1, 0, 0, 1, 1, ...
                  0, 0, 0, 1, 0, 1, ...
                  1, 0, 0, 0, 1, 1, ...
                  0, 0, 1, 1, 0, 1];       % let 0 be P residue and let 1 be H residue
xvals = 1:length(hydrophobicity);          % initial x-values of every residue in polypeptide
yvals = zeros(1, length(hydrophobicity));  % initial y-values of every residue in polypeptide

% folding variables
temperature = 1;                           % temperature
maxiter = 10^7;                            % total no. of iterations of Metropolis algorithm
energyArray = zeros(1, maxiter+1);         % energy at the end of all trial moves

% start the timer: takes <300 seconds (<5 minutes) with 10^7 iterations
tic;
[xvals, yvals, energyArray] = fold(hydrophobicity, xvals, yvals, energyArray, temperature, maxiter);
toc;
