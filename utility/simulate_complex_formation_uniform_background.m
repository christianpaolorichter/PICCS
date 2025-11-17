function [listA,listB,PrcntComplex] = ...
    simulate_complex_formation_uniform_background(Area,MolDens,PrcntComplex,Sigma,CorrLimit)
% SIMULATE_UNIFORM Generates synthetic localization microscopy data for two
% species (A and B) that are primarily distributed uniformly but contain a
% specified percentage of binary complexes (A-B pairs).
%
% The positions are simulated with localization noise (Sigma) and the
% boundary effects are mitigated by removing localizations near the edge.
%
%   [listA,listB,PrcntComplex] = simulate_uniform(Area,MolDens,PrcntComplex,Sigma,CorrLimit)
%
%   Inputs:
%       Area:           Scalar, the total area of the simulation region (in µm^2).
%       MolDens:        1x2 vector, the molecular density [Density_A, Density_B] (molecules per µm^2).
%       PrcntComplex:   Scalar, the target percentage of species B that
%       forms a complex with A. [x100%]
%       Sigma:          1x2 vector, the localization precision [Sigma_A, Sigma_B] (in µm).
%       CorrLimit:      Scalar, the maximum correlation length used for boundary exclusion.
%
%   Outputs:
%       listA:          N_A x 2 matrix of final simulated localizations for species A.
%       listB:          N_B x 2 matrix of final simulated localizations for species B (after boundary exclusion).
%       PrcntComplex:   The actual, discrete percentage of B molecules that formed a complex.
%
%   Copyright (c) 2025 Christian Paolo Richter
%   University of Osnabrueck

% --- 1. Determine Total Number of Molecules ---

% Determine the total number of molecules A & B
nTotal = ceil(MolDens*Area); % Calculate total count: [N_A, N_B]

% Calculate the number of binary complexes AB/B
% The number of complexes is based on the target percentage of the less abundant species (B).
nAB = ceil(PrcntComplex*nTotal(2));

% Calculate the number of free molecules A_ & _B
% Free molecules A_ (A not complexed) and _B (B not complexed).
nA_ = nTotal(1)-nAB;
nB_ = nTotal(2)-nAB;

% Calculate the true (discrete) percentage of B molecules that are complexed.
PrcntComplex = nAB/nTotal(2); %true percantage complex due to discrete number of particles

% --- 2. Draw Start Positions (Ideal Uniform Distribution) ---

% Draw Uniform distribution of the starting positions for A_ & _B & AB
% Initial positions (r0) are drawn uniformly in the area sqrt(Area) x sqrt(Area).
r0A_ = rand([nA_,2])*sqrt(Area); %[um] - Free A positions
r0B_ = rand([nB_,2])*sqrt(Area); %[um] - Free B positions
r0AB = rand([nAB,2])*sqrt(Area); %[um] - Complex core positions (common to A and B in the complex)

% --- 3. Apply Localization Precision (Simulating Measurement Noise) ---

% Calculate molecule positions considering the
% shift due to finite localization precision of A & B
% Positions are shifted from the ideal center by a random amount drawn from a normal distribution (mu=0, sigma=Sigma).

% Apply noise to free A and the A component of complexes
r0A_ = r0A_ + randn([nA_,2])*Sigma(1); % Free A localizations
r0A0 = r0AB + randn([nAB,2])*Sigma(1); % Complex A localizations (shifted from complex core r0AB)
listA = [r0A_;r0A0]; % Combine free and complex A localizations

% Apply noise to free B and the B component of complexes
r0B_ = r0B_ + randn([nB_,2])*Sigma(2); % Free B localizations
r0B0 = r0AB + randn([nAB,2])*Sigma(2); % Complex B localizations (shifted from complex core r0AB)
listB = [r0B_;r0B0]; % Combine free and complex B localizations

% --- 4. Boundary Exclusion ---

% Reduce src area according to the max. correlation length to evade boundary
% effects (=infinite boundaries)
% Filter B localizations that are too close to the boundary (x or y < CorrLimit or x or y > sqrt(Area) - CorrLimit).
% This prevents artificial correlation/clustering caused by edge truncation (common method for simulating infinite boundaries).
good = all(listB > CorrLimit,2) & all(listB < sqrt(Area)-CorrLimit,2);
listB = listB(good,:);
end %fun