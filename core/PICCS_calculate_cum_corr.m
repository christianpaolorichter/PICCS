function Ccum = PICCS_calculate_cum_corr(A,B,r)
% PICCS_CALCULATE_CUM_CORR Calculates the cumulative cross-correlation function (Ccum)
% between two sets of spatial localizations, A and B, up to a maximum radius r.
%
% This function is often used in localization microscopy (e.g., PALM/STORM)
% to quantify the degree of spatial co-localization or clustering between
% two distinct molecular species (A and B).
%
%   Ccum = calculate_cum_corr(A,B,r)
%
%   Inputs:
%       A:      N_A x 2 matrix of localizations for species A.
%       B:      N_B x 2 matrix of localizations for species B.
%       r:      M x 1 vector of radii (distance thresholds) at which to calculate the correlation.
%
%   Output:
%       Ccum:   M x 1 vector where Ccum(i) is the number of localizations of
%               species A that are within radius r(i) of *any* localization
%               of species B, normalized by the total number of B localizations.
%
%   Copyright (c) 2025 Christian Paolo Richter
%   University of Osnabrueck

% Check for empty input data sets.
if isempty(A) || isempty(B)
    % If either set is empty, the cumulative correlation is zero for all radii.
    Ccum = zeros(size(r));
else
    % --- 1. Calculate All-Pairs Distances ---

    % Calculate the Euclidean distances between every localization in A and every localization in B.
    % distMat will be an N_A x N_B matrix.
    distMat = pdist2(A,B); % euclidean distances between localizations of [A] and [B]

    % --- 2. Calculate Correlation and Normalization ---

    % The complex expression calculates the number of A localizations that are
    % within a certain radius of B localization, normalized.
    Ccum = sum(cumsum(histc(distMat,[0 rowvec(r)]+eps,1),1),2)/size(B,1); % [AB]/[B] & [A]

    % Remove the last element, which corresponds to the overflow bin (distances > r(end)).
    Ccum(end) = [];
end %if
end %fun