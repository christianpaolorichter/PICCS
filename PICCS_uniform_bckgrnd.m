function [paramEst,paramSE] = PICCS_uniform_bckgrnd(Ccum,r,verbose)
% PICCS_UNIFORM_BCKGRND Fits a theoretical model to the
% cumulative cross-correlation function (Ccum) to estimate spatial interaction
% parameters (clustering strength and cluster size) and background density.
%
% The model assumes the total correlation is the sum of a term for random
% (non-interactive) co-localization and an exponential term for clustered
% (interactive) co-localization.
%
%   [paramEst,paramSE] = evaluate_cross_corr_uniform_bckgrnd(Ccum,r,verbose)
%
%   Inputs:
%       Ccum:       The calculated cumulative cross-correlation values.
%       r:          The radial distances (radii) corresponding to Ccum.
%       verbose: Boolean flag (1 or 0) to display the fitting plots and results.
%
%   Outputs:
%       paramEst:   Vector of estimated parameters: [rho,alpha,epsilon]
%                   where: rho = Effective Background Density (per square micron),
%                          alpha = Interaction Strength/Fraction (unitless),
%                          epsilon = Localization Precision/Cluster size (microns).
%       paramSE:    Vector of standard errors for the estimated parameters.
%
%   Copyright (c) 2025 Christian Paolo Richter
%   University of Osnabrueck

% Ensure input vectors are row vectors for consistent use in polyfit and model evaluation.
r = rowvec(r);
Ccum = rowvec(Ccum);

%% define model (Cumulative Cross-Correlation Function, Ccum)

% Interactive/Clustering component (C_int): Exponential term modeling the
% correlation caused by molecules clustering together (Gaussian interaction volume).
% param(2) is the Interaction Strength (alpha)
% param(3) is the Cluster Size/Localization Precision (epsilon)
interactiveCorr = @(param,r)param(2)*(1-exp(-r.^2/(4*param(3)^2)));

% Non-Interactive/Background component (C_bckg): Proportional to the search area
% (pi * r^2), representing random co-localization across a uniform background.
% param(1) is the Effective Background Density (rho)
nonInteractiveCorr = @(param,r)param(1)*pi*r.^2;

% Combined Model: Ccum is the sum of the interactive and non-interactive terms.
model = @(param,r)feval(interactiveCorr,param,r) + ...
    feval(nonInteractiveCorr,param,r);

%% initial guess (x0)

% Select data points at larger radii where the random background term dominates the correlation.
take = r > r(floor(numel(r)/2));

% Use linear regression (polyfit) on Ccum versus pi*r^2 for the background-dominated region.
% The slope of this fit estimates the effective background density term, rho (param(1)).
x0 = polyfit(pi*r(take).^2,Ccum(take),1); % Initial estimates [slope (rho), intercept (alpha)]

% Assemble the initial guess vector x0 = [rho,alpha,epsilon]:
x0 = [x0 0.02]; %[rho alpha epsilon]
% Define bounds for optimization parameters [rho,alpha,epsilon]
lb = [0.001 0 0.001]; % Lower bounds
ub = [10 1 1];       % Upper bounds

%% Non-linear Least Squares Fitting (lsqcurvefit)

% Set optimization tolerances and limits.
fitOptions = ...
    optimset(...
    'TolFun', 10^-12,...    % Termination tolerance on the function value
    'TolX', 10^-12,...      % Termination tolerance on the parameter values
    'MaxIter', 1000,...     % Maximum number of iterations
    'MaxFunEvals', 3000,... % Maximum number of function evaluations
    'Display', 'off');  

% Perform the non-linear least squares fit using the defined model.
[paramEst,~,resid,~,~,~,J] = ...
    lsqcurvefit(model,x0,r,Ccum,lb,ub,fitOptions);

% Calculate the standard errors (SE) for the estimated parameters.
% The 95% Non-linear Parameter Confidence Interval (CI) width is divided by 3.92 (2 * 1.96).
paramSE = rowvec(diff(nlparci(paramEst,resid,'jacobian',J),[],2)/3.92);

% --- Calculate Goodness-of-Fit (R-squared) ---
SST = sum(bsxfun(@minus,Ccum,mean(Ccum)).^2);     % Total Sum of Squares
SSR = sum((Ccum-model(paramEst,r)).^2);           % Sum of Squares of Residuals
rSquare = 1-SSR./SST;                             % R-squared value

%% Display Results (Optional)

if verbose
    figure('Color',[1 1 1]);
    ax(1) = axes(...
        'Position',[0.15 0.4 0.8 0.45],...
        'XTickLabel','',...             % Hide X-ticks on top plot
        'NextPlot','add');

    % Plot raw data points (Ccum vs r^2)
    plot(ax(1),r.^2,Ccum,'ko')
    % Plot the total fitted model
    plot(ax(1),r.^2,feval(model,paramEst,r),'y-','Linewidth',3)
    % Plot the Interactive/Clustering component (C_int)
    plot(ax(1),r.^2,feval(interactiveCorr,paramEst,r),'r--','Linewidth',2);
    % Plot the Non-Interactive/Background component (C_bckg)
    plot(ax(1),r.^2,feval(nonInteractiveCorr,paramEst,r),'b--','Linewidth',2);
    grid on
    box on

    % Title showing the estimated parameters and R-squared.
    title(ax(1),sprintf('Interaction [x100%%] = %.2e\\pm%.2e;\nParticle Density [µm^{-2}] = %.2e\\pm%.2e (R^2 = %.2f)\nLoc. Precision [µm] = %.1e\\pm%.1e',...
        paramEst(2), paramSE(2), paramEst(1), paramSE(1), rSquare, paramEst(3), paramSE(3)))
    ylabel(ax(1),'C_{cum}')

    % --- Residuals Plot (Bottom axes) ---
    ax(2) = axes(...
        'Position',[0.15 0.15 0.8 0.2],...
        'NextPlot','add');
    % Plot the residuals (Data - Fit)
    stem(ax(2),r.^2,resid,'k')
    ylabel(ax(2),'Residuals')
    xlabel(ax(2),'r^2 [µm^2]')
    grid on
    box on

    % Finalize plot appearance: adjust limits and link X-axes.
    axis(ax,'tight')
    linkaxes(ax,'x')
end %if
end %fun