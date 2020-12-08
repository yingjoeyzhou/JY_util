function W = trainfun_pimBH(cfg, X, phi)
% INPUT:
%   cfg.numC        : number of hypothetical orientation channels
%   cfg.tuning_curve: function handle
%   cfg.gamma       : regularization parameter.
%   X               : sensor-by-trial matrix.
%   phi             : 1-by-trial orientation labels.
%
% OUTPUT:
%   W: an nSensor-by-nOriChannel matrix.
% 
% JY

numF = size(X, 1);      % Number of features
numN = size(X, 2);      % Number of trials
numC = cfg.numC;        % Number of channels

% Construct hypothetical channel responses
C = zeros(numC, numN);

for in = 1:numN
    C(:, in) = cfg.tuning_curve((0:(numC-1))' * (pi/numC) - phi(in)*(pi/180));
end

C = C - repmat(mean(C, 2), [1, numN]);

% Iterate over channels
W = zeros(numF, numC);

for ic = 1:numC
    % Estimate spatial pattern
    l = (X*C(ic, :)')/(C(ic, :)*C(ic, :)');
    
    % Calculate noise
    eta = X - l*C(ic, :);

    % Estimate noise covariance
    S = eta*eta' / (numN-1);

    % Regularize
    if (isfield(cfg, 'gamma'))
        nu = trace(S)/numF;
      
        S = (1-cfg.gamma)*S + cfg.gamma*nu*eye(numF);
    end
        
    % Calculate filter
    W(:, ic) = S\l;
    W(:, ic) = W(:, ic) / (l' * W(:, ic));
end    
    
end    
    