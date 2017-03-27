function NuOut = PhiLin(NuIn, Params)
%
% NuOut = PhiLin(NuIn, Params)
%
% Params.TauARP % Absolute refractory period.
% Params.Crec   % average total synaptic coupling: MuRec = Crec*NuIn.
% Params.Mu0    % rheobase current
% Params.Sigma0 % current fluctuation size.
%
%   Copyright 2010 Maurizio Mattia @ Ist. Super. Sanità, Rome - Italy
%   Version: 1.0 - Sep. 24, 2010
%


Mu = Params.Crec * NuIn + Params.Mu0;

Csi = Mu/Params.Sigma0;

NuOut = 1./(1./(2*Csi.*Mu).*(exp(-2.*Csi)-exp(-2*0.*Csi)) + (1-0)./Mu + Params.TauARP);