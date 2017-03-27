function NuOut = PhiExp(NuIn, Params)
%
% NuOut = PhiExp(NuIn, Params)
%
% Params.TauARP % Absolute refractory period.
% Params.TauV % Membrane potential decay constant.
% Params.Crec % average total synaptic coupling: MuRec = Crec*NuIn.
% Params.Mu0  % rheobase current
Params.Sigma0  % current fluctuation size.
%x²
%   Copyright 2010 Maurizio Mattia @ Ist. Super. Sanit�, Rome - Italy
%   Version: 1.0 - Sep. 23, 2010
%


Mu = Params.Crec * NuIn + Params.Mu0;
a = (0 - Mu * Params.TauV) / (Params.Sigma0 * sqrt(Params.TauV));
b = (1 - Mu * Params.TauV) / (Params.Sigma0 * sqrt(Params.TauV));

ndxPos = find(b<=5);

NuOut = zeros(size(NuIn));

for n = ndxPos
%    disp(NuIn(n));
   NuOut(n) = 1 / (Params.TauV * sqrt(pi) * quad(@auxPhiExp, a(n), b(n)) + Params.TauARP);
end


function out = auxPhiExp(w)

for i = 1:length(w)
   if w(i)>=0
      out(i) = exp(w(i)^2)*(1 + erf(w(i)));
   else
      str = ['exp(-x.^2).*x./sqrt((' num2str(w(i)) ').^2 + x.^2)'];
      out(i) = 2/sqrt(pi)*quad(str, 0, 10);
   end
end