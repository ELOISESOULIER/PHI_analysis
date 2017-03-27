function [NuIn, NuOut] = searchAndPlotGainFunction(TauARP,Crec,Sigma0,Nu0)
%
%   [NuIn, NuOut] = searchAndPlotGainFunction(TauARP,Crec,Sigma0,Nu0)
%
%   Copyright 2010 Maurizio Mattia @ Ist. Super. Sanit�, Rome - Italy
%   Version: 1.0 - Sep. 24, 2010
%

NU_RANGE = [0 45];

% ORANGE... Alpha1 = 0.97
% Params.TauARP = 1/50; % Absolute refractory period.
% Params.Crec = 7.15; % average total synaptic coupling: MuRec = Crec*NuIn.
% Params.Mu0 = - Params.Crec * 30;  % rheobase current
% Params.Sigma0 = 83.0;

% BLUE... Alpha1 = 0.97
% Params.TauARP = 1/40; % Absolute refractory period.
% Params.Crec = 91; % average total synaptic coupling: MuRec = Crec*NuIn.
% Params.Mu0 = - Params.Crec * 70;  % rheobase current
% Params.Sigma0 = 1300;

% ... Alpha1 = 0.80
% Params.TauARP = 1/40; % Absolute refractory period.
% Params.Crec = 5.87; % average total synaptic coupling: MuRec = Crec*NuIn.
% Params.Mu0 = - Params.Crec * 30;  % rheobase current
% Params.Sigma0 = 67;

% Mavi's example...
%
% 1/100,1.5,6,10
% Params.TauARP = 1/100; % Absolute refractory period.
% Params.Crec = 1.5; % average total synaptic coupling: MuRec = Crec*NuIn.
% Params.Mu0 = - Params.Crec * 10;  % rheobase current
% Params.Sigma0 = 6;

Params.TauARP = TauARP; % Absolute refractory period.
Params.Crec = Crec; % average total synaptic coupling: MuRec = Crec*NuIn.
Params.Mu0 = - Params.Crec * Nu0;  % rheobase current
Params.Sigma0 = Sigma0;

NuIn = linspace(NU_RANGE(1), NU_RANGE(2), diff(NU_RANGE)/2+1);

NuOut = PhiLin(NuIn, Params);

% figure
% hold on
% 
% plot(NuIn, NuOut, '.-');
% hold on
% xlabel('\nu');
% ylabel('\Phi(\nu)');
% set(gca, 'XLim', NU_RANGE, 'YLim', NU_RANGE);

csPhi = csape(NuIn, NuOut);
% fnplt(csPhi, 'r')

k = 0;
for Nu0 = [15.9 22.5]
   k = k + 1;
   
   zrs = fnzeros(csape(NuIn, NuOut-Nu0),NU_RANGE);

   Nu(k) = zrs(1);
   fprintf('Phi(%.3g) = %.3g', Nu(k), fnval(csPhi,Nu(k)));
   DPhiDNu(k) = fnval(fnder(csPhi), Nu(k));
   fprintf('Phi''(%.3g) = %.3g', Nu(k), DPhiDNu(k));
end

disp(['Phi''_1 / Phi''_2 = ' num2str(DPhiDNu(1)/DPhiDNu(2),5)])

% figure
% fnplt(fnder(csPhi))
% xlabel('\nu');
% ylabel('\Phi''(\nu)')


PLOT_SAMPLES = 20;
% DAMPING_FACTOR = 0.05;
% DAMPING_FACTOR = zeros(Net.P,1) + 0.05;
% DAMPING_FACTOR(Net.ndxI) = 0.25;
SEARCH_TOLERANCE = 0.01;
FINITE_SIZE_NOISE = 1;

%
% Fits the effective rate-to-rate transfer function...
%
p = polyfit(NuIn, NuOut, 5);
X = linspace(NuIn(1), NuIn(end), 5*PLOT_SAMPLES);
Y = polyval(p, X);

%
% Search for fixed points...
%
q = p;
q(end-1) = q(end-1) - 1;
s = roots(q)';
r = [];
for k = 1:length(s)
   if isreal(s(k))
      if s(k) >= NuIn(1) & s(k) <= NuIn(end)
         r = [r s(k)];
      end
   end
end

%
% Computes an effective potential energy...
%
qint = q./ (length(q):-1:1);
qint = [qint -(NuIn(1)-NuOut(1))*diff(NuIn(1:2))];
Yint = -polyval(qint, X);


%
% Plots the effective rate-to-rate transfer function and related
% functions...
%
figure


subplot(2,2,1);
hold on

plot(NuIn, NuOut, 'r.', X, Y, 'r', NuIn, NuIn, 'k--');
plot(r, r, 'bo', 'MarkerFaceColor', 'w');

for k = 1:length(r)
   text(r(k), r(k), ['  ' num2str(r(k),3) ' Hz']);
end

set(gca, 'XLim', [NuIn(1) NuIn(end)]);
xlabel('\nu_{in} (Hz)');
ylabel('\nu_{out} (Hz)');


subplot(2,2,3);
hold on

plot(NuIn, NuOut - NuIn, 'r.', X, Y-X, 'r', NuIn, NuIn*0, 'k--');
plot(r, r*0, 'bo', 'MarkerFaceColor', 'w');

% set(gca, 'XLim', [NuIn(1) NuIn(end)], 'YLim', [-4 2]);
set(gca, 'XLim', [NuIn(1) NuIn(end)], 'YLim', [-2 1]);
xlabel('\nu_{in} (Hz)');
ylabel('\nu_{out} - \nu_{in} (Hz)');


subplot(2,2,2);

% plot(NuIn, cumsum(NuIn - NuOut) * diff(NuIn(1:2)), 'r.', ...
%    X, Yint, 'r');
plot(X, Yint - min(Yint), 'r');

% set(gca, 'XLim', [NuIn(1) NuIn(end)], 'YLim', [0 2]);
set(gca, 'XLim', [NuIn(1) NuIn(end)], 'YLim', [0 4]);
xlabel('\nu_{in} (Hz)');
ylabel('\int d\nu_{in} (\nu_{out} - \nu_{in}) (Hz^2)');


subplot(2,2,4);

% plot(NuIn, cumsum(NuIn - NuOut) * diff(NuIn(1:2)), 'r.', ...
%    X, Yint, 'r');
NormFactor = 1/ (sum(exp(-Yint/FINITE_SIZE_NOISE)) * diff(X(1:2)));
plot(X, NormFactor * exp(-Yint/FINITE_SIZE_NOISE), 'r');

% set(gca, 'XLim', [NuIn(1) NuIn(end)], 'YLim', [0 0.21]);
set(gca, 'XLim', [NuIn(1) NuIn(end)], 'YLim', [0 0.1]);
xlabel('\nu_{in} (Hz)');
ylabel('Prob. density func. (a.u.)');

set(gcf, 'PaperUnits', 'inch', 'PaperPosition', [1 2.5 6 6]);
print('-deps2c', 'EffectiveGainFunction.eps');
