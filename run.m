% ORANGE... Alpha1 = 0.97
TauARP = 1/50; % Absolute refractory period.
Crec = 7.15; % average total synaptic coupling: MuRec = Crec*NuIn.
Mu0 = - Crec * 30;  % rheobase current
Sigma0 = 83.0;
Nu0 = 30;

sampling = 50;

load('/home/eloise/Stages/ISS/code/moritzAugustin/adex_comparison/phi0.mat')

mu_sampled = mu_tot(1:sampling:end);
sigma_sampled = sigma_tot(1:sampling:end);
TF = zeros(length(mu_sampled),length(sigma_sampled));


for i=1:length(mu_tot)
    for j =1:length(sigma_tot)
        if  mod((i-1)*length(sigma_tot) + j,sampling) == 0 
            phi0 = eval(['phi' num2str( (i-1)*length(sigma_tot) + j ) ]);
            der = diff(phi0);
            phi_ms = -0.5*sigma_tot(j)^2*der(end);
            i,j
            TF(i/sampling,j/sampling) = phi_ms;
            figure(1)
            plot(phi0)
            hold on
                
        end
    end
end


%searchAndPlotGainFunction(TauARP,Crec,Sigma0,Nu0)