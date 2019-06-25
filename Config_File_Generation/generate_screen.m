function generate_screen()
% Basic parameters
parameters.R=2.5;
parameters.rho=0.004;
parameters.gradient_min=0.4;
parameters.gradient_max=0.6;
parameters.nReceptor=10000;
parameters.Ngprotein=2500;
parameters.Dg=0.002;
parameters.Dr=0;
parameters.eps=0.1;
parameters.tf=600;
parameters.dt=1e-4;
parameters.sampling_rate=100000;
parameters.kinact=0.2734; % From calibration to the N=10000 dataset
parameters.densityGradientMode='align';
parameters.density_gradient_min=0.25;
parameters.density_gradient_max=0.75;
parameters.activity_gradient_min=0.4;
parameters.activity_gradient_max=0.6;

% Screen through variable receptor density
activity_mode={'ratiometric','nonratiometric'};
densityGradient_modes = {'uniform','align','oppose','orthogonal','diagonal'};

parameters.Dr = 0;% 0.0005;

n_realizations = 50;

% Bash file setup is needed to launch all jobs simultaneously.
% As below, set to run from the UNC Longleaf compute cluster.
% Smoldyn-2.56 was compiled and installed to run from 
% /nas/longleaf/home/mikepab/smoldyn-2.56/cmake/smoldyn

fid=fopen('run_cfg.sh','w');
fprintf(fid,'#!/bin/bash\n\n');
for j=1:numel(activity_mode)
    for k=1:n_realizations
        for dg=1:numel(densityGradient_modes)
            parameters.densityGradientMode=densityGradient_modes{dg};
            parameters.modelType=activity_mode{j};
            parameters.random_seed=k;

            filebase=sprintf('20190328_R10000_%s_seed%i_denGrad%s',activity_mode{j},k,parameters.densityGradientMode);
            parameters.output_file=filebase;
            create_cfg([filebase '.cfg'],parameters);
            fprintf(fid,'sbatch -p general -N 1 -J gpcr -t 80:00:00 --wrap="/nas/longleaf/home/mikepab/smoldyn-2.56/cmake/smoldyn %s"\n',[filebase '.cfg']);
        end
    end
end
fclose(fid);
end