function create_cfg(filename,parameters)
% This function writes a Smoldyn configuration file to filename, using the information
% contained within the parameters variable.
%
% INPUTS: filename: filename for the Smoldyn configuration file. I typically named these .cfg files.
%         parameters: struct containing the following fields...
%
%         R, radius of simulated spherical cell (um)
%         rho, interaction radius (um)
%         nReceptor, number receptor molecules
%         Ngprotein, number gprotein molecules
%         densityGradientMode = {'uniform',... % Uniform receptor density distribution
%                                'align',... % Increase with activity gradient (0 degree)
%                                'oppose',... % Decrease with activity gradient (180 degree)
%                                'orthogonal',... % Decrease along 90 degree angle to activity gradient
%                                'diagonal'} % Increase along 145 degree angle to activity gradient
%         density_gradient_min, gradient min for receptor density
%         density_gradient_max, gradient max for receptor density
%         activity_gradient_min, gradient min for activity density
%         activity_gradient_max, gradient max for activity density
%         modelType : {'ratiometric','nonratiometric'}
%         --> if modelType is 'nonratiometric', also need field kinact
%         Dg, diffusion coefficient of G protein (um2/s)
%         Dr, diffusion coefficient of receptor (um2/s)
%         eps, metaparameter for adjusting smoldyn simulation boundaries. Typically used 0.1 um.
%         tf, final timepoint for simulation (s).
%         dt, simulation timestep (s)
%         sampling_rate, # of timesteps between writing simulation data to file.
%         output_file, name of the smoldyn output file
%         random_seed, random seed to control random number generation
%
% ====================================================================================================

% Some basic parameter checking; additional errors are thrown during later checks.
assert(parameters.gradient_min>0 && parameters.gradient_min<parameters.gradient_max,'Gradient error');
assert(parameters.gradient_max>0 && parameters.gradient_max<1,'Gradient error')
assert(sqrt((parameters.Dr+parameters.Dg)*parameters.dt)/parameters.rho<0.15,'RMS:rho too large')
assert(parameters.dt*parameters.sampling_rate<parameters.tf,'check timestep or sampling rate');


fid=fopen(filename,'w');
[XYZa,XYZi]=generate_receptor_distribution(parameters);
fprintf(fid,'random_seed %i\n',parameters.random_seed);
fprintf(fid,'dim 3\n');
simbounds = [0-parameters.eps,parameters.R*2+parameters.eps];
fprintf(fid,'boundaries x %g %g\n',simbounds(1),simbounds(2));
fprintf(fid,'boundaries y %g %g\n',simbounds(1),simbounds(2));
fprintf(fid,'boundaries z %g %g\n',simbounds(1),simbounds(2));
fprintf(fid,'species Ri Ra Gi Ga RaGi_transient RiGa_transient\n');
fprintf(fid,'difc Ri(all) %g\n',parameters.Dr);
fprintf(fid,'difc Ra(all) %g\n',parameters.Dr);
fprintf(fid,'difc Gi(all) %g\n',parameters.Dg);
fprintf(fid,'difc Ga(all) %g\n',parameters.Dg);

fprintf(fid,'start_surface membrane\n');
fprintf(fid,'polygon both edge\n');
fprintf(fid,'panel sphere %g %g %g %g 50 50 panel_memb\n',parameters.R,...
                                                          parameters.R,...
                                                          parameters.R,...
                                                         -parameters.R);
fprintf(fid,'end_surface\n');

simbounds_surface = [0-parameters.eps,parameters.R*2+2*parameters.eps];
fprintf(fid,'start_surface sim_boundary\n');
fprintf(fid,'action both all reflect\n');
fprintf(fid,'polygon both none\n');
fprintf(fid,'panel rect +x %g %g %g %g %g\n',-simbounds_surface(1),...
                                             -simbounds_surface(1),...
                                             -simbounds_surface(1),...
                                              simbounds_surface(2),...
                                              simbounds_surface(2));
fprintf(fid,'panel rect -x %g %g %g %g %g\n', simbounds_surface(2),...
                                             -simbounds_surface(1),...
                                             -simbounds_surface(1),...
                                              simbounds_surface(2),...
                                              simbounds_surface(2));

fprintf(fid,'panel rect +y %g %g %g %g %g\n',-simbounds_surface(1),...
                                             -simbounds_surface(1),...
                                             -simbounds_surface(1),...
                                              simbounds_surface(2),...
                                              simbounds_surface(2));
fprintf(fid,'panel rect -y %g %g %g %g %g\n',-simbounds_surface(1),...
                                             -simbounds_surface(2),...
                                             -simbounds_surface(1),...
                                              simbounds_surface(2),...
                                              simbounds_surface(2));
                                          
 
fprintf(fid,'panel rect +z %g %g %g %g %g\n',-simbounds_surface(1),...
                                             -simbounds_surface(1),...
                                             -simbounds_surface(1),...
                                              simbounds_surface(2),...
                                              simbounds_surface(2));
fprintf(fid,'panel rect -z %g %g %g %g %g\n',-simbounds_surface(1),...
                                             -simbounds_surface(1),...
                                             -simbounds_surface(2),...
                                              simbounds_surface(2),...
                                              simbounds_surface(2));     
fprintf(fid,'end_surface\n');
fprintf(fid,'reaction RaGi_enc Ra(down) + Gi(down) -> RaGi_transient(down)\n');
fprintf(fid,'reaction_probability RaGi_enc 1\n');
fprintf(fid,'binding_radius RaGi_enc %g\n',parameters.rho);

fprintf(fid,'reaction RaGi_cat RaGi_transient(down) -> Ra(down) + Ga(down)\n');
fprintf(fid,'reaction_probability RaGi_cat 1\n');

if strcmpi(parameters.modelType,'ratiometric')
    fprintf(fid,'reaction RiGa_enc Ri(down) + Ga(down) -> RiGa_transient(down)\n');
    fprintf(fid,'reaction_probability RiGa_enc 1\n');
    fprintf(fid,'binding_radius RiGa_enc %g\n',parameters.rho);
    fprintf(fid,'reaction RiGa_cat RiGa_transient(down) -> Ri(down) + Gi(down)\n');
    fprintf(fid,'reaction_probability RiGa_cat 1\n');
    
elseif strcmpi(parameters.modelType,'nonratiometric')
    fprintf(fid,'reaction Ga2Gi Ga(down) -> Gi(down) %g\n',parameters.kinact);
else
    error('invalid modelType specified')
end

fprintf(fid,'start_compartment cell_domain\n');
fprintf(fid,'surface membrane\n');
fprintf(fid,'point %g %g %g\n',parameters.R,parameters.R,parameters.R);
fprintf(fid,'end_compartment\n');

fprintf(fid,'time_start 0\n');
fprintf(fid,'time_stop %g\n',parameters.tf);
fprintf(fid,'time_step %g\n',parameters.dt);

fprintf(fid,'output_files %s.xyz\n',parameters.output_file);
fprintf(fid,'cmd B molpos Ga(all) %s.xyz\n',parameters.output_file);
fprintf(fid,'cmd B molpos Gi(all) %s.xyz\n',parameters.output_file);
fprintf(fid,'cmd B molpos Ra(all) %s.xyz\n',parameters.output_file);
fprintf(fid,'cmd B molpos Ri(all) %s.xyz\n',parameters.output_file);

fprintf(fid,'cmd N %i molpos Ga(all) %s.xyz\n',parameters.sampling_rate,...
                                               parameters.output_file);
fprintf(fid,'cmd N %i molpos Gi(all) %s.xyz\n',parameters.sampling_rate,...
                                               parameters.output_file);
fprintf(fid,'cmd N %i molpos Ra(all) %s.xyz\n',parameters.sampling_rate,...
                                               parameters.output_file);
fprintf(fid,'cmd N %i molpos Ri(all) %s.xyz\n',parameters.sampling_rate,...
                                               parameters.output_file);
                                           
fprintf(fid,'surface_mol %i Ga(down) membrane sphere panel_memb\n',parameters.Ngprotein);

for i=1:size(XYZa,1)
   fprintf(fid,'surface_mol 1 Ra(down) membrane sphere panel_memb %g %g %g\n',XYZa(i,1),XYZa(i,2),XYZa(i,3));
end
for i=1:size(XYZi,1)
   fprintf(fid,'surface_mol 1 Ri(down) membrane sphere panel_memb %g %g %g\n',XYZi(i,1),XYZi(i,2),XYZi(i,3));
end
fclose(fid);
end

function [XYZa,XYZi] = generate_receptor_distribution(parameters)
% This function creates two lists of (x,y,z) coordinate matrices for
% active receptor and inactive receptor, according to the parameters.
% Smoldyn 2.56 does not have built-in functions set up such gradients.
%
% Required inputs:
% parameters.nReceptor
% parameters.densityGradientMode
% parameters.density_gradient_min
% parameters.density_gradient_max
% parameters.activity_gradient_min
% parameters.activity_gradient_max
% parameters.R
% ================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate coordinates; do rejection sampling if density gradient needed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
assert(strcmpi(parameters.densityGradientMode,'uniform') || ...
       strcmpi(parameters.densityGradientMode,'align') || ...
       strcmpi(parameters.densityGradientMode,'oppose') || ...
       strcmpi(parameters.densityGradientMode,'orthogonal') || ...
       strcmpi(parameters.densityGradientMode,'diagonal'),...
       'Invalid densityGradientMode specified');
       
N=parameters.nReceptor;
TH=2*pi*rand(1,N);
PH=asin(-1+2*rand(1,N));
R=parameters.R;
[X,Y,Z]=sph2cart(TH,PH,R);
X=X+R;
Y=Y+R;
Z=Z+R;

if strcmpi(parameters.densityGradientMode,'uniform')
    % No need to do rejection sampling; continue to activity gradient setup
    XYZ = [X',Y',Z'];
else    
    % Set up probability gradient for accepting a coordinate along x
    if strcmpi(parameters.densityGradientMode,'align') || strcmpi(parameters.densityGradientMode,'oppose')
        Pr_setDensity = (parameters.density_gradient_max-parameters.density_gradient_min) / ...
                        (2*R)*X +parameters.density_gradient_min;

        if strcmpi(parameters.densityGradientMode,'align')
            % Accept some proportion based on probability sampling
            toAccept = rand(size(X))<Pr_setDensity;
        elseif strcmpi(parameters.densityGradientMode,'oppose')
            toAccept = ~(rand(size(X))<Pr_setDensity);
        end

        XYZ=[X(toAccept)',Y(toAccept)',Z(toAccept)'];

        % While total accepted number < target number
        remaining_N = N - size(XYZ,1); 
        while remaining_N > 0
            % Re-draw more random numbers
            TH=2*pi*rand(1,remaining_N);
            PH=asin(-1+2*rand(1,remaining_N));
            [X,Y,Z]=sph2cart(TH,PH,R);
            X=X+R; Y=Y+R; Z=Z+R;

            Pr_setDensity = (parameters.density_gradient_max-parameters.density_gradient_min) ...
                            / (2*R)*X +parameters.density_gradient_min;

            if strcmpi(parameters.densityGradientMode,'align')
                % Accept some proportion based on probability sampling
                toAccept =   rand(size(X))<Pr_setDensity;
            elseif strcmpi(parameters.densityGradientMode,'oppose')
                toAccept = ~(rand(size(X))<Pr_setDensity);
            end

            newXYZ=[X(toAccept)',Y(toAccept)',Z(toAccept)'];

            % Re-accept more random numbers, up to remaining_N
            if size(newXYZ,1) > remaining_N
                newXYZ = newXYZ(1:remaining_N,:);
            end

            XYZ = [XYZ; newXYZ];
            remaining_N = N - size(XYZ,1); 
        end
    elseif strcmpi(parameters.densityGradientMode,'orthogonal')% Set up probability gradient along y [orthogonal]
        Pr_setDensity = (parameters.density_gradient_max-parameters.density_gradient_min) / ...
                        (2*R)*Y + parameters.density_gradient_min;
        % Accept some proportion based on probability sampling
        toAccept = rand(size(X))<Pr_setDensity;
        XYZ=[X(toAccept)',Y(toAccept)',Z(toAccept)'];
        % While total accepted number < target number
        remaining_N = N - size(XYZ,1); 
        while remaining_N > 0
            % Re-draw more random numbers
            TH=2*pi*rand(1,remaining_N);
            PH=asin(-1+2*rand(1,remaining_N));
            [X,Y,Z]=sph2cart(TH,PH,R);
            X=X+R; Y=Y+R; Z=Z+R;

            Pr_setDensity = (parameters.density_gradient_max-parameters.density_gradient_min) ...
                            / (2*R)*Y +parameters.density_gradient_min;

            % Accept some proportion based on probability sampling
            toAccept = rand(size(X))<Pr_setDensity;

            newXYZ=[X(toAccept)',Y(toAccept)',Z(toAccept)'];

            % Re-accept more random numbers, up to remaining_N
            if size(newXYZ,1) > remaining_N
                newXYZ = newXYZ(1:remaining_N,:);
            end

            XYZ = [XYZ; newXYZ];
            remaining_N = N - size(XYZ,1); 
        end
    elseif strcmpi(parameters.densityGradientMode,'diagonal') % set up probability gradient along (-1,1,0);
        diagonal_projection = zeros(1,numel(X));
        for i=1:numel(X)
            diagonal_projection(i) = dot([X(i)-R,Y(i)-R,Z(i)-R],[-1/sqrt(2),1/sqrt(2),0]/norm([-1/sqrt(2),1/sqrt(2),0]));
        end
        diagonal_projection = (diagonal_projection+R)/(2*R);
        Pr_setDensity = (parameters.density_gradient_max-parameters.density_gradient_min) ...
                *diagonal_projection + parameters.density_gradient_min;
        % Accept some proportion based on probability sampling
        toAccept = rand(size(X))<Pr_setDensity;
        XYZ=[X(toAccept)',Y(toAccept)',Z(toAccept)'];
        % While total accepted number < target number
        remaining_N = N - size(XYZ,1); 
        while remaining_N > 0
            % Re-draw more random numbers
            TH=2*pi*rand(1,remaining_N);
            PH=asin(-1+2*rand(1,remaining_N));
            [X,Y,Z]=sph2cart(TH,PH,R);
            X=X+R; Y=Y+R; Z=Z+R;


            diagonal_projection = zeros(1,remaining_N);
            for i=1:remaining_N
                diagonal_projection(i) = dot([X(i)-R,Y(i)-R,Z(i)-R],[-1/sqrt(2),1/sqrt(2),0]/norm([-1/sqrt(2),1/sqrt(2),0]));
            end
            diagonal_projection = (diagonal_projection+R)/(2*R);
            Pr_setDensity = (parameters.density_gradient_max-parameters.density_gradient_min) ...
                    * diagonal_projection + parameters.density_gradient_min;
            
            % Accept some proportion based on probability sampling
            toAccept = rand(size(X))<Pr_setDensity;

            newXYZ=[X(toAccept)',Y(toAccept)',Z(toAccept)'];

            % Re-accept more random numbers, up to remaining_N
            if size(newXYZ,1) > remaining_N
                newXYZ = newXYZ(1:remaining_N,:);
            end

            XYZ = [XYZ; newXYZ];
            remaining_N = N - size(XYZ,1); 
        end
    else
        error('invalid density setting specified');
    end
end % Generated XYZ coordinates


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Apply activity gradient
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pr_setActivity = (parameters.activity_gradient_max-parameters.activity_gradient_min) / ...
     (2*R)*XYZ(:,1) + parameters.activity_gradient_min;
 
isActive=rand(size(XYZ(:,1)))<Pr_setActivity;
XYZa = XYZ(isActive,:);
XYZi = XYZ(~isActive,:);
end