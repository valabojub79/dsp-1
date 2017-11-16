%% Unified Particle Swarm Optimisation method [1]
% 
% UPSO - performs a metaheuristic procedure to solve a simple-constrained
%        optimisation problem, such as:
% 
%                min FOBJ(x),    where X = [X_1, X_2, ..., X_Nd]',  
%                 X           
%                 subject to     X_lb <= X <= X_ub, 
%
% XOPT = UPSO(FOBJ,BOUNDARIES) searchs a simple-constrained minimal XOPT of
% the fitness function (or objective function) FOBJ. The feasible region is
% drawn by the simple constraints given by BOUNDARIES. BOUNDARIES is a 
% matrix of Nd-by-2, where Nd is the number of dimensions or variables. 
% Elements of the first column of BOUNDARIES are lower boundaries of each
% design variable, and the elements of the second one are upper boundaries
% of each design variable.
%
% XOPT = UPSO(FOBJ,BOUNDARIES,PARAMETERS) searchs a simple-constrained 
% minimal XOPT of the FOBJ function using a different set of tune parameters 
% of the method. PARAMETERS contains tune parameters for UPSO and it is a 
% struct defined as follows:
% 
%                                      (default values)
%     PARAMETERS = struct(...
%                           'NAG'       ,   25   , ...  % population size
%                           'CP'        ,   2.0  , ...  % particular coef.
%                           'CG'        ,   2.5  , ...  % global coef.
%                           'CHI'       ,   0.7  , ...  % constriction
%                           'U'         ,   0.5  , ...  % unification
%                           'EPS1'      ,   0.5  , ...  % historical eps.
%                           'EPS2'      ,   1e-3 , ...  % population eps.
%                           'MSAT'      ,   10   , ...  % max. saturation
%                           'MITE'      ,   1e12 , ...  % max. iteration
%                           'UNCONST'   ,   false  ...  % unconstrained
%                         );
%
% [XOPT,FOPT,DETAILS] = UPSO(FOBJ,BOUNDARIES,PARAMETERS) performs the above
% process but in this case the final value of fitness function (FOPT) and
% some additional information (DETAILS) are returned. DETAILS is also a 
% struct which contains:
%
%       DETAILS.time    - time required (in seconds) to find XOPT
%              .fevs    - evaluations of fitness function done
%              .steps   - number of steps or iterations performed
%              .outmsg  - flag of convergence (1 is convergence)
%
% --------------------------------------------------------------------------
% Reference:
% [1] PARSOPOULOS, K.E. and VRAHATIS, M.N. 2004. UPSO: A unified particle 
%     swarm optimization scheme. Proceedings of the International Conference 
%     of "Computational Methods in Sciences and Engineering. Zeist, The 
%     Netherlands: VSP International Science Publishers, pp. 868�873.
% --------------------------------------------------------------------------
%
% {Copyright (c) Jorge M. Cruz-Duarte. All rights reserved.}
%
% Departamento de Ingenier�a El�ctrica,
% Divisi�n de Ingenier�as, Campus Irapuato-Salamanca,
% Universidad de Guanajuato, Salamanca, Guanajuato, M�xico.
%
% Grupo de Investigaci�n CEMOS,
% Escuela de Ingenier�as El�ctrica, Electr�nica y de Telecomunicaciones,
% Universidad Industrial de Santander, Bucaramanga, Santander, Colombia.
%
% Modifications:
%                   2015-oct ver :: v0
%                   2014-apr ver :: v1
%                   2016-apr ver :: v2
%
% Contact us: jorge.cruz@ugto.mx

function [Pg,fPg,details] = UPSO(fObj,bnd,parameters)

% Read parameters
if nargin < 3,
    Na      = 40;
    cp      = 2*1.494;
    cg      = 2*1.494;
    chi     = 0.729;
    u       = 0.5;
    
    eps1    = 1e-1;
    eps2    = 1e-3;
    M       = 1e3;
    msat    = 10;
    
    unconst = false;
else
    Na      = parameters.NAG;
    cp      = parameters.CP;
    cg      = parameters.CG;
    chi     = parameters.CHI;
    u       = parameters.U;
    
    eps1    = parameters.EPS1;
    eps2    = parameters.EPS2;
    M       = parameters.MITE;
    msat    = parameters.MSAT;
    
    unconst = parameters.UNCONST;
end

% Read dimensions
Nd          = size(bnd,1);

% Define a quick function to get Pi values
get_Pi      = @(condition,current_position,best_curren_position) ...
    repmat(condition,1,Nd).*current_position + ...
    repmat(~condition,1,Nd).*best_curren_position;

% Define a quick function to obtain f(X)
    function the_results = evaluate_function(the_function,the_positions)
        the_results   = nan(Na,1);
        for s = 1 : Na,
            the_results(s) = the_function(the_positions(s,:));
        end
    end

% Define the ring topology
SSw         = diag(ones(Na,1)) + diag(ones(Na - 1,1),1) + ...
    diag(ones(Na - 1,1),-1) + diag(ones(1),Na - 1) + ...
    diag(ones(1),(1 - Na));

% Define the boundaries for each dimension
bnd         = [min(bnd,[],2) max(bnd,[],2)];
bnd_1       = repmat(bnd(:,1)',Na,1);
bnd_2       = repmat(bnd(:,2)',Na,1);

% Calculate the initial positions
X           = bnd_1 + rand(Na,Nd).*(bnd_2 - bnd_1);
U           = zeros(Na,Nd);
Pi          = X;

% Calculate the initial fitness values
fX          = evaluate_function(fObj,X);
fPi         = fX;

% Found the neighbourhood best position (initial step)
[~,gi]      = min(repmat(fPi,1,Na).*(SSw == 1) + ~(SSw == 1)*100);
Pgi         = Pi(gi,:);

% Found the swarm best position (initial step)
[fPg,g]     = min(fPi); Pg = X(g,:);

% Set auxiliar variables
steps       = 0;
%msatc       = 0;
%sumAVG      = 0;
%sumSD       = 0;

% Plot objective function value evolution                       [toDelete]
% topl        = []; topl1       = []; figure('Name','Enjambre','Color','White','MenuBar','None','Units','Normalized','Position',[0.025 0.025 0.95 0.925]);

%% Main process
tic,
while steps <= M %&& msatc < msat,
    % Update step
    steps    = steps + 1;
    
    % Update the velocity for each particle
    U      = u*(chi*U + cp*rand(Na,Nd).*(Pi - X) + cg*rand(Na,Nd).*...
        (repmat(Pg,Na,1) - X)) + (1 - u)*(chi*U + cp*rand(Na,Nd).*(Pi - X)+ ...
        cg*rand(Na,Nd).*(Pgi - X));
    
    % Update the position for each particle
    X       = X + U;
    
    % Check if the particle is inside the search space
    if unconst == false,
        check = X < bnd_1; X = ~check.*X + check.*bnd_1;
        check = X > bnd_2; X = ~check.*X + check.*bnd_2;
    end
    
    % Evaluate objective function in the new position for each particle
    fX      = evaluate_function(fObj,X);
    
    % Found the best position for each particle
    Pi      = get_Pi(fX < fPi,X,Pi); fPi     = min(fX,fPi);
    
    % Found the neighbourhood's best position
    [~,gi]  = min(repmat(fPi,1,Na).*(SSw == 1) + ~(SSw == 1)*100);
    Pgi     = Pi(gi,:);
    
    fPg_ = fPg;
    % Found the swarm best position
    [fPg,g] = min(fPi); %Pg = X(g,:),
    if fPg < fPg_, Pg = X(g,:);
    else fPg = fPg_; end
    
    % Statistics related to the historical evolution of the fitness
    %sumAVG      = sumAVG + fPg;
    %sumSD       = sumSD + fPg^2;
    %currAVG     = sumAVG/steps;
    %currSD      = eps1*sqrt(sumSD/steps - currAVG^2);
    
    % Statistics related to dispersion of the particles
    %radii       = sqrt(sum((repmat(Pg,Na,1) - X).^2,2));
    %radious     = max(radii);
    
    % Stop criteria
    %if fPg < currAVG + currSD && radious < eps2,
    %    msatc = msatc + 1;
    %else
    %    msatc = 0;
    %end
    
    %fprintf('iteration = %5d, nu = %5.5g, sigma = %5.5g, lse = %5.5g\n',...
    %    steps,Pg(1),Pg(2),fPg);
    
    % Plot objective function value evolution                    [toDelete]
  %   topl        = [topl; [steps,fPg,currAVG,currSD]];
 %    topl1       = [topl1; [steps,fPg,radious,radious]];
    
    % Plot objective function value evolution
%    plot_things(topl,topl1,X,Pg,[1 steps+1;1 steps+1;bnd(1,:)],[nan nan;0 eps2*1.5;bnd(2,:)],eps2,fPg,steps);
    
    
end
t       = toc;

% Plot objective function value evolution                       [toDelete]
% plot_things(topl,topl1,X,Pg,[1 steps+1;1 steps+1;bnd(1,:)],[nan nan;0 eps2*1.5;bnd(2,:)],eps2,fBest,steps);
   
% Save details
details = struct('time',t,'fevs',(Na + 1)*(steps-1),'steps',steps-1);%,...
%     'historical',f_hist);
end
