function returnData=ldacpso(fitfuncHandle,nDim,varargin)
% PSO minimizer template code for LDAC
% S=PSO(Fhandle,N)
% Runs PSO on the fitness function with handle Fhandle. If Fname is the name
% of the function, Fhandle = @(x) <Fname>(x, FP). The optional
% parameters to supply to the fitness function should be given in the structure FP.
% N is the dimensionality of the  fitness function. If N is empty, then N=2. 
% The output is returned in the structure S. The field of S are:
%  bestLocation : Best location found (in standardized coordinates)
%  bestSNR: Best fitness value found
%  totalFuncEvals: Total number of fitness function evaluations
%
%S=PTAPSO(Fhandle,N,P)
%overrides the default PSO parameters with those provided in structure P.
%Set any of the fields of P to [] in order to invoke the corresponding
%default value. The fields of P are as follows.
%     'popSize': Number of PSO particles
%     'maxSteps': Number of iterations for termination
%     'c1','c2': acceleration constant
%     'maxInitialVelocity': maximum value for each initial velocity
%                           component
%     'maxVelocity': maximum value for each velocity component for all
%                     subsequent iterations
%     'startInertia': Starting value of inertia weight
%     'endInertia': End value of inertia at iteration# <maxSteps>
%     'boundaryCond': Set to '' for invisible walls boundary conditions.
%                     Any other value is passed onto the fitness function
%                     to handle. 
%     'fminMaxFunEvals': maximum number of function evaluations for
%                        fminsearch around new gbest.
%Setting P to [] will invoke the default values for all pso parameters.
%NOTE: The optional P argument is normally to be used only for testing and
%debugging (e.g., reduce the number of particles and/or iterations for a
%quick run). Stable values should be hardcoded in the function and used as
%the default.
%
%Example:
%nDim = 5;
%fitFuncParams = struct('rmin',-5*ones(1,nDim),...
%                           'rmax',5*ones(1,nDim));
%fitFuncHandle = @(x) ldacpsotestfunc(x,fitFuncParams);
%psoOut = ldacpso(fitFuncHandle,5);

% Authors
% Soumya D. Mohanty, Dec 2017
% Light weight version of a PSO code used in PTA data analysis.

popsize=40;
maxSteps= 2000;
c1=2;
c2=2;
max_initial_velocity = 0.5;
max_velocity = 0.2;
dcLaw_a = 0.9;
dcLaw_b = 0.4;
dcLaw_c = maxSteps;
dcLaw_d = 0.2;
bndryCond = '';%invisible walls boundary condition
nmOptions = optimset('fminsearch');
nmOptions.MaxFunEvals=200;
nmOptions.Display = 'off';
%Override defaults if needed
nreqArgs = 2;
if nargin-nreqArgs
    for largs = 1:(nargin-nreqArgs)
        if isempty(varargin{largs})
        else
            switch largs
                case 1
                    psoParams = varargin{largs};
                    %pso parameters
                    psoParfieldNames = fieldnames(psoParams);
                    for lpfields = 1:length(psoParfieldNames)
                        fieldVal = getfield(psoParams,psoParfieldNames{lpfields});
                        if isempty(fieldVal)
                            %do nothing and retain default value
                        else
                            switch psoParfieldNames{lpfields}
                                case 'popSize'
                                    popsize = fieldVal;
                                case 'maxSteps'
                                    maxSteps = fieldVal;
                                case  'c1'
                                    c1 = fieldVal;
                                case  'c2'
                                    c2 = fieldVal;
                                case 'maxInitialVelocity'
                                    max_initial_velocity = fieldVal;
                                case 'maxVelocity'
                                    max_velocity = fieldVal;
                                case 'startInertia'
                                    dcLaw_a = fieldVal;
                                case 'endInertia'
                                    dcLaw_b = fieldVal;
                                case 'boundaryCond'
                                    bndryCond = fieldVal;
                                case 'fminMaxFunEvals'
                                    nmOptions.MaxFunEvals = fieldVal;
                            end
                        end
                    end
            end
        end
    end
end

%Information about each particle stored as a row of a matrix ('pop').
%Specify which column stores what information.
%(The fitness function for matched filtering is SNR, hence the use of 'snr'
%below.)
partCoordCols = 1:nDim; %Particle location
partVelCols = (nDim+1):2*nDim; %Particle velocity
partPbestCols = (2*nDim+1):3*nDim; %Particle pbest
partSnrPbestCols = 3*nDim+1; %Fitness value at pbest
partSnrCurrCols = partSnrPbestCols+1; %Fitness value at current iteration
partSnrLbestCols = partSnrCurrCols+1; %Fitness value at local best location
partInertiaCols = partSnrLbestCols+1; %Inertia weight
partLocalBestCols = (partInertiaCols+1):(partInertiaCols+nDim); %Particles local best location
partFlagFitEvalCols = partLocalBestCols(end)+1; %Flag whether fitness should be computed or not
partFitEvalsCols = partFlagFitEvalCols+1; %Number of fitness evaluations

nColsPop = length([partCoordCols,partVelCols,partPbestCols,partSnrPbestCols,...
                   partSnrCurrCols,partSnrLbestCols,partInertiaCols,partLocalBestCols,...
                   partFlagFitEvalCols,partFitEvalsCols]);
pop=zeros(popsize,nColsPop);

% Best value found by the swarm over its history
best_in_history = inf;  
% Location of the best value found by the swarm over its history
gbest_loc_in_history = 2*ones(1,nDim);
% Best value found by the swarm at the current iteration
best_fitness = inf;
pop(:,partCoordCols)=rand(popsize,nDim);
pop(:,partVelCols)= -pop(:,partCoordCols)+rand(popsize,nDim);
pop(:,partPbestCols)=pop(:,partCoordCols);
pop(:,partSnrPbestCols)= inf;
pop(:,partSnrCurrCols)=0;
pop(:,partSnrLbestCols)= inf;
pop(:,partLocalBestCols) = 0;
pop(:,partFlagFitEvalCols)=1;
pop(:,partInertiaCols)=0;
pop(:,partFitEvalsCols)=0;
bestLocHistPart = 1;

fitnessValues = zeros(popsize,1);
for lpc_steps=2:maxSteps  
    if isempty(bndryCond)
        %Invisible wall boundary condition
        fitnessValues = fitfuncHandle(pop(:,partCoordCols));
    else
        %Special boundary condition handled by fitness function
        [fitnessValues,~,pop(:,partCoordCols)] = fitfuncHandle(pop(:,partCoordCols));
    end
    for k=1:popsize
        startCoord = pop(k,partCoordCols);
        endFitnessVal = pop(k,partSnrCurrCols);
        computeOK = pop(k,partFlagFitEvalCols);
        funcCount = 0;
        if computeOK
            endFitnessVal = fitnessValues(k);
            endCoord = startCoord;
            funcCount = 1;
        else
            endCoord = startCoord;
            funcCount = 0;
        end
        pop(k,partCoordCols)=endCoord;
        pop(k,partSnrCurrCols)=endFitnessVal;
        pop(k,partFitEvalsCols)=pop(k,partFitEvalsCols)+funcCount;
        if pop(k,partSnrPbestCols) > pop(k,partSnrCurrCols)
            pop(k,partSnrPbestCols) = pop(k,partSnrCurrCols);
            pop(k,partPbestCols) = pop(k,partCoordCols);
        end
    end
    [best_fitness,bestfitParticle]=min(pop(:,partSnrCurrCols)); 
    if best_in_history > best_fitness
        best_in_history     =   best_fitness;
        gbest_loc_in_history  =  pop(bestfitParticle,partCoordCols);
        startCoord = gbest_loc_in_history;
        [endCoord,endFitnessVal,~,nmInfo] = fminsearch(fitfuncHandle,...
            startCoord,...
            nmOptions);
        funcCount = nmInfo.funcCount;
        pop(bestfitParticle,partCoordCols)=endCoord;
        pop(bestfitParticle,partSnrCurrCols)=endFitnessVal;
        pop(bestfitParticle,partFitEvalsCols)=pop(bestfitParticle,partFitEvalsCols)+funcCount;
        gbest_loc_in_history = endCoord;
        best_in_history = endFitnessVal;
    else
    end
    for k = 1:popsize
           switch k
               case 1
                   ringNbrs = [popsize,k,k+1];
               case popsize
                   ringNbrs = [k-1,k,1];
               otherwise
                   ringNbrs = [k-1,k,k+1];
           end
           [~,lbestPart] = min(pop(ringNbrs,partSnrCurrCols));
           lbestTruIndx = ringNbrs(lbestPart);
           lbestCurrSnr = pop(lbestTruIndx,partSnrCurrCols);
           if lbestCurrSnr < pop(k,partSnrLbestCols)
               pop(k,partSnrLbestCols) = lbestCurrSnr;
               pop(k,partLocalBestCols) = pop(lbestTruIndx,partCoordCols);
           end
    end
     inertiaWt = max(dcLaw_a-(dcLaw_b/dcLaw_c)*(lpc_steps-1),dcLaw_d);
    for k=1:popsize
        pop(k,partInertiaCols)=inertiaWt;
        partInertia = pop(k,partInertiaCols);
        chi1 = diag(rand(1,nDim));
        chi2 = diag(rand(1,nDim));
        pop(k,partVelCols)=partInertia*pop(k,partVelCols)+...
            c1*(pop(k,partPbestCols)-pop(k,partCoordCols))*chi1+...
            c2*(pop(k,partLocalBestCols)-pop(k,partCoordCols))*chi2;
        maxvBustCompPos = find(pop(k,partVelCols) > max_velocity);
        maxvBustCompNeg = find(pop(k,partVelCols) < -max_velocity);
        if ~isempty(maxvBustCompPos)
            pop(k,partVelCols(maxvBustCompPos))= max_velocity;
        end
        if ~isempty(maxvBustCompNeg)
            pop(k,partVelCols(maxvBustCompNeg))= -max_velocity(1);
        end
        pop(k,partCoordCols)=pop(k,partCoordCols)+pop(k,partVelCols);        
        if any(pop(k,partCoordCols)> 1 | ...
                pop(k,partCoordCols)< 0)
                    pop(k,partSnrCurrCols)= inf;
                    pop(k,partFlagFitEvalCols)= 0;
        else
            pop(k,partFlagFitEvalCols)=1;
        end
    end
end

stepsToEnd=lpc_steps-1;
actualEvaluations = sum(pop(:,partFitEvalsCols));
returnData = struct('totalSteps',stepsToEnd,...
    'totalFuncEvals',actualEvaluations,...
    'bestLocation',gbest_loc_in_history,...
    'bestSNR',best_in_history);





