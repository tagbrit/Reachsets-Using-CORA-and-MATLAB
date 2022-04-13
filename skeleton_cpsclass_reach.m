function skeleton_cpsclass_reach(sysnb, timeStep)
% sysnb selects which system to run: 1 or 2
% if no argument is given, runs system 1 by default
if nargin ==0
    sysnb = 1;
end
%------------- System 1 --------------
if (find(sysnb==1))
    % dimension of the system
    dim=2;
    % bound on input norm
    mu = 0.05;
    
    options.tStart=0; %start time
    %%%%%%%%%%%%%%%%%%
    options.tFinal=2; %final time
    %%%%%%%%%%%%%%%%%%
    options.x0=[1; 0]; %initial state for simulation
    options.R0=zonotope([options.x0,[0.1 0.1; 0 0]]); %initial state for reachability analysis
    %%%%%%%%%%%%%%%%%%
    options.timeStep=0.02; %time step size for reachable set computation
    %%%%%%%%%%%%%%%%%%
    options.taylorTerms=4; %number of taylor terms for reachable sets
    %%%%%%%%%%%%%%%%%%
    options.zonotopeOrder=10; %zonotope order
    %%%%%%%%%%%%%%%%%%
    options.originContained=0;
    options.reductionTechnique='girard';
    options.linAlg = 1;
    
    % center of input set
    options.uTrans=[0; 0]; 
    % input set is a zonotope with above center, and generator matrix
    % mu*identity matrix
    options.U=zonotope([zeros(dim,1),mu*eye(dim)]); %input for reachability analysis
    
    %specify continuous dynamics
    A=[-1 -4; 
        4 -1];
    B=1;
    twoDimSys=linearSys('twoDimSys',A,B); %initialize system
    
    
    %compute reachable set using zonotopes
    Rcont = reach(twoDimSys, options);
    
    % safety check
    % obstacle set is a zonotope with given center and generators
    obstacle_center = [0;0];
    obstacle_generators = 0.1*[1 0.5; 0.5 1];
    % The intersection method, and, requires two zonotope bundles
    Obs  = zonotopeBundle([obstacle_center, obstacle_generators ]);
    % This variable is needed for plotting later.
    zObs = zonotope([obstacle_center, obstacle_generators ]);
    % Get the interesection and convert back to a zonotope object for
    % determining if it's empty or not. 
    Zint = zonotope(and(Obs, Rcont));
    if Zint.isempty
        disp('Safe!')
    else
        disp('Unsafe!')
    end
        
    plotreachsets(Rcont, 'sys1');
    hold on
    plot(zObs);
end
%------------- System 2 --------------
if (find(sysnb==2))
    dim=5;
    mu = 0.01;
    
    options.tStart=0; %start time
    options.tFinal=1; %final time
    % The initial set is a zonotope centered on [1;1;1;1;1], and with
    % generators 0.1*identity matrix
    options.x0=ones(dim,1) ; %initial state for simulation
    options.R0=zonotope([options.x0,0.1*eye(dim)]); %initial state for reachability analysis
    %%%%%%%%%%%%%%%%%%%%
    if nargin < 2
        timeStep = 0.005;
        options.timeStep=timeStep; %time step size for reachable set computation
    else
        options.timeStep = timeStep;
    end
    %%%%%%%%%%%%%%%%%%%%
    options.taylorTerms=4; %number of taylor terms for reachable sets
    %%%%%%%%%%%%%%%%%%%%
    options.zonotopeOrder=40; %zonotope order
    %%%%%%%%%%%%%%%%%%%%
    options.originContained=0;
    options.reductionTechnique='girard';
    options.linAlg = 1;
    
    
    options.uTrans= zeros(dim,1); %input for simulation
    options.U=zonotope([zeros(dim,1),mu*eye(dim)]); %input for reachability analysis
    
    %specify continuous dynamics
    A=[-1 -4 0 0 0;
        4 -1 0 0 0;
        0 0 -3 1 0;
        0 0 -1 -3 0;
        0 0 0 0 -2];
    B=1;
    fiveDimSys=linearSys('fiveDimSys',A,B); %initialize system
    
    
    %compute reachable set using zonotopes
    Rcont = reach(fiveDimSys, options);
    

    plotreachsets(Rcont, 'sys2');
    title(['sys2 step size ', num2str(timeStep)])

end
end
function plotreachsets(reachset, figtitle)
% Plot a collection of zonotopes
figure;
Z=reachset{1};
center = Z.center;
plot(reachset{1});
plot(center(1), center(2),'ro');
hold on
for ii=2:length(reachset)
    plot(reachset{ii});
end
title([figtitle]);
end