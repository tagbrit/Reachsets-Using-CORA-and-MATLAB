function cpsclass_nnv_skeleton(sysnb)
% cpsclass_nnv()
%   run all systems
% cpsclass_nnv(sysnb)
%   sysnb is a vector containing at most the numbers 1,2,3
supported_sys = [1,2,3];
if nargin ==0
    sysnb = supported_sys;
else
    temp = sort(unique(sysnb));
    for ii=1:length(sysnb)
        assert(~isempty(find(supported_sys==temp(ii),1)), ['testnb must be one of ', num2str(supported_sys)])
    end
end
% If the parallel toolbox on your machine works fine, go for it and set
% this to 4 or 6. Otherwise, keep it at 1.
numcores = 1;
%__________________________________________________________________________________
% This first system is an example which demonstrates the basic use of the
% toolbox. No code modifications are required, you can run it.
if (find(sysnb==1))
    
    %% test 1: FFNNS falsify
    % Weights of NN layer. You can tell it takes in 2 inputs and has two
    % neurons
    W = [1 1; 0 1];
    % Bias of NN layer
    b = [0; 0.5];
    L = LayerS(W, b, 'poslin'); % poslin is ReLU
    Layers = [L];
    
    % Create the forward NN
    F = FFNNS(Layers);
    
    % Input star is simply the box given by the bounds below
    lb = [-1; -1];
    ub = [1; 1];    
    I = Star(lb, ub);
    
    % Reach!
    [R, ~] = F.reach(I, 'exact-star');
    
    %This represents the unsafe set of states to be avoided
    % U = {x | Gx <= g}
    % So this is -x1 <= -1.5, i.e. x1 >= 1.5
    G = [-1 0];
    g = [-1.5];        
    U = HalfSpace(G, g);
    
    % Plot the reachset and unsafe set
    R.plots(R);
    hold on
    U.plot();
    title('Reachset and Unsafe set')
    
    % Intersection with unsafe set?
    %figure
    nintersections = 0;
    for s=1:length(R) % in general, it's a union of stars
        SI = R(s).intersectHalfSpace(U.G, U.g);
        if (isa(SI, 'Star') && ~SI.isEmptySet())
            fprintf('Star %i intersects unsafe set\n', s)
            nintersections = nintersections+1;
        end
    end
    if nintersections > 0
        disp('Unsafe!! AAAAAAA!!!')
    end
    
    % Try to find points in the input set I which witness violation, i.e.
    % whose corresponding NN outputs are in the unsafe set
    n_samples = 1000;
    
    % Method falsify does a naive uniform random sampling. As such there is
    % no guarantee that it will find falsifying inputs even if they exist
    counter_inputs = F.falsify(I, U, n_samples);
    counter_outputs = F.sample(counter_inputs);
    
    figure;
    subplot(1, 2, 1);
    I.plot;
    hold on;
    plot(counter_inputs(1, :), counter_inputs(2, :), 'o');
    title('Input Set and counter inputs');
    subplot(1, 2, 2);
    Star.plots(R);
    hold on;
    plot(counter_outputs(1, :), counter_outputs(2, :), 'o');
    title('Output set and counter outputs');
    
end
%__________________________________________________________________________________

%tests below originally taken from test_FFNNS_verify_DFS.m
if (find(sysnb==2))    
    
    % Layer 1 has the following parameters:
    % x1 = z1 - z2 - 1
    % x2 = 0.5z1 + 2z2 + 0.5
    % x3 = -z1 + z2
    % Activation function is poslin (i.e., ReLU)
    W1 = [1 -1; 0.5 2; -1 1];
    b1 = [-1; 0.5; 0];
    L1 = LayerS(W1, b1, 'poslin');
    
    % Layer 2 has following parameters:
    % y1 = -2x1 + x2 + x3 - 0.5
    % y2 = 0.5x1 + x2 + x3 -0.5
    % Activation function is purelin
    W2 = [-2 1 1; 0.5 1 1];
    b2 = [-0.5; -0.5];
    L2 = LayerS(W2, b2, 'purelin');
    
    F = FFNNS([L1 L2]); % construct Feedforward neural network
    
    % Input set is the box -2 <= z1 <= 2, -1 <= z2 <= 2
    lb = [-2; -1];
    ub = [2; 2];
    I = Star(lb, ub); % construct input set
    
    % Compute the reachset using the exact-star method
    [R, t] = F.reach(I, 'exact-star', numcores, []); % compute the exact reachable set
    
    % Also compute it using the approx-star method
    [Rapx, t] = F.reach(I, 'approx-star', numcores, []);
        
    % plot reachable set
    fig = figure;
    subplot(1, 2, 1);
    I.plot;
    title('Input Set', 'FontSize', 20);
    xlabel('x_1', 'FontSize', 16);
    ylabel('x_2', 'FontSize', 16);
    
    subplot(1, 2, 2)
    Star.plots(Rapx, 'g')
    hold on
    Star.plots(R)
    
    title('Output Set', 'FontSize', 20);
    xlabel('y_1', 'FontSize', 16);
    ylabel('y_2', 'FontSize', 16);
        
    % unsafe region: y[1] >= 5
    G = [-1 0];
    g = [-5];
    U = HalfSpace(G, g);
    subplot(1,2,2);
    hold on
    U.plot();
    
    % Check whether the reachset intersects the unsafe set
    
    nintersections = 0;
    for s=1:length(R) % in general, it's a union of stars
        SI = R(s).intersectHalfSpace(U.G, U.g);
        if (isa(SI, 'Star') && ~SI.isEmptySet())
            fprintf('Star %i intersects unsafe set\n', s)
            nintersections = nintersections+1;
        end
    if nintersections > 0
        disp('Exact Unsafe!! AAAAAAA!!!')
    end    
    end
    
    nintersections = 0;
    for s=1:length(Rapx) % in general, it's a union of stars
        SI = Rapx(s).intersectHalfSpace(U.G, U.g);
        if (isa(SI, 'Star') && ~SI.isEmptySet())
            fprintf('Star %i intersects unsafe set\n', s)
            nintersections = nintersections+1;
        end
    end
    if nintersections > 0
        disp('Approx Unsafe!! AAAAAAA!!!')
    end
end
% %__________________________________________________________________________________

if (find(sysnb==3))
    %% Test 6: ACAS Xu
    load ACASXU_run2a_2_1_batch_2000.mat;
    Layers = [];
    n = length(b);
    % Instantiate the first n-1 layers
    for i=1:n - 1
        bi = cell2mat(b(i));
        Wi = cell2mat(W(i));
        Li = LayerS(Wi, bi, 'poslin');
        Layers = [Layers Li];
    end
    % Instantiate the last layer, which has the purelin activation
    bn = cell2mat(b(n));
    Wn = cell2mat(W(n));
    Ln = LayerS(Wn, bn, 'purelin');
    
    Layers = [Layers Ln];
    F = FFNNS(Layers);
    
    % Input Constraints
    % 1500 <= i1(\rho) <= 1800,
    % -0.06 <= i2 (\theta) <= 0.06,
    % 3.1 <= i3 (\shi) <= 3.14
    % 980 <= i4 (\v_own) <= 1200,
    % 960 <= i5 (\v_in) <= 1000
    
    lb = [1500; -0.06; 3.1; 980; 960];
    ub = [1800; 0.06; 3.14; 1200; 1000];    
    % normalize input
    for i=1:5
        lb(i) = (lb(i) - means_for_scaling(i))/range_for_scaling(i);
        ub(i) = (ub(i) - means_for_scaling(i))/range_for_scaling(i);
    end    
    I = Star(lb, ub);
    
    % Output: [x1 = COC; x2 = Weak Left; x3 = Weak Right; x4 = Strong Left; x5 = Strong Right]
    
    % [R, ~] = F.reach(I, 'exact-star', numCores); % exact reach set using polyhdedron
    % F.print('F_exact_star.info'); % print all information to a file
    
    % Compute the reach-set using approx-star and the specified numcores
    [R, ~] = F.reach(I, 'approx-star', numcores, []);
    F.print('F_approx_star.info'); % print all information to a file
    %R2.plot('r', [1 0 0 0 0 ; 0 1 0 0 0]);
    %title('Projection of reachset on first 2 dimensions')
    
    % unsafe region: COC is the minimal score: x1 <= x2; x1 <= x3; x1 <= x4, x1<= x5    
    G = [1 -1 0 0 0;
        1 0 -1 0 0;
        1 0 0 -1 0;
        1 0 0 0 -1;
        0 0 0 0 0];
    g = [0; 0; 0; 0; 0];
    
    
    t = tic;
    fprintf('\nVerifying reach set...\n');
    unsafe = 0;
    n = length(R);
    for i=1:n
        S = R(i).intersectHalfSpace(G, g);
        if ~isempty(S)            
            fprintf('\nThe %d^th star output set reaches the unsafe region', i);
            unsafe = unsafe + 1;
        end
    end
   
end