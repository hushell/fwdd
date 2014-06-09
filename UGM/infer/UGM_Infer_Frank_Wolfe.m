function [nodeBel, edgeBel, logZ] = UGM_Infer_Frank_Wolfe(nodePot, edgePot, edgeStruct, ...
    nodePot_T1, edgePot_T1, edgeStruct_T1, nodePot_T2, edgePot_T2, edgeStruct_T2)
%UGM_Infer_Frank_Wolfe
% INPUT
% Original Problem:
% nodePot(node,class)
% edgePot(class,class,edge) where e is referenced by V,E (must be the same
% between feature engine and inference engine)
%
% Subproblems T1, T2: same except with _T1 or _T2
% Subproblems must be trees that summed together get back the original
%
% OUTPUT
% nodeBel(node,class) - marginal beliefs
% edgeBel(class,class,e) - pairwise beliefs
% logZ - negative of free energy

%% ***** Frank Wolfe Algorithm *****

%% parameters
MAX_ITERS = 5;

%% initialization
iter = 0;
EXIT_CONDITION = 0;

omega.t1.nodePot = nodePot_T1;
omega.t1.edgePot = edgePot_T1;
omega.t1.edgeStruct = edgeStruct_T1;
omega.t2.nodePot = nodePot_T2;
omega.t2.edgePot = edgePot_T2;
omega.t2.edgeStruct = edgeStruct_T2;

%% main loop
while ~EXIT_CONDITION && iter < MAX_ITERS
    fprintf('iter %d...\n', iter);
    %% direction-finding subproblem
    
    % solve subproblem
    tic
    fprintf('\tsolving subproblem...');
    [nu, omega_full] = subproblem(omega, edgeStruct);
    fprintf('done (%.2fs)\n', toc);
    nu_vec = [nu.t1.nodes(:); nu.t1.edges(:); nu.t2.nodes(:); nu.t2.edges(:)];
    n = length(nu_vec);
    n1 = length(nu.t1.nodes(:));
    n2 = length(nu.t1.edges(:));
    n3 = length(nu.t2.nodes(:));
    s = zeros(n, 1);
    
    % solve optimization problem
    tic
    fprintf('\tsolving direction...');
    cvx_begin quiet
        variable s(n);
        minimize( s' * nu_vec );
        subject to
            s_t1_nodes = s(1:n1);
            s_t1_edges = s(1+n1:n1+n2);
            s_t2_nodes = s(1+n1+n2:n1+n2+n3);
            s_t2_edges = s(1+n1+n2+n3:end);
            
            s_t1_nodes + s_t2_nodes == nodePot;
            s_t1_edges + s_t2_edges == edgePot;
    cvx_end
    fprintf('done (%.2fs)\n', toc);
    
    % unpack s variable
    s_t1_nodes = s(1:n1);
    s_t1_edges = s(1+n1:n1+n2);
    s_t2_nodes = s(1+n1+n2:n1+n2+n3);
    s_t2_edges = s(1+n1+n2+n3:end);
    
    %% step size determination
    stepSize = 2 / (iter + 2);
    
    %% update step
    omega.t1.nodePot = omega.t1.nodePot + stepSize * (s_t1_nodes - omega.t1.nodePot);
    omega.t1.edgePot = omega.t1.edgePot + stepSize * (s_t1_edges - omega.t1.edgePot);
    omega.t2.nodePot = omega.t1.nodePot + stepSize * (s_t2_nodes - omega.t1.nodePot);
    omega.t2.edgePot = omega.t1.edgePot + stepSize * (s_t2_edges - omega.t1.edgePot);
    
    iter = iter + 1;
    
    %% check exit condition
    EXIT_CONDITION = 0; % TODO
end

if iter >= MAX_ITERS
    fprintf('Warning: reached maximum number of iterations!\n');
end

%% set nodeBel, edgeBel, logZ based on omega
[nodeBel,edgeBel,logZ] = UGM_Infer_Tree(omega.t1.nodePot,omega.t1.edgePot,omega.t1.edgeStruct);

end
