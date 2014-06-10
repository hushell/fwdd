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
MAX_ITERS = 10;

%% initialization
iter = 0;
EXIT_CONDITION = 0;

% nS = edgeStruct.nStates(1);
nV = edgeStruct.nNodes;
nE = edgeStruct.nEdges;
[nNodes,maxState] = size(nodePot);
nEdges = size(edgePot,3);
edgeEnds = edgeStruct.edgeEnds;
nStates = edgeStruct.nStates;
V = edgeStruct.V;
E = edgeStruct.E;
nS = maxState;

[~, t1_id] = intersect(edgeStruct.edgeEnds, edgeStruct_T1.edgeEnds, 'rows');
[~, t2_id] = intersect(edgeStruct.edgeEnds, edgeStruct_T2.edgeEnds, 'rows');

edges_A_T1 = zeros(nS,nS,edgeStruct.nEdges);
edges_A_T2 = zeros(nS,nS,edgeStruct.nEdges);

edges_A_T1(:,:,t1_id) = edgePot_T1;
edges_A_T2(:,:,t2_id) = edgePot_T2;

omega.t1.nodePot = nodePot_T1;
omega.t1.edgePot = edges_A_T1;
omega.t1.edgeStruct = edgeStruct_T1;
omega.t2.nodePot = nodePot_T2;
omega.t2.edgePot = edges_A_T2;
omega.t2.edgeStruct = edgeStruct_T2;


%% main loop
while ~EXIT_CONDITION && iter < MAX_ITERS
    fprintf('iter %d...\n', iter);
    %% direction-finding subproblem
    
    % solve subproblem
    tic
    fprintf('\tsolving subproblem...\n');
    nu = subproblem(omega, edgeStruct, t1_id, t2_id);
    fprintf('done (%.2fs)\n', toc);
    nu_vec = [nu.t1.nodes(:); nu.t1.edges(:); nu.t2.nodes(:); nu.t2.edges(:)];
    
    nu1 = nu.t1.nodes(:); 
    nu2 = nu.t1.edges(:); 
    nu3 = nu.t2.nodes(:); 
    nu4 = nu.t2.edges(:);
    
    n = length(nu_vec);
    n1 = length(nu.t1.nodes(:));
    n2 = length(nu.t1.edges(:));
    n3 = length(nu.t2.nodes(:));
    n4 = length(nu.t2.edges(:));
    s = zeros(n, 1);
    
    % solve optimization problem
    tic
    fprintf('\tsolving direction...\n');
%     cvx_begin %quiet
%           variables s1(n1) s2(n2) s3(n3) s4(n4)
%           minimize(  s1'*nu1 + s2'*nu2 + s3'*nu3 + s4'*nu4 )
%           subject to
%             s1 + s3 == nodePot(:);
%             s2 + s4 == edgePot(:);
%         variable s(n);
%         minimize( s' * nu_vec );
%         subject to
%             s_t1_nodes = s(1:n1);
%             s_t1_edges = s(1+n1:n1+n2);
%             s_t2_nodes = s(1+n1+n2:n1+n2+n3);
%             s_t2_edges = s(1+n1+n2+n3:end);
%             s_t1_nodes + s_t2_nodes == nodePot(:);
%             s_t1_edges + s_t2_edges == edgePot(:);
%             s(1:n1) + s(1+n1+n2:n1+n2+n3) == nodePot(:);
%             s(1+n1:n1+n2) + s(1+n1+n2+n3:end) == edgePot(:);
%     cvx_end
    
    t1_lid = zeros(1,edgeStruct.nEdges);
    t1_lid(t1_id) = 1;
    t2_lid = zeros(1,edgeStruct.nEdges);
    t2_lid(t2_id) = 1;
    
    Aeq_t1 = zeros(1,n2);
    Aeq_t1(1:4:end) = t1_lid;
    Aeq_t1(2:4:end) = t1_lid;
    Aeq_t1(3:4:end) = t1_lid;
    Aeq_t1(4:4:end) = t1_lid;
    
    Aeq_t2 = zeros(1,n4);
    Aeq_t2(1:4:end) = t2_lid;
    Aeq_t2(2:4:end) = t2_lid;
    Aeq_t2(3:4:end) = t2_lid;
    Aeq_t2(4:4:end) = t2_lid;
    
    Aeq = [sparse(diag(1:n1)), sparse(zeros(n1,n2)), sparse(diag(1:n3)), sparse(zeros(n1,n4)); 
        sparse(zeros(n2,n1)), sparse(diag(Aeq_t1)), sparse(zeros(n2,n3)), sparse(diag(Aeq_t2))];
    beq = [nodePot(:); edgePot(:)];
    
    [s,fval,exitflag] = linprog(nu_vec, [], [], Aeq, beq, zeros(n,1));
    fprintf('done (%.2fs)\n', toc);
    assert(exitflag == 1);
    
    % unpack s variable
    s_t1_nodes = s(1:n1);
    s_t1_edges = s(1+n1:n1+n2);
    s_t2_nodes = s(1+n1+n2:n1+n2+n3);
    s_t2_edges = s(1+n1+n2+n3:end);
    
    s_t1_nodes = reshape(s_t1_nodes, [nV,nS]);
    s_t1_edges = reshape(s_t1_edges, [nS,nS,nE]);
    s_t2_nodes = reshape(s_t2_nodes, [nV,nS]);
    s_t2_edges = reshape(s_t2_edges, [nS,nS,nE]);
    
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
omega_edge_T1 = omega.t1.edgePot(:,:,t1_id);
omega_edge_T2 = omega.t2.edgePot(:,:,t2_id);

[nodes_T1,edges_T1, logZ_T1] = UGM_Infer_Tree(omega.t1.nodePot,omega_edge_T1,omega.t1.edgeStruct);
[nodes_T2,edges_T2, logZ_T2] = UGM_Infer_Tree(omega.t2.nodePot,omega_edge_T2,omega.t2.edgeStruct);
% [nodeBel,edgeBel,logZ] = UGM_Infer_Tree(omega.t1.nodePot,omega.t1.edgePot,omega.t1.edgeStruct);

nodeBel = nodes_T1;

edgeBel = zeros(nS,nS,edgeStruct.nEdges);
edgeBel(:,:,t1_id) = edges_T1;
[~, loc] = ismember(setdiff(1:edgeStruct.nEdges, t1_id), t2_id);
edgeBel(:,:,t2_id(loc)) = edges_T2(:,:,loc);

if nargout > 2
   % Compute Bethe free energy 
   % (Z could also be computed as normalizing constant for any node in the tree
   %    if unnormalized messages are used)
   Energy1 = 0; Energy2 = 0; Entropy1 = 0; Entropy2 = 0;
   nodeBel = nodeBel+eps;
   edgeBel = edgeBel+eps;
   for n = 1:nNodes
      edges = E(V(n):V(n+1)-1);
      nNbrs = length(edges);

      % Node Entropy (can get divide by zero if beliefs at 0)
      Entropy1 = Entropy1 + (nNbrs-1)*sum(nodeBel(n,1:nStates(n)).*log(nodeBel(n,1:nStates(n))));

      % Node Energy
      Energy1 = Energy1 - sum(nodeBel(n,1:nStates(n)).*log(nodePot(n,1:nStates(n))));
   end
   for e = 1:nEdges
      n1 = edgeEnds(e,1);
      n2 = edgeEnds(e,2);

      % Pairwise Entropy (can get divide by zero if beliefs at 0)
      eb = edgeBel(1:nStates(n1),1:nStates(n2),e);
      Entropy2 = Entropy2 - sum(eb(:).*log(eb(:)));

      % Pairwise Energy
      ep = edgePot(1:nStates(n1),1:nStates(n2),e);
      Energy2 = Energy2 - sum(eb(:).*log(ep(:)));
   end
   F = (Energy1+Energy2) - (Entropy1+Entropy2);
   logZ = -F;
end

end
