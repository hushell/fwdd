function [nodeBel, edgeBel, logZ, gaps] = UGM_Infer_FWDD(nodePot, edgePot, edgeStruct, ...
    nodePot_T1, edgePot_T1, edgeStruct_T1, nodePot_T2, edgePot_T2, edgeStruct_T2)
    %nodePotTrees, edgePotTrees, edgeStructTrees)
%% parameters
MAX_ITERS = 100;
gap_threshold = 1e+2;

%% initialization
iter = 1;
EXIT_CONDITION = 0;

nS = edgeStruct.nStates(1);
nV = edgeStruct.nNodes;
nE = edgeStruct.nEdges;
% [nNodes,maxState] = size(nodePot);
% nEdges = size(edgePot,3);
% edgeEnds = edgeStruct.edgeEnds;
% nStates = edgeStruct.nStates;
% nS = maxState;

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
omega_vec = [omega.t1.nodePot(:);omega.t1.edgePot(:);omega.t2.nodePot(:);omega.t2.edgePot(:)];


%% main loop
gaps = zeros(1,MAX_ITERS);
while ~EXIT_CONDITION && iter < MAX_ITERS
    fprintf('iter %d...\n', iter);
    %% direction-finding subproblem
    
    % solve subproblem
    tic
    fprintf('\tsolving subproblem...\n');
    [nu, nodeBel, edgeBel] = subproblem(omega, edgeStruct, t1_id, t2_id);
    fprintf('done (%.2fs)\n', toc);
    nu_vec = [nu.t1.nodes(:); nu.t1.edges(:); nu.t2.nodes(:); nu.t2.edges(:)];
    
    logZ_t = bethe_free_energy(nodeBel, edgeBel, nodePot, edgePot, edgeStruct);
    fprintf('>>> current logZ: %g\n', logZ_t);
    
    n = length(nu_vec);
    n1 = length(nu.t1.nodes(:));
    n2 = length(nu.t1.edges(:));
    n3 = length(nu.t2.nodes(:));
    n4 = length(nu.t2.edges(:));
    s = zeros(n,1);
    
    % solve LB optimization problem
    tic
    fprintf('\tsolving lower bound...\n');

    cvx_begin quiet
        variable s(n)
        minimize( nu_vec'* s )
        subject to
            s(1:n1) + s(1+n1+n2:n1+n2+n3) == nodePot(:);
            s(1+n1:n1+n2) + s(1+n1+n2+n3:end) == edgePot(:);
            s >= 0;
    cvx_end
    
    fprintf('done (%.2fs)\n', toc);
    
    %% check exit condition
    duality_gap = (omega_vec - s)'*nu_vec;
    gaps(iter) = duality_gap;
    fprintf('*** current gap: %g, lower bound: %g\n', duality_gap, logZ_t - duality_gap);
    if duality_gap < gap_threshold
        EXIT_CONDITION = 1; % TODO gap isn't like what I expected
        break;
    end
    
    %% step size determination
    stepSize = 2 / (iter + 2);
    
    %% update step    
    omega_vec = omega_vec + stepSize * (s - omega_vec);
    
    omega.t1.nodePot = reshape(omega_vec(1:n1), [nV,nS]);
    omega.t1.edgePot = reshape(omega_vec(1+n1:n1+n2), [nS,nS,nE]);
    omega.t2.nodePot = reshape(omega_vec(1+n1+n2:n1+n2+n3), [nV,nS]);
    omega.t2.edgePot = reshape(omega_vec(1+n1+n2+n3:end), [nS,nS,nE]);
    
    iter = iter + 1;
end

if iter >= MAX_ITERS
    fprintf('Warning: reached maximum number of iterations!\n');
end

%% set nodeBel, edgeBel, logZ based on omega
[~, nodeBel, edgeBel] = subproblem(omega, edgeStruct, t1_id, t2_id);

if nargout > 2
    logZ = bethe_free_energy(nodeBel, edgeBel, nodePot, edgePot, edgeStruct);
end

end

function logZ = bethe_free_energy(nodeBel, edgeBel, nodePot, edgePot, edgeStruct)
   % Compute Bethe free energy 
   % (Z could also be computed as normalizing constant for any node in the tree
   %    if unnormalized messages are used)
   V = edgeStruct.V;
   E = edgeStruct.E;
   nStates = edgeStruct.nStates;
   nNodes = edgeStruct.nNodes;
   nEdges = edgeStruct.nEdges;
   edgeEnds = edgeStruct.edgeEnds;
   
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
