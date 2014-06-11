function [nu, nodeBel, edgeBel] = subproblem(omega, edgeStruct, t1_id, t2_id)

% [~, t1_id] = intersect(edgeStruct.edgeEnds, omega.t1.edgeStruct.edgeEnds, 'rows');
% [~, t2_id] = intersect(edgeStruct.edgeEnds, omega.t2.edgeStruct.edgeEnds, 'rows');

omega_edgePot_T1 = omega.t1.edgePot(:,:,t1_id);
omega_edgePot_T2 = omega.t2.edgePot(:,:,t2_id);

[nodeBel_T1,edgeBel_T1,logZ_T1] = UGM_Infer_Tree(omega.t1.nodePot,omega_edgePot_T1,omega.t1.edgeStruct);
[nodeBel_T2,edgeBel_T2,logZ_T2] = UGM_Infer_Tree(omega.t2.nodePot,omega_edgePot_T2,omega.t2.edgeStruct);


% TODO: edgeStruct.nStates is not scalar
nS = edgeStruct.nStates(1);
nu_edge_T1 = zeros(nS,nS,edgeStruct.nEdges);
nu_edge_T2 = zeros(nS,nS,edgeStruct.nEdges);

nu_edge_T1(:,:,t1_id) = edgeBel_T1;
nu_edge_T2(:,:,t2_id) = edgeBel_T2;

nu.t1.nodes = nodeBel_T1;
nu.t1.edges = nu_edge_T1;
nu.t2.nodes = nodeBel_T2;
nu.t2.edges = nu_edge_T2;
nu.t1.logZ = logZ_T1;
nu.t2.logZ = logZ_T2;

if nargout > 1
%     nodeBel = nodeBel_T1;
% 
%     edgeBel = zeros(nS,nS,edgeStruct.nEdges);
%     edgeBel(:,:,t1_id) = edgeBel_T1;
%     [~, loc] = ismember(setdiff(1:edgeStruct.nEdges, t1_id), t2_id);
%     edgeBel(:,:,t2_id(loc)) = edgeBel_T2(:,:,loc);
    nodeBel = nodeBel_T1 + nodeBel_T2;
    nodeBel = bsxfun(@rdivide, nodeBel, sum(nodeBel,2));
    
    edgeBel = zeros(nS,nS,edgeStruct.nEdges);
    edgeBel(:,:,t1_id) = edgeBel_T1;
    edgeBel(:,:,t2_id) = edgeBel(:,:,t2_id) + edgeBel_T2;
    edgeBel = bsxfun(@rdivide, edgeBel, sum(sum(edgeBel)));
end

end