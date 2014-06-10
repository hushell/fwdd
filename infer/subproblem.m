function [nu] = subproblem(omega, edgeStruct, t1_id, t2_id)

% [~, t1_id] = intersect(edgeStruct.edgeEnds, omega.t1.edgeStruct.edgeEnds, 'rows');
% [~, t2_id] = intersect(edgeStruct.edgeEnds, omega.t2.edgeStruct.edgeEnds, 'rows');

omega_edge_T1 = omega.t1.edgePot(:,:,t1_id);
omega_edge_T2 = omega.t2.edgePot(:,:,t2_id);

[nodes_T1,edges_T1] = UGM_Infer_Tree(omega.t1.nodePot,omega_edge_T1,omega.t1.edgeStruct);
[nodes_T2,edges_T2] = UGM_Infer_Tree(omega.t2.nodePot,omega_edge_T2,omega.t2.edgeStruct);


% TODO: edgeStruct.nStates is not scalar
nS = edgeStruct.nStates(1);
nu_edge_T1 = zeros(nS,nS,edgeStruct.nEdges);
nu_edge_T2 = zeros(nS,nS,edgeStruct.nEdges);

nu_edge_T1(:,:,t1_id) = edges_T1;
nu_edge_T2(:,:,t2_id) = edges_T2;

nu.t1.nodes = nodes_T1;
nu.t1.edges = nu_edge_T1;
nu.t2.nodes = nodes_T2;
nu.t2.edges = nu_edge_T2;

end