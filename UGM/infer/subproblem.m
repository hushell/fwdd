function [nu] = subproblem(omega, edgeStruct)

%[shared_edges, sh_id] = intersect(omega.t1.edgeStruct.edgeEnds, omega.t1.edgeStruct.edgeEnds, 'rows');

[nodes_T1,edges_T1] = UGM_Infer_Tree(omega.t1.nodePot,omega.t1.edgePot,omega.t1.edgeStruct);
[nodes_T2,edges_T2] = UGM_Infer_Tree(omega.t2.nodePot,omega.t2.edgePot,omega.t2.edgeStruct);

[~, nt1_id] = setdiff(edgeStruct.edgeEnds, omega.t1.edgeStruct.edgeEnds, 'rows');
[~, nt2_id] = setdiff(edgeStruct.edgeEnds, omega.t2.edgeStruct.edgeEnds, 'rows');

% TODO: edgeStruct.nStates is not scalar
edges_N_T1 = zeros(edgeStruct.nStates,edgeStruct.nStates,length(nt1_id));
edges_N_T2 = zeros(edgeStruct.nStates,edgeStruct.nStates,length(nt2_id));

edges_A_T1 = zeros(edgeStruct.nStates,edgeStruct.nStates,edgeStruct.nEdges);
edges_A_T2 = zeros(edgeStruct.nStates,edgeStruct.nStates,edgeStruct.nEdges);

edges_A_T1(:,:,nt1_id) = edges_N_T1;
edges_A_T2(:,:,nt2_id) = edges_N_T2;

edges_A_T1(:,:,setdiff(1:edgeStruct.nEdges, nt1_id)) = edges_T1;
edges_A_T2(:,:,setdiff(1:edgeStruct.nEdges, nt2_id)) = edges_T2;

nu.t1.nodes = nodes_T1;
nu.t1.edges = edges_A_T1;
nu.t2.nodes = nodes_T2;
nu.t2.edges = edges_A_T2;

end