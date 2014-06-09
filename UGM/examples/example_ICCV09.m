clear all

[nInstances,nNodes] = size(y);

%% Make edgeStruct
nStates = max(y);
adj = 
edgeStruct = UGM_makeEdgeStruct(adj,nStates);
nEdges = edgeStruct.nEdges;
maxState = max(nStates); 
