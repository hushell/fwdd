% test_FWDD
clear 
close all

getNoisyX_DD

nTrees = 2;

edgeStruct = UGM_makeEdgeStruct(adj,nStates);
edgeStruct_T1 = UGM_makeEdgeStruct(adj_t1,nStates);
edgeStruct_T2 = UGM_makeEdgeStruct(adj_t2,nStates);

% initial nodePot and edgePot
if 0 % Hand-picked parameters
else % Learned optimal sub-modular parameters
    Xstd = UGM_standardizeCols(reshape(X,[1 1 nNodes]),1);
    nodePot = zeros(nNodes,nStates);
    nodePot(:,1) = exp(-1-2.5*Xstd(:));
    nodePot(:,2) = 1;

    edgePot = zeros(nStates,nStates,edgeStruct.nEdges);
    for e = 1:edgeStruct.nEdges
        n1 = edgeStruct.edgeEnds(e,1);
        n2 = edgeStruct.edgeEnds(e,2);

        pot_same = exp(1.8 + .3*1/(1+abs(Xstd(n1)-Xstd(n2))));
        edgePot(:,:,e) = [pot_same 1;1 pot_same];
    end
    
    Xstd = UGM_standardizeCols(reshape(X,[1 1 nNodes]),1);
    nodePot_T1 = zeros(nNodes,nStates);
    nodePot_T1(:,1) = exp(-1-2.5*Xstd(:))/nTrees;
    nodePot_T1(:,2) = 1/nTrees;
    nodePot_T2 = zeros(nNodes,nStates);
    nodePot_T2(:,1) = exp(-1-2.5*Xstd(:))/nTrees;
    nodePot_T2(:,2) = 1/nTrees;

    edgePot_T1 = zeros(nStates,nStates,edgeStruct_T1.nEdges);
    edgePot_T2 = zeros(nStates,nStates,edgeStruct_T2.nEdges);
    
    [shared_edges, sh_id] = intersect(edgeStruct_T1.edgeEnds, edgeStruct_T2.edgeEnds, 'rows');
    
    % shared edges
    for ei = 1:size(shared_edges,1)
        n1 = shared_edges(ei,1);
        n2 = shared_edges(ei,2);

        pot_same = exp(1.8 + .3*1/(1+abs(Xstd(n1)-Xstd(n2))))/nTrees;
        edgePot_T1(:,:,sh_id(ei)) = [pot_same 1;1 pot_same];
        edgePot_T2(:,:,sh_id(ei)) = [pot_same 1;1 pot_same];
    end
    
    % unique edges in T1
    for e = setdiff(1:edgeStruct_T1.nEdges, sh_id)
        n1 = edgeStruct_T1.edgeEnds(e,1);
        n2 = edgeStruct_T1.edgeEnds(e,2);

        pot_same = exp(1.8 + .3*1/(1+abs(Xstd(n1)-Xstd(n2))));
        edgePot_T1(:,:,e) = [pot_same 1;1 pot_same];
    end
    
    % unique edges in T2
    for e = setdiff(1:edgeStruct_T2.nEdges, sh_id)
        n1 = edgeStruct_T2.edgeEnds(e,1);
        n2 = edgeStruct_T2.edgeEnds(e,2);

        pot_same = exp(1.8 + .3*1/(1+abs(Xstd(n1)-Xstd(n2))));
        edgePot_T2(:,:,e) = [pot_same 1;1 pot_same];
    end
end


% [nodeBel_T1,edgeBel_T1,logZ_T1] = UGM_Infer_Tree(nodePot_T1,edgePot_T1,edgeStruct_T1);
% [nodeBel_T2,edgeBel_T2,logZ_T2] = UGM_Infer_Tree(nodePot_T2,edgePot_T2,edgeStruct_T2);
% 
% figure;
% imagesc(reshape(nodeBel_T1(:,2),nRows,nCols));
% colormap gray
% title('Tree1 Belief Propagation Estimates of Marginals');
% 
% figure;
% imagesc(reshape(nodeBel_T2(:,2),nRows,nCols));
% colormap gray
% title('Tree2 Belief Propagation Estimates of Marginals');

%% MICHAEL
[nodeBel,edgeBel,logZ] = UGM_Infer_Frank_Wolfe(nodePot,edgePot,edgeStruct,...
    nodePot_T1,edgePot_T1,edgeStruct_T1, nodePot_T2,edgePot_T2,edgeStruct_T2);

figure;
imagesc(reshape(nodeBel(:,2),nRows,nCols));
colormap gray
title('Full Belief Propagation Estimates of Marginals');