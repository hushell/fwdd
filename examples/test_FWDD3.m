% test_FWDD
clear 
close all

% run('/home/hushell/working/softwares/cvx/cvx_setup')
% addpath(genpath('.'))

%% noisy data
load X.mat

% X = X(1:16,1:16);
% figure;
% imagesc(X);
% colormap gray
% title('Original X');

figure;
X = X + randn(size(X))/2;
imagesc(X);
colormap gray
title('Noisy X');

[nRows,nCols] = size(X);
nNodes = nRows*nCols;
nStates = 2;
nTrees = 2;
RANDOM_TREES = 0;

% Original grid graph
adj = sparse(nNodes,nNodes);

% Add Down Edges
ind = 1:nNodes;
exclude = sub2ind([nRows nCols],repmat(nRows,[1 nCols]),1:nCols); % No Down edge for last row
ind = setdiff(ind,exclude);
adj(sub2ind([nNodes nNodes],ind,ind+1)) = 1;

% Add Right Edges
ind = 1:nNodes;
exclude = sub2ind([nRows nCols],1:nRows,repmat(nCols,[1 nRows])); % No right edge for last column
ind = setdiff(ind,exclude);
adj(sub2ind([nNodes nNodes],ind,ind+nRows)) = 1;

% Add Up/Left Edges
adj = adj+adj';
edgeStruct = UGM_makeEdgeStruct(adj,nStates);

%% get spanning trees
% edgeAppears, spanningTrees, adj_trees, edgeStructTrees
if RANDOM_TREES % random spanning trees
    edgeEnds = edgeStruct.edgeEnds;
    nEdges = edgeStruct.nEdges;
    i = 0;
    edgeAppears = zeros(nEdges,1);
    spanningTrees = {};
    while 1
        i = i+1;
        tree = minSpan(nNodes,[edgeEnds rand(nEdges,1)]);
        edgeAppears = edgeAppears + tree;
        spanningTrees{end+1} = logical(tree);
        if all(edgeAppears > 0)
            break;
        end
    end
    nTrees = i;
    assert(nTrees == numel(spanningTrees));
	mu = edgeAppears/i;
    
    adj_trees = cell(1,numel(spanningTrees));
    edgeStructTrees = cell(1,numel(spanningTrees));
    
    for t = 1:nTrees
        disp(t);
        adj_trees{t} = sparse(nNodes,nNodes);
        edgeEnds_t = edgeEnds(spanningTrees{t},:);
        adj_trees{t}(sub2ind([nNodes nNodes],edgeEnds_t(:,1),edgeEnds_t(:,2))) = 1;
        adj_trees{t} = adj_trees{t} + adj_trees{t}';
        edgeStructTrees{t} = UGM_makeEdgeStruct(adj_trees{t},nStates);
        % DEBUG
        draw_adj(adj_trees{t});
        pause
    end
else % 2 trees
    adj_t1 = sparse(nNodes,nNodes);
    
    % Add Down Edges
    %ind = (nCols-1)*nRows+1:nRows*nCols;
    ind = 1:nRows;
    exclude = sub2ind([nRows nCols],repmat(nRows,[1 nCols]),1:nCols); % No Down edge for last row
    ind = setdiff(ind,exclude);
    adj_t1(sub2ind([nNodes nNodes],ind,ind+1)) = 1;

    % Add Right Edges
    ind = 1:nNodes;
    exclude = sub2ind([nRows nCols],1:nRows,repmat(nCols,[1 nRows])); % No right edge for last column
    ind = setdiff(ind,exclude);
    adj_t1(sub2ind([nNodes nNodes],ind,ind+nRows)) = 1;

    % Add Up/Left Edges
    %assert(tree_graph(adj_t1));
    adj_t1 = adj_t1+adj_t1';
      
    adj_t2 = sparse(nNodes,nNodes);
    
    % Add Down Edges
    ind = 1:nNodes;
    exclude = sub2ind([nRows nCols],repmat(nRows,[1 nCols]),1:nCols); % No Down edge for last row
    ind = setdiff(ind,exclude);
    adj_t2(sub2ind([nNodes nNodes],ind,ind+1)) = 1;

    % Add Right Edges
    ind = 1:nRows:nNodes;
    exclude = sub2ind([nRows nCols],1:nRows,repmat(nCols,[1 nRows])); % No right edge for last column
    ind = setdiff(ind,exclude);
    adj_t2(sub2ind([nNodes nNodes],ind,ind+nRows)) = 1;

    % Add Up/Left Edges
    %assert(tree_graph(full(adj_t2)));
    adj_t2 = adj_t2+adj_t2';
       
    % DEBUG
    adj2 = adj_t1 + adj_t2;
    adj2 = adj2 >= 1;
    assert(all(all(adj2 == adj)));
    
    edgeStruct_T1 = UGM_makeEdgeStruct(adj_t1,nStates);
    edgeStruct_T2 = UGM_makeEdgeStruct(adj_t2,nStates);
end


%% initial nodePotTrees and edgePotTrees
if RANDOM_TREES % Hand-picked parameters
    % original 
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
    
    % spanning trees
    nodePotTrees = cell(1,numel(edgeStructTrees));
    edgePotTrees = cell(1,numel(edgeStructTrees));
    for t = 1:nTrees
        Xstd = UGM_standardizeCols(reshape(X,[1 1 nNodes]),1);
        nodePotTrees{t} = zeros(nNodes,nStates);
        nodePotTrees{t}(:,1) = exp(-1-2.5*Xstd(:))/nTrees;
        nodePotTrees{t}(:,2) = 1/nTrees;
        
        [shared_edges, sh_id] = intersect(edgeStruct.edgeEnds, edgeStructTrees{t}.edgeEnds, 'rows');
        assert(length(sh_id) == edgeStructTrees{t}.nEdges);
        
        edgePotTrees{t} = zeros(nStates,nStates,edgeStructTrees{t}.nEdges);
        for e = 1:size(shared_edges,1)
            n1 = shared_edges(e,1);
            n2 = shared_edges(e,2);

            pot_same = exp(1.8 + .3*1/(1+abs(Xstd(n1)-Xstd(n2))));
            edgePotTrees{t}(:,:,e) = [pot_same 1;1 pot_same] ./ edgeAppears(sh_id(e));
        end
    end
    
else % Learned optimal sub-modular parameters
    % original 
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
    
    % spanning trees
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

        pot_same = exp(1.8 + .3*1/(1+abs(Xstd(n1)-Xstd(n2))));
        edgePot_T1(:,:,sh_id(ei)) = [pot_same 1;1 pot_same] ./nTrees;
        edgePot_T2(:,:,sh_id(ei)) = [pot_same 1;1 pot_same] ./nTrees;
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

%% baseline
% nodePot_T1 = nodePotTrees{1};
% nodePot_T2 = nodePotTrees{2};
% edgePot_T1 = edgePotTrees{1};
% edgePot_T2 = edgePotTrees{2};
% edgeStruct_T1 = edgeStructTrees{1};
% edgeStruct_T2 = edgeStructTrees{2};
[nodeBel_T1,edgeBel_T1,logZ_T1] = UGM_Infer_Tree(nodePot_T1,edgePot_T1,edgeStruct_T1);
[nodeBel_T2,edgeBel_T2,logZ_T2] = UGM_Infer_Tree(nodePot_T2,edgePot_T2,edgeStruct_T2);

figure;
imagesc(reshape(nodeBel_T1(:,2),nRows,nCols));
colormap gray
title('Tree1 Belief Propagation Estimates of Marginals');

figure;
imagesc(reshape(nodeBel_T2(:,2),nRows,nCols));
colormap gray
title('Tree2 Belief Propagation Estimates of Marginals');

%% FWDD
[nodeBel,edgeBel,logZ,gaps] = UGM_Infer_FWDD(nodePot,edgePot,edgeStruct,...
    nodePot_T1,edgePot_T1,edgeStruct_T1, nodePot_T2,edgePot_T2,edgeStruct_T2);

figure;
imagesc(reshape(nodeBel(:,2),nRows,nCols));
colormap gray
title('FWDD Estimates of Marginals');

[~, nodeLabels] = max(nodeBel,[],2);
figure;
imagesc(reshape(nodeLabels,nRows,nCols));
colormap gray
title('Max of FWDD Marginals');
