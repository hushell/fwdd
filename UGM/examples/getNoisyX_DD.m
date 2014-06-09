clear
close all

load X.mat

[nRows,nCols] = size(X);
nNodes = nRows*nCols;
nStates = 2;

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
% edgeStruct = UGM_makeEdgeStruct(adj,nStates);

% get spanning trees
if 0 % random spanning trees
    tree = minSpan(nNodes,[edgeStruct.edgeEnds rand(edgeStruct.nEdges,1)]);
else % 2 trees
    adj_t1 = sparse(nNodes,nNodes);
    
    % Add Down Edges
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
    assert(tree_graph(adj_t1));
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
    assert(tree_graph(full(adj_t2)));
    adj_t2 = adj_t2+adj_t2';
    
    
    % DEBUG
    adj2 = adj_t1 + adj_t2;
    adj2 = adj2 >= 1;
    assert(all(all(adj2 == adj)));
    
    %edgeStruct = UGM_makeEdgeStruct(adj_t1,nStates);
end


