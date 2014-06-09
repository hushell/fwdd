%% Make noisy X
getNoisyX_DD

%% Tree-Reweighted Belief Propagation

% fprintf('Running tree-reweighted belief propagation for inference...\n');
% [nodeBelTRBP,edgeBelTRBP,logZTRBP] = UGM_Infer_TRBP(nodePot,edgePot,edgeStruct);
% 
% figure;
% imagesc(reshape(nodeBelTRBP(:,2),nRows,nCols));
% colormap gray
% title('Tree-Reweighted Belief Propagation Estimates of Marginals');
% % fprintf('(paused)\n');
% % pause
% 
% fprintf('Running tree-reweighted belief propagation and computing max of marginals\n');
% maxOfMarginalsTRBPdecode = UGM_Decode_MaxOfMarginals(nodePot,edgePot,edgeStruct,@UGM_Infer_TRBP);
% 
% figure;
% imagesc(reshape(maxOfMarginalsTRBPdecode,nRows,nCols));
% colormap gray
% title('Max of Tree-Reweighted Belief Propagation Marginals');

%%

fprintf('Running tree2 lief propagation for inference...\n');
[nodeBelTRBP,edgeBelTRBP,logZTRBP] = UGM_Infer_TRBP(nodePot,edgePot,edgeStruct);

figure;
imagesc(reshape(nodeBelTRBP(:,2),nRows,nCols));
colormap gray
title('tree2 Belief Propagation Estimates of Marginals');
% fprintf('(paused)\n');
% pause

fprintf('Running tree2 belief propagation and computing max of marginals\n');
maxOfMarginalsTRBPdecode = UGM_Decode_MaxOfMarginals(nodePot,edgePot,edgeStruct,@UGM_Infer_TRBP);

figure;
imagesc(reshape(maxOfMarginalsTRBPdecode,nRows,nCols));
colormap gray
title('Max of tree2 Belief Propagation Marginals');

