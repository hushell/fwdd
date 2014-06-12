function draw_adj(adj)

R = sqrt( size(adj,1) );
[xloc,yloc] = meshgrid(1:R,1:R);

gplot(adj,[xloc(:) yloc(:)],'-*');
set(gca,'ydir','reverse');
set(gca,'xlim',[0 R]+0.5);
set(gca,'ylim',[0 R]+0.5);