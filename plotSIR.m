% function plotSIR
% agents gives the number of agents at each ROI (or vertex) at any time
% e.g. Rnor_all, 20000

function plotSIR(agents, time)
lh_verts = load('surface_data.mat').lh_verts;
lh_faces = load('surface_data.mat').lh_faces;
lh_aparc = load('surface_data.mat').lh_aparc;
figure;
plotSurfaceROIBoundary( struct('vertices', lh_verts, 'faces', lh_faces), lh_aparc, agents((1:34),time), 'faces', parula(100), 2);
axis equal; axis off; colorbar;
end
