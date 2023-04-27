% function plotSIR

figure;
plotSurfaceROIBoundary( struct('vertices', lh_verts, 'faces', lh_faces), lh_aparc, lh_aparc, 'faces', parula(100), 2);
axis equal; axis off; colorbar;
