function plot3Dspheres(cdInit,numParticle,phys_cond,L)
for j=1:numParticle
    Plot_Spheres(cdInit(j,1),cdInit(j,2),cdInit(j,3),phys_cond.R*50);
    hold on
end
xlim([-L L]);
ylim([-L L]);
zlim([-L L]);
end