function Plot_Spheres(X,Y,Z,R)
[x,y,z]=sphere;
colour=zeros(size(x));% 获得矩阵大小和x相同
for i=1:1:length(colour(1,:))
    for j=1:1:length(colour(:,1))
        colour(i,j,1)=1;
        colour(i,j,2)=0.85;
        colour(i,j,3)=0;   % 金色RGB 1 0.85 0
    end
end
mesh(R*x+X,R*y+Y,R*z+Z,colour)
hold on
grid on
end