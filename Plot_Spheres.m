function Plot_Spheres(X,Y,Z,R)
[x,y,z]=sphere;
colour(1:size(x),1:size(x),1)=1;
colour(1:size(x),1:size(x),2)=0.85;
colour(1:size(x),1:size(x),3)=0;
mesh(R*x+X,R*y+Y,R*z+Z,colour)
hold on
grid on
end