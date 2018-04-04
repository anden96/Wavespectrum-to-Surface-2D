
x = ()linspace(0,2*pi,64);
y = (1:64);%0.25.*sin(4*x);
z = repmat(y,length(x),1);
%z = 1;
[x, y]=meshgrid(x,y);
plot(x,y)

surf(z)
SurfaceSave = [x y z];
dlmwrite('mySinus.txt',SurfaceSave,'delimiter',' ');

