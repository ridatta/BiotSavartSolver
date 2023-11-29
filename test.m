R = 6e-2;
z_sep = 6e-2;
c = z_sep / (2 * pi);


t = linspace(0,4*2*pi,4*100);

x = R * cos(t);
y = R * sin(t);
z = c * t; 


figure
plot3(x*100,y*100,z*100);
axis equal