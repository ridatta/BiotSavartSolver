clc; close all; clear;
% Biot-Savart Code to calculate magnetic field generated by Coils
% R. Datta, June 2022

% SET the Total current 
I0 = 10e3; % [A]

% (1) Create a Helmholtz coil
% This will create a coil with 2 sections centered at +/- L/2
% Each section has a length (in z direction) of w and N turns
% R is the radius of the coils

N = 3; % number of turns in each section
w = 7e-3; % [m], coil width
R = 30e-3; % [m] coil radius
L = 66e-3; % [m] coil separation
num_el = 50; % Number of elements for descretization for each coil, total # elemets is 2 * N * num_el
[elmc, elm, I, dl] = createHelmholtzCoil(N,R,L,w,I0,num_el);% this function generates the lements

% (2) Show the coil
figure
plot3(elmc(:,1)*1e3,elmc(:,2)*1e3,elmc(:,3)*1e3,'mo'); hold on; % plot element centers
plot3(elm(:,1)*1e3,elm(:,2)*1e3,elm(:,3)*1e3,'g-','linewidth',2); hold on; % plot element centers
axis equal
xlim([-1,1]*1.25*R*1e3);ylim([-1,1]*1.25*R*1e3);
xlabel('x [mm]'); 
ylabel('y [mm]'); 
zlabel('z [mm]'); 

% (3) Caluclating the field 
pts = [0,0,0]; % [m] where to calculate field
plot3(pts(1)*1e3,pts(2)*1e3,pts(3)*1e3,'xr'); % show point on figure

% view(2)
B = getB(pts,elmc,I,dl); % this function calculates the field

% (4) Output
disp('Calculated Field [T] = ');
display(B);
load physicalConstants-SI.mat mu0
Bth = sqrt(64/125) * mu0 * N * I0 / R; % compare to analytical relation for a Helmholtz coil for which R = L
fprintf('Theoretical Field for an IDEAL Helmholtz Coil [T] = %0.3f\n', Bth);

% For multiple points use this code
% x = linspace(-R*1.5,R*1.5,10); y = x; z = linspace(0,l*1.5,10);
% [xx,yy,zz] = meshgrid(x,y,z);
% P = [xx(:), yy(:), zz(:)]; B = zeros(length(xx(:)),3);
% f = waitbar(0,'Please wait...');
% for ii =1:length(xx(:))
%     pts = P(ii,:);
%     B(ii,:) = getB(pts,elmc,I,dl);
%     waitbar(ii/length(xx(:)),f,'Calculating....');
% end
% 
% quiver3(xx(:)*1e3,yy(:)*1e3,zz(:)*1e3,B(:,1),B(:,2),B(:,3),'Color','b','AutoScaleFactor',20);

% FUNCTIONS

function [elmc, elm, I, dl] = createHelmholtzCoil(N,R,L,w,I0,num_el)
% N = # turns
% R = coil radius [m]
% L = coil separation [m]
% w = coil width [m]
% I0 = total current [A]
% num_el = # of discretization elements for each single coil
    zpos = linspace(-w/2,w/2,N); zpos = [zpos+L/2, zpos-L/2];
    elmc = []; elm = []; I = [];
    for ii =1:length(zpos)
        [xc,e,dl,i] = createCoil(num_el,R,zpos(ii),I0);
        elmc = [elmc; xc];
        elm = [elm; e];
        I = [I; i];
    end
end

function [elmc, elm, I, dl] = createNcoil(N,l,R,num_el,I0)
    zpos = linspace(0,l,N);
    elmc = []; elm = []; I = [];
    for ii =1:length(zpos)
        [xc,e,dl,i] = createCoil(num_el,R,zpos(ii),I0);
        elmc = [elmc; xc];
        elm = [elm; e];
        I = [I; i];
    end
end

function out = getB(pts,elmc,I,dl)
    % Do biot savart
    load physicalConstants-SI.mat mu0
    B = 0;
    for ii = 1:size(elmc,1) % Cycle each elem
        r = pts - elmc(ii,:); % vector from element dl to point P
        h1 = quiver3(elmc(ii,1)*1e3,elmc(ii,2)*1e3,elmc(ii,3)*1e3,r(1)*1e3,r(2)*1e3,r(3)*1e3,'Color','r');
        h2 = quiver3(elmc(ii,1)*1e3,elmc(ii,2)*1e3,elmc(ii,3)*1e3,...
            I(ii,1)/norm(I(ii,:),2),I(ii,2)/norm(I(ii,:),2),I(ii,3)/norm(I(ii,:),2),'Color','b','AutoScaleFactor',0.001);
        % Magnetic field
        dB = mu0 / (4 * pi) * cross(I(ii,:),r) / norm(r,2)^3 * dl;
        B = B + dB;
        pause(0.01);
        delete(h1);
        delete(h2);
    end
    out = B; % [1 x 3] vector
end



function [elmc,elm,dl,I] = createCoil(num_el,R,zpos,I0)
    % Creates a coil with radius R,  num_el elements at z position z-pos
    % Number of elements
    % R = coil radius [m]
    % zpos = z position [m]
    % I0 = current in coil [A]
    th = linspace(0,360-360/num_el,num_el);
    xc = R .* cosd(th); yc = R .* sind(th); zc = zpos* ones(size(xc)); %[m]
    elmc = [xc', yc', zc']; % element centers, [m]
    dl = 2 * R * tan(pi / num_el); 
    I = I0 * [-1*sind(th'), cosd(th'), 0*th'];
    
    th = [th, 360];
    xc = R .* cosd(th + 0.5*360/num_el); yc = R .* sind(th + 0.5*360/num_el); zc = zpos* ones(size(xc)); % plot elements
    elm = [xc', yc', zc'];
    
end