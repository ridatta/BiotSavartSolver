clc; close all; clear;

% Total current 
I0 = 1e3; % [A]

% create a coil
N =10;
l = 5e-2;
num_el = 50 
R = 1e-2; % [m] coil width
[elmc, elm, I, dl] = createNcoil(N,l,R,num_el,I0);


figure
plot3(elmc(:,1),elmc(:,2),elmc(:,3),'mo'); hold on; % plot element centers
plot3(elm(:,1),elm(:,2),elm(:,3),'g-','linewidth',2); hold on; % plot element centers
axis equal
xlim([-1,1]*1.25*R);ylim([-1,1]*1.25*R);

pts = [0,0,l/2]; 

plot3(pts(1),pts(2),pts(3),'xr');

% view(2)
B = getB(pts,elmc,I,dl);
disp('Calculated Field [T} = ');
display(B);
load physicalConstants-SI.mat mu0
fprintf('B_solenoid[T] = %0.3f\n', mu0 * N/l * I0);
fprintf('B_coil[T] = %0.3f\n', mu0 * N * I0 / (2 * R));

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
% quiver3(xx(:),yy(:),zz(:),B(:,1),B(:,2),B(:,3),'Color','b','AutoScaleFactor',20);

% FUNCTIONS

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
        h1 = quiver3(elmc(ii,1),elmc(ii,2),elmc(ii,3),r(1),r(2),r(3),'Color','r');
        h2 = quiver3(elmc(ii,1),elmc(ii,2),elmc(ii,3),...
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