clc; close all; clear;
% Biot-Savart Code to calculate magnetic field generated by Helical Coils
% R. Datta, June 2022

% Total current 
I0 = 0.7e6; % [A]

% (1) Create a Helmholtz coil
% This will create a helical coil 

R = 2e-3; % [m]
z_sep = 1e-3; % [m] vertical separation between turns
c = z_sep / (2 * pi);
num_el = 100;
N = 6; % number of turns
[elmc, elm, dl, I] = createHelicalCoil(num_el,R,z_sep,N,I0);% this function generates the lements

% (2) Show the coil
figure
plot3(elmc(:,1)*1e3,elmc(:,2)*1e3,elmc(:,3)*1e3,'mo'); hold on; % plot element centers
% plot3(elm(:,1)*1e3,elm(:,2)*1e3,elm(:,3)*1e3,'g-','linewidth',2); hold on; % plot element centers
axis equal
xlim([-1,1]*1.25*R*1e3);ylim([-1,1]*1.25*R*1e3);
xlabel('x [mm]'); 
ylabel('y [mm]'); 
zlabel('z [mm]'); 

% (3) Caluclating the field 

% For a single point
pts = [0,0,z_sep * N /2]; % [m] where to calculate field
plot3(pts(1)*1e3,pts(2)*1e3,pts(3)*1e3,'xr'); % show point on figure

% view(2)
B = getB(pts,elmc,I,dl); % this function calculates the field

% (4) Output
disp('Calculated Field [T] = ');
display(B);
load physicalConstants-SI.mat mu0
display('Analytical:')
fprintf('B_solenoid[T] = %0.3f\n', mu0 * N/(z_sep*N) * I0);
fprintf('B_coil[T] = %0.3f\n', mu0 * N * I0 / (2 * R));

% For multiple points use this code
% x = linspace(-R*2,R*2,10); y = x; z = N * z_sep / 2 * ones(size(x));
% [xx,yy,zz] = meshgrid(x,y,z);
y = linspace(R,R+10e-3,50);
x = 0 * y;
z = N * z_sep / 2 * ones(size(x));
P = [x', y', z']; B = zeros(length(x),3);
% f = waitbar(0,'Please wait...');
for ii =1:length(x)
    pts = P(ii,:);
    B(ii,:) = getB(pts,elmc,I,dl);
%     waitbar(ii/length(x),f,'Calculating....');
end

% quiver3(x'*1e3,y'*1e3,z'*1e3,B(:,1),B(:,2),B(:,3),'Color','b','AutoScaleFactor',20);

% straight cathode
B_str = mu0 * I0 ./ (2 * pi * y); 



figure
% plot(y*1e3,B(:,2),DisplayName='By',linewidth=2); hold on;
plot(y*1e3,B(:,3),DisplayName='|B_z|',linewidth=3); hold on;
plot(y*1e3,abs(B(:,1)),DisplayName='|B_\theta|',linewidth=3,linestyle='--'); hold on;
plot(y*1e3,sqrt(B(:,1).^2 + B(:,3).^2),DisplayName='|B|',linewidth=3,linestyle='-',color='r'); hold on;
plot(y*1e3,abs(B_str),DisplayName='|B (straight K)|',linewidth=3,linestyle='-',color=[0.1 0.5 0.1]); hold on;
xlabel('Distance from coil center [mm]'); ylabel('B (T)')
ylim([0,50]); 
xlim([R*1e3+1,R*1e3+10])
formatPlots(600);
grid on;
title({['Cathode R [mm] = ', num2str(R*1e3), ', I [MA] = ', num2str(I0/1e6)],...
    ['\Delta z [mm] = ' num2str(z_sep*1e3) ]})
set(gcf,'Position',[0   0   524   317]*2);
saveas(gcf,['figures/solenoidal_cathode' num2str(randi(500)), '.png']);



% FUNCTIONS


function out = getB(pts,elmc,I,dl)
    % Do biot savart
    load physicalConstants-SI.mat mu0
    B = 0;
    for ii = 1:size(elmc,1) % Cycle each elem
        r = pts - elmc(ii,:); % vector from element dl to point P
%         h1 = quiver3(elmc(ii,1)*1e3,elmc(ii,2)*1e3,elmc(ii,3)*1e3,r(1)*1e3,r(2)*1e3,r(3)*1e3,'Color','r');
%         h2 = quiver3(elmc(ii,1)*1e3,elmc(ii,2)*1e3,elmc(ii,3)*1e3,...
%             I(ii,1)/norm(I(ii,:),2),I(ii,2)/norm(I(ii,:),2),I(ii,3)/norm(I(ii,:),2),'Color','b','AutoScaleFactor',10);
        % Magnetic field
        dB = mu0 / (4 * pi) * cross(I(ii,:),r) / norm(r,2)^3 * dl;
        B = B + dB;
%         pause(0.01);
%         delete(h1);
%         delete(h2);
    end
    out = B; % [1 x 3] vector
end



function [elmc,elm,dl,I] = createHelicalCoil(num_el,R,zsep,nturns,I0)
    % Creates a coil with radius R,  num_el elements at z position z-pos
    % Number of elements per turn
    % R = coil radius [m]
    % zsep = z separation between turns [m]
    % nturns = number of turns
    % I0 = current in coil [A]
    th = linspace(0,nturns*2*pi,num_el*nturns);
    c = zsep / (2 * pi);
    xc = R .* sin(th); yc = R .* cos(th); zc = c * th; %[m]
    elmc = [xc', yc', zc']; % element centers, [m]
    L = nturns * sqrt( (2 * pi * R).^2 + zsep.^2);  % helix length
    dl = L / (num_el * nturns); 
    I = I0 * [R*cos(th'), -R*sin(th'), c*ones(size(th'))] / sqrt(R^2 + c^2); % current tangent to curve; check norm(I(ii,:),2) should be I0
    
%     th = [th, 360];
%     xc = R .* cosd(th + 0.5*360/num_el); yc = R .* sind(th + 0.5*360/num_el); zc = zpos* ones(size(xc)); % plot elements
    elm = [xc', yc', zc'];
    
end