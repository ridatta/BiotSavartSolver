classdef ArrayBiotSavart
    % Use this class of functions to visualze the magnetic field in z-pinch
    % arrays
    properties
        I0 % Total current, MA
        wire_pos % Array containing wire positions
        wire_pos_cat % Array containing cathode wire positions
    end
    methods
        function obj = ArrayBiotSavart(I0)
            obj.I0 = I0; % current, MA
            obj.wire_pos = [];
            obj.wire_pos_cat = [];
        end
        function obj = createPlanarArray(obj,nwire,dx)
            % nwire = number of wires
            % dx = wire-separation
            x = 0:dx:(nwire-1)*dx; x = x - dx*(nwire-1) / 2;
            obj.wire_pos = [x', 0*x', 0*x']; % wire positions in array
        end
        function obj = createExplodingArray(obj,nwire,R0,nwire_cat,R0_cat)
            % nwire = number of wires
            % R0 = Radius
            % Wires
            th = 2*pi / nwire / 2 +  linspace(0,2*pi-2*pi/nwire,nwire)';
            x = R0 * cos(th); y = R0 * sin(th);
            obj.wire_pos = [x, y, 0*x];
            % Cathode
            th_cat = linspace(0,2*pi-2*pi/nwire_cat,nwire_cat)';
            x = R0_cat * cos(th_cat); y = R0_cat * sin(th_cat);
            obj.wire_pos_cat = [x, y, 0*x];
        end
        function obj = createImplodingArray(obj,nwire,R0)
            % nwire = number of wires
            % R0 = Radius
            th = linspace(0,2*pi-2*pi/nwire,nwire)';
            x = R0 * cos(th); y = R0 * sin(th);
            obj.wire_pos = [x, y, 0*x];
        end
        function obj = addCathode(obj,nwire_cat,dx_cat,del_y)
            % nwire_cat = number of wires cathode
            % dx_ca = wire separation cathode
            % del_y = chathode position from wires
            x = 0:dx_cat:nwire_cat*dx_cat; x = x - dx_cat*nwire_cat / 2;
            obj.wire_pos_cat = [x', 0*x' - del_y, 0*x'];
        end
        function out = getB(obj,N,pts)
            mu0 = 1.2566e-06;
            I = obj.I0 * [0, 0 ,1] / size(obj.wire_pos,1); % Current per wire, A
            dl = 1 / N; % [m]
            B = 0;
            for jj = 1:size(obj.wire_pos,1) % Cycle each wire
                wire_now = obj.wire_pos(jj,:);
                for ii = 1:N
                    cen = wire_now' + [0; 0; (ii-1)*dl]; % element center
                    r = pts' - cen; % vector from element dl to point P
                    % (5) Magnetic field
                    dB = mu0 / (4 * pi) * cross(I,r) / norm(r,2)^3 * dl;
                    B = B + dB;
                end
            end
            % (5b)Cathode Current
            I = -obj.I0 * [0, 0 ,1] / size(obj.wire_pos_cat,1); % Current per wire on cathode, A
            for jj = 1:size(obj.wire_pos_cat,1) % Cycle each wire on cathode
                wire_now_cat = obj.wire_pos_cat(jj,:);
                for ii = 1:N % cycle each differential current element
                    cen = wire_now_cat' + [0; 0; (ii-1)*dl]; % element center
                    r = pts' - cen; % vector from element dl to point P
                    % (5) Magnetic field
                    dB = mu0 / (4 * pi) * cross(I,r) / norm(r,2)^3 * dl;
                    B = B + dB;
                end
            end
            out = B; % T
        end
        function showArray(obj)
            for jj = 1:size(obj.wire_pos,1) % Cycle each wire
                wire_now = obj.wire_pos(jj,:);
                plot3([1,1]*wire_now(1),[1,1]*wire_now(2),[wire_now(3),1],...
                    'm','LineWidth',3); hold on;
            end
            for jj = 1:size(obj.wire_pos_cat,1) % Cycle each wire
                wire_now = obj.wire_pos_cat(jj,:);
                plot3([1,1]*wire_now(1),[1,1]*wire_now(2),[wire_now(3),1],...
                    'g','LineWidth',3); hold on;
            end
            zlim([0,1]*1); xlim([-10e-3,10e-3]); ylim([-10e-3,10e-3]);
            xlabel('x'); ylabel('y'); zlabel('z'); grid on; view([120 20]);
        end
        
    end
end
