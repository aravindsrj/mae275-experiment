plume4occupancy;

% Sample waypoints
    % First Argument:
        % 1  ---> out of plume (static)
        % 2  ---> in plume (static)
        % 3  ---> out of plume (dynamic)
        % 4  ---> in plume (dynamic)
        % 5  ---> in plume, then out plume
        % 6  ---> out plume, then in plume
     
    % Second Argument -> Pertaining to different environment and plumes
        % 1  ---> 2000x1000 environment with stability type 'F'
        % 2  ---> 500x400 environment with stability type 'A'
    
% Pro tip: Do NOT change the second argument here -> Won't work without changing some other things    
P_uav = waypoints(2,2);

%% Plume and Simulation characteristics
Duration   = length(P_uav); %800;
dt         = 1;                    % Time step
N          = length(0:dt:Duration);
Vwind      = [props.U+0.01,0.01];  % Constant wind vector
Wind       = zeros(N,2);           % Wind vector
sx         = 10;                   % Standard deviations may have to be modified
sy         = 5;                   % Standard deviations may have to be modified  
mu         = 0.9;                  % Sensor accuracy
lambda     = 25;                   % (meters); can be removed; a parameter to reduce horizon time so that plume estimate is not touching boundary
%tk_max0    = (gridMap.xlims(2)-gridMap.xlims(1) - lambda)/Vwind(1); % Max horizon time
L          = 1;                    % lower time bound
K          = 0;                    % Upper time bound; or current time step
T          = zeros(1,N);
plume_start = 400; % Time between the start of the plume and the start of the mission
gamma_wt   = 10;   % A weighting variable for gamma so that its significance is increased

pos        = 0;    % UAV position pointer

%% Environment dimensions
global Lx Ly m
Lx         = 40;    % Length of probability cells
Ly         = 20;    % Length of probability cells
m          = floor((gridMap.xlims(2)-gridMap.xlims(1))/Lx) + 1;    % number of cells in x-direction 
n          = floor((gridMap.ylims(2)-gridMap.ylims(1))/Ly) + 1;    % number of cells in y-direction    
M          = m*n;   % Total cells

x          = gridMap.xlims(1):Lx:gridMap.xlims(2);
y          = gridMap.ylims(1):Ly:gridMap.ylims(2);
[Y,X]      = meshgrid(y,x);
xcell      = X(:); % 
ycell      = Y(:); %

%% Prob map Initialization
alpha      = (1/M)*ones(M,N+plume_start);   % Initializing probability map
Sij        = zeros(M,1);
beta       = zeros(M,M);          % Detection map
gamma      = (1/M)*ones(M,M);     % Non-detection map
figure(2)
view(0,90);
surf(X,Y,reshape(alpha(:,K+1),[m,n]))

%%

for Time = dt:dt:Duration+plume_start
    
    K = K + 1;
    index = find_index(P_uav(pos+1,1),P_uav(pos+1,2), gridMap.ylims(1));
    T(K+1) = Time;
%     if Time > tk_max0
%         L = L + 1;
%     end
    tk_max = find_tkmax(P_uav(pos+1,1), gridMap.xlims, lambda, Vwind(1)); % Currently this always assumes wind is from west (180 degree)
    L = max(1,K-floor(tk_max));
       
    % =============== Append uav position ===========================
%     if K == 1
%         continue
%     else
%         % Code
%     end
    % ===============================================================
    
    % =============== Wind data =====================================
    Wind(K,:) = Vwind + [normrnd(0,0.5),normrnd(0,0.5)];
            windsum = sum(Wind(L:K,:),1);
            Vx = windsum(1)*dt; Vy = windsum(2)*dt;
    % ===============================================================
    
    if Time <= plume_start
        alpha(:,K+1) = alpha(:,K);
        continue
    end
    pos = pos + 1;
    
    
    detection = plume.conc(P_uav(pos,1),P_uav(pos,2)) > plume.threshold; % Have to change this to binary method
    beta(:,index) = zeros(M,1);
    for tl = L:K
        for i = 1:M % Calculation of Sij
            
            deltax = P_uav(pos,1) - xcell(i) - Vx;
            deltay = P_uav(pos,2) - ycell(i) - Vy;
            deviation_x = sqrt(T(K+1)-T(tl))*sx;
            deviation_y = sqrt(T(K+1)-T(tl))*sy;
            
            %if abs(deltax) < 10*deviation_x && ...
            %       abs(deltay) < 10*deviation_y
            % Equation (19)
            
            Sij(i) = Lx*Ly * exp((-deltax^2)/(2*deviation_x^2)) * ...
                exp((-deltay^2)/(2*deviation_y^2)) / ...
                (2*pi*deviation_x*deviation_y);
            
            %else
            %   Sij(i) = 0;
            %end
        end % // end cell traversal loop -> Calculation of Sij
        if all(Sij == 0)
            keyboard % If it reaches here, probably would have to adjust sx, sy
        end
        
        Sij = Sij./sum(Sij); % Equation 20
        
        if detection
            beta(:,index) = beta(:,index) + Sij;
        else
            gamma(:,index) = gamma(:,index) .* (1 - mu*Sij);
            %             Sij(Sij < 1e-10) = 0;
            %             gamma = gamma .* (1 - mu*Sij);
        end
        
    end % // end loop L:K
    
    if detection
        beta(:,index) = beta(:,index)./length(L:K);   %  Equation (11)
        alpha(:,K+1) = M * beta * alpha(:,K) ;
    else
        alpha(:,K+1) = (M/sum(sum(gamma))) * gamma * alpha(:,K);
%         alpha(:,K+1) = (M/sum(sum(gamma))) * alpha(:,K) .* gamma ;
    end
    pause(0.5);
    alpha(:,K+1) = alpha(:,K+1)./sum(alpha(:,K+1)); % <-------- Not part of the paper
    surf(X,Y,reshape(alpha(:,K+1),[m,n]));
    colorbar
    view(0,90);
    
end

%% 
function ind = find_index(xpos,ypos, ylimit)
        
    global m Lx Ly
    if mod(xpos,Lx) ~= 0
        xbar = floor(xpos/Lx);
        if xbar >= Lx/2
            xpos = xbar*Lx + Lx;
        else
            xpos = xbar*Lx;
        end
    end
    if mod(ypos,Ly) ~= 0
        ybar = floor(ypos/Ly);
        if ybar >= Ly/2
            ypos = ybar*Ly + Ly;
        else
            ypos = ybar*Ly;
        end
    end

    a = round(xpos/Lx + 1);
    b = round((ypos-ylimit)/Ly + 1);
    ind = a + (b-1)*m;
end

function tk_max = find_tkmax(pos, xlims, lambda, windspeed)
    tk_max = (pos(1)-xlims(1)-lambda)/windspeed;

end

function P_uav = waypoints(c, Set)
    % 1  ---> out of plume (static)
    % 2  ---> in plume (static)
    % 3  ---> out of plume (dynamic)
    % 4  ---> in plume (dynamic)
    % 5  ---> in plume, then out plume
    % 6  ---> out plume, then in plume
    switch Set
        case 1
            switch c
                case 1
                    P_uav = [
                        500 -200
                        500 -200
                        500 -200
                        500 -200
                        500 -200
                        500 -200
                        500 -200
                        500 -200
                        500 -200
                        500 -200
                        ];
                case 2
                    P_uav = [
                        1320 80
                        1320 80
                        1320 80
                        1320 80
                        1320 80
                        1320 80
                        1320 80
                        1320 80
                        1320 80
                        1320 80
                        1320 80
                        1320 80
                        1320 80
                        1320 80
                        1320 80
                        1320 80
                        1320 80
                        1320 80
                        1320 80
                        1320 80
                        1320 80
                        1320 80
                        1320 80
                        1320 80
                        1320 80
                        1320 80
                        1320 80
                        1320 80
                        1320 80
                        1320 80
                        1320 80
                        ];
                case 3
                    P_uav = [
                        1320 200
                        1320 190
                        1320 170
                        1320 150
                        1320 120
                        1320 100
                        1320 200
                        1320 190
                        1320 170
                        1320 150
                        1320 120
                        1320 100
                        1320 200
                        1320 190
                        1320 170
                        1320 150
                        1320 120
                        1320 100
                        ];
                case 4
                    P_uav = [
                        1320 0
                        1300 0
                        1280 0
                        1260 0
                        1200 0
                        1190 0
                        1160 0
                        1140 0
                        ];
                case 5
                    P_uav = [
                        1320 200
                        1320 150
                        1320 100
                        1320 80
                        1320 60
                        1320 40
                        ];
                    P_uav = flip(P_uav);
                    
                case 6
                    P_uav = [
                        1320 200
                        1320 150
                        1320 100
                        1320 80
                        1320 60
                        1320 40
                        ];
                otherwise
                    disp('No Values');
                    P_uav =[];
            end
            
        case 2
            switch c
                case 1
                    P_uav = [
                        400 140
                        400 140
                        400 140
                        400 140
                        400 140
                        400 140
                        400 140
                        400 140
                        ];
                case 2
                    P_uav = [
                        400 50
                        400 50
                        400 50
                        400 50
                        400 50
                        400 50
                        400 50
                        400 50
                        400 50
                        ]; 
                case 4
                    P_uav = [
                        400 10
                        390 10
                        380 10
                        370 10
                        360 10
                        350 10
                        340 10
                        330 10
                        320 10
                        310 10
                        300 10
                        ];
                case 5
                    P_uav = [
                        460 80
                        450 80
                        440 80
                        430 80
                        420 80
                        410 80
                        400 80
                        390 80
                        380 80
                        380 60
                        380 40
                        380 20
                        380 0
                        ];
                    P_uav = flip(P_uav);
                    
                case 6
                    P_uav = [
                        460 80
                        450 80
                        440 80
                        430 80
                        420 80
                        410 80
                        400 80
                        390 80
                        380 80
                        380 60
                        380 40
                        380 20
                        380 0
                        ];
                otherwise
                    disp('No Values');
                    P_uav =[];
            end
        otherwise
            disp('No Values');
            P_uav =[];     
    end

end

    