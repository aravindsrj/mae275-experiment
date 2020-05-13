plume4occupancy;

% Sample waypoints
    % 1  ---> out of plume (static)
    % 2  ---> in plume (static)
    % 3  ---> out of plume (dynamic)
    % 4  ---> in plume (dynamic)
    % 5  ---> in plume, then out plume
    % 6  ---> out plume, then in plume
P_uav = waypoints(2);

Duration   = length(P_uav); %800;
dt         = 1;             % Time step
N          = length(0:dt:Duration);
Vwind      = [props.U+0.01,0.01];  % Constant wind vector
Wind       = zeros(N,2);  % Wind vector


global Lx Ly m
Lx         = 40;    % Length of probability cells
Ly         = 20;    % Length of probability cells
m          = 51;    % number of cells in x-direction 
n          = 51;    % number of cells in y-direction     
M          = m*n;   % Total cells

x          = gridMap.xlims(1):Lx:gridMap.xlims(2);
y          = gridMap.ylims(1):Ly:gridMap.ylims(2);
[Y,X]      = meshgrid(y,x);
xcell      = X(:); % 
ycell      = Y(:); %

sx         = 20;   % Standard deviations may have to be modified
sy         = 10;   % Standard deviations may have to be modified  
mu         = 0.9;  % Sensor accuracy

L          = 1;     % lower time bound
K          = 0;     % Upper time bound; or current time step
T          = zeros(1,N);

tk_max     = (gridMap.xlims(2)-gridMap.xlims(1)-1500)/Vwind(1); % Max horizon time
plume_start = 400; % Time between the start of the plume and the start of the mission

alpha      = (1/M)*ones(M,N+plume_start);   % Initializing probability map
Sij        = zeros(M,1);
beta       = zeros(M,M);     % Detection map
gamma      = ones(M,M);     % Non-detection map

figure(2)
view(0,90);
surf(X,Y,reshape(alpha(:,K+1),[51,51]))

for Time = dt:dt:Duration+plume_start
    
    K = K + 1;
    T(K+1) = Time;
    if Time > tk_max
        L = L + 1;
    end
       
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

    detection = plume.conc(P_uav(K-plume_start,1),P_uav(K-plume_start,2)) > plume.threshold; % Have to change this to binary method
    view(0,90);
    for tl = L:K
        for i = 1:M

            deltax = P_uav(K-plume_start,1) - xcell(i) - Vx;
            deltay = P_uav(K-plume_start,2) - ycell(i) - Vy;
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
        end
       if all(Sij == 0)
           keyboard
       end
       Sij = Sij./sum(Sij); % Equation 20
        index = find_index(P_uav(K-plume_start,1),P_uav(K-plume_start,2));
        if detection
            beta(:,index) = beta(:,index) + Sij;
        else
            gamma(:,index) = gamma(:,index) .* (1 - mu*Sij);
%             Sij(Sij < 1e-10) = 0;
%             gamma = gamma .* (1 - mu*Sij);
        end               
    end
    
    if detection
        beta(:,index) = beta(:,index)./length(L:K);   %  Equation (11)
        alpha(:,K+1) = M * beta * alpha(:,K) ;
    else
        alpha(:,K+1) = (M/sum(sum(gamma))) * gamma * alpha(:,K)  ;
%         alpha(:,K+1) = (M/sum(sum(gamma))) * alpha(:,K) .* gamma ;
    end
    pause(0.5);
    alpha(:,K+1) = alpha(:,K+1)./sum(alpha(:,K+1)); % <-------- Not part of the paper
    surf(X,Y,reshape(alpha(:,K+1),[51,51]));
    colorbar
    view(0,90);
    
end

function ind = find_index(xpos,ypos)
    global m Lx Ly
    a = round(xpos/Lx + 1);
    b = round((ypos+500)/Ly + 1);
    ind = a + (b-1)*m;
end

function P_uav = waypoints(c)
    % 1  ---> out of plume (static)
    % 2  ---> in plume (static)
    % 3  ---> out of plume (dynamic)
    % 4  ---> in plume (dynamic)
    % 5  ---> in plume, then out plume
    % 6  ---> out plume, then in plume
    
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

end

    