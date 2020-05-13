plume4occupancy;

% Sample waypoints
P_uav =[
%     500 -200
%     500 -200
%     500 -200
%     500 -200
%     500 -200
%     500 -200
%     500 -200
%     500 -200
%     500 -200
%     500 -200
%    200 -300
%     500 -200
%     500 -200
%     500 -200
%     500 -200
%     500 -200
%     500 -200
%     500 -200
%     500 -200
%     500 -200
%     500 -200
%     500 -200
%     500 -200
%     500 -200
%     500 -200
%     500 -200
%     500 -200
%     500 -200
%     500 -200
%     500 -200
%     500 -200
%    
%     500 0
%     500 10
%     500 20
%     500 30
%     500 40
%     500 50
%     500 100
%     500  120
%     500 0
%     490 0
%     480 0
%     490 0
%     480 0
%     480 0
%     480 0
%     480 0

      1320 80
      1320 70
      1320 60
      1320 50
      1320 40
      1320 30
      1320 20
      1320 80
      1320 80
      1320 80
      1320 80
      1320 80
      10 80
      10 80
      10 80
      10 80
      10 80
      10 80
      10 80
      10 80
      10 80
      10 80
      10 80
      10 80
      
%       1000 80
%       1000 80
%       1000 80
%       1000 80
%       1000 80
%       1000 80
%       1000 80
%       1000 80
%       1000 80
%       1000 80
%       1000 80
%       1000 80
%       1000 80
%       1000 80
%       1000 80
%       1000 80
%       1000 80
%       1000 80
      
    
    ];    

Duration   = length(P_uav); %800;
dt         = 1;             % Time step
N          = length(0:dt:Duration);
%P_uav      = zeros(N,2);  % position of uav
Vwind      = [10*props.U+0.01,0.01];  % Constant wind vector
Wind       = zeros(N,2);  % Wind vector

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

sx         = 10;   % Standard deviations may have to be modified
sy         = 10;   % Standard deviations may have to be modified  
mu         = 0.9;  % Sensor accuracy

L          = 1;     % lower time bound
K          = 0;     % Upper time bound
T          = zeros(1,N);

tk_max     = (gridMap.xlims(2)-gridMap.xlims(1)-1000)/Vwind(1); % Max horizon time



figure(2)
view(0,90);
surf(X,Y,reshape(alpha(:,K+1),[51,51]))

A = 400; % 0 or 400 -> To use for the time look-back
KK = K + A;

alpha      = (1/M)*ones(M,N+A);   % Initializing probability map
Sij        = zeros(M,1);

for Time = dt:dt:Duration+A
    beta       = zeros(M,1);     % Detection map
    gamma      = ones(M,1);     % Non-detection map

    K = K + 1;
    KK = KK + 1;
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
    
    if Time <= A
        alpha(:,K+1) = alpha(:,K);
        continue
    end
%     
    
    detection = plume.conc(P_uav(K-A,1),P_uav(K-A,2)) > plume.threshold; % Have to change this to binary method
    view(0,90);
    for tl = L:K
        for i = 1:M

            deltax = P_uav(K-A,1) - xcell(i) - Vx;
            deltay = P_uav(K-A,2) - ycell(i) - Vy;
            deviation_x = sqrt(T(K+1)-T(tl))*sx;
            deviation_y = sqrt(T(K+1)-T(tl))*sy;
            
            %if abs(deltax) < 10*deviation_x && ...
             %       abs(deltay) < 10*deviation_y
                % Equation (19)
                Sij(i) = Lx*Ly * exp((-deltax^2)/(2*deviation_x^2)) * ...
                        exp((-deltay^2)/(2*deviation_y^2)) / ...
                        (2*pi*deviation_x*deviation_y);   
%                     if i == 1513 && tl == K
%                         keyboard
%                     end
            %else
             %   Sij(i) = 0;
            %end 
            if isnan(Sij(i))
                keyboard
            end
            
        end

%        Sij = Sij./sum(Sij); % Equation 20
%         if tl == K && ((K-A) == length(P_uav))
%             keyboard
%         end
        if detection
            beta = beta + Sij;
%             beta = beta .* Sij;
        else
            gamma = gamma .* (1 - mu*Sij);
        end               
    end
    
    if detection
        beta = beta./length(L:K);   %  Equation (11)
        alpha(:,K+1) = M * alpha(:,K) .* beta ;
%           alpha(:,K+1) = (M/sum(sum(beta))) * alpha(:,K)' .* beta ;
    else
        alpha(:,K+1) = (M/sum(sum(gamma))) * alpha(:,K) .* gamma ;
    end
    pause(0.5);
    alpha(:,K+1) = alpha(:,K+1)./sum(alpha(:,K+1)); % <-------- Not part of the paper
    surf(X,Y,reshape(alpha(:,K+1),[51,51]));
    colorbar
    view(0,90);
    
end

    