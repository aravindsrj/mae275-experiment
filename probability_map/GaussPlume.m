classdef GaussPlume
    properties
        props
        stability
        concentration
        gridMap
        theta
        Xqr
        Yqr
        threshold
    end
    methods
        function self = GaussPlume(props, stability, gridMap, theta)
            self.props = props;
            self.stability = stability;
            self.gridMap = gridMap;
            self.theta = theta;
            self.concentration = self.calculateConcentration();
%             [self.Xqr,self.Yqr] = self.rot();
            self.threshold = 1e-5;
        end
        
        function Conc = calculateConcentration(self)
            Conc = zeros(size(self.gridMap.xMesh));
            x = self.props.source(1);
            y = self.props.source(2);

            a = self.props.rate/(2 * pi * self.props.U .*self.stability.sy(self.gridMap.xMesh - x)...
                .* self.stability.sy(self.gridMap.xMesh - x));
            b = exp(-((self.gridMap.yMesh-y).^2)./(2*self.stability.sy(self.gridMap.xMesh - x).^2));
            c = exp(-(self.gridMap.zMesh-self.props.H).^2./(2*self.stability.sz(self.gridMap.xMesh - x).^2))...
                +exp(-(self.gridMap.zMesh+self.props.H).^2./(2*self.stability.sz(self.gridMap.xMesh - x).^2));
            Conc = squeeze(Conc + a.*b.*c);
         end

        function C = conc(self,x,y)
            if abs(x) < 1e-10
                x = 0;
            end
            if abs(y) < 1e-10
                y = 0;
            end
            Y = (y - self.gridMap.ylims(1))/self.gridMap.yStep + 1;
            X = (x - self.gridMap.xlims(1))/self.gridMap.xStep + 1;
            if mod(y,self.gridMap.yStep) ~=0 || mod(x,self.gridMap.xStep) ~=0
                if mod(x,self.gridMap.xStep) == 0
                    y1 = [floor(Y); ceil(Y)];
                    C2 = [self.concentration(y1(1),X); self.concentration(y1(2),X)];
                    C = interp1(y1,C2,Y);
                elseif mod(y,self.gridMap.yStep) == 0
                    x1 = [floor(X); ceil(X)];
                    C2 = [self.concentration(Y,x1(1)); self.concentration(Y,x1(2))];
                    C = interp1(x1,C2,X);
                else
                    x1 = [floor(X); ceil(X)]; y1 = [floor(Y); ceil(Y)];
                    [Y1,X1] = meshgrid(y1,x1);
                    C2 = [
                        self.concentration(y1(1),x1(1)) self.concentration(y1(1),x1(2));
                        self.concentration(y1(2),x1(1)) self.concentration(y1(2),x1(2))];
                    C = interp2(Y1,X1,C2,Y,X);
                    
                end
                %error('Input coordinates should be multiples of step size');
                %                 x = round(x); y = round(y);
            else
                C = self.concentration(Y,X);
            end

        end
        function [Xqr, Yqr] = rot(self)
           R=[cosd(self.theta) -sind(self.theta); sind(self.theta) cosd(self.theta)]; 
           XY = [self.gridMap.Xq(:) self.gridMap.Yq(:)];
           rotXY = XY*R';
           Xqr = reshape(rotXY(:,1),size(self.gridMap.Xq,1),[]);
           Yqr = reshape(rotXY(:,2),size(self.gridMap.Yq,1),[]);
        end
    end
    methods(Static)
        function [Xqr,Yqr] = rotation(Xq,Yq,theta)
           theta = -theta;
           R=[cosd(theta) -sind(theta); sind(theta) cosd(theta)]; 
           XY = [Xq(:) Yq(:)];
           rotXY = XY*R';
           Xqr = reshape(rotXY(:,1),size(Xq,1),[]);
           Yqr = reshape(rotXY(:,2),size(Yq,1),[]);            
        end
    end
end

        