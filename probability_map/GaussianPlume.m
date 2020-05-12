classdef GaussianPlume
    properties
        props
        stability
        concentration
        threshold
        epsilon
    end
    methods
        function self = GaussianPlume(props, stability, epsilon)
            self.props = props;
            self.stability = stability;
            self.epsilon = epsilon;
        end
        
        function [Conc,self] = calculateConcentration(self)
            Conc = zeros(size(self.props.xMesh));
            x = self.props.source(1);
            y = self.props.source(2);
            
            a = self.props.rate/(2 * pi * self.props.U .*self.stability.sy(self.props.xMesh - x)...
                .* self.stability.sy(self.props.xMesh - x));
            b = exp(-((self.props.yMesh-y).^2)./(2*self.stability.sy(self.props.xMesh - x).^2));
            c = exp(-(self.props.zMesh-self.props.H).^2./(2*self.stability.sz(self.props.xMesh - x).^2))...
                +exp(-(self.props.zMesh+self.props.H).^2./(2*self.stability.sz(self.props.xMesh - x).^2));
            Conc = Conc + a.*b.*c;
            self.concentration = Conc;
            self.threshold = self.epsilon * max(max(Conc));
        end
        
        function C = conc(self,x,y)
            if abs(x) < 1e-10
                x = 0;
            end
            if abs(y) < 1e-10
                y = 0;
            end
            Y = (y - self.props.ymin)/self.props.ystep + 1;
            X = (x - self.props.xmin)/self.props.xstep + 1;
            if mod(y,self.props.ystep) ~=0 || mod(x,self.props.xstep) ~=0
                if mod(x,self.props.xstep) == 0
                    y1 = [floor(Y); ceil(Y)];
                    C2 = [self.concentration(y1(1),X); self.concentration(y1(2),X)];
                    C = interp1(y1,C2,Y);                    
                elseif mod(y,self.props.ystep) == 0
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
    end
end