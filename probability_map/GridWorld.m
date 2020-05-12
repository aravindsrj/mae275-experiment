classdef GridWorld
    properties
        xlims
        ylims
        totalLen
        xStep
        yStep
        Xq
        Yq
        xMesh
        yMesh
        zMesh
    end
    methods
        function self = GridWorld(varargin)
            Defaults = {101, [0 2000], [-500 500], 10};
            optargs(1:length(Defaults)) = Defaults;
            if nargin > 0
                optargs(1:nargin) = varargin;
            end
            [self.totalLen, self.xlims, self.ylims, z] = optargs{:};
            x = linspace(self.xlims(1),self.xlims(2),self.totalLen);
            y = linspace(self.ylims(1),self.ylims(2),self.totalLen);
            self.xStep = (self.xlims(2) - self.xlims(1))/(self.totalLen-1);
            self.yStep = (self.ylims(2) - self.ylims(1))/(self.totalLen-1);
            [self.Xq,self.Yq] = meshgrid(x,y);
            [self.zMesh, self.yMesh, self.xMesh] = meshgrid(z,y,x);
        end
    end
end