classdef stabilityClass
    properties
        Iy
        Iz
        Jy
        Jz
        Ky
        Kz
    end
    methods
        function self = stabilityClass(letter)
            if letter == 'A'
                self.Iy = -1.104;
                self.Jy = 0.9878;
                self.Ky = -0.0076;
                
                self.Iz = 4.679;
                self.Jz = -1.7172;
                self.Kz = 0.2770;
            
            elseif letter == 'B'
                self.Iy = -1.634;
                self.Jy = 1.0350;
                self.Ky = -0.0096;
                
                self.Iz = -1.999;
                self.Jz = 0.8752;
                self.Kz = 0.0136;
                
            elseif letter == 'C'
                self.Iy = -2.054;
                self.Jy = 1.0231;
                self.Ky = -0.0076;
                
                self.Iz = -2.341;
                self.Jz = 0.9477;
                self.Kz = -0.0020;
            
            elseif letter == 'D'
                self.Iy = -2.555;
                self.Jy = 1.0423;
                self.Ky = -0.0087;
                
                self.Iz = -3.186;
                self.Jz = 1.1737;
                self.Kz = -0.0316;
                
            elseif letter == 'E'
                self.Iy = -2.754;
                self.Jy = 1.0106;
                self.Ky = -0.0064;
                
                self.Iz = -3.783;
                self.Jz = 1.3010;
                self.Kz = -0.0450;
                
            elseif letter == 'F'
                self.Iy = -3.143;
                self.Jy = 1.0148;
                self.Ky = -0.0070;
                
                self.Iz = -4.490;
                self.Jz = 1.4024;
                self.Kz = -0.0540;       
            end
        end
            
        function val = sy(self, dist)
            val = exp(self.Iy + self.Jy*log(dist) + self.Ky*log(dist).^2);
        end
        function val = sz(self, dist)
            val = exp(self.Iz + self.Jz*log(dist) + self.Kz*log(dist).^2);
        end
    end
end
        