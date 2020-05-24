%clear; clc;
close all
gridMap = GridWorld(10, 8, [0 500], [-200 200]);
% gridMap = GridWorld();
props = struct('rate',2,'U',5,'H',0,'source',[0,0,0]);

% Stability Class: Input must be a letter from A to F
stab_class = 'A';
stability = stabilityClass(stab_class);
theta = 0;
plume = GaussPlume(props,stability, gridMap, theta);


% plume.conc(100,10); % Calculates plume concentration at particular coordinates

% plume.concentration(plume.Xqr<gridMap.xlims(1) | plume.Xqr>gridMap.xlims(2)) = nan;
% plume.concentration(plume.Yqr<gridMap.ylims(1) | plume.Yqr>gridMap.ylims(2)) = nan;

% theta = [10:-1:-10 -10:1:10];
% theta = [theta theta];
[Xqr,Yqr] = plume.rotation(gridMap.Xq, gridMap.Yq,0);
c = zeros(size(Xqr,1),size(Xqr,2),1);
c(:,:,1) = interp2(gridMap.Xq,gridMap.Yq,plume.concentration, Xqr,Yqr);
c(isnan(c(:,:,1))) = 0;
% s = surf(gridMap.Xq,gridMap.Yq,c(:,:,1));
% colormap jet
% colorbar
% view(0,90)
inplume = double(c > plume.threshold);
figure(1)
set(gcf, 'Position', [465.8000 41.8000 628.8000 740.8000]);
subplot(2,1,1)
t = surf(gridMap.Xq,gridMap.Yq,inplume(:,:,1));
colorbar
view(0,90)
hold on