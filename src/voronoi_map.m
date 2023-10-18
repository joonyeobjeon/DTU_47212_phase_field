

%% Poisson point field 
% create map with nuclei
function [AA] = voronoi_map(nmax, NP)
randseed = 05062019;    % historic reason

rng(randseed);

AA = zeros(nmax);       % creates field


xP = floor(rand(NP,2)*nmax)+1; % list with coordinates of NP points 
for iP=1:NP
    AA(xP(iP,2),xP(iP,1))=iP;  % mark point with index
end

Nnuclei = nnz(AA);
Missednuclei = NP - Nnuclei;
% figure; 
% im = imagesc(AA);
% axis equal;
% axis off;

%% Voronoi tessellation (based on Poisson point field)
% (simplest case of tessellation)
% Calculate Eucledian distancess from all nuclei
% and give point same index as nearest nucleus

for j=1:nmax                    % rows
    for i=1:nmax                % columns
        dy = abs(xP(:,2)-j);
        dy = min(dy,nmax-dy);   % periodic boundary condition
        dx = abs(xP(:,1)-i);
        dx = min(dx,nmax-dx);   % periodic boundary condition
        dist2 = dx.^2+dy.^2;    % Eucledian distance      
        [mindist2, nearestpoint] = min(dist2);
        AA(j,i)=nearestpoint;
    end
end

AAstart = AA;
return 
% 
% figure; 
% imagesc(AAstart);
% axis equal;
% axis off;
