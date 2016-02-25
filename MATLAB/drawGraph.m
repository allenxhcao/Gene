function [bg, bgfig] = drawGraph(A,ids)

% transpose the interaction matrix A
% because in original A, Aij means from j to i, but we need from i to j
% here

n = size(A,1); % number of nodes
selfInteraction = diag(A);  % self interaction
otherInteraction = A - eye(n).*A;   % except self interaction

edgeIdx = find(otherInteraction~=0);
edgeSign = sign(otherInteraction(edgeIdx));
edgeMag = sqrt(abs(otherInteraction(edgeIdx)));

nodeSign = sign(selfInteraction);
% nodeMag = sqrt(abs(selfInteraction));
nodeMag = abs(selfInteraction);

if nargin < 2
    bg = biograph(otherInteraction');
else
    bg = biograph(otherInteraction',ids);
end

for k = 1:n
    
    if nodeSign(k) == 0
        bg.nodes(k).LineColor = [0,0,0];
    elseif nodeSign(k) > 0
        bg.nodes(k).LineWidth = nodeMag(k);
        bg.nodes(k).LineColor = [0,1,0];
    else
        bg.nodes(k).LineWidth = nodeMag(k);
        bg.nodes(k).LineColor = [1,0,0];
    end
    bg.nodes(k).Shape = 'ellipse';
    bg.nodes(k).Size = [100 50];
    
end

for k = 1:length(edgeIdx)
    bg.edges(k).LineWidth = edgeMag(k);
    if edgeSign(k) >0
        bg.edges(k).LineColor = [0,1,0];
    else
        bg.edges(k).LineColor = [1,0,0];
    end
end

view(bg)
set(0,'ShowHiddenHandles','on');
bgfig = gcf;
close all
    
    
    
