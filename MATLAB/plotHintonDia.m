function out = plotHintonDia(A,maxV)

n = size(A,1);
maxV = ones(1,n)*maxV; 

h = plotwb([],[A;maxV],[]);

h = get(h,'Children');


ylim([0.51,n+0.5])

set(h(end-1),'XLabel',[])
set(h(end-1),'YLabel',[])
set(h(end-1),'XTickLabel',[])
set(h(end-1),'YTickLabel',[])

set(h(end),'XLabel',[])
set(h(end),'YLabel',[])
set(h(end),'XTickLabel',[])
set(h(end),'YTickLabel',[])
out = gcf;

