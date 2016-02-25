% Input:
%   tData: a matrix with time series data. The first row is the time stamps
%       each row after represents a variable
%   wSize: window size
%   incrPerc: increment percentage when the window shifts
%   startAndEnd: (optional) starting and ending time stamps
%   lambda: (optional) regularization coefficient for smoothness
%   dt: (optional) time interval

% output:
%   cA: a cell array of interaction matrix A

function cA = LDSSolver (tData,wSize,incrPerc,startAndEnd,lambda,dt)

if nargin<4
    startAndEnd = [tData(1,1) tData(1,end)];
    lambda = 1e-5;
    dt = min(diff(tData(1,:)));
elseif nargin < 5
    lambda = 1e-5;
    dt = min(diff(tData(1,:)));
elseif nargin < 6
    dt = min(diff(tData(1,:)));
end


tStart = startAndEnd(1);
tEnd = startAndEnd(2);

tData = tData(:,tData(1,:)>=tStart&tData(1,:)<=tEnd);

ts = tData(1,:);
param.lambda = lambda;
param.dt = dt;

firstTs = tData(1,1);
lastTs = tData(1,end);

incr = floor(wSize*incrPerc);

tsStart = firstTs:incr*dt:lastTs-wSize*dt;
tsEnd = tsStart + wSize*dt - dt;

nw = length(tsStart);
cA = cell(nw,1);
p = size(tData,1)-1;
Aest = -.1*eye(p);

n = ceil(sqrt(nw));

for iw = 1:nw
    
    t1 = tic;
    inputID = (ts - tsStart(iw)).*(ts - tsEnd(iw))<=0;
    inputTData = tData(:,inputID);
    [Aest, fval] = optDepMtx(inputTData,Aest,param);
    cA{iw} = Aest;
    t2 = toc(t1);
    disp(['Window #' num2str(iw) ' takes ' num2str(t2) ' fVal = ' num2str(fval)])
    
    estSeries = dlds(Aest,inputTData(2:end,1),tsEnd(iw)-tsStart(iw),dt);
    estSeries(1,:) = estSeries(1,:) + tsStart(iw);
    
    subplot(n,n,iw)
    plot(inputTData(2:end,:)','o');
    hold on
    plot(estSeries(2:end,:)')
    hold off
    ylim([0 1])
    title(num2str(iw))

    
end