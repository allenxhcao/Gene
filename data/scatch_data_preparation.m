clear
close all

load('gse54514_top10');

%% re-format the data and noralized the data
[nTp,nFea,nSub] = size(data);
ts = (0:nTp-1)';

ts_temp = repmat(ts,1,1,nSub);

stackData = cat(2,ts_temp,data);
nFea = nFea + 1;

temp = reshape(permute(stackData,[1,3,2]),[],nFea);

temp_normalized = nan(size(temp));
temp_normalized(:,1) = temp(:,1);



for k = 2:size(temp,2)
    mea = temp(:,k);
    nanInd = isnan(mea);
    mea(nanInd) = [];
    [normalizedMea, setting] = normalize_lr(mea,ones(length(mea),1));
    temp_normalized(:,k) = normalize_lr(temp(:,k),ones(length(temp(:,k)),1),setting);
end

normalizedStackData = reshape(temp_normalized,nTp,nSub,[]);
normalizedStackData = permute(normalizedStackData,[1,3,2]);

%% choose subjects that have 3 or more data points
temp = squeeze(stackData(:,2,:));
chosen_ind = sum(~isnan(temp))>=3;
save('gse54514_top10_processed','stackData','normalizedStackData','chosen_ind','label');