

load sepsisGeneExp.mat

load dataForXi.mat


SubjN=length(ind);
TimeP=5;
FeatD=size(data,1);

dataMatrix=NaN(FeatD,TimeP,SubjN)

for i = 1:SubjN
    tempS=ind(i);
    temp=info(info(:,4)==tempS,:);
    for Tpoints=1:size(temp,1)
        time=temp(Tpoints,3);
        index=temp(Tpoints,1);
        dataMatrix(:,time,i)=data(:,index);
    end
end

