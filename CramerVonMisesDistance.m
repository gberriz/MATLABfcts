function [dist CDF1 CDF2 BinCenters] = CramerVonMisesDistance(x1,x2)

BinCenters    =  [min([x1;x2])-1 ; sort([x1;x2]) ; max([x1;x2])+1];

binCounts1  =  hist (x1 , BinCenters);
binCounts2  =  hist (x2 , BinCenters);

CDF1  =  cumsum(binCounts1)'./sum(binCounts1);
CDF2  =  cumsum(binCounts2)'./sum(binCounts2);

dist = sum( ((CDF1(1:end-1)-CDF2(1:end-1)).^2).*(BinCenters(2:end)-BinCenters(1:end-1)));
