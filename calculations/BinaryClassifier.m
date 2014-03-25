function [MCC, ConfMx] = BinaryClassifier(x,thres,group)
% [MCC, ConfMx] = BinaryClassifier(x,thres,group)
%

classes = unique(group);
assert(length(classes)==2, 'Only two classes (true/false) are accepted')

ConfMx = zeros(2,2);
for i=1:2
    ConfMx(i,:) = [mean(x(group==classes(i))>thres) ...
        mean(x(group==classes(i))<thres)];
end

if trace(ConfMx)<trace(ConfMx(:,[2 1]))
    ConfMx = ConfMx(:,[2 1]);
end

MCC = MCCscoring(ConfMx);
