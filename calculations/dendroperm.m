function [perm,T] = dendroperm(Z)
%DENDROPERM Generate dendrogram premuation.
%   DENDROPERM(Z) generates the permutations for the dendrogram plot of the
%   hierarchical binary cluster tree represented by Z without actually
%   plotting the dendrogram.  Z is an (M-1)-by-3 matrix, generated by the
%   LINKAGE function, where M is the number of objects in the original
%   dataset.
%
%   perm = DENDROGRAM(...) returns the permuations.
%
%   [perm,T] = DENDROGRAM(...) generates a dendrogram and returns T, a vector
%   of size M that contains the leaf node number for each object in the
%   original dataset.  T is useful when P is less than the total number of
%   objects, so some leaf nodes in the display correspond to multiple
%   objects.  For example, to find out which objects are contained in leaf
%   node k of the dendrogram, use find(T==k). When there are fewer than P
%   objects in the original data, all objects are displayed in the
%   dendrogram.  In this case, T is the identity map, i.e., T = (1:M)',
%   where each node contains only a single object.
%
%
%   Example:
%      X = rand(100,2);
%      Y = pdist(X,'cityblock');
%      Z = linkage(Y,'average');
%      [H, T] = dendrogram(Z);
%   
%      rng('default')
%      X=rand(10,3);
%      Z=linkage(X,'ave');
%      % Draw the dendrogram using the default setting
%      perm = dendrogram(Z);
%
%   See also LINKAGE, PDIST, CLUSTER, CLUSTERDATA, COPHENET, INCONSISTENT,
%   SILHOUETTE.

%   Copyright 1993-2011 The MathWorks, Inc.


numLeaves = size(Z,1)+1; %the number of observations

p = 0;
    
color = false;
orientation = 't'; %default top
obslabels = [];
threshold = 0.7 * max(Z(:,3));
leafOrder = [];
horz = false;


% For each node currently labeled numLeaves+k, replace its index by
% min(i,j) where i and j are the nodes under node numLeaves+k.
Z = transz(Z);
T = (1:numLeaves)'; 



A = zeros(4,numLeaves-1);
B = A;
X = 1:numLeaves; %the initial points for observation 1:n
Y = zeros(numLeaves,1);

    r = Y;
    % arrange Z into W so that there will be no crossing in the dendrogram.
    W = zeros(size(Z));
    W(1,:) = Z(1,:);
    
    nsw = zeros(numLeaves,1); rsw = nsw;
    nsw(Z(1,1:2)) = 1; rsw(1) = 1;
    k = 2; s = 2;
    
    while (k < numLeaves)
        i = s;
        while rsw(i) || ~any(nsw(Z(i,1:2)))
            if rsw(i) && i == s
                s = s+1;
            end
            i = i+1;
        end
        
        W(k,:) = Z(i,:);
        nsw(Z(i,1:2)) = 1;
        rsw(i) = 1;
        if s == i
            s = s+1;
        end
        k = k+1;
    end
    
    % initialize X based on W
    g = 1;
    for k = 1:numLeaves-1
        i = W(k,1); %the left node in W(k,:)
        if ~r(i),
            X(i) = g;
            g = g+1;
            r(i) = 1;
        end
        i = W(k,2); %the right node in W(k,:)
        if ~r(i),
            X(i) = g;
            g = g+1;
            r(i) = 1;
        end
    end
    perm(X) = 1:numLeaves;




% ---------------------------------------
function Z = transz(Z)
%TRANSZ Translate output of LINKAGE into another format.
%   This is a helper function used by DENDROGRAM and COPHENET.

%   In LINKAGE, when a new cluster is formed from cluster i & j, it is
%   easier for the latter computation to name the newly formed cluster
%   min(i,j). However, this definition makes it hard to understand
%   the linkage information. We choose to give the newly formed
%   cluster a cluster index M+k, where M is the number of original
%   observation, and k means that this new cluster is the kth cluster
%   to be formed. This helper function converts the M+k indexing into
%   min(i,j) indexing.

numLeaves = size(Z,1)+1;

for i = 1:(numLeaves-1)
    if Z(i,1) > numLeaves
        Z(i,1) = traceback(Z,Z(i,1));
    end
    if Z(i,2) > numLeaves
        Z(i,2) = traceback(Z,Z(i,2));
    end
    if Z(i,1) > Z(i,2)
        Z(i,1:2) = Z(i,[2 1]);
    end
end


function a = traceback(Z,b)

numLeaves = size(Z,1)+1;

if Z(b-numLeaves,1) > numLeaves
    a = traceback(Z,Z(b-numLeaves,1));
else
    a = Z(b-numLeaves,1);
end
if Z(b-numLeaves,2) > numLeaves
    c = traceback(Z,Z(b-numLeaves,2));
else
    c = Z(b-numLeaves,2);
end

a = min(a,c);



