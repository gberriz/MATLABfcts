function drugs_struct = DesignCombo_PlateLayout(PrimaryDrugs, PrimaryDoses, PrimarySingleDoses, ...
    SecondaryDrugs, SecondaryDoses, SecondarySingleDoses, randomize, edge_ctrl)
% function drugs_struct = plate_design_combo(PrimaryDrugs, PrimaryDoses, PrimarySingleDoses, ...
%    SecondaryDrugs, SecondaryDoses, SecondarySingleDoses, randomize, edge_ctrl)
%

if ~exist('edge_ctrl','var')
    edge_ctrl = true;
end

% well volume
well_volume = 6e-5;

% conc in mM

drugs_struct = struct('name', [PrimaryDrugs;SecondaryDrugs], 'conc', 10, ...
    'layout',zeros(16,24));


total_cnt = length(PrimaryDrugs)*length(PrimarySingleDoses) + ...
    length(SecondaryDrugs)*length(SecondarySingleDoses) + ...
    numel(PrimaryDoses)*numel(SecondaryDoses);


for iPD = 1:length(PrimaryDrugs)
    drugs_struct(iPD).Primary = true;
    drugs_struct(iPD).Doses = PrimaryDoses(iPD,:);
    drugs_struct(iPD).SingleDoses = PrimarySingleDoses(iPD,:);
    assert(all(ismember(PrimaryDoses(iPD,:), PrimarySingleDoses(iPD,:))),...
        'Dose in the combo not found in the titration for %s',PrimaryDrugs{iPD})
end
for iSD = 1:length(SecondaryDrugs)
    drugs_struct(length(PrimaryDrugs)+iSD).Primary = false;
    drugs_struct(length(PrimaryDrugs)+iSD).SingleDoses = SecondarySingleDoses(iSD,:);
    drugs_struct(length(PrimaryDrugs)+iSD).Doses = SecondaryDoses(iSD,:);
    if size(SecondarySingleDoses,2)>0
        assert(all(ismember(SecondaryDoses(iSD,:), SecondarySingleDoses(iSD,:))),...
            'Dose in the combo not found in the titration for %s',SecondaryDrugs{iSD})
    end
end

        
fprintf('Number of test wells: %i\n',total_cnt);
assert(total_cnt<384)

ctrl_cnt = 384 - total_cnt;

if edge_ctrl
    assert(ctrl_cnt>=12, 'For controls on the edge, need at least 12 controls')
    ctrlpos = zeros(16,24);
        
    if ctrl_cnt<14 % only 6 controls on the edge
        disp('Only 6 control on the edge, would be better with at least 14 controls')
        ctrlpos([1 4 12 16],[7 18])=1;
        ctrlpos(8,[1 7 18 24])=1;
    else % 2 controls on each edge, 6 regularly spread in the middle.
        ctrlpos([1 4 12 16],[7 18])=1;
        ctrlpos([5 11],[1 12 24])=1;
    end
    
    if ctrl_cnt>=20 % put the corners as control (should be then discarded)
        ctrlpos([1 end], [1 end]) = 1;
    end
    
    % complete the number of control with randomized positions
    temp = [ 384*ones(1,24); ...
        384*ones(14,1) reshape(randperm(14*22),14,22) 384*ones(14,1);
        384*ones(1,24)];
    temp(ctrlpos==1) = 384;
    
    cutoff = sort(temp(:));
    cutoff = cutoff(ctrl_cnt-sum(ctrlpos(:)));
    ctrlpos(temp<=cutoff) = 1;
    ctrlidx = find(ctrlpos);
    
    assert(ctrl_cnt==length(ctrlidx));
else
    ctrlidx = find(reshape(randperm(16*24),16,24)<=(ctrl_cnt));
end

if ~exist('randomize','var') || randomize>0
    trtidx = randperm(384);
    trtidx = trtidx(~ismember(trtidx,ctrlidx));
else
    trtidx = setdiff(1:384,ctrlidx);
end

letter_label = 'A':'P';

fprintf('Control wells (%i):\n', ctrl_cnt);
for i=1:length(ctrlidx)
    fprintf('\t%s%i', letter_label(mod(ctrlidx(i)-1,16)+1), ceil(ctrlidx(i)/16));
end
disp(' ')



cnt = 1;

for iPD = 1:length(PrimaryDrugs)
    
    % titration
    for iPSDo = 1:length(drugs_struct(iPD).SingleDoses)
        drugs_struct(iPD).layout(trtidx(cnt)) = drugs_struct(iPD).SingleDoses(iPSDo);
        cnt = cnt+1;
    end
    
    for iSD = 1:length(SecondaryDrugs)
        
        % combos
        for iPDo = 1:length(drugs_struct(iPD).Doses)
            for iSDo = 1:length(drugs_struct(length(PrimaryDrugs)+iSD).Doses);
                
                drugs_struct(iPD).layout(trtidx(cnt)) = drugs_struct(iPD).Doses(iPDo);
                
                drugs_struct(length(PrimaryDrugs)+iSD).layout(trtidx(cnt)) = ...
                    drugs_struct(length(PrimaryDrugs)+iSD).Doses(iSDo);
                
                cnt = cnt+1;
            end
        end
        
    end
end
        

for iSD = 1:length(SecondaryDrugs)    
    % titration
    for iSSDo = 1:length(drugs_struct(length(PrimaryDrugs)+iSD).SingleDoses)
        drugs_struct(length(PrimaryDrugs)+iSD).layout(trtidx(cnt)) = ...
            drugs_struct(length(PrimaryDrugs)+iSD).SingleDoses(iSSDo);
        cnt = cnt+1;
    end
    
end



assert( cnt==length(trtidx)+1)

for iPD = 1:length(PrimaryDrugs)
    assert(sum(drugs_struct(iPD).layout(:)~=0) == (length(drugs_struct(iPD).Doses)*...
        sum((~[drugs_struct.Primary]).*cellfun(@length, {drugs_struct.Doses}))+...
        length(drugs_struct(iPD).SingleDoses)))
    % volume in nl from stock (conc in mM)
    drugs_struct(iPD).volume = 1e9*well_volume*sum(drugs_struct(iPD).layout(:))/...
        (1e3*drugs_struct(iPD).conc);
end

for iSD = 1:length(SecondaryDrugs)
    assert(sum(drugs_struct(length(PrimaryDrugs)+iSD).layout(:)~=0)==...
        (length(drugs_struct(length(PrimaryDrugs)+iSD).Doses)*...
        sum([drugs_struct.Primary].*cellfun(@length, {drugs_struct.Doses}))+...
        length(drugs_struct(length(PrimaryDrugs)+iSD).SingleDoses)))
    drugs_struct(length(PrimaryDrugs)+iSD).volume = ...
        1e9*well_volume*sum(drugs_struct(length(PrimaryDrugs)+iSD).layout(:))/...
        (1e3*drugs_struct(length(PrimaryDrugs)+iSD).conc);
end