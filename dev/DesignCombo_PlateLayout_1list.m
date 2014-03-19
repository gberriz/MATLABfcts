function drugs_struct = DesignCombo_PlateLayout_1list(Drugs, Doses, SingleDoses, ...
    randomize, edge_ctrl, nominal_conc)
% function drugs_struct = plate_design_combo(Drugs, Doses, SingleDoses, ...
%     randomize, edge_ctrl, nominal_conc)
%
%   Doses are in uM
%   nominal_conc in mM (stock)
%   default volume is 60 ul
%   minimal dispensed volume is 20 pl, step of 10pl
%   maximal dispensed volume is 120 nl to have a dilution of 500-fold minimum

if ~exist('edge_ctrl','var') || isempty(edge_ctrl)
    edge_ctrl = true;
end

% well volume
well_volume = 6e-5;
min_volumedrop = 1e-11; % minimum step of the drop is 10pl for the 
max_volumedrop = 1.2e-7; % need a dilution of 500-fold minimum

% stock conc in mM
if ~exist('nominal_conc','var')
    nominal_conc = 10;
elseif isvector(nominal_conc)
    nominal_conc = num2cell(ToColumn(nominal_conc));
end



drugs_struct = struct('name', Drugs, 'nominal_conc', nominal_conc, ...
    'layout',zeros(16,24), 'well_volume', well_volume);


total_cnt = length(Drugs)*length(SingleDoses) + ...
    length(Drugs)*(length(Drugs)-1)*(length(Doses)^2)/2;


for iD = 1:length(Drugs)
    
    dose_step = 1e3*drugs_struct(iD).nominal_conc *min_volumedrop/well_volume; 
    dose_max = 1e3*drugs_struct(iD).nominal_conc *max_volumedrop/well_volume; 
    temp_d = Doses(iD,:);
    if any(temp_d<2*dose_step)
        disp(['!! Doses below minimal conc (' num2str(2*dose_step) ') for ' Drugs{iD}])
        disp(['       at doses : ' num2str(temp_d(temp_d<2*dose_step))])
        temp_d = max(temp_d, 2*dose_step);
    end
    if any(temp_d<50*dose_step)
        idx = temp_d<50*dose_step;
        disp([' - Rouding doses for ' Drugs{iD}])
        disp([' - former doses: ' num2str(temp_d(idx))])
        temp_d(idx) = dose_step*round(temp_d(idx)/dose_step);
        disp(['        now as: ' num2str(temp_d(idx))])
    end
    if any(temp_d>dose_max)
        disp(['!! Doses above maximal conc (' num2str(dose_max) ') for ' Drugs{iD}])
        disp(['       at doses : ' num2str(temp_d(temp_d>dose_max))])
        temp_d = min(temp_d, dose_max);
    end
    Doses(iD,:) = temp_d;
    
    temp_d = SingleDoses(iD,:);
    if any(temp_d<2*dose_step)
        disp(['!! Doses below minimal conc (' num2str(2*dose_step) ') for ' Drugs{iD}])
        disp(['       at doses : ' num2str(temp_d(temp_d<2*dose_step))])
        temp_d = max(temp_d, 2*dose_step);
    end
    if any(temp_d<50*dose_step)
        idx = temp_d<50*dose_step;
        disp([' - Rouding doses for ' Drugs{iD}])
        disp([' - former doses: ' num2str(temp_d(idx))])
        temp_d(idx) = dose_step*round(temp_d(idx)/dose_step);
        disp(['        now as: ' num2str(temp_d(idx))])
    end
    if any(temp_d>dose_max)
        disp(['!! Doses above maximal conc (' num2str(dose_max) ') for ' Drugs{iD}])
        disp(['       at doses : ' num2str(temp_d(temp_d>dose_max))])
        temp_d = min(temp_d, dose_max);
    end
    SingleDoses(iD,:) = temp_d;
    
    drugs_struct(iD).Doses = Doses(iD,:);
    drugs_struct(iD).SingleDoses = SingleDoses(iD,:);
    assert(all(ismember(Doses(iD,:), SingleDoses(iD,:))),...
        'Dose in the combo not found in the titration for %s',Drugs{iD})
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
    
    if ctrl_cnt>sum(ctrlpos(:))        
        % complete the number of control with randomized positions
        temp = reshape(randperm(16*24),16,24);
        temp(ctrlpos==1) = 384;
        cutoff = sort(temp(:));
        cutoff = cutoff(ctrl_cnt-sum(ctrlpos(:)));
        ctrlpos(temp<=cutoff) = 1;
    end
    
    ctrlidx = find(ctrlpos);    
    assert(ctrl_cnt==length(ctrlidx));
else
    ctrlidx = find(reshape(randperm(16*24),16,24)<=(ctrl_cnt));
end

if ~exist('randomize','var') || isempty(randomize) || randomize
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

for iD = 1:length(Drugs)
    
    % titration
    for iPSDo = 1:length(drugs_struct(iD).SingleDoses)
        drugs_struct(iD).layout(trtidx(cnt)) = drugs_struct(iD).SingleDoses(iPSDo);
        cnt = cnt+1;
    end
    
    for iD2 = (iD+1):length(Drugs)        
        % combos
        for iDo = 1:length(drugs_struct(iD).Doses)
            for iDo2 = 1:length(drugs_struct(iD2).Doses);                
                drugs_struct(iD).layout(trtidx(cnt)) = drugs_struct(iD).Doses(iDo);                
                drugs_struct(iD2).layout(trtidx(cnt)) = ...
                    drugs_struct(iD2).Doses(iDo2);                
                cnt = cnt+1;
            end
        end        
    end
end
        
assert( cnt==length(trtidx)+1)

for iD = 1:length(Drugs)
    % assert that the number of well with drug is equal to the number of
    % combo and titration
    assert(sum(drugs_struct(iD).layout(:)~=0) == (length(drugs_struct(iD).Doses)*...
        sum(cellfun(@length, {drugs_struct(~ismember(1:end,iD)).Doses}))+...
        length(drugs_struct(iD).SingleDoses)))
    % volume in nl from stock (conc in mM)
    drugs_struct(iD).volume = 1e9*well_volume*sum(drugs_struct(iD).layout(:))/...
        (1e3*drugs_struct(iD).nominal_conc);
end
