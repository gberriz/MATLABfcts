function Designs = TreatmentDesign(DrugNames, HMSLids, SingleDoses, nReps, varargin)
% Designs = TreatmentDesign(DrugNames, HMSLid, SingleDoses, Nreps, varargin)
%
%   Produce randomized treatment desgins based on the specifications
%
%   DrugNames       label used to identify the drug
%   HMSLid          corresponding HMSL id 
%   SingleDoses     cell array of doses vectors, order corresponds to
%                       DrugNames, dispensed as a single agent
%   nReps           Number of repeats (randomized replicates)
%
%
%   Optional arguments (name-value pairs):
%
%   ComboDoses      cell array of doses vectors, order corresponds to
%                       DrugNames
%   DrugPairs       [Nx2] list of drug combinations to test (with
%                       ComboDoses)
%   Seed            initial seed, set to 0 for no randomization (default 1).
%                       For each repeat the seed is 'Seed+i-1'
%   edge_ctrl       default is true
%   stock_conc      in uM, default is 1e4 (= 10mM)
%   well_volume     in uL, default is 60 (for 384-well plates)
%   plate_size      either [N of rows ,  N of columns] or [N of wells],
%                       default is [16 24]
%   max_DMSOpct     maximum percent of DMSO dispensed (default 0.2%)
%   min_volume      step volume dispensed (in uL, default is 2e-5)
%   step_volume     step volume dispensed (in uL, default is 1.3e-5)
%
%   Output: array od Design structure matching the DrugResponseSuite
%
%


%% inputs and parser
p = inputParser;

addParameter(p,'ComboDoses',[], @(x) all(cellfun(@isnumeric,x)));
addParameter(p,'DrugPairs',zeros(0,2), @isnumeric);
addParameter(p,'Seed',1,  @(x) isscalar(x) && isnumeric(x));
addParameter(p,'edge_ctrl',true, @(x) islogical(x) & isscalar(x));
addParameter(p,'stock_conc',1e4, @(x) isvector(x) && isnumeric(x));     % in uM
addParameter(p,'well_volume',60, @(x) isscalar(x) && isnumeric(x));       % in uL
addParameter(p,'plate_size',[16 24], @(x) isvector(x) && isnumeric(x));

% based on the specifications of the D300
addParameter(p,'min_volume', 1.3e-5, @(x) isscalar(x) && isnumeric(x)); % minimum volume of 13pl (in uL)
addParameter(p,'step_volume', 2e-5, @(x) isscalar(x) && isnumeric(x)); % minimum step of 20pl (in uL)
addParameter(p,'max_DMSOpct', .2, @(x) isscalar(x) && isnumeric(x)); % maximum 0.2% percent of DMSO

parse(p, varargin{:})
p = p.Results;
for i = {'ComboDoses' 'DrugPairs' 'Seed' 'edge_ctrl' 'stock_conc' ...
        'well_volume' 'plate_size' 'min_volume' 'step_volume' 'max_DMSOpct'}
    eval([i{:} ' = p.' i{:} ';'])
end
% avoid too high level of DMSO (max is 0.2%)
max_volume = well_volume*max_DMSOpct/100;

%% checks
DrugNames = ToColumn(DrugNames);
assert((length(DrugNames)==length(HMSLids)) || isempty(HMSLids))
assert(length(DrugNames)==length(SingleDoses))

assert(all(DrugPairs(:)<=length(DrugNames)) || isempty(DrugPairs))
assert((length(DrugNames)==length(ComboDoses)) || isempty(ComboDoses))
assert(size(DrugPairs,2)==2)

assert((length(DrugNames)==length(stock_conc)) || (length(stock_conc)==1))
assert(all(stock_conc>20 & stock_conc<5e4), 'Stock concentration should be in uM')
if length(stock_conc)==1
    stock_conc = stock_conc*ones(length(DrugNames),1);
end
stock_conc = num2cell(ToColumn(stock_conc));

if length(plate_size)==1
    plate_size = sqrt(plate_size*[1/1.5 1.5]);
end
assert(round(log2(plate_size(1)))==log2(plate_size(1)))

nWells = prod(plate_size);
if well_volume<1e2 && nWells<100
    warnprintf('Default volume is %.0f set for 384-well plate; check if correct', ...
        well_volume)
end
assert(well_volume>1 && well_volume<1e4, 'Well volume should be given in uL')

if Seed==0 && nReps>1
    warnprintf('Seed is set to 0 (i.e. no randomization) and Nreps>1')
    warnprintf('This means all replicates will have the same layout!')
    warnprintf('Confirm to proceed')
    pause
end


%% initialization
s = RandStream('mt19937ar','Seed',Seed);
RandStream.setGlobalStream(s);


Drugs = struct('DrugName', DrugNames, 'HMSLid', ToColumn(HMSLids),...
    'stock_conc', stock_conc, 'layout', zeros(plate_size));
Designs = struct('plate_dims', repmat({plate_size}, nReps, 1), ...
    'treated_wells', repmat({true(plate_size)}, nReps, 1), ...
    'well_volume', repmat({well_volume}, nReps, 1), ...
    'Drugs', repmat({Drugs}, nReps, 1), 'Seed', num2cell(Seed-1+(1:nReps)'));

SingleDoses = cellfun2(@(x,y,z) round_Doses(x,y,z,'Single',min_volume, ...
    step_volume, max_volume, well_volume), SingleDoses, stock_conc, DrugNames);

if ~isempty(ComboDoses)
    ComboDoses = cellfun2(@(x,y,z) round_Doses(x,y,z,'Combo',min_volume, ...
        step_volume, max_volume, well_volume), ComboDoses, stock_conc, DrugNames);
    if any(cellfun(@(x,y) any(~ismember(x,y)), ComboDoses, SingleDoses))
        fprintf('\n')
        warnprintf('Some doses for the combo are not part of the single doses for %s', ...
            strjoin(DrugNames(cellfun(@(x,y) any(~ismember(x,y)), ComboDoses, SingleDoses)),', '))
    end
    
    nTreatments = sum(cellfun(@length,SingleDoses)) + ...
        sum(cellfun(@(x,y) length(x)*length(y), ComboDoses(DrugPairs(:,1)), ...
        ComboDoses(DrugPairs(:,2))));
else
    nTreatments = sum(cellfun(@length,SingleDoses));
end



%% define the position of the fixed controls
fprintf('Number of test wells: %i\n',nTreatments);
ctrl_cnt = nWells - nTreatments;
assert(ctrl_cnt>5, 'Too many well used (%i out of %i), need at least 6 control wells',...
    nTreatments, nWells)


fixedctrlpos = false(plate_size);
if edge_ctrl
    % specific position of the controls on the edges
    
    assert(ctrl_cnt>=12, 'For controls on the edge, need at least 12 controls')
    
    if ctrl_cnt<14 % only 6 controls on the edge
        disp('Only 6 controls on the edge, would be better with at least 14 controls in total')
        fixedctrlpos(round(1+((plate_size(1)-1)*(0:3))/3), round(1+((plate_size(2)-1)*(1:2))/3)) = true;
        fixedctrlpos(round(1+(plate_size(1)-1)/2),round(1+(plate_size(2)*(0:3))/3)) = true;
    else % 2 controls on each edge, 6 regularly spread in the middle.
        fixedctrlpos(round(1+((plate_size(1)-1)*(0:3))/3), round(1+((plate_size(2)-1)*(1:2))/3)) = true;
        fixedctrlpos(round(1+((plate_size(1)-1)*(1:2))/3), round(1+((plate_size(2)-1)*(0:3))/3)) = true;
    end
    
    if ctrl_cnt>=20 % put the corners as control (should be then discarded)
        fixedctrlpos([1 end], [1 end]) = true;
    elseif ctrl_cnt > (prod(plate_size-1)+6) % remove all edges --> not treated
        fixedctrlpos([1 end], :) = true;
        fixedctrlpos(:, [1 end]) = true;
        for iR=1:nReps
            Designs(iR).treated_wells([1 end], :) = true;
            Designs(iR).treated_wells(:, [1 end]) = true;
        end
    end
end
assert(all(size(fixedctrlpos)==plate_size))
ctrlidx = find(fixedctrlpos);


%% define all possible treatments
allTreatments = zeros(length(DrugNames), nTreatments);
cnt = 0;
for iD = 1:length(DrugNames)
    allTreatments(iD, cnt+(1:length(SingleDoses{iD}))) = SingleDoses{iD};
    cnt = cnt + length(SingleDoses{iD});
end
for iCo = 1:size(DrugPairs,1)
    for iD1 = 1:length(ComboDoses{DrugPairs(iCo,1)})
        allTreatments(DrugPairs(iCo,1), cnt+(1:length(ComboDoses{DrugPairs(iCo,2)}))) = ...
            ComboDoses{DrugPairs(iCo,1)}(iD1);
        allTreatments(DrugPairs(iCo,2), cnt+(1:length(ComboDoses{DrugPairs(iCo,2)}))) = ...
            ComboDoses{DrugPairs(iCo,2)};
        cnt = cnt + length(ComboDoses{DrugPairs(iCo,2)});
    end
end
assert(cnt==nTreatments)
assert(all(any(allTreatments>0)))


%% assign the design for each replicate
for iR = 1:nReps
    
    s = RandStream('mt19937ar','Seed',Seed+iR-1);
    RandStream.setGlobalStream(s);
    if Seed>0
        idx = randperm(nWells); % randomize the order
    else
        idx = 1:nWells;
    end
    idx(ctrlidx) = nWells+1; % put the 'fixed control' as control
    order = sortidx(idx,'ascend'); % find the order for the treatment on the plate
    
    for iD = 1:length(DrugNames)
        Designs(iR).Drugs(iD).layout(order(1:nTreatments)) = allTreatments(iD,:);
    end
    
    allDrugs = reshape([Designs(iR).Drugs.layout], [plate_size length(DrugNames)]);
    nDrugs = sum(allDrugs>0,3);
    assert(all(squeeze(sum(sum((allDrugs>0).*repmat(nDrugs==1,1,1,length(DrugNames)),2),1))==...
        cellfun(@length,SingleDoses)))
    for iCo = 1:size(DrugPairs,1)
        assert(sum(sum( all(allDrugs(:,:,DrugPairs(iCo,:))>0,3)))== ...
            length(ComboDoses{DrugPairs(iCo,1)})*length(ComboDoses{DrugPairs(iCo,2)}))
    end
    
    Treated = any(allDrugs>0,3);
    assert(~any(Treated(ctrlidx)))
    assert(sum(~Treated(:))==ctrl_cnt)
    
end


end


function new_Doses = round_Doses(old_Doses, nominal_conc, DName, tag, ...
    min_volume, step_volume, max_volume, well_volume)
    % adjust doses to match the specifications of the D300

min_dose = nominal_conc *min_volume/well_volume;
step_dose = nominal_conc *step_volume/well_volume;
max_dose = nominal_conc *max_volume/well_volume;
new_Doses = old_Doses;
fprintf('Processing %s (%s doses):\n', DName, tag)
if any(new_Doses<min_dose)
    warnprintf('%i dose(s) below minimal conc (set to %.2g)', ...
        sum(new_Doses<min_dose), min_dose);
    new_Doses = max(new_Doses, min_dose);
end
if any(new_Doses>max_dose)
    warnprintf('%i dose(s) above max conc (set to %.2g)', ...
        sum(new_Doses>max_dose), max_dose);
    new_Doses = min(new_Doses, max_dose);
end
if any( ((mod(new_Doses,step_dose)./new_Doses)>.02) & (new_Doses>1.02*min_dose) )
    idx = ((mod(new_Doses,step_dose)./new_Doses)>.02) & (new_Doses>1.02*min_dose);
    fprintf('\tNote: %i doses round down because more than 2%% difference\n', ...
        sum(idx));
    % case of the lowest doses
    idx = (new_Doses<mean([min_dose, step_dose])) & (new_Doses>1.02*min_dose);
    new_Doses(idx) = min_dose;
    idx = (new_Doses>mean([min_dose, step_dose])) & (new_Doses<step_dose/1.02);
    new_Doses(idx) = step_dose;    
    % other cases
    idx = ((mod(new_Doses,step_dose)./new_Doses)>.02) & (new_Doses>1.02*step_dose);
    new_Doses(idx) = step_dose*round(new_Doses(idx)/step_dose);
end
if any(new_Doses~=old_Doses)
    fprintf('\tChanged doses:')
    for iDoses = ToRow(find(new_Doses~=old_Doses))
        fprintf('\t%-7.2g', old_Doses(iDoses))
    end
    fprintf('\n\t    (%3i)  -->', sum(new_Doses~=old_Doses))
    for iDoses = ToRow(find(new_Doses~=old_Doses))
        fprintf('\t%-7.2g', new_Doses(iDoses))
    end
    fprintf('\n')
end

end

