function Designs = TreatmentDesign(DrugNames, HMSLids, SingleDoses, nReps, varargin)
% Designs = TreatmentDesign(DrugNames, HMSLids, SingleDoses, nReps, varargin)
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
%   plate_dims      either [N of rows ,  N of columns] or [N of wells],
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
addParameter(p,'ComboLists',[], @(x) all(cellfun(@isnumeric,x(:))));
addParameter(p,'DrugPairs',zeros(0,2), @isnumeric);
addParameter(p,'Seed',1,  @(x) isscalar(x) && isnumeric(x));
addParameter(p,'edge_ctrl',true, @(x) islogical(x) & isscalar(x));
addParameter(p,'stock_conc',1e4, @(x) isvector(x) && isnumeric(x));     % in uM
addParameter(p,'well_volume',60, @(x) isscalar(x) && isnumeric(x));       % in uL
addParameter(p,'plate_dims',[16 24], @(x) isvector(x) && isnumeric(x));

% based on the specifications of the D300
addParameter(p,'min_volume', 1.3e-5, @(x) isscalar(x) && isnumeric(x)); % minimum volume of 13pl (in uL)
addParameter(p,'step_volume', 2e-5, @(x) isscalar(x) && isnumeric(x)); % minimum step of 20pl (in uL)
addParameter(p,'max_DMSOpct', .2, @(x) isscalar(x) && isnumeric(x)); % maximum 0.2% percent of DMSO

parse(p, varargin{:})
p = p.Results;
for i = {'ComboDoses' 'ComboLists' 'DrugPairs' 'Seed' 'edge_ctrl' 'stock_conc' ...
        'well_volume' 'plate_dims' 'min_volume' 'step_volume' 'max_DMSOpct'}
    eval([i{:} ' = p.' i{:} ';'])
end
% avoid too high level of DMSO (max is 0.2%)
max_volume = well_volume*max_DMSOpct/100;

%% checks
DrugNames = ToColumn(DrugNames);
assert((length(DrugNames)==length(HMSLids)) || isempty(HMSLids))
assert(length(DrugNames)==length(SingleDoses))

assert(isempty(DrugPairs) || (~isempty(ComboDoses) || ~isempty(ComboLists)))
assert(all(DrugPairs(:)<=length(DrugNames)) || isempty(DrugPairs))
assert((length(DrugNames)==length(ComboDoses)) || isempty(ComboDoses))
assert(size(DrugPairs,2)==2)
assert(all(size(DrugPairs)==size(ComboLists)) || isempty(ComboLists))

assert((length(DrugNames)==length(stock_conc)) || (length(stock_conc)==1))
assert(all(stock_conc>20 & stock_conc<5e4), 'Stock concentration should be in uM')
if length(stock_conc)==1
    stock_conc = stock_conc*ones(length(DrugNames),1);
end
stock_conc = num2cell(ToColumn(stock_conc));

if length(plate_dims)==1
    plate_dims = sqrt(plate_dims*[1/1.5 1.5]);
end
assert(round(log2(plate_dims(1)))==log2(plate_dims(1)))

nWells = prod(plate_dims);
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


Drugs = struct('DrugName', DrugNames, 'HMSLid', ToColumn(HMSLids),...
    'stock_conc', stock_conc, 'layout', zeros(plate_dims));
Designs = struct('plate_dims', repmat({plate_dims}, nReps, 1), ...
    'treated_wells', repmat({true(plate_dims)}, nReps, 1), ...
    'well_volume', repmat({well_volume}, nReps, 1), ...
    'Drugs', repmat({Drugs}, nReps, 1), 'Seed', num2cell(Seed-1+(1:nReps)'));

SingleDoses = cellfun2(@(x,y,z) round_Doses(x,y,z,'Single',min_volume, ...
    step_volume, max_volume, well_volume), SingleDoses, stock_conc, DrugNames);

% if isempty(ComboLists) && ~isempty(ComboDoses)
% ComboLists = [ComboDoses(DrugPairs(:,1)) ComboDoses(DrugPairs(:,2))];
% end

if ~isempty(ComboDoses)
    ComboDoses = cellfun2(@(x,y,z) round_Doses(x,y,z,'Combo',min_volume, ...
        step_volume, max_volume, well_volume), ComboDoses, stock_conc, DrugNames);
    if any(cellfun(@(x,y) any(~ismember(x,y) & ~isempty(y)), ComboDoses, SingleDoses))
        fprintf('\n')
        warnprintf('Some doses for the combo are not part of the single doses for %s', ...
            strjoin(ToRow(DrugNames(cellfun(@(x,y) any(~ismember(x,y) & ~isempty(y)), ComboDoses, SingleDoses))),', '))
    end
    
    nTreatments = sum(cellfun(@length,SingleDoses)) + ...
        sum(cellfun(@(x,y) length(x)*length(y), ComboDoses(DrugPairs(:,1)), ...
        ComboDoses(DrugPairs(:,2))));
    
elseif ~isempty(ComboLists)
    ComboLists(:) = cellfun2(@(x,y,z) round_Doses(x,y,z,'Combo',min_volume, ...
        step_volume, max_volume, well_volume), ComboLists(:), stock_conc, DrugNames);
    if any(cellfun(@(x,y) any(~ismember(x,y) & ~isempty(y)), ComboLists(:), SingleDoses))
        fprintf('\n')
        warnprintf('Some doses for the combo are not part of the single doses for %s', ...
            strjoin(ToRow(DrugNames(cellfun(@(x,y) any(~ismember(x,y) & ~isempty(y)), ComboLists(:), SingleDoses))),', '))
    end
    
    nTreatments = sum(cellfun(@length,SingleDoses)) + ...
        sum(cellfun(@(x,y) length(x)*length(y), ComboLists(:,1), ...
        ComboLists(:,2)));
else
    nTreatments = sum(cellfun(@length,SingleDoses));
end



%% define the position of the fixed controls
fprintf('Number of test wells: %i\n',nTreatments);
ctrl_cnt = nWells - nTreatments;
assert(ctrl_cnt>5, 'Too many well used (%i out of %i), need at least 6 control wells',...
    nTreatments, nWells)


if edge_ctrl
    [ctrlidx, treated_wells] = DefineFixedControlPositions(plate_dims, ctrl_cnt);
    for iR=1:nReps
        Designs(iR).treated_wells = treated_wells;
    end
else
    ctrlidx = [];
end

%% define all possible treatments
allTreatments = zeros(length(DrugNames), nTreatments);
cnt = 0;
for iD = 1:length(DrugNames)
    allTreatments(iD, cnt+(1:length(SingleDoses{iD}))) = SingleDoses{iD};
    cnt = cnt + length(SingleDoses{iD});
end
for iCo = 1:size(DrugPairs,1)
    if ~isempty(ComboDoses)
        for iD1 = 1:length(ComboDoses{DrugPairs(iCo,1)})
            allTreatments(DrugPairs(iCo,1), cnt+(1:length(ComboDoses{DrugPairs(iCo,2)}))) = ...
                ComboDoses{DrugPairs(iCo,1)}(iD1);
            allTreatments(DrugPairs(iCo,2), cnt+(1:length(ComboDoses{DrugPairs(iCo,2)}))) = ...
                ComboDoses{DrugPairs(iCo,2)};
            cnt = cnt + length(ComboDoses{DrugPairs(iCo,2)});
        end
    else % using DrugList
        for iD1 = 1:length(ComboLists{iCo,1})
            allTreatments(DrugPairs(iCo,1), cnt+(1:length(ComboLists{iCo,2}))) = ...
                ComboLists{iCo,1}(iD1);
            allTreatments(DrugPairs(iCo,2), cnt+(1:length(ComboLists{iCo,2}))) = ...
                ComboLists{iCo,2};
            cnt = cnt + length(ComboLists{iCo,2});
        end
        
    end
    
end

assert(cnt==nTreatments)
assert(all(any(allTreatments>0)))


%% assign the design for each replicate
for iR = 1:nReps
    
    Designs(iR).Drugs = RandomizePlatePositions(Designs(iR).Drugs, ...
        allTreatments, nWells, plate_dims, ctrlidx, ctrl_cnt, Seed+iR-1);
    
    allDrugs = reshape([Designs(iR).Drugs.layout], [plate_dims length(DrugNames)]);
    nDrugs = sum(allDrugs>0,3);
    assert(all(squeeze(sum(sum((allDrugs>0).*repmat(nDrugs==1,1,1,length(DrugNames)),2),1))==...
        cellfun(@length,SingleDoses)))
%     for iCo = 1:size(DrugPairs,1)
%         assert(sum(sum( all(allDrugs(:,:,DrugPairs(iCo,:))>0,3)))== ...
%             length(ComboDoses{DrugPairs(iCo,1)})*length(ComboDoses{DrugPairs(iCo,2)}))
%     end
    
end


end



