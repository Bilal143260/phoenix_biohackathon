% clear; clc;

% RAVEN path
addpath(genpath('/Users/aordl/Library/Application Support/MathWorks/MATLAB Add-Ons/Collections/RAVEN Toolbox'));

% (Optional) solver for feasibility checks / later FBA (GLPK is fine)
try
    setRavenSolver('glpk');
catch
end

% Common inputs
hmmSet = 'euk90_kegg105';
baseOutDir = fullfile(pwd, 'raven_batch_rxn');
if ~exist(baseOutDir,'dir'); mkdir(baseOutDir); end

jobs = struct( ...
  'orgID', { ...
    'GCF_000182925.2', ...
    'GCA_002083745.1', ...
    'GCA_003344705.1', ...
    'GCF_000184455.2'  ...
  }, ...
  'faaFile', { ...
    '/Users/aordl/Desktop/biohackathon26/GCF_000182925.2/protein.faa', ...
    '/Users/aordl/Desktop/biohackathon26/GCA_002083745.1/protein.faa', ...
    '/Users/aordl/Desktop/biohackathon26/GCA_003344705.1/protein.faa', ...
    '/Users/aordl/Desktop/biohackathon26/GCF_000184455.2/protein.faa'  ...
  } ...
);

for i = 1:numel(jobs)
    orgID_raw = jobs(i).orgID;
    faaFile   = jobs(i).faaFile;

    orgID = strrep(orgID_raw, '.', '_');  % SBML-safe id

    outDir = fullfile(baseOutDir, ['raven_' orgID]);
    if ~exist(outDir,'dir'); mkdir(outDir); end

    fprintf('\n==== [%d/%d] %s ====\n', i, numel(jobs), orgID_raw);

    try
        % 1) Build draft model (KEGG + HMM KO mapping)
        model = getKEGGModelForOrganism(orgID, faaFile, hmmSet, outDir, ...
            false,false,false,false, 1e-30,0.8,0.3,-1);

        % 2) Draft robustness: remove forced positive fluxes (common infeasibility cause)
        model.lb(model.lb > 0) = 0;

        % 3) Add exchanges (so you can define medium by opening uptakes)
        %    This is MVP-style: create EX_ for dead-end metabolites.
        model = addDeadEndExchanges(model, 800);  % tune 200–2000 depending on size

        % 4) Ensure there is a solvable objective:
        %    - use biomass if present
        %    - else add ATPM using KEGG compounds and maximize it (proxy objective)
        model = addATPM_KEGG(model);
        [model, objLabel] = setObjectiveBiomassElseATPM(model);
        fprintf('Objective: %s\n', objLabel);

        try
            % If your model uses KEGG-style H+, it may be C00080 (often safer than 'H+')
            % Try one; if it fails, try the other.
            try
                model = removeBadRxns(model, 1, {'C00080'}, true);
            catch
                model = removeBadRxns(model, 1, {'H+'}, true);
            end
        catch ME
            fprintf(2, '[%s] removeBadRxns FAILED. Continuing without it.\n', orgID_raw);
            fprintf(2, 'Reason: %s\n', ME.message);
        
            % Optional: write full debug report to file
            fid = fopen(fullfile(outDir, 'removeBadRxns_error.txt'), 'w');
            if fid ~= -1
                fprintf(fid, '%s\n', getReport(ME, 'extended'));
                fclose(fid);
            end
        end
        % 6) Export
        xmlFile   = fullfile(outDir, [orgID '.xml']);
        matFile   = fullfile(outDir, [orgID '.mat']);
        excelFile = fullfile(outDir, [orgID '.xlsx']);

        model = padModelForExport(model);
        nRxns = numel(model.rxns);

        if isfield(model,'eccodes')
            % normalize to cell column
            if isstring(model.eccodes); model.eccodes = cellstr(model.eccodes); end
            if ischar(model.eccodes);   model.eccodes = cellstr(model.eccodes); end
            model.eccodes = model.eccodes(:);
        
            % pad
            if numel(model.eccodes) < nRxns
                model.eccodes(end+1:nRxns,1) = {''};
            end
        else
            model.eccodes = repmat({''}, nRxns, 1);
        end
        save(matFile, 'model');

        exportModel(model, xmlFile);

        exportToExcelFormat(model, excelFile);

        fprintf('Exported:\n  %s\n  %s\n  %s\n', xmlFile, matFile, excelFile);

    catch ME
        fprintf(2, 'FAILED for %s\n%s\n', orgID_raw, getReport(ME, 'extended'));
    end
end

%% ---------- Local helper functions ----------

function ok = isFeasible(model)
    tmp = model;
    tmp.c = zeros(numel(tmp.rxns),1);
    sol = solveLP(tmp);
    ok = isfield(sol,'stat') && sol.stat == 1;
end

function model = addATPM_KEGG(model)
% Adds ATPM: ATP + H2O -> ADP + Pi + H+
% Using KEGG compound IDs:
% ATP C00002, ADP C00008, Pi C00009, H2O C00001, H+ C00080

    if any(strcmpi(model.rxns,'ATPM'))
        return
    end

    atp = find(contains(model.mets,'C00002'),1);
    adp = find(contains(model.mets,'C00008'),1);
    pi_ = find(contains(model.mets,'C00009'),1);
    h2o = find(contains(model.mets,'C00001'),1);
    h   = find(contains(model.mets,'C00080'),1);

    if any([isempty(atp), isempty(adp), isempty(pi_), isempty(h2o), isempty(h)])
        % If your mets aren’t KEGG-style, we skip silently (or you can error)
        warning('ATPM not added (missing KEGG ATP/ADP/Pi/H2O/H+ metabolites).');
        return
    end

    newCol = sparse(numel(model.mets),1);
    newCol(atp) = -1;
    newCol(h2o) = -1;
    newCol(adp) =  1;
    newCol(pi_) =  1;
    newCol(h)   =  1;

    model.S(:,end+1)        = newCol;
    model.rxns{end+1,1}     = 'ATPM';
    model.rxnNames{end+1,1} = 'ATP maintenance';
    model.lb(end+1,1)       = 0;
    model.ub(end+1,1)       = 1000;
    if isfield(model,'rev'); model.rev(end+1,1) = 0; end
end

function model = addDeadEndExchanges(model, maxAdded)
% Adds exchange reactions to dead-end metabolites (MVP approach).
% This creates EX_<metID> with stoich -1*met.

    if nargin < 2; maxAdded = inf; end

    S = model.S;
    nM = numel(model.mets);

    consOnly = false(nM,1);
    prodOnly = false(nM,1);

    for m = 1:nM
        coeffs = full(S(m, S(m,:)~=0));
        if isempty(coeffs); continue; end
        if all(coeffs < 0), consOnly(m) = true; end
        if all(coeffs > 0), prodOnly(m) = true; end
    end

    deadEnds = find(consOnly | prodOnly);
    deadEnds = deadEnds(1:min(numel(deadEnds), maxAdded));

    for k = 1:numel(deadEnds)
        m = deadEnds(k);

        rxnID = ['EX_' regexprep(model.mets{m},'[^a-zA-Z0-9_]','_')];
        if any(strcmpi(model.rxns, rxnID))
            continue
        end

        newCol = sparse(nM,1);
        newCol(m) = -1;

        model.S(:,end+1)        = newCol;
        model.rxns{end+1,1}     = rxnID;
        model.rxnNames{end+1,1} = ['Exchange: ' model.metNames{m}];

        model.lb(end+1,1) = 0;      % uptake closed by default
        model.ub(end+1,1) = 1000;   % secretion allowed
        if isfield(model,'rev'); model.rev(end+1,1) = 1; end
    end
end

function [model, label] = setObjectiveBiomassElseATPM(model)
    bioIdx = find(contains(lower(model.rxnNames),"biomass") | contains(lower(model.rxns),"biomass"), 1);

    model.c = zeros(numel(model.rxns),1);
    model.osenseStr = 'max';

    if ~isempty(bioIdx)
        model.c(bioIdx) = 1;
        label = "biomass";
        return
    end

    atpmIdx = find(strcmpi(model.rxns,'ATPM'), 1);
    if ~isempty(atpmIdx)
        model.c(atpmIdx) = 1;
        label = "ATPM (proxy objective)";
        return
    end

    % last resort: placeholder (not useful for media optimisation)
    model.c(1) = 1;
    label = "placeholder objective (rxn1)";
    warning('No biomass/ATPM objective found. Using rxn1 placeholder (not meaningful).');
end

function model = padModelForExport(model)
% Pad common RAVEN fields so exportModel doesn't index past array lengths
% after you manually add reactions/metabolites.

nRxns  = numel(model.rxns);
nMets  = numel(model.mets);
nGenes = numel(model.genes);

%% ---------- Reaction-level fields ----------
% numeric vectors
model = padNumeric(model, 'lb', nRxns, 0);
model = padNumeric(model, 'ub', nRxns, 1000);
model = padNumeric(model, 'c',  nRxns, 0);
model = padNumeric(model, 'rev', nRxns, 0);
model = padNumeric(model, 'rxnConfidenceScores', nRxns, NaN);

% cell vectors
model = padCell(model, 'rxnNames', nRxns, '');
model = padCell(model, 'grRules',  nRxns, '');
model = padCell(model, 'rules',    nRxns, '');
model = padCell(model, 'rxnNotes', nRxns, '');
model = padCell(model, 'rxnReferences', nRxns, '');

% SBML/MIRIAM (cell-of-cells usually, but empty [] is fine)
model = padCell(model, 'rxnMiriams', nRxns, []);

% GPR matrix
if isfield(model,'rxnGeneMat')
    if size(model.rxnGeneMat,1) < nRxns
        model.rxnGeneMat(end+1:nRxns, 1:nGenes) = 0;
    end
end

%% ---------- Metabolite-level fields ----------
model = padCell(model, 'metNames', nMets, '');
model = padCell(model, 'metNotes', nMets, '');
model = padCell(model, 'metReferences', nMets, '');
model = padCell(model, 'metMiriams', nMets, []);
model = padCell(model, 'metFormulas', nMets, '');
model = padNumeric(model, 'metCharges', nMets, NaN);

%% ---------- Gene-level fields (sometimes used in export) ----------
model = padCell(model, 'geneShortNames', nGenes, '');
model = padCell(model, 'geneMiriams', nGenes, []);
end

%% ===== helpers =====
function model = padNumeric(model, field, n, fillVal)
    if isfield(model,field)
        v = model.(field);
        v = v(:);
        if numel(v) < n
            v(end+1:n,1) = fillVal;
        end
        model.(field) = v;
    else
        model.(field) = repmat(fillVal, n, 1);
    end
end

function model = padCell(model, field, n, fillVal)
    if isfield(model,field)
        c = model.(field);
        c = c(:);
        if numel(c) < n
            c(end+1:n,1) = {fillVal};
        end
        model.(field) = c;
    else
        model.(field) = repmat({fillVal}, n, 1);
    end
end
function model = padModelForExcel(model)
nRxns = numel(model.rxns);

% subSystems is used by exportToExcelFormat for the reaction sheet
if isfield(model,'subSystems') && ~isempty(model.subSystems)
    ss = model.subSystems;

    % normalize to cell column
    if isstring(ss); ss = cellstr(ss); end
    if ischar(ss);   ss = cellstr(ss); end
    if ~iscell(ss);  ss = {}; end
    ss = ss(:);

    % flatten cell-of-cells to a single string per reaction
    for i = 1:numel(ss)
        if iscell(ss{i})
            ss{i} = strjoin(string(ss{i}), '; ');
        elseif isstring(ss{i})
            ss{i} = char(ss{i});
        end
        if isempty(ss{i}); ss{i} = ''; end
    end

    % pad/truncate to nRxns
    if numel(ss) < nRxns
        ss(end+1:nRxns,1) = {''};
    elseif numel(ss) > nRxns
        ss = ss(1:nRxns);
    end
else
    ss = repmat({''}, nRxns, 1);
end

model.subSystems = ss;
end
