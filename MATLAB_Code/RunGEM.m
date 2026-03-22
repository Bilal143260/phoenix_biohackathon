clear; clc;

% RAVEN path
addpath(genpath('/Users/aordl/Library/Application Support/MathWorks/MATLAB Add-Ons/Collections/RAVEN Toolbox'));

% Common inputs
hmmSet = 'euk90_kegg105';
baseOutDir = fullfile(pwd, 'raven_batch_rxn');
if ~exist(baseOutDir,'dir'); mkdir(baseOutDir); end

% --- Provide ONE FAA per genome/assembly ---
% Update these paths to the correct protein FASTA for each ID
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

    % Normalize ID for folders/filenames (and optionally for RAVEN call)
    orgID = strrep(orgID_raw, '.', '_');   % e.g., GCF_000182925_2

    outDir = fullfile(baseOutDir, ['raven_' orgID]);
    if ~exist(outDir,'dir'); mkdir(outDir); end

    fprintf('\n==== [%d/%d] %s ====\n', i, numel(jobs), orgID_raw);

    try
        % 1) Build draft model
        model = getKEGGModelForOrganism(orgID, faaFile, hmmSet, outDir, ...
            false,false,false,false, 1e-30,0.8,0.3,-1);

        % 2) Optional: remove "something from nothing" reactions
        % model = removeBadRxns(model, 1, {'H+'}, true);
        % 2) Optional cleanup: remove "something from nothing" reactions
    % Wrap in try/catch so batch doesn't fail
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


        % 3) Set an objective (helps SBML-FBC / general solver setup)
        bioIdx = find(contains(lower(model.rxnNames),"biomass") | contains(lower(model.rxns),"biomass"), 1);
        model.c = zeros(numel(model.rxns),1);
        if ~isempty(bioIdx)
            model.c(bioIdx) = 1;
        else
            model.c(1) = 1; % placeholder
        end
        model.osenseStr = 'max';

        % 4) Export
        xmlFile   = fullfile(outDir, [orgID '.xml']);
        matFile   = fullfile(outDir, [orgID '.mat']);
        excelFile = fullfile(outDir, [orgID '.xlsx']);

        exportModel(model, xmlFile);            % SBML
        exportToExcelFormat(model, excelFile);  % Excel
        save(matFile, 'model');                 % MAT

        fprintf('Exported:\n  %s\n  %s\n  %s\n', xmlFile, matFile, excelFile);

    catch ME
        fprintf(2, 'FAILED for %s\n%s\n', orgID_raw, getReport(ME, 'extended'));
    end
end
