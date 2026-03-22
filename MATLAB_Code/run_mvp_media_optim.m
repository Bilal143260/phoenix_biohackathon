function results = run_mvp_media_optim(matModelPath)
    addpath(genpath('/Users/aordl/Library/Application Support/MathWorks/MATLAB Add-Ons/Collections/RAVEN Toolbox'));

    pickSolver(); % <-- won’t pick gurobi unless it exists

    S = load(matModelPath,'model');
    model = S.model;

    % Draft-model robustness: remove forced positive fluxes (often breaks feasibility)
    model.lb(model.lb > 0) = 0;

    % --- Find exchanges/boundary rxns robustly (works even if not named EX_*) ---
    exIdx = find(sum(model.S ~= 0, 1) == 1);  % reactions involving exactly 1 metabolite

    Tb = listBoundaryRxns(model, exIdx);
    disp(Tb(1:min(50,height(Tb)),:));   % inspect first 50


    % --- Define a starting medium: close all uptakes ---
    % Uptake convention: negative flux = uptake -> enforce lb = 0 to block uptake
    model.lb(exIdx) = 0;

    % Optional: open "all" uptakes weakly to test feasibility, then refine
    % model.lb(exIdx) = -1;

    % --- Objective: biomass if available (otherwise we stop with guidance) ---
    % bioIdx = find(contains(lower(model.rxnNames),"biomass") | contains(lower(model.rxns),"biomass"), 1);
    [model,objIdx,objLabel] = setObjectiveBiomassElseATPM(model);
    fprintf("Objective set to: %s (%s)\n", model.rxns{objIdx}, objLabel);


    % if isempty(bioIdx)
    %     warning(['No biomass reaction found (common for pure KEGG drafts). ' ...
    %              'You can still do feasibility checks, but growth/media optimisation ' ...
    %              'needs a biomass equation. See section 3 below.']);
    % 
    %     % Feasibility solve (objective = 0)
    %     tmp = model;
    %     tmp.c = zeros(numel(tmp.rxns),1);
    %     sol = solveLP(tmp);
    % 
    %     results = struct();
    %     results.isFeasible = isfield(sol,'stat') && sol.stat == 1;
    %     results.note = 'No biomass objective available, growth optimisation skipped.';
    %     fprintf('Feasible under current bounds? %d\n', results.isFeasible);
    %     return;
    % end

    model.lb(exIdx) = -1;   % allow small uptake of anything (discovery)

    model.c = zeros(numel(model.rxns),1);
    model.c(bioIdx) = 1;
    model.osenseStr = 'max';


    if muRef < 1e-9
        error("Reference objective is ~0. Open a carbon/energy source uptake before essentiality scan.");
    end

    % --- Solve FBA (max growth) ---
    sol = solveLP(model);
    assert(sol.stat == 1, 'FBA failed: infeasible/unbounded under current medium.');

    v  = sol.x;
    mu = sol.f;
    fprintf('Max growth = %.6f\n', mu);

    % --- Outputs for non-specialists ---
    % 1) What nutrients are consumed (uptake fluxes among boundary reactions)
    uptakeIdx = exIdx(v(exIdx) < -1e-9);
    T_uptake = table(model.rxns(uptakeIdx), model.rxnNames(uptakeIdx), v(uptakeIdx), ...
        'VariableNames', {'rxn','name','flux'});
    T_uptake = sortrows(T_uptake,'flux'); % most negative at top

    figure('Name','Uptakes used at optimum');
    n = min(15,height(T_uptake));
    bar(T_uptake.flux(1:n));
    set(gca,'XTick',1:n,'XTickLabel',T_uptake.rxn(1:n));
    xtickangle(45);
    ylabel('Flux (negative = uptake)');
    title('Top consumed media components (model optimum)');

    % 2) Essential nutrient screen (close one uptake at a time)
    T_ess = essentialUptakeScan(model, exIdx, mu);

    figure('Name','Essential nutrients');
    ess = T_ess(T_ess.isEssential,:);
    if ~isempty(ess)
        bar(categorical(ess.rxn), ess.growthRatio);
        ylim([0 1]);
        ylabel('Growth ratio after removing uptake');
        title('Essential nutrient uptakes (growth collapses when removed)');
    else
        text(0.1,0.5,'No essential uptakes detected under this setup'); axis off;
    end

    % 3) Dose-response for the strongest uptake (helps “optimal ranges”)
    if height(T_uptake) >= 1
        rxn = T_uptake.rxn{1};
        figure('Name',['Dose-response: ' rxn]);
        doseResponseCurve(model, rxn, 0:2:40);
    end

    % Return results
    results = struct();
    results.growth = mu;
    results.uptakeTable = T_uptake;
    results.essentialTable = T_ess;
end

function T = essentialUptakeScan(model, exIdx, muRef)
    rxns = model.rxns(exIdx);
    growthRatio = nan(numel(exIdx),1);

    for k = 1:numel(exIdx)
        m2 = model;
        m2.lb(exIdx(k)) = 0;  % close uptake for that boundary rxn
        sol = solveLP(m2);
        if isfield(sol,'stat') && sol.stat == 1
            growthRatio(k) = sol.f / muRef;
        else
            growthRatio(k) = 0;
        end
    end

    T = table(rxns, growthRatio, growthRatio < 0.01, ...
        'VariableNames', {'rxn','growthRatio','isEssential'});
    T = sortrows(T,'growthRatio','ascend');
end

function doseResponseCurve(model, rxnID, uptakeCaps)
    idx = find(strcmpi(model.rxns, rxnID), 1);
    if isempty(idx), title('Rxn not found'); return; end

    mu = nan(numel(uptakeCaps),1);
    for i = 1:numel(uptakeCaps)
        m2 = model;
        m2.lb(idx) = -abs(uptakeCaps(i));  % allow more uptake
        sol = solveLP(m2);
        if isfield(sol,'stat') && sol.stat == 1
            mu(i) = sol.f;
        end
    end

    plot(uptakeCaps, mu, '-o');
    xlabel(['Max uptake capacity for ' rxnID]);
    ylabel('Max growth');
    grid on;
end

% ----------------- Helpers -----------------
function solver = pickSolver()
    % Prefer commercial solvers if installed; otherwise fall back.
    if exist('gurobi','file') == 2
        solver = 'gurobi';
    elseif exist('cplexlp','file') == 2 || exist('cplex','file') == 2
        solver = 'cplex';
    elseif exist('glpkmex','file') == 3 || exist('glpk','file') == 2
        solver = 'glpk';
    else
        error(['No LP solver found. Install one of: ' ...
               'Gurobi, CPLEX, or GLPK (glpkmex).']);
    end

    setRavenSolver(solver);
    fprintf('Using solver: %s\n', solver);
end


function model = setBiomassObjective(model)
    bioIdx = find(contains(lower(model.rxnNames),"biomass") | contains(lower(model.rxns),"biomass"), 1);
    model.c = zeros(numel(model.rxns),1);
    model.osenseStr = 'max';
    if ~isempty(bioIdx)
        model.c(bioIdx) = 1;
    else
        model.c(1) = 1; % fallback so the LP is well-defined
        warning('No biomass found; using rxn 1 as placeholder objective.');
    end
end

function exIdx = findExchangeRxnsHeuristic(model)
    % Heuristic only. Adjust once you know your model’s exchange naming.
    exIdx = find(startsWith(model.rxns,"EX_") | ...
                 contains(lower(model.rxnNames),"exchange") | ...
                 startsWith(model.rxnNames,"EX_"));
end

function model = closeAllUptakes(model)
    exIdx = findExchangeRxnsHeuristic(model);
    % close uptake: lb = 0 (so flux cannot go negative)
    model.lb(exIdx) = 0;
end

function model = openUptakeIfExists(model, rxnID, maxUptake)
    idx = find(strcmpi(model.rxns, rxnID), 1);
    if ~isempty(idx)
        model.lb(idx) = -abs(maxUptake);
    end
end

function Tb = listBoundaryRxns(model, exIdx)
    metIdx = zeros(numel(exIdx),1);
    for k = 1:numel(exIdx)
        j = exIdx(k);
        metIdx(k) = find(model.S(:,j)~=0, 1, 'first');
    end

    Tb = table( ...
        model.rxns(exIdx), model.rxnNames(exIdx), ...
        model.mets(metIdx), model.metNames(metIdx), ...
        model.lb(exIdx), model.ub(exIdx), ...
        'VariableNames', {'rxn','rxnName','met','metName','lb','ub'} ...
    );
end

function [model,objIdx,objLabel] = setObjectiveBiomassElseATPM(model)
    % Try biomass first
    bioIdx = find(contains(lower(model.rxnNames),"biomass") | contains(lower(model.rxns),"biomass"), 1);
    if ~isempty(bioIdx)
        model.c = zeros(numel(model.rxns),1);
        model.c(bioIdx) = 1;
        model.osenseStr = 'max';
        objIdx = bioIdx;
        objLabel = "biomass";
        return;
    end

    % Else try ATP maintenance
    atpmIdx = find(strcmpi(model.rxns,'ATPM') | ...
                   contains(lower(model.rxnNames),'atp maintenance') | ...
                   contains(lower(model.rxnNames),'maintenance'), 1);

    if isempty(atpmIdx)
        error("No biomass and no ATP maintenance (ATPM) found. Need to add an objective reaction.");
    end

    model.c = zeros(numel(model.rxns),1);
    model.c(atpmIdx) = 1;
    model.osenseStr = 'max';
    objIdx = atpmIdx;
    objLabel = "ATPM_proxy";
end
