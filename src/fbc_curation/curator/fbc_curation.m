function fbc_curation(fileName)
% Performs FBC curation + FROG analysis on a given model
%
% USAGE:
%
%    fbc_curation(fileName)
%
% INPUT:
%    fileName:                  COBRA model file
%
% OUTPUTS:
%    01_objective.tsv:          Computed growth rate ratio between deletion strain and wild type
%    02_fva.tsv:                Deletion strain growth rates (1/h)
%    03_gene_deletion.tsv:      Wild type growth rate (1/h)
%    04_reaction_deletion.tsv:  Does a reaction deletion affect anything
%
% .. Authors:
%       - Karthik Raman 2020/11/04

t = tic; fprintf('Loading model from %s... ', fileName);
model = readCbModel(fileName);
toc(t);
nRxns = numel(model.rxns);

dir_name = replace(model.description, '.xml', '');

[success,msg,~] = mkdir(dir_name);

if (success)
    fprintf('Created directory `%s` successfully.\n', dir_name);
else
    fprintf(2,[msg '\n']);
    return
end


%% [01] FBA
fname_obj = sprintf('%s/%s', dir_name, '01_objective.tsv');
fid = fopen(fname_obj,'w');

sol = optimizeCbModel(model);

fprintf(fid, 'model\tobjective\tstatus\tvalue\n');
if (nnz(model.c) > 1)
    error('Model does not have a single objective reaction.');
end
fprintf(fid, '%s\t%s\t%s\t%f\n', fileName, model.rxnNames{model.c~=0}, sol.origStat, sol.f);
fprintf('[01] Wrote FBA objective results to %s.\n', fname_obj);
fclose(fid);

%% [02] FVA
fname_fva = sprintf('%s/%s', dir_name, '02_fva.tsv');
fid = fopen(fname_fva,'w');
optPercentage = 100;
[minFlux, maxFlux] = fluxVariability(model, optPercentage);

fprintf(fid, 'model\tobjective\treaction\tflux\tstatus\tminimum\tmaximum (optPercentage = %d)\n', optPercentage);
for k = 1:nRxns
    fprintf(fid, '%s\t%s\t%s\t%f\t%s\t%f\t%f\n', fileName, 'obj', model.rxns{k}, sol.x(k), 'optimal', minFlux(k), maxFlux(k));
end
fprintf('[02] Wrote FVA results (optPercentage = %d) to %s.\n', optPercentage, fname_fva);
fclose(fid);

%% [03] Reaction deletion results 
fname_genedel = sprintf('%s/%s', dir_name, '03_gene_deletion.tsv');
fid = fopen(fname_genedel,'w');
[grRatio, grRateKO, grRateWT, hasEffect] = singleGeneDeletion(model);
nGenes = numel(model.genes);

fprintf(fid, 'model\tobjective\tgene\tstatus\tvalue\n');
for k = 1:nGenes
    if (~isnan(grRateKO(k)))
        fprintf(fid, '%s\t%s\tG_%s\t%s\t%f\n', fileName, 'obj', model.genes{k}, 'optimal', grRateKO(k));
    else
        fprintf(fid, '%s\t%s\tG_%s\t%s\t%f\n', fileName, 'obj', model.genes{k}, 'infeasible', grRateKO(k));
    end
end
fprintf('[03] Wrote gene deletion results to %s.\n', fname_genedel);
fclose(fid);

%% [04] Gene deletion results 
fname_rxndel = sprintf('%s/%s', dir_name, '04_reaction_deletion.tsv');
fid = fopen(fname_rxndel,'w');
[grRatio, grRateKO, grRateWT, hasEffect] = singleRxnDeletion(model);

fprintf(fid, 'model\tobjective\treaction\tstatus\tvalue\n');
for k = 1:nRxns
    if (~isnan(grRateKO(k)))
        fprintf(fid, '%s\t%s\tR_%s\t%s\t%f\n', fileName, 'obj', model.rxns{k}, 'optimal', grRateKO(k));
    else
        fprintf(fid, '%s\t%s\tR_%s\t%s\t%f\n', fileName, 'obj', model.rxns{k}, 'infeasible', grRateKO(k));
    end
end
fclose(fid);
fprintf('[04] Wrote gene deletion results to %s.\n', fname_rxndel);
fprintf('Total '); toc(t);