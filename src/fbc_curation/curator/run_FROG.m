function run_FROG(fileName, curator_name, objective_name)
% Performs FBC curation + FROG analysis on a given model
%
% USAGE:
%
%    fbc_curation(fileName, curator_name)
%
% INPUT:
%    fileName:                  COBRA model file
%    curator_name:              optional, character array
%    objective_name:            optional, name given to the objective (default: 'obj')
%
% OUTPUTS:
% Creates a folder with model name as folder name, containing the below
% files
%    metadata.json:          File containing the metadata
%    01_objective.tsv:          Computed growth rate ratio between deletion strain and wild type
%    02_fva.tsv:                Deletion strain growth rates (1/h)
%    03_gene_deletion.tsv:      Wild type growth rate (1/h)
%    04_reaction_deletion.tsv:  Does a reaction deletion affect anything
%
% Also creates a COMBINE archive file (.zip) which is a zipped version of
% above folder with an additional manifest.xml file
%    manifest.xml               manifest file will be at the root of the
%                               archive. It is an xml file containing the
%                               list of names of all files in the archive.
%                               It also describes each file's type and
%                               location inside the archive
% .. Authors:
%       - Karthik Raman 2020/11/04

problemTypeParams = parseSolverParameters('LP');
solver = problemTypeParams.solver;
if (isempty(solver))
    error('No solver found! Please run initCobraToolbox before running fbc_curation!');
else
    fprintf('Using %s solver.\n', solver);
end

t = tic; fprintf('Loading model from %s... ', fileName);

%Getting file name from the path
exact_file_name = split(fileName, '\');

%load the file
if (extractAfter(fileName,".") == 'xml')
        model = readCbModel(fileName);
elseif (extractAfter(fileName,".") == 'mat')
        updated_filename = replace(fileName, "'", '');
        data = load(updated_filename);
        varname = fieldnames(data);
        model = data.(varname{1});
end

toc(t);

dir_name = 'FROG';

[success,msg,~] = mkdir(dir_name);

if (success)
    fprintf('Created directory ./%s successfully.\n', dir_name);
else
    fprintf(2,[msg '\n']);
    return
end

%% [00] METADATA FILE
fname_meta = sprintf('%s/%s', dir_name, 'metadata.json');
fid = fopen(fname_meta,'w');
FROG_software_name='FBC-Curation-Matlab';
FROG_software_version = '1.0';
FROG_software_url='https://github.com/RamanLab/fbc_curation_matlab/';
software_name='CobraToolbox';
software_version = '3.0';
software_url='https://opencobra.github.io/cobratoolbox/stable/';
FID = fopen(fileName, 'r');
S = fread(FID, inf, 'uchar=>char');
fclose(FID);
Model_MD5 = GetMD5(S, '8bit');
solverName=solver;

if ~exist('curator_name','var')
    curator_name = '';
end

fprintf(fid, '{\n');
fprintf(fid, '  "software": {\n');
fprintf(fid, '    "frog": {\n');
fprintf(fid, '      "name":\t"%s",\n', FROG_software_name);
fprintf(fid, '      "version":\t"%s",\n', FROG_software_version);
fprintf(fid, '      "url":\t"%s"\n', FROG_software_url);
fprintf(fid, '    },\n');
fprintf(fid, '    "toolbox": {\n');
fprintf(fid, '      "name":\t"%s",\n', software_name);
fprintf(fid, '      "version":\t"%s",\n', software_version);
fprintf(fid, '      "url":\t"%s"\n', software_url);
fprintf(fid, '    },\n');
fprintf(fid, '    "solver":{\n');
fprintf(fid, '      "name":\t"%s",\n', solverName);
fprintf(fid, '      "url":\t"Null"\n');
fprintf(fid, '    }\n');
fprintf(fid, '  },\n');
fprintf(fid, '  "model_filename":\t"%s",\n', model.description);
fprintf(fid, '  "frog_date":\t"%s",\n', date);
fprintf(fid, '  "frog_version":\t"%s",\n', '1.0');
fprintf(fid, '  "frog_curators":\t"%s",\n', curator_name);
fprintf(fid, '  "model_md5":\t"%s",\n', Model_MD5);
fprintf(fid, '  "environment":\t"%s, %s"\n', getenv('OS'), system_dependent('getwinsys'));
fprintf(fid, '}');
fprintf('[00] Wrote Metadata details to %s.\n', fname_meta);
fclose(fid);

%% [01] FBA
fname_obj = sprintf('%s/%s', dir_name, '01_objective.tsv');
fid = fopen(fname_obj,'w');

sol = optimizeCbModel(model);

% Finding the status of the solution
switch sol.stat
    case 0
        status = 'infeasible problem';
    case 1
        status = 'optimal';
    case 2
        status = 'unbounded solution';
    case 3
        status = 'almost optimal solution';
    case -1
        status = 'some other problem';
end

fprintf(fid, 'model\tobjective\tstatus\tvalue\n');
if (nnz(model.c) == 0)
    error('Model does not have an objective reaction.');
end

if ~exist('objective_name','var')
    objective_name = 'obj';
end

fprintf(fid, '%s\t%s\t%s\t%f\n', exact_file_name{end}, objective_name, status, sol.f);
fprintf('[01] Wrote FBA objective results to %s.\n', fname_obj);
fclose(fid);

%% [02] FVA
fname_fva = sprintf('%s/%s', dir_name, '02_fva.tsv');
fid = fopen(fname_fva,'w');
optPercentage = 100;
[minFlux, maxFlux] = fluxVariability(model, optPercentage);

fprintf(fid, 'model\tobjective\treaction\tflux\tstatus\tminimum\tmaximum\tfraction_optimum\n');
nRxns = numel(model.rxns);
for k = 1:nRxns
    fprintf(fid, '%s\t%s\t%s\t%f\t%s\t%f\t%f\t%f\n', exact_file_name{end}, objective_name, model.rxns{k}, sol.x(k), 'optimal', minFlux(k), maxFlux(k), optPercentage/100);
end
fprintf('[02] Wrote FVA results (optPercentage = %d) to %s.\n', optPercentage, fname_fva);
fclose(fid);

%% [03] Gene deletion results 
fname_genedel = sprintf('%s/%s', dir_name, '03_gene_deletion.tsv');
fid = fopen(fname_genedel,'w');
[grRatio, grRateKO, grRateWT, hasEffect] = singleGeneDeletion(model);
nGenes = numel(model.genes);

fprintf(fid, 'model\tobjective\tgene\tstatus\tvalue\n');
for k = 1:nGenes
    if (~isnan(grRateKO(k)))
        fprintf(fid, '%s\t%s\t%s\t%s\t%f\n', exact_file_name{end}, objective_name, model.genes{k}, 'optimal', grRateKO(k));
    else
        fprintf(fid, '%s\t%s\t%s\t%s\t%f\n', exact_file_name{end}, objective_name, model.genes{k}, 'infeasible', grRateKO(k));
    end
end
fprintf('[03] Wrote gene deletion results to %s.\n', fname_genedel);
fclose(fid);

%% [04] Reaction deletion results 
fname_rxndel = sprintf('%s/%s', dir_name, '04_reaction_deletion.tsv');
fid = fopen(fname_rxndel,'w');
[grRatio, grRateKO, grRateWT, hasEffect] = singleRxnDeletion(model);

fprintf(fid, 'model\tobjective\treaction\tstatus\tvalue\n');
for k = 1:nRxns
    if (~isnan(grRateKO(k)))
        fprintf(fid, '%s\t%s\t%s\t%s\t%f\n', exact_file_name{end}, objective_name, model.rxns{k}, 'optimal', grRateKO(k));
    else
        fprintf(fid, '%s\t%s\t%s\t%s\t%f\n', exact_file_name{end}, objective_name, model.rxns{k}, 'infeasible', grRateKO(k));
    end
end
fclose(fid);
fprintf('[04] Wrote gene deletion results to %s.\n', fname_rxndel);
fprintf('Total '); toc(t);
%% Writing Manifest file for COMBINE archive
docNode = com.mathworks.xml.XMLUtils.createDocument('manifest');
manifest = docNode.getDocumentElement;
product = docNode.createElement('omexManifest');
product.setAttribute('xmlns','https://identifiers.org/combine.specifications/omex-manifest');
curr_node = docNode.createElement('content');
curr_node.setAttribute('location', '.')
curr_node.setAttribute('format', 'https://identifiers.org/combine.specifications:omex')
product.appendChild(curr_node);
curr_node = docNode.createElement('content');
curr_node.setAttribute('location', ['./' fileName])
curr_node.setAttribute('format', 'https://identifiers.org/combine.specifications:sbml')
curr_node.setAttribute('master', 'True')
product.appendChild(curr_node);
curr_node = docNode.createElement('content');
curr_node.setAttribute('location', './manifest.xml')
curr_node.setAttribute('format', 'https://identifiers.org/combine.specifications:omex-manifest')
product.appendChild(curr_node);
files = {'metadata.json', '01_objective.tsv', '02_fva.tsv', '03_gene_deletion.tsv', '04_reaction_deletion.tsv'};
formats = { ...
    'https://identifiers.org/combine.specifications:frog-metadata-version-1' ...
    'https://identifiers.org/combine.specifications:frog-objective-version-1' ...
    'https://identifiers.org/combine.specifications:frog-fva-version-1' ...
    'https://identifiers.org/combine.specifications:frog-genedeletion-version-1' ...
    'https://identifiers.org/combine.specifications:frog-reactiondeletion-version-1'};
for idx = 1:numel(files)
    curr_node = docNode.createElement('content');
    curr_node.setAttribute('location', ['./' dir_name '/' files{idx}])
    curr_node.setAttribute('format', formats{idx})
    product.appendChild(curr_node);
end
manifest.appendChild(product);
xmlwrite('manifest.xml',docNode);
%% Zipping the files
zip_file_name = replace(model.description, '.xml', '');
zip(zip_file_name,{dir_name,'manifest.xml',fileName});
fprintf('Created COMBINE archive file ./%s.zip successfully.\n', zip_file_name);

%% Renaming the zip archive to omex archive
[rename_status,msg,~] = movefile([zip_file_name, '.zip'], [zip_file_name, '.omex']);

if (rename_status)
    fprintf('%s file renamed to ./%s successfully.\n', [zip_file_name, '.zip'], [zip_file_name, '.omex']);
else
    fprintf(2,[msg '\n']);
    return
end
%% Deleting the temporary files
rmdir(dir_name, 's')
delete('manifest.xml')
