function run_FROG(fileName)
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
% Creates a folder with model name as folder name, containing the below
% files
%    00_metadata.json:          File containing the metadata
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
fname_meta = sprintf('%s/%s', dir_name, '00_metadata.json');
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

prompt = 'FROG curator name:';
curator = input(prompt,'s');

fprintf(fid, '{\n');
fprintf(fid, '\t"FROG_date":\t"%s",\n', date);
fprintf(fid, '\t"FROG_version":\t"%s",\n', '1.0');
fprintf(fid, '\t"FROG_curators":\t"%s",\n', curator);
fprintf(fid, '"FROG_software":{\n');
fprintf(fid, '\t"name":\t"%s",\n', FROG_software_name);
fprintf(fid, '\t"version":\t"%s",\n', FROG_software_version);
fprintf(fid, '\t"url":\t"%s",\n', FROG_software_url);
fprintf(fid, '},\n');
fprintf(fid, '"software":{\n');
fprintf(fid, '\t"name":\t"%s",\n', software_name);
fprintf(fid, '\t"version":\t"%s",\n', software_version);
fprintf(fid, '\t"url":\t"%s",\n', software_url);
fprintf(fid, '},\n');
fprintf(fid, '"solver":{\n');
fprintf(fid, '\t"name":\t"%s",\n', solverName);
fprintf(fid, '\t"url":\t"Null"\n');
fprintf(fid, '},\n');
fprintf(fid, '\t"model.filename":\t"%s",\n', model.description);
fprintf(fid, '\t"model.md5":\t"%s",\n', Model_MD5);
fprintf(fid, '\t"environment":\t"%s, %s",\n', getenv('OS'), system_dependent('getwinsys'));
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
        status = 'Infeasible problem';
    case 1
        status = 'Optimal';
    case 2
        status = 'Unbounded solution';
    case 3
        status = 'Almost optimal solution';
    case -1
        status = 'Some other problem';
end

fprintf(fid, 'model\tobjective\tstatus\tvalue\n');
if (nnz(model.c) > 1)
    error('Model does not have a single objective reaction.');
end
fprintf(fid, '%s\t%s\t%s\t%f\n', fileName, model.rxnNames{model.c~=0}, status, sol.f);
fprintf('[01] Wrote FBA objective results to %s.\n', fname_obj);
fclose(fid);

%% [02] FVA
fname_fva = sprintf('%s/%s', dir_name, '02_fva.tsv');
fid = fopen(fname_fva,'w');
optPercentage = 100;
[minFlux, maxFlux] = fluxVariability(model, optPercentage);

fprintf(fid, 'model\tobjective\treaction\tflux\tstatus\tminimum\tmaximum\tFracton_optimum\n');
nRxns = numel(model.rxns);
for k = 1:nRxns
    fprintf(fid, '%s\t%s\t%s\t%f\t%s\t%f\t%f\t%f\n', fileName, 'obj', model.rxns{k}, sol.x(k), 'optimal', minFlux(k), maxFlux(k), optPercentage/100);
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
curr_node.setAttribute('location', './manifest.xml')
curr_node.setAttribute('format', 'https://identifiers.org/combine.specifications:omex-manifest')
product.appendChild(curr_node);
files = {'00_metadata.json', '01_objective.tsv', '02_fva.tsv', '03_gene_deletion.tsv', '04_reaction_deletion.tsv'};
formats = {'https://identifiers.org/combine.specifications:frog-metadata-version-1' ...
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
zip(zip_file_name,{dir_name,'manifest.xml'});
