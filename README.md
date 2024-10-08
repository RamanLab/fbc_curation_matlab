<img title="FROG logo" src="https://github.com/RamanLab/fbc_curation_matlab/blob/main/FROG_analysis_white_BG_outline.svg" height="100">


# fbc_curation_matlab
fbc_curation_matlab contains MATLAB/COBRA helpers for checking reproducibility of fbc models.

## Pre-requisites
* [CobraToolBox](https://opencobra.github.io/cobratoolbox/stable/installation.html)

## Components in the repository
* The main file for running fbc_curation in MATLAB is present in [`run_FROG.m` from the `src/curator` directory](https://github.com/RamanLab/fbc_curation_matlab/blob/main/src/fbc_curation/curator/run_FROG.m)
* An additional function [GetMD5](https://in.mathworks.com/matlabcentral/fileexchange/25921-getmd5) (Copyright (c) 2017-2019, Jan Simon). This is required to generate MD5 signature for model files.

## Installation
1. Clone the repository to your system.
2. Extract the folders.
3. Add the folders to you MATLAB path using [`addpath(genpath(path to fbc_curation_matlab))`](https://www.mathworks.com/help/matlab/ref/addpath.html)
4. Installation is complete!

## Usage
Step 1: Initialize CobraToolBox using the command 
```
initCobraToolbox
```
fbc_curation_matlab uses functions from CobraToolBox. So it has to be initialized every time MATLAB is started. Running this command once every session is enough.

Step 2: Running FBC-curation
```
run_FROG(path to model)
```
Step 2 with FBC Curator name
```
run_FROG(path to model, curator_name)
```
This command will start the fbc curation and FROG report files will be available in the current directory as [COMBINE archive](https://co.mbine.org/documents/archive) format. It is basically a zip file with all the report files in a specific file structure.

## Verify installation
To make sure if the tool is installed properly, run the following command.
```
run_FROG(path to fbc_curation_matlab/src/fbc_curation/examples/models/e_coli_core.xml)
```
If the command runs and give the following output without any error. The installation is successfull.
```
Using gurobi solver.
Loading model from e_coli_core.xml... Elapsed time is 3.178258 seconds.
Created directory ./FROG successfully.
[00] Wrote Metadata details to FROG/metadata.json.
[01] Wrote FBA objective results to FROG/01_objective.tsv.
[02] Wrote FVA results (optPercentage = 100) to FROG/02_fva.tsv.
Single gene deletion analysis in progress ...
100%    [........................................]
[03] Wrote gene deletion results to FROG/03_gene_deletion.tsv.
Single reaction deletion analysis in progress ...
100%    [........................................]
[04] Wrote gene deletion results to FROG/04_reaction_deletion.tsv.
Total Elapsed time is 4.759139 seconds.
Created COMBINE archive file ./e_coli_core.zip successfully.
e_coli_core.zip file renamed to ./e_coli_core.omex successfully.
```

## Example Output
```
>> run_FROG('fbc_curation_matlab-main/src/fbc_curation/examples/models/iJR904.xml', 'Karthik Raman');
Using glpk solver.
Loading model from fbc_curation_matlab-main/src/fbc_curation/examples/models/iJR904.xml... Elapsed time is 39.618113 seconds.
Created directory ./FROG successfully.
FROG curator name:Karthik Raman
[00] Wrote Metadata details to FROG/00_metadata.json.
[01] Wrote FBA objective results to FROG/01_objective.tsv.
[02] Wrote FVA results (optPercentage = 100) to FROG/02_fva.tsv.
Single gene deletion analysis in progress ...
100%    [........................................]
[03] Wrote gene deletion results to FROG/03_gene_deletion.tsv.
Single reaction deletion analysis in progress ...
100%    [........................................]
[04] Wrote gene deletion results to FROG/04_reaction_deletion.tsv.
Total Elapsed time is 111.763478 seconds.
Created COMBINE archive file ./iJR904.zip successfully.
iJR904.zip file renamed to ./iJR904.omex successfully.
>> 
```
See the [files in `fbc_curation/examples/results/`](https://github.com/RamanLab/fbc_curation_matlab/tree/main/src/fbc_curation/examples/results) for how the output files look like.

## Possible minor issues and troubleshooting
* In macos systems, you might encounter some security prompts saying 'The system cannot open GetMD5.mexmaci64 file as the devoluper is unknown'. This can be solved by going into System preferences > Security & Privacy > General tab > allow the GetMD5.mexmaci64 file in 'allow apps downloaded from' section. Now, running the `run_FROG` will generate the output successfully.

## License

* [LGPLv3](http://opensource.org/licenses/LGPL-3.0)

``fbc_curation_matlab`` source is released under both the GPL and LGPL licenses version 2 or
later. You may choose which license you choose to use the software under.

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License or the GNU Lesser General Public
License as published by the Free Software Foundation, either version 2 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.


## Acknowledgements

* The author thanks the entire SBML community, particularly the [SBML Team](https://github.com/sbmlteam), the SBML 2020 HARMONY & COMBINE participants for useful discussions. Specials thanks to [Matthias König](https://github.com/matthiaskoenig) and his excellent [`fbc-curation` Python Package](https://github.com/matthiaskoenig/fbc_curation), [Brett Olivier](https://github.com/bgoli) and [Rahuman Sheriff](https://www.ebi.ac.uk/about/people/rahuman-sheriff).
* [Initiative for Biological Systems Engineering](https://ibse.iitm.ac.in/)
* [Robert Bosch Centre for Data Science and Artificial Intelligence (RBCDSAI)](https://rbcdsai.iitm.ac.in/)

<img title="IBSE logo" src="https://github.com/RBC-DSAI-IITM/rbc-dsai-iitm.github.io/blob/master/images/IBSE_logo.png" height="100"><img title="RBC-DSAI logo" src="https://github.com/RBC-DSAI-IITM/rbc-dsai-iitm.github.io/blob/master/images/logo.jpg" height="100">

© 2020 Karthik Raman
