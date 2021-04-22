# fbc_curation_matlab
fbc_curation contains MATLAB/COBRA helpers for reproducibility of fbc models.

## Installation
Just download the single script [`fbc_curation.m` from the `src/curator` directory](https://github.com/RamanLab/fbc_curation_matlab/blob/main/src/fbc_curation/curator/fbc_curation.m).

Additional addon: GetMD5
(https://in.mathworks.com/matlabcentral/fileexchange/25921-getmd5).

## Example Output
```
>> fbc_curation('iJR904.xml');
Using glpk solver.
Loading model from fbc-curation-project/iJR904.xml... Elapsed time is 41.577884 seconds.
Created directory `iJR904` successfully.
[00] Wrote Metadata details to iJR904/00_metadata.json.
[01] Wrote FBA objective results to iJR904/01_objective.tsv.
[02] Wrote FVA results (optPercentage = 100) to iJR904/02_fva.tsv.
Single gene deletion analysis in progress ...
100%    [........................................]
[03] Wrote gene deletion results to iJR904/03_gene_deletion.tsv.
Single reaction deletion analysis in progress ...
100%    [........................................]
[04] Wrote rxn deletion results to iJR904/04_reaction_deletion.tsv.
Total Elapsed time is 100.817721 seconds.
>> 
```
See the [files in `fbc_curation/examples/results/`](https://github.com/RamanLab/fbc_curation_matlab/tree/main/src/fbc_curation/examples/results) for how the output files look like.

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
