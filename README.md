# fbc_curation_matlab
fbc_curation contains MATLAB/COBRA helpers for reproducibility of fbc models.

## Installation
Just download the single script [`fbc_curation.m` from the `src/curator` directory](https://github.com/RamanLab/fbc_curation_matlab/blob/main/src/fbc_curation/curator/fbc_curation.m).

## Example Output
```
>> fbc_curation('e_coli_core.xml');
Loading model from e_coli_core.xml... Elapsed time is 2.558526 seconds.
Created directory `e_coli_core` successfully.
[01] Wrote FBA objective results to e_coli_core/01_objective.tsv.
[02] Wrote FVA results (optPercentage = 100) to e_coli_core/02_fva.tsv.
Single gene deletion analysis in progress ...
100%    [........................................]
[03] Wrote gene deletion results to e_coli_core/03_gene_deletion.tsv.
Single reaction deletion analysis in progress ...
100%    [........................................]
[04] Wrote gene deletion results to e_coli_core/04_reaction_deletion.tsv.
Total Elapsed time is 4.136678 seconds.
>> 
```
See the [files in `fbc_curation/examples/results/`](https://github.com/RamanLab/fbc_curation_matlab/tree/main/src/fbc_curation/examples/results) for how the output files look like.

## Acknowledgements

* The author thanks the entire SBML community, particularly the [SBML Team](https://github.com/sbmlteam), the SBML 2020 HARMONY & COMBINE participants for useful discussions. Specials thanks to [Matthias König](https://github.com/matthiaskoenig) and his excellent [`fbc-curation` Python Package](https://github.com/matthiaskoenig/fbc_curation), [Brett Olivier](https://github.com/bgoli) and [Rahuman Sheriff](https://www.ebi.ac.uk/about/people/rahuman-sheriff).
* [Initiative for Biological Systems Engineering](https://ibse.iitm.ac.in/)
* [Robert Bosch Centre for Data Science and Artificial Intelligence (RBCDSAI)](https://rbcdsai.iitm.ac.in/)

<img title="IBSE logo" src="https://github.com/RBC-DSAI-IITM/rbc-dsai-iitm.github.io/blob/master/images/IBSE_logo.png" height="100"><img title="RBC-DSAI logo" src="https://github.com/RBC-DSAI-IITM/rbc-dsai-iitm.github.io/blob/master/images/logo.jpg" height="100">
