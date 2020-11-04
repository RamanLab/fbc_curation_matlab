# fbc_curation_matlab
fbc_curation contains MATLAB/COBRA helpers for reproducibility of fbc models.

## Installation
Just download the single script `fbc_curation.m` from the `src/curator` directory.

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

## Acknowledgements

The author thanks the entire SBML community, particularly @sbmlteam, the SBML 2020 HARMONY & COMBINE participants for usual discussions. Specials thanks to @matthiaskoenig, @bgoli and Rahuman Sheriff.
