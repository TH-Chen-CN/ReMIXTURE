![](images/rm_logo.png)

# ReMIXTURE

Implementation of the ReMIXTURE method to quantify how (genetic or other) diversity is distributed and shared between groups, and to create pretty and intuitive plots showing this on a globe.

![](images/rmDemoVitis.png)

- Outer Circle=Total diversity
- Inner Circle=Region-unique diversity
- Line width=Overlapped diversity

ReMIXTURE is an R package, and to run it requires:

1) A symmetrical numeric matrix of pairwise sample-to-sample distances (any distance metric will do in principle), whose rownames and colnames give the region to which the sample is assigned.
2) A data.table or data.frame with column names `region`, `lon`, and `lat`, providing the position (numeric) on the globe given to each region (character).

Alternatively, ReMIXTURE can construct the distance matrix directly from a `.vcf` or `.vcf.gz` file using `SNPRelate`, provided you also supply a `sample_table` with character columns named `sample` and `region`. Direct `.bcf` input is not supported by this lightweight path.

# History

The ReMIXTURE concept was first attempted in Tripodi & Rabanus-Wallace, et al. (2021) _Global range expansion history of pepper (*Capsicum spp.*) revealed by over 10,000 genebank accessions_. PNAS. Newer versions have very significant improvements. The algorithm currently in use is not published.

# Using ReMIXTURE.

Install using `devtools::install_github("https://github.com/mtrw/ReMIXTURE")`, and follow the tutorial in ?ReMIXTURE (in the examples section). If you want to build a distance matrix from VCF input, install `SNPRelate` from Bioconductor first.

For composite map figures, first inspect H behaviour with the existing ReMIXTURE workflow, then use the selected run for display:

```r
rm$run(...)
rm$plot_h_optimisation()
rm$plot_results_grid()
rm$plot_distance_densities(HdistFromRun = selected_run)

rm$plot_maps_composite(run = selected_run)
```

# Future

Any questions, suggestions, feedback please email me! tim.rabanuswallace@unimelb.edu.au.
