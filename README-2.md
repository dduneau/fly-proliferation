# Fly Proliferation

**Version 1.0.** A web app that compares the early proliferation of bacteria in two groups of flies, as described in the protocol (section "Early proliferation"). It fits a generalized additive model of log2 load against time, in which both curves start at the inoculum measured by injecting into PBS. From it, the app gives:

- where the groups differ, as a simultaneous band of the difference curve;
- five parameters of each curve and their differences, with 95% intervals: growth, lag, onset of control, time to control and maximum;
- the protocol's figure, a methods paragraph and the R code of the analysis.

The analysis is R itself (package `mgcv`). The app is built with shinylive, so R runs in the browser: there is no server, and the data are never uploaded. The first visit downloads R and the packages, which takes a few seconds.

## The data

The **per-fly table saved by [CFU Drop Counter](https://dduneau.github.io/cfu-drop-counter/)** (Load step, log2), with one row per fly and per inoculum tube. The loads are estimated there, once, so that every app uses the same estimate. The layout file given to CFU Drop Counter must therefore carry the columns this analysis needs: the group (sex, genotype …), the time post-injection and, if there is more than one, the replicate experiment.

- **The inoculum tubes of the Time 0 plate are part of the table, at time 0.** They are the common start of both curves. The dose was measured in PBS and belongs to neither group, so the app enters it once, under the reference group (protocol, step 2), whatever the label of these rows. The app stops if no row is at time 0. If some time 0 rows are not the injected dose, **Time 0 rows used as the inoculum** keeps only those with one label. The app warns when the time 0 rows differ by label by more than twofold.
- **Times are in hours.** Use only times before any mortality (protocol, CRITICAL note in "Quantification of Early Bacterial Proliferation"). Times after 12 h start unselected, since the protocol stops there for *P. rettgeri*.
- **Only flies with the status `Detected` enter the model.** Flies without a load are listed in the Flies tab and counted in a warning, and are shown at their detection limit in the plot of flies.

Any table with one row per fly and a log2 (or log10) load works as well, separated by commas, semicolons or tabs, with decimal points or decimal commas. Column names are free: you tell the app which column plays which role.

### Example data set

`example_proliferation_CS_sexdim.csv` is the per-fly table saved by CFU Drop Counter for the experiment of the protocol's figure. Canton-S females and males were injected with *P. rettgeri* (OD 0.1) and sampled hourly to 16 h, 8 per sex and time point. The loads were estimated in the Load step with the protocol's settings: censored drops above a ceiling of 30, a homogenate of 250 µL, log2. They agree with the lab's R estimator to 10⁻¹³.

This experiment predates the protocol's plate design. Each dilution was spotted on its own lid, so in CFU Drop Counter each lid was given its own plate name (plate01_dil5 …) and each fly one layout row per lid. Its time 0 rows are the 8 inoculum tubes, labelled `Time 0`; a second set of 8 tubes was not used as the inoculum and is left out of the file. The app opens the example with seed 99, the seed of the original analysis, so it gives the figure of the protocol.

## Using the app

1. **Data.** Load the per-fly table, or try the example.
2. **Comparison.** Choose the column of the **groups to compare**; it must hold exactly two groups among the flies kept. Then set their **order**: the first is the reference, which carries the inoculum, and differences are the second minus the first.
   - **Keep only flies where** restricts the analysis, for example to one genotype when the comparison is between the sexes. The inoculum tubes are always kept.
   - **Replicate experiment** (optional): each replicate is fitted separately, with its own inoculum, and the differences are compared between replicates rather than pooled (protocol, Note on replicated experiments).
3. **Model.**
   - **Time points used**: leave some out to see how much the parameters depend on the sampling density; time 0 is always kept.
   - **Basis dimension k**: start with 6, and increase it if the basis check shows too little flexibility.
   - **Random seed**: reported with the results.

   The 10,000 draws and the definition of time to control (growth below 10% of its maximum) are fixed, as in the protocol.
4. **Equivalence margins** (optional). Set them before looking at the results. A parameter is called equivalent when the 90% interval of its difference lies within the margin; otherwise the comparison is inconclusive. Without a margin, no equivalence is declared.

## What it gives

- **Curves**, the protocol's figure, in three panels. Above them are the difference for each of the five parameters, with its 95% interval and reading. PDF, PNG or SVG.
  - Fitted log2 load, with the end of the lag (dotted tangent at the steepest point) and the maximum.
  - Growth rate, with the fastest growth (onset of control) and the time to control.
  - The difference between the groups, with its simultaneous 95% band and the time from which the band excludes zero.
- **Flies over time**: every fly and the inoculum, with medians from the common start and the detection limit.
- **Results**:
  - the test of the difference smooth;
  - the basis check;
  - the stretches where the simultaneous band excludes zero;
  - the five parameters of each group and their differences, with 95% intervals, the 90% interval, the margin and the reading;
  - the share of draws in which each parameter exists in both groups;
  - with replicates, the differences of every replicate side by side.

  "Differ" means the 95% interval of the difference excludes zero. Five parameters are examined, and the protocol asks for the pattern to be weighed rather than one verdict.
- **Flies**: every row, with its role (fly or inoculum) and whether it entered the model.
- **Methods text**: a paragraph with the model, the settings and the results.
- **R code**: the functions the app runs (`prolif_core.R`), then the analysis with the settings chosen. It runs on the table from the Flies tab and reproduces the app's results exactly.

## Repository layout and updating

```
app_proliferation/   the app: app.R, prolif_core.R (the analysis) and the example data set (edit here)
docs/                the published site, made by shinylive from app_proliferation/ (do not edit)
```

To update the app, edit `app_proliferation/app.R` or `prolif_core.R`, then in R, from the repository folder:

```r
shinylive::export("app_proliferation", "docs")
file.create("docs/.nojekyll")
httpuv::runStaticServer("docs")   # check it locally
```

The empty file `docs/.nojekyll` lets GitHub Pages serve all of shinylive's files. Then commit and push. GitHub Pages serves `docs/` from the `main` branch.

To run the app in R instead: `shiny::runApp("app_proliferation")`.

## Citation

If you use this app, please cite the protocol: [authors, title, journal, year, DOI].
