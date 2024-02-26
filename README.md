# StackbarExtended: R Package for Visualizing and Analyzing Compositional Data


<br /> <br />
<a href="https://github.com/ThibaultCuisiniere/stackbarextended/issues/new?assignees=&labels=bug&template=01_BUG_REPORT.md&title=bug%3A+">Report
a Bug</a> ·
<a href="https://github.com/ThibaultCuisiniere/stackbarextended/issues/new?assignees=&labels=enhancement&template=02_FEATURE_REQUEST.md&title=feat%3A+">Request
a Feature</a> .
<a href="https://github.com/ThibaultCuisiniere/stackbarextended/issues/new?assignees=&labels=question&template=04_SUPPORT_QUESTION.md&title=support%3A+">Ask
a Question</a>

------------------------------------------------------------------------

## Overview

 StackbarExtended is an R package including methods to automatically plot and perform differential abundance analysis of compositional data.
A major advantage of StackbarExtended is that it integrates graphical representation of the samples' composition and differential abundance results into a single, publication-ready figure, while allowing a high level of personalization and being easy to use.

The main function `plot_gut_microbiota()` function uses in combination the [`DESeq2`](https://github.com/thelovelab/DESeq2) and [`ggplot2`](https://github.com/tidyverse/ggplot2) packages. 

<details>

<summary>Screenshot</summary>

<br>


|                               Graphical representation output                                |
|:----------------------------------:|
| <img src="docs/images/Example_graph.png" title="Graphical output" width="100%"/> | 



</details>

## Getting Started


### Installation

StackbarExtended package can be installed via Github. Users should have R installed on their computer before installing StackbarExtended. R can be downloaded here: http://www.r-project.org/. To install the latest version of StackbarExtended package via Github, run following commands in R:

```
if (!require("devtools"))
install.packages("devtools")
devtools::install_github("ThibaultCuisiniere/StackbarExtended")

```

## Usage

```
data(GlobalPatterns)
ps_unfiltered <- GlobalPatterns

my_plot <- plot_gut_microbiota(
ps_object = ps_unfiltered,
exp_group = 'SampleType',
subset_group = c("Feces", "Skin"),
sample_name = 'X.SampleID',
hues = c("Greens", "Blues", "Oranges", "Purples"),
differential_analysis = T,
sig_lab = T,
fdr_threshold = 0.05
)

print(my_plot$plot)
print(my_plot$significant_table)

```


## Roadmap

See the [open
issues](https://github.com/ThibaultCuisiniere/stackbarextended/issues)
for a list of proposed features (and known issues).

-   [Top Feature
    Requests](https://github.com/ThibaultCuisiniere/stackbarextended/issues?q=label%3Aenhancement+is%3Aopen+sort%3Areactions-%2B1-desc)
    (Add your votes using the 👍 reaction)
-   [Top
    Bugs](https://github.com/ThibaultCuisiniere/stackbarextended/issues?q=is%3Aissue+is%3Aopen+label%3Abug+sort%3Areactions-%2B1-desc)
    (Add your votes using the 👍 reaction)
-   [Newest
    Bugs](https://github.com/ThibaultCuisiniere/stackbarextended/issues?q=is%3Aopen+is%3Aissue+label%3Abug)



Reach out to the maintainer at one of the following places:

-   [GitHub
    issues](https://github.com/ThibaultCuisiniere/stackbarextended/issues/new?assignees=&labels=question&template=04_SUPPORT_QUESTION.md&title=support%3A+)
-   Contact options listed on [this GitHub
    profile](https://github.com/ThibaultCuisiniere)

## Contributing

First off, thanks for taking the time to contribute! Contributions are
what make the open-source community such an amazing place to learn,
inspire, and create. Any contributions you make will benefit everybody
else and are **greatly appreciated**.


## License

This project is licensed under the **GNU General Public License v3**.

See [LICENSE](https://github.com/ThibaultCuisiniere/StackbarExtended/blob/master/LICENSE.txt) for more information.


## Citation


## Applications of StackbarExtended

