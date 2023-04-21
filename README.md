# mathml
### Translate R expressions to MathML

`mathml` allows rendering R terms in pretty mathematical equations
bridging the gap between computational needs, presentation of results, and their
reproducibility. 

Researchers or teachers can already use RMarkdown to conduct analyses and show
results, `mathml` smoothes this process and allows for integrated
calculations and output. The package `mathml` can help in fact to improve data analyses and 
statistical reports from an aesthetical perspective, as well as regarding 
reproducibility of research, by allowing also for a better detection of possible
mistakes in R programs. 

The package supports both MathML and Latex/MathJax for use in R Markdown documents, 
presentations and ShinyApp webpages.

## License

This R package is distributed under a BSD-2 simplified license (see the file LICENSE).

## Installation

1. Download and install a recent R from https://www.r-project.org/

2. Download and install a recent RStudio from https://www.rstudio.com/

3. R> `install.packages("mathml")`

# Example

````
library(mathml)
term <- quote(a^b + c*d - a^2*(a+d))
mathout(term)
````

$(a^b + c*(d+3) - a^2*(a+d))$
