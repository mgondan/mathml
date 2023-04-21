# mathml
### Translate R expressions to MathML

mathml allows R\ to render its terms in pretty mathematical equations
bridging the gap between computational needs, presentation of results, and their
reproducibility. 

Researchers or teachers can already use RMarkdown to conduct analyses and show
results, and mathml smoothes this process and allows for integrated
calculations and output. mathml can help to improve data analyses and 
statistical reports from an aesthetical perspective, as well as regarding 
reproducibility of research, by allowing also for a better detection of possible
mistakes in R\ programs. 

The package supports both MathML and Latex/MathJax for use in R\ Markdown documents, 
presentations and ShinyApp webpages.

## License

This R package is distributed under a 

## Installation

1. Download and install a recent R from https://www.r-project.org/

2. Download and install a recent RStudio from https://www.rstudio.com/

3. R> `install.packages("mathml")`


# Examples

## Example 1

````
library(mathml)
term <- quote(a - ((b + c)) - d*e + f*(g + h) + i/j + k^(l + m) + (n*o)^{p + q})
mathout(term)
````

$(a - ((b + c)) - d*e + f*(g + h) + i/j + k^(l + m) + (n*o)^{p + q})$


