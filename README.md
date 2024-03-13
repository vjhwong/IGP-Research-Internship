# IGP Research Internship
As a Research Intern at the Department of Immunology, Genetics, and Pathology, I contributed to a project exploring the application of deconvolution techniques to RNA-sequencing data.

This repository contains the data analysis work I did to research the utility of the SIA metric described in this article [The Lancet - Article](https://www.thelancet.com/journals/ebiom/article/PIIS2352-3964(23)00017-8/fulltext#:~:text=By%20combining%20the%20prognostic%20information,%2Dof%2Dart%20immune%20score.)

# SIA 
Signal of Immune Activation, SIA, is a metric that can be used to predict cancer patients' response to certain immunotherapy treatments. The metrics can be calculated in two ways. The first is the gene expression-based SIA which is calculated by the following fraction:
$\frac{CD8^+}{C1QX}$
where C1QX is one of C1QA, C1QB or C1QC. In the notebook, I compared the gene expression-based SIA score with the cell type-based SIA score. The genes selected in the gene expression-based SIA are meant to be proxies for the cell types CD8+ T cells and CD68+CD163+ Macrophages. I used [CIBERSORTx](https://cibersortx.stanford.edu/), a deconvolution tool developed by Stanford to impute cell fractions of the dataset I worked with which was a melanoma dataset. 

# Plot from the notebook
![plot](https://i.postimg.cc/V67DGZn6/fig-5-plot.png)
These strip plots display the four different SIA scores divided into the three categories of responders, partial-responders, and non-responders. The text box in each strip plot shows the number of non-responders and partial-responders whose values are below the lowest responder's values. This is to give a clear indication of the distribution of scores within each category and how they compare to each other. Additionally, the color scheme used in the plots helps differentiate between the different categories for easy visual interpretation. These strip plots provide valuable insights into the distribution of SIA scores across different responder groups, aiding in the analysis and understanding of the data.
