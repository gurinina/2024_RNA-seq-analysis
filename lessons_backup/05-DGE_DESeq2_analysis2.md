# Model and hypothesis testing

## Generalized Linear Model fit for each gene

The final step in the DESeq2 workflow is fitting the Negative Binomial model for each gene and performing differential expression testing.

![Alt text](img/deseq2_workflow_separate_2018.png){ width=400 }


The first step in performing differential expression analysis is to model the RNA-seq counts. DESeq2 fits a negative binomial generalized linear model (GLM) to estimatae the counts for each gene. A generalized linear model (GLM) is a statistical method that allows you to model relationships between different variables. It's similar to a linear regression model, which fits data to a line, but it can handle data that is not normally distributed, such as data that is binary or count data. The two parameters required to fit the negative binomial GLM are the size factor, and the dispersion estimate. The model incorporates the estimated dispersion and the design matrix, which specifies the experimental conditions and any covariates of interest. The negative binomial model is shown below:

![Alt text](img/NB_model_formula.png){ width=600 }

where $K_{ij}$ are the raw counts for gene i in sample j. 

After the model is fit, coefficients are estimated for each sample group along with their standard error. These coefficients are the estimates for the log2 fold changes using the following model:

![Alt text](img/NB_model_formula_betas.png){ width=600 }

An easier way of writing this is:

$$\LARGE y = \beta_0 + \beta_1x_{1} + \beta_2x_{2} + ...$$

Here, $y$ is the log2 fitted counts for each gene taking dispersion into account, $x_{i}$ is the design factor i (e.g control = 0, case = 1)  and the $\beta$ coefficients are the Log2 Fold Changes (LFC) for each sample.

In our dataset the design has three levels, and can be written like this:

$$\LARGE y = \beta_0 + \beta_1 x_{mov10kd} + \beta_2 x_{mov10oe}$$

Where:

* $y$ is the log2 of the fitted counts for each gene

* $\beta_0$ is the log2 average of the reference group (in this example “control”)

* $\beta_1$ is the log2 fold difference between “Mov10_kd” and the reference group “control”

* $\beta_2$ is the log2 fold difference between “Mov10_oe” and the reference group “control”

* $x_{mov10kd}$ indicates whether a sample belongs to group Mov10_kd or not.

* $x_{mov10oe}$ indicates whether a sample belongs to group Mov10_oe or not.

The $\beta$ coefficents are the estimates for the **log2 fold changes** for each sample group. However, **log2 fold changes** are inherently noisier when counts are low due to the large dispersion we observe with low read counts. To avoid this, the **log2 fold changes calculated by the model need to be adjusted**. 

## Shrunken log2 foldchanges (LFC)

To generate more accurate LFC estimates, DESeq2 allows for the **shrinkage of the LFC estimates toward zero** when the information for a gene is low, which could include:

- Low counts
- High dispersion values

As with the shrinkage of dispersion estimates, LFC shrinkage uses **information from all genes** to generate more accurate estimates. Specifically, the distribution of LFC estimates for all genes is used (as a prior) to shrink the LFC estimates of genes with little information or high dispersion toward more likely (lower) LFC estimates. 

![Alt text](img/deseq2_shrunken_lfc.png){ width=500 }

*Illustration taken from the [DESeq2 paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8).*

For example, in the figure above, the green gene and purple gene have the same mean values for the two sample groups (C57BL/6J and DBA/2J), but the green gene has little variation while the purple gene has high levels of variation. For the green gene with low variation, the **unshrunken LFC estimate** (vertex of the green **solid line**) is very similar to the shrunken LFC estimate (vertex of the green dotted line), but the LFC estimates for the purple gene are quite different due to the high dispersion. So even though two genes can have similar normalized count values, they can have differing degrees of LFC shrinkage. Notice the **LFC estimates are shrunken toward the prior (black solid line)**.

To generate the shrunken log2 fold change estimates, you have to run an additional step on your results object (that we will create below) with the function `lfcShrink()`.

## Statistical test for LFC estimates: Wald test

In DESeq2, the Wald test is the default used for hypothesis testing when comparing two groups. The Wald test is a test usually performed on the LFC estimates.

DESeq2 implements the Wald test by:

* Taking the LFC and dividing it by its standard error, resulting in a z-statistic
* The z-statistic is compared to a standard normal distribution, and a p-value is computed reporting the probability that a z-statistic at least as extreme as the observed value would be selected at random
* If the p-value is small we reject the null hypothesis (LFC = 0) and state that there is evidence against the null (i.e. the gene is differentially expressed).

### Overall DESeq2 workflow:

Estimate Size Factor

Estimate Dispersions

Fit negative binomial GLM

Hypothesis Testing

### Creating contrasts

To indicate to DESeq2 the two groups we want to compare, we can use **contrasts**. Contrasts are then provided to DESeq2 to perform differential expression testing using the Wald test. Contrasts can be provided to DESeq2 a couple of different ways:

1. Do nothing. Automatically DESeq2 will use the base factor level of the condition of interest as the base for statistical testing. The base level is chosen based on alphabetical order of the levels.
2. In the `results()` function you can specify the comparison of interest, and the levels to compare. The level given last is the base level for the comparison. The syntax is given below:
	
	```r
	
	# DO NOT RUN!
	contrast <- c("condition", "level_to_compare", "base_level")
	results(dds, contrast = contrast, alpha = alpha_threshold)
	
	```
The alpha threshold is the significance cutoff used for optimizing the independent filtering. After testing your genes and getting p-values, DESeq2 performs a step which is called "independent filtering" and the goal is to remove genes with very low counts because these usually have low power to be detected as significant in the first place, because of high dispersion (you remove them even without looking at their p-values). The alpha in threshold in results controls how this procedure is made, and it should be set to the threshold that you want to use in your adjusted p-values (e.g. 0.05).

#### MOV10 DE analysis: contrasts and Wald tests

We have three sample classes so we can make three possible pairwise comparisons:

1. Control vs. Mov10 overexpression
2. Control vs. Mov10 knockdown
3. Mov10 knockdown vs. Mov10 overexpression

**We are really only interested in #1 and #2 from above**. Using the design formula we provided `~ sampletype`, indicating that this is our main factor of interest.

### Building the results table

To build our results table we will use the `results()` function. To tell DESeq2 which groups we wish to compare, we supply the contrasts we would like to make using the`contrast` argument. For this example we will save the unshrunken and shrunken versions of results to separate variables. Additionally, we are including the `alpha` argument and setting it to 0.05. This is the significance cutoff used for optimizing the independent filtering (by default it is set to 0.1). If the adjusted p-value cutoff (FDR) will be a value other than 0.1 (for our final list of significant genes), `alpha` should be set to that value.

```r
## Define contrasts, extract results table, and shrink the log2 fold changes

contrast_oe <- c("sampletype", "MOV10_overexpression", "control")

res_tableOE_unshrunken <- results(dds, contrast=contrast_oe, alpha = 0.05)

resultsNames(dds)

res_tableOE <- lfcShrink(dds = dds,coef="sampletype_MOV10_overexpression_vs_control",res=res_tableOE_unshrunken)
```

**The order of the names determines the direction of fold change that is reported.** The name provided in the second element is the level that is used as baseline. So for example, if we observe a log2 fold change of -2 this would mean the gene expression is lower in Mov10_oe relative to the control. 

### MA Plot

A plot that can be useful to exploring our results is the MA plot. The MA plot shows the mean of the normalized counts versus the log2 foldchanges for all genes tested. The genes that are significantly DE are colored to be easily identified. This is also a great way to illustrate the effect of LFC shrinkage. The DESeq2 package offers a simple function to generate an MA plot. 

**Let's start with the unshrunken results:**

```r
plotMA(res_tableOE_unshrunken, ylim=c(-2,2))
```

![Alt text](img/maplot_unshrunken.png){ width=600 }

**And now the shrunken results:**

```r
plotMA(res_tableOE, ylim=c(-2,2))
```

![Alt text](img/MA_plot.png){ width=600 }

In addition to the comparison described above, this plot allows us to evaluate the magnitude of fold changes and how they are distributed relative to mean expression. Generally, we would expect to see significant genes across the full range of expression levels. 

### MOV10 DE analysis: results exploration

The results table looks very much like a dataframe and in many ways it can be treated like one (i.e when accessing/subsetting data). However, it is important to recognize that it is actually stored in a `DESeqResults` object. When we start visualizing our data, this information will be helpful. 

```r
class(res_tableOE)
```

Let's go through some of the columns in the results table to get a better idea of what we are looking at. To extract information regarding the meaning of each column we can use `mcols()`:

```r
mcols(res_tableOE, use.names=T)
```

* `baseMean`: mean of normalized counts for all samples
* `log2FoldChange`: log2 fold change
* `lfcSE`: standard error
* `stat`: Wald statistic
* `pvalue`: Wald test p-value
* `padj`: BH adjusted p-values
 

Now let's take a look at what information is stored in the results:

```r
res_tableOE %>% data.frame() %>% head()
```

```
log2 fold change (MAP): sampletype MOV10_overexpression vs control 
Wald test p-value: sampletype MOV10_overexpression vs control 
DataFrame with 6 rows and 6 columns
               baseMean log2FoldChange      lfcSE       stat    pvalue       padj
              <numeric>      <numeric>  <numeric>  <numeric> <numeric>  <numeric>
1/2-SBSRNA4  45.6520399     0.26976764 0.18775752  1.4367874 0.1507784 0.25242910
A1BG         61.0931017     0.20999700 0.17315013  1.2128030 0.2252051 0.34444163
A1BG-AS1    175.6658069    -0.05197768 0.12366259 -0.4203185 0.6742528 0.77216278
A1CF          0.2376919     0.02237286 0.04577046  0.4888056 0.6249793         NA
A2LD1        89.6179845     0.34598540 0.15901426  2.1758136 0.0295692 0.06725157
A2M           5.8600841    -0.27850841 0.18051805 -1.5428286 0.1228724 0.21489067
```

> **NOTE: on p-values set to NA**
> > 
> 1. If within a row, all samples have zero counts, the baseMean column will be zero, and the log2 fold change estimates, p-value and adjusted p-value will all be set to NA.
> 2. If a row contains a sample with an extreme count outlier then the p-value and adjusted p-value will be set to NA. These outlier counts are detected by Cook’s distance. 
> 3. If a row is filtered by automatic independent filtering, for having a low mean normalized count, then only the adjusted p-value will be set to NA. 


### Multiple test correction

Note that we have pvalues and p-adjusted values in the output. Which should we use to identify significantly differentially expressed genes?

If we used the `p-value` directly from the Wald test with a significance cut-off of p < 0.05, that means there is a 5% chance it is a false positives. Each p-value is the result of a single test (single gene). The more genes we test, the more we inflate the false positive rate. **This is the multiple testing problem.** For example, if we test 20,000 genes for differential expression, at p < 0.05 we would expect to find 1,000 genes by chance. If we found 3000 genes to be differentially expressed total, 150 of our genes are false positives. We would not want to sift through our "significant" genes to identify which ones are true positives.

DESeq2 helps reduce the number of genes tested by removing those genes unlikely to be significantly DE prior to testing, such as those with low number of counts and outlier samples (gene-level QC). However, we still need to correct for multiple testing to reduce the number of false positives, and there are a few common approaches:

- **Bonferroni:** The adjusted p-value is calculated by: p-value * m (m = total number of tests). **This is a very conservative approach with a high probability of false negatives**, so is generally not recommended.
- **FDR/Benjamini-Hochberg:** Benjamini and Hochberg (1995) defined the concept of FDR and created an algorithm to control the expected FDR below a specified level given a list of independent p-values. **An interpretation of the BH method for controlling the FDR is implemented in DESeq2 in which we rank the genes by p-value, then multiply each ranked p-value by m/rank where m is the number of tests.**

In DESeq2, the p-values attained by the Wald test are corrected for multiple testing using the Benjamini and Hochberg method by default. There are options to use other methods in the `results()` function. The p-adjusted values should be used to determine significant genes. The significant genes can be output for visualization and/or functional analysis.

### MOV10 DE analysis: Control versus Knockdown

Now that we have results for the overexpression results, let's do the same for the **Control vs. Knockdown samples**. Use contrasts in the `results()` to extract a results table and store that to a variable called `res_tableKD`.  

```r
## Define contrasts, extract results table and shrink log2 fold changes
contrast_kd <-  c("sampletype", "MOV10_knockdown", "control")

res_tableKD_unshrunken <- results(dds, contrast=contrast_kd, alpha = 0.05)

resultsNames(dds)

res_tableKD <- lfcShrink(dds, coef = 2, res=res_tableKD_unshrunken)

# We will save these results for later use in the data directory using the following command:

saveRDS(res_tableKD, file = "data/res_tableKD.rds")

# just for security, save all our objects in ".RData" again:
save.image()
```

Take a quick peek at the results table containing Wald test statistics for the Control-Knockdown comparison we are interested in and make sure that format is similar to what we observed with the OE.

## Summarizing results

To summarize the results table, a handy function in DESeq2 is `summary()`. Confusingly it has the same name as the function used to inspect data frames. This function when called with a DESeq results table as input, will summarize the results using the alpha threshold: FDR < 0.05 (padj/FDR is used even though the output says `p-value < 0.05`). Let's start with the OE vs control results:

```r
## Summarize results
summary(res_tableOE)
```

In addition to the number of genes up- and down-regulated at the default threshold, **the function also reports the number of genes that were tested (genes with non-zero total read count), and the number of genes not included in multiple test correction due to a low mean count**.


### Extracting significant differentially expressed genes

What we noticed is that the FDR threshold on it's own doesn't appear to be reducing the number of significant genes. With large significant gene lists it can be hard to extract meaningful biological relevance. To help increase stringency, one can also **add a fold change threshold**. _The `summary()` function doesn't have an argument for fold change threshold._

Let's first create variables that contain our threshold criteria:

```r
### Set thresholds
padj.cutoff <- 0.05
lfc.cutoff <- 0.58
```

The `lfc.cutoff` is set to 0.58; remember that we are working with log2 fold changes so this translates to an actual fold change of 1.5 which is pretty reasonable. 

We can easily subset the results table to only include those that are significant using the `filter()` function, but first we will convert the results table into a tibble:

```r
res_tableOE_tb <- res_tableOE %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()
```

Now we can subset that table to only keep the significant genes using our pre-defined thresholds:

```r
sigOE <- res_tableOE_tb %>%
        filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)
```

**How many genes are differentially expressed in the Overexpression compared to Control, given our criteria specified above? Does this reduce our results?**
	
Using the same thresholds as above (`padj.cutoff < 0.05` and `lfc.cutoff = 0.58`), subset `res_tableKD` to report the number of genes that are up- and down-regulated in Mov10_knockdown compared to control.

```r

res_tableKD_tb <- res_tableKD %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()
  
sigKD <- res_tableKD_tb %>%
        filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)
```

**How many genes are differentially expressed in the Knockdown compared to Control?** 

Now that we have subsetted our data, we are ready for visualization!

***