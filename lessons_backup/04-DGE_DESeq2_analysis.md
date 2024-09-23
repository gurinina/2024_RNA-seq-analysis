# DGE analysis workflow

Differential expression analysis with DESeq2 involves multiple steps as displayed in the flowchart below in blue. Briefly, DESeq2 will:

![Alt text](img/DESeq2_workflow_2018.png){ width=600 }

* normalize the raw counts using size factors to account for differences in library depth

* estimate the gene-wise dispersions

* shrink these estimates to generate more accurate dispersion estimates

* fit the negative binomial (NG) model and perform hypothesis testing using the Wald test.

This final step in the differential expression analysis workflow of fitting the raw counts to the NB model and performing the statistical test for differentially expressed genes, is the step we care about. This is the step that determines whether the mean expression levels of different sample groups are significantly different.


![Alt text](img/de_theory.png){ width=600 }

The [DESeq2 paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8) was published in 2014, but the package is continually updated and available for use in R through Bioconductor. 

## Running DESeq2

Prior to performing the differential expression analysis, it is a good idea to know what **sources of variation** are present in your data, either by exploration during the QC and/or prior knowledge. Once you know the major sources of variation, you can remove them prior to analysis or control for them in the statistical model by including them in your **design formula**.

### Design formula

A design formula tells the statistical software the known sources of variation to control for, as well as, the factor of interest to test for during differential expression testing. For example, if you know that sex is a significant source of variation in your data, then sex should be included in your model. The design formula should have all of the factors in your metadata that account for major sources of variation in your data. The last factor entered in the formula should be the condition of interest.

For example, suppose you have the following metadata:

![Alt text](img/meta_example.png){ width=300 }

If you want to examine the expression differences between treatments, and you know that major sources of variation include sex and age, then your design formula would be:

`design <- ~ sex + age + treatment`

The tilde (`~`) should always proceed your factors and tells DESeq2 to model the counts using the following formula. Note the **factors included in the design formula need to match the column names in the metadata**. 

***
**Exercises**

1. Suppose you wanted to study the expression differences between the two age groups in the metadata shown above, and major sources of variation were `sex` and `treatment`, how would the design formula be written?

2. Based on our **Mov10** `metadata` dataframe, which factors could we include in our design formula?

3. What would you do if you wanted to include a factor in your design formula that is not in your metadata? 

***

### MOV10 DE analysis 

Now that we know how to specify the model to DESeq2, we can run the differential expression pipeline on the **raw counts**. 

**To get our differential expression results from our raw count data, we only need to run 2 lines of code!**

First we create a DESeqDataSet as we did in the ['Count normalization'](https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html#2-create-deseq2-object) lesson and specify the location of our raw counts and metadata, and input our design formula:

```r
## Create DESeq object
dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ sampletype)
```

Then, to run the actual differential expression analysis, we use a single call to the function `DESeq()`. 

```r
## Run analysis
dds <- DESeq(dds)
```

By re-assigning the results of the function back to the same variable name (`dds`), we can fill in the `slots` of our `DESeqDataSet` object.

![deseq1](img/deseq_obj2.png)

**Everything from normalization to linear modeling was carried out by the use of a single function!** This function will print out a message for the various steps it performs: 

```
estimating size factors
estimating dispersions
gene-wise dispersion estimates
mean-dispersion relationship
final dispersion estimates
fitting model and testing
``` 

## DESeq2 differential gene expression analysis workflow

With the 2 lines of code above, we just completed the workflow for the differential gene expression analysis with DESeq2. The steps in the analysis are output below:

![Alt text](img/deseq2_workflow_separate.png){ width=200 }

We will be taking a detailed look at each of these steps to better understand how DESeq2 is performing the statistical analysis and what metrics we should examine to explore the quality of our analysis.

### Step 1: Estimate size factors

The first step in the differential expression analysis is to estimate the size factors, which is exactly what we already did to normalize the raw counts. 

![Alt text](img/deseq2_workflow_separate_sf.png){ width=200 }

DESeq2 will automatically estimate the size factors when performing the differential expression analysis. However, if you have already generated the size factors using `estimateSizeFactors()`, as we did earlier, then DESeq2 will use these values.

To normalize the count data, DESeq2 calculates size factors for each sample using the *median of ratios method* discussed previously in the Count normalization (Chapter \@ref(count-normalization)) lesson.

#### MOV10 DE analysis: examining the size factors

Let's take a quick look at size factor values we have for each sample:

```
## Check the size factors
sizeFactors(dds)

Mov10_kd_2 Mov10_kd_3 Mov10_oe_1 Mov10_oe_2 Mov10_oe_3 Irrel_kd_1 Irrel_kd_2 Irrel_kd_3 
 1.5646728  0.9351760  1.2016082  1.1205912  0.6534987  1.1224020  0.9625632  0.7477715  
```
 
Take a look at the total number of reads for each sample:

```r
## Total number of raw counts per sample
colSums(counts(dds))
```

*How do the numbers correlate with the size factor?*


### Step 2: Estimate gene-wise dispersion

The next step in the differential expression analysis is the estimation of gene-wise dispersions. Before we get into the details, we should have a good idea about what dispersion is referring to in DESeq2.

![Alt text](img/deseq2_workflow_separate_dis.png){ width=200 }

**What is dispersion?** 

Dispersion is a measure of spread or variability in the data. Variance, standard deviation, IQR, among other measures, can all be used to measure dispersion. However, DESeq2 uses a specific measure of dispersion (α) related to the mean (μ) and variance of the data: `Var = μ + α*μ^2`.  For genes with moderate to high count values, the square root of dispersion will be equal to the coefficient of variation (`σ / μ`). So 0.01 dispersion means 10% variation around the mean expected across biological replicates. 

**What does the DESeq2 dispersion represent?** 

The DESeq2 dispersion estimates are **inversely related to the mean** and **directly related to variance**. **Based on this relationship, the dispersion is higher for small mean counts and lower for large mean counts.** The dispersion estimates for genes with the same mean will differ only based on their variance. **Therefore, the dispersion estimates reflect the variance in gene expression for a given mean value.** 

The plot of mean versus variance in count data below shows the variance in gene expression increases with the mean expression (each black dot is a gene). Notice that the relationship between mean and variance is linear on the log scale, and for higher means, we could predict the variance relatively accurately given the mean. However, **for low mean counts, the variance estimates have a much larger spread; therefore, the dispersion estimates will differ much more between genes with small means**. 

![Alt text](img/deseq_mean_vs_variance.png){ width=600 }

**How does the dispersion relate to our model?** 

To accurately model sequencing counts, we need to generate accurate estimates of within-group variation (variation between replicates of the same sample group) for each gene. With only a few (3-6) replicates per group, the **estimates of variation for each gene are often unreliable** (due to the large differences in dispersion for genes with similar means). 

To address this problem, DESeq2 **shares information across genes** to generate more accurate estimates of variation based on the mean expression level of the gene using a method called 'shrinkage'. **DESeq2 assumes that genes with similar expression levels have similar dispersion.** 

**Estimating the dispersion for each gene separately:**

The first step in estimating dispersion is to get the dispersion estimates for each gene. The dispersion for each gene is estimated using maximum likelihood estimation (MLE). MLE can be thought about as **given the count values of the replicates, the most likely estimate of dispersion is calculated**.

### Step 3: Fit curve to gene-wise dispersion estimates

The next step in the workflow is to fit a curve to the dispersion estimates for each gene. This curve represents the expected value of the dispersion given the expression values of the genes. The idea behind fitting a curve to the data is that different genes will have different scales of biological variability, but, over all genes, there will be a distribution of reasonable estimates of dispersion.

![Alt text](img/deseq2_workflow_separate_fit.png){ width=200 }

This curve is displayed as a red line in the figure below, which plots the estimate for the **expected dispersion value for genes of a given expression strength**. Each black dot is a gene with an associated mean expression level and maximum likelihood estimation (MLE) of the dispersion (Step 1).

![Alt text](img/deseq_dispersion1.png){ width=400 }

### Step 4: Shrink gene-wise dispersion estimates toward the values predicted by the curve

The next step in the workflow is to shrink the gene-wise dispersion estimates toward the expected dispersion values to get the final dispersion estimates. 

![Alt text](img/deseq2_workflow_separate_shr.png){ width=200 }

The curve allows for more accurate identification of differentially expressed genes when sample sizes are small, and the strength of the shrinkage for each gene depends on :
	
- how close gene dispersions are from the curve 
- sample size (more samples = less shrinkage)

**This shrinkage method is particularly important to reduce false positives in the differential expression analysis.** Genes with low dispersion estimates are shrunken towards the curve, and the more accurate, higher shrunken values are output for fitting of the model and differential expression testing.

Dispersion estimates that are slightly above the curve are also shrunk toward the curve for better dispersion estimation; however, genes with extremely high dispersion values are not. This is due to the likelihood that the gene does not follow the modeling assumptions and has higher variability than others for biological or technical reasons and therefore is likely to be an outlier. Shrinking the values toward the curve could result in false positives, so these values are not shrunken. These dispersion outliers are shown surrounded by blue circles below. 

![Alt text](img/deseq_dispersion2.png){ width=600 }

**This is a good plot to examine to ensure your data is a good fit for the DESeq2 model.** You expect your data to generally scatter around the curve, with the dispersion decreasing with increasing mean expression levels. If you see a cloud or different shapes, then you might want to explore your data more to see if you have contamination (mitochondrial, etc.) or outlier samples.  Note how much shrinkage you get across the whole range of means in the `plotDispEsts()` plot for any experiment with low degrees of freedom.

Examples of **worrisome dispersion plots** are shown below:

![Alt text](img/bad_dispersion1.png){ width=600 }

![Alt text](img/bad_dispersion2.png){ width=600 }


#### MOV10 DE analysis: exploring the dispersion estimates and assessing model fit

Let's take a look at the dispersion estimates for our MOV10 data:

```r
## Plot dispersion estimates
plotDispEsts(dds)
```

<img src="img/plotDispersion.png">
 
**Since we have a small sample size, for many genes we see quite a bit of shrinkage. Do you think our data are a good fit for the model?**

***







