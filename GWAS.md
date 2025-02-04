- Written by Lavinia Paternoster; adapted by Charikleia Chatzigeorgiou (fb23004@bristol.ac.uk)
# Practical : Genome-wide association study of BMI 

## Objectives
1. Get familiar working with Plink
2. Clean a genetic dataset
3. Conduct a GWAS analysis
4. Produce Manhattan, QQ and LocusZoom plots
5. Interpret GWAS results

In this practical you will run a GWAS for BMI on some simulated data with the use of Plink.

## Set up files and directories

### Log In
Log into bluecrystal using PuTTY (for windows users) or with Terminal (for Mac users).

Run the following command to access a compute node:

```bash
srun --account=SSCM034564  --partition=teach_cpu  --time=3:00:00  --pty bash -i
```

### Directory set-up
Copy the GWAS directory to your home space on bluecrystal

 ```bash
 cp -R /mnt/storage/private/mrcieu/training/mol_epi/gwas ~
```

Check this has worked: 

```bash
tree gwas
```


### Set up an alias for the data and scripts directories:

To avoid us having to keep typing full paths, we will set up shorter alias names for the key directories.

```bash
datadir="/mnt/storage/private/mrcieu/training/genetic_epidemiology/GWAS/data"
scriptsdir="$HOME/gwas/scripts"
outdir="$HOME/gwas/output"
```



## Introduction to PLINK
PLINK is a widely used open-source software tool for whole-genome association studies (GWAS) and other genetic analyses. The software is primarily used for managing, analysing, and manipulating genotype and phenotype data, as well as performing various statistical analyses in genetics. 

`https://www.cog-genomics.org/plink/2.0/`

### PLINK data format

Conceptually, genetic data is stored in matrix form – e.g. rows for individuals, columns for SNPs.
In practice, this can take many different shapes, styles and conventions. We will use PLINK
format (and mostly binary plink format). You can find information about it here:

`https://www.cog-genomics.org/plink/2.0/formats#bed`

1. The .fam file - this has information about the individuals
2. The .bim file - this has information about the SNPs
3. The .bed file - this is the matrix of genotypes for all SNPs and individuals. Note that this
is not human readable.


### Loading Plink
We will load Plink using the following code (already included in the script provided):  

```bash
module add apps/plink2
```


## Exercise 1a

Open a linux command line, navigate to the data directory and use appropriate linux
commands to look at the files:

Look at the structure and size of the geno_raw.fam file


```bash
head ${datadir}/geno_unclean.fam
```

You should see that the .fam contains 6 columns:

1. Family ID
2. Individual ID
3. Father ID
4. Mother ID
5. Sex (1=male, 2=female, 0=unknown)
6. Phenotype (-9=missing)

```bash
wc -l ${datadir}/geno_unclean.fam
```

You should see that the .fam file contains 8237 lines that correspond to 8237 individuals

## Exercise 1b
Look at the structure and the size of the geno_raw.bim file and see what each column
contains:

```bash
head ${datadir}/geno_unclean.bim
```

The .bim file also contains 6 columns:

1. Chromosome
2. SNP ID
3. Genetic position
4. Physical position
5. Allele 1 (normally the minor allele)
6. Allele 2

```bash
wc -l ${datadir}/geno_unclean.bim
```

You should see that the .bim file contains 463,080 lines that correspond to 463,080 SNPs


## Exercise 2 - Phenotypic information of the particupants

Using the Linux commands from above section:

Take a look at the phenotype data:

```bash
head ${datadir}/phen.txt
```

The phen.txt file contains 7 columns:

1. Family ID
2. Individual ID
3. BMI - Body Mass Index
4. DBP - Diastolic Blood Pressure
5. SBP - Systolic Blood Pressure
6. CRP - C-reactive protein
7. HT - Hypertension 

phen.txt contains phenotypic information for 5 outcomes 


## Exercise 3 - Quality Control of data
There are many steps to a good QC procedure (see Weale 2010. Quality control for genome-wide association studies. Methods in Molecular Biology 628:341-372). Here, we assume that related individuals and non-white Europeans have already been removed to proceed with the final steps of the QC.

We want to filter the data to the following parameters:
1. SNPs with more than 5% missing values removed
2. Individuals with more than 5% missing values removed
3. SNPs with allele frequency < 0.01 removed
4. SNPs with Hardy Weinberg disequilibrium p value < 1e-6 removed


Look at the `qc.sh` script to see the suggested exclusions for this dataset. 

```bash
head -20 ${scriptsdir}/qc.sh
```

**Run the `qc.sh` script to generate new ‘cleaned’ data files `geno_qc.bed` `.bim` `.fam`**

You can run this command in two different ways, either copy and paste the code below: 

```bash
plink2 \
        --bfile ${datadir}/geno_unclean \
        --maf 0.01 \
        --hwe 1e-6 \
        --geno 0.05 \
        --mind 0.05 \
        --make-bed \
        --out ${outdir}/geno_qc
```

OR the more common way to do this is to run the existing script. Make sure the script is executable and then run it using these commands:

 ```bash
chmod +x ${scriptsdir}/qc.sh
${scriptsdir}/qc.sh
```

All the details for every analysis that is conducted with the use of Plink can be found at the corresponding .log file. Take a look at this file

```bash
head -45 ${outdir}/geno_qc.log
```

**_Question:_**

> For each of the parameters below, What numbers have been removed from the data?


>     Minor allele frequency?
>     Hardy Weinberg p-value?
>     SNP missingness?
>     Individual missingness?



## Exercise 4 - Running the GWAS

Now you are ready to run the clean GWAS this is saved in this script file:`${scriptsdir}/clean_gwas.sh`, you can run this script or copy the script from below:


```bash
plink2 \
        --bfile ${datadir}/geno_qc \
        --chr 1-22 \
        --linear \
        --pheno ${datadir}/phen.txt \
        --pheno-name BMI \
        --covar ${datadir}/covs.txt \
        --covar-name age sex, PC1-PC5 \
        --covar-variance-standardize \
        --ci 0.95 \
        --out ${outdir}/bmi       

```

```bash
chmod +x ${scriptsdir}/clean_gwas.sh
${scriptsdir}/clean_gwas.sh
```

We'll come back and look at the output in part 2.

## End of GWAS practical part 1


# GWAS practical and interpretation part 2

For this part of the practical we will use the results from analysis that we did at part 1 and we will use R in order to plot and further explore the output of this analysis.

## Exercise 5 -Make QQ and Manhattan plots.

In order to generate those plots we will use "QCGWAS" package in R.

Navigate to the `outputs` directory.

```bash
cd ${outdir}
```


```bash
export R_LIBS="~/R_libs"
mkdir ~/R_libs
module add languages/R/4.1.2
Rscript ${scriptsdir}/gwas_graphs.R /mnt/storage/private/mrcieu/training/genetic_epidemiology/GWAS/output/bmi.assoc.linear.add ${outdir}/bmi_clean
```

## Tranfering data from bluecrystal to your local disk in order to view plots 
### Transferring files to a windows machine
Open up WinSCP (from the Start menu), and connect using the same credentials as you have used in PuTTY. Once connected you should be able to navigate to the folder ${outdir}. Local files are shown on the left hand panel. BC4 files are shown in the right hand panel. You should be able to find your files and images that were generated during Practical part 1. Try to download those by dragging it to the local folder location you want to move it to.
> Open WINSCP (in windows) and open the graphs you’ve created.

### Transferring files to a mac machine
Go to Applications and open Filezilla. Connect to bluecrystal using the same credentials as to login and set port to 22. On the left site select your local site where you want to copy the file to your local folder. On the right site, navigate to the folder ${outdir}. At the bottom right part you should be able to find your files and images. Try to download those.


For the rest of the questions, go to the `results/` directory where you can see the pre-made results file `bmi_clean16.BMI.glm.linear.add`, this contains the clean gwas results for chromosome 16 only.

```bash
cd /mnt/storage/private/mrcieu/training/mol_epi/gwas/results
```

**_Questions (use clean **gwas results for chromosome 16 only provided in the results folder**) :_**

> How many individuals in are the final analysis?

```bash
head bmi_clean16.BMI.glm.linear.add
```

> How many genome-wide significant (p<5x10-8) signals do you have?

```bash
awk '$12<0.00000005' bmi_clean16.BMI.glm.linear.add | wc -l
```

```bash
awk '{if(NR==1 || $12<0.00000005) print $0}' bmi_clean16.BMI.glm.linear.add
```

> Are these likely to all be independent?

```bash
grep -v NA bmi_clean16.BMI.glm.linear.add | sort -g -k 12 | head
```

> What is the top signal?

> How might we go about confirming this finding?

> Go to OpenTargets website `https://genetics.opentargets.org/`, enter this SNP name and see what you can find out about this variant



## Exercise 6 -Optional - Generating a LocusZoom Plot
We will now generate a LocusZoom plot of the signal on chromosome 16. You will find the file `bmi_clean16.BMI.glm.linear.add` in your `results/` directory. 

This file was generated by extracting only the chr16 results from `bmi_clean16.BMI.glm.linear.add`. You will need to transfer this file to your local computer in order to upload to the locuszoom website.

Go to the locuszoom website (`http://locuszoom.org/`) and click 'Single plot - Your Data'.

Upload the regional association results already prepared for you, `bmi_clean16.BMI.glm.linear.add`.

Select `Plink data` format and specify the `rs8050136` SNP.

Specify the "Marker Column Name" as `ID`

Check the genome build - is the default correct for your data?

Now, create your plot!

Save this plot to your directory, as you will need to look at this again during the imputation session.

**_Questions:_**
> Are all associated SNPs in high LD?

> What genes are nearby?

