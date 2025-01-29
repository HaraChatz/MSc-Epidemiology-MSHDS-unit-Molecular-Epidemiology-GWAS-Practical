- Charikleia Chatzigeorgiou
# Practical : Genome-wide association study of BMI 

## Objectives
1. Introduction to Plink
2. Clean a dataset
3. Conduct a GWAS analysis
4. Produce Manhattan, QQ and LocusZoom plots
5. Interpret GWAS results

In this practical you will run a GWAS for BMI on some simulated data with the use of Plink.

### PLINK data format

Conceptually, genetic data is stored in matrix form – e.g. rows for individuals, columns for SNPs.
In practice, this can take many different shapes, styles and conventions. We will use PLINK
format (and mostly binary plink format). You can find information about it here:

`https://www.cog-genomics.org/plink/1.9/formats#bed`

1. The .fam file - this has information about the individuals
2. The .bim file - this has information about the SNPs
3. The .bed file - this is the matrix of genotypes for all SNPs and individuals. Note that this
is not human readable.

### Exercise 1a

Open a linux command line, navigate to the data directory and use appropriate linux
commands to look at the files:

Look at the size and structure of the geno_raw.fam file

`head geno_raw.fam`
`wc -l geno_raw.fam`

### Exercise 1b

Look at the size and structure of the geno_raw.bim file and see what each column
contains:

`wc -l geno_raw.bim`
`head geno_raw.bim`

How many individuals and SNPs are in the dataset?



### Log In
Log into bluecrystal using PuTTY.

Run the following command to access a compute node:

`srun --nodes=1 --ntasks-per-node=1 --time=01:00:00 --reservation=SSCM026902 --account=SSCM026902 --pty bash -i`

### Data
You should already have most of the files for this practical cloned from GitHub to your **home scratch** directory. 
Type the following to move into the right directory:

`cd ~/scratch/genetic-epidemiology-practicals/GWAS`

Set up an alias to the complete data directory:

`datadir="/mnt/storage/private/mrcieu/training/genetic_epidemiology/GWAS/data"`

You should save your output to `/GWAS/output/`, make this (local) directory, using the `mkdir output`.
It might be easier to also set up an alias to your local outcome directory (modify path accordingly, $HOME will point to your home directory)

`outdir="$HOME/scratch/genetic-epidemiology-practicals/GWAS/output"`

Scripts (files containing commands) should be saved in `/GWAS/scripts/`.

If you get really stuck, example scripts and ready-made output are available in `/GWAS/results/` (no peaking unless you have to!).

### Loading Plink
We will load Plink using the following code (already included in the script provided):  

`module add apps/plink/2.00`

### Exercise 2 - GWAS
First, we will run a GWAS of BMI dataset.

Using the Linux commands from above section:

Take a look at the genetic data `${datadir}/geno_clean.bim` & `.fam` using the "less" and "head" commands.

Similarly, take a look at the phenotype data `${datadir}/phen.txt`.

**_Questions:_**
> (1) How many individuals are there?

> (2) How many SNPs are there?


### Exercise 3 - Cleaning the GWAS data
There are many steps to a good QC procedure (see Weale 2010. Quality control for genome-wide association studies. Methods in Molecular Biology 628:341-372). Here, we assume that related individuals and non-white Europeans have already been removed to proceed with the final steps of the QC.

**_Question:_**
> (4) For each of the parameters below, which SNPs or individuals do you think should be removed?

>     Minor allele frequency?
>     Hardy Weinberg p-value?
>     SNP missingness?
>     Individual missingness?

Navigate to the `scripts` directory. Now look at the `qc.sh` script to see the suggested exclusions for this dataset. 

**Run the `qc.sh` script to generate new ‘cleaned’ data files `geno_qc.bed` `.bim` `.fam` (again, you may need to run the `chmod +x qc.sh` script to make sure the script is executable).**


### Exercise 4 - Running the GWAS

We will now run a ‘clean’ GWAS for BMI. This requires you to edit some scripts in order to apply Qc filters.

To do this, copy the `unclean_gwas.sh` to a new script using `cp unclean_gwas.sh {newfilename}`.

Edit the file (using nano) to use the new ‘clean’ GWAS data (you will need to be careful here, as the clean GWAS data should be in your local data directory not the ${datadir} we set up). Then, modify the covar-name option so it includes all principle component covariates in the model (as well as age and sex) and output to a new file named `clean`.

**Now you are ready to run the clean GWAS `./clean_gwas.sh`**

Ypu can also make new QQ and Manhattan plots.

Go back to `./scripts`

Copy the `unclean_gwas_graphs.sh` file to a new file named `clean_gwas_graphs.sh`.

Edit (using nano) to use the new clean results and output to new files. Run this new script with the `./clean_gwas_graphs.sh` script.

Open WINSCP (in windows) and open the graphs you’ve just created.

For the rest of the questions, go to the `results/` directory where you can see the pre-made results file `clean16.BMI.glm.linear.add`, this contains the clean gwas results for chromosome 16 only.

`cd ~/scratch/genetic-epidemiology-practicals/GWAS/results`

**_Questions (use clean **gwas results for chromosome 16 only provided in the results folder**) :_**

> (5) How many individuals in are the final analysis?

`head clean16.BMI.glm.linear.add`

> (6) How many genome-wide significant (p<5x10-8) signals do you have?

`awk '$12<0.00000005' clean16.BMI.glm.linear.add | wc -l`

`awk '{if(NR==1 || $12<0.00000005) print $0}' clean16.BMI.glm.linear.add`

> (7) Are these likely to all be independent?

`grep -v NA clean16.BMI.glm.linear.add | sort -g -k 12 | head`

> (8) What is the top signal?

> (9) How might we go about confirming this finding?



### Exercise 5 -Optional - Generating a LocusZoom Plot
We will now generate a LocusZoom plot of the signal on chromosome 16. You will find the file `clean16.BMI.glm.linear.add` in the `results/` directory.
This file was generated by extracting only the chr16 results from `clean16.BMI.glm.linear.add`.

Go to the locuszoom website (`http://locuszoom.org/`) and click 'Single plot - Your Data'.

Upload the regional association results already prepared for you, `clean16.BMI.glm.linear.add`.

Select `Plink data` format and specify the `rs8050136` SNP.

Check the genome build - is the default correct for your data?

Now, create your plot!

Save this plot to your directory, as you will need to look at this again during the imputation session.

**_Questions:_**
> (10) Are all associated SNPs in high LD?

> (11) What genes are nearby?

