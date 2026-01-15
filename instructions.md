# Application Instructions

This guide provides instructions on how to use the different tools available in this Shiny application. Each section below details the functionality of a tool, the required data format, and steps for usage.

---

## 1. GWAS Phenotype Statistical Analysis

### Functionality
This tool performs statistical analysis on Genome-Wide Association Study (GWAS) phenotype data. It allows you to compare different genotypes for various phenotypic traits and visualize the results as violin plots. Supported statistical methods include T-tests, ANOVA with Tukey HSD, and Kruskal-Wallis with Dunn's Test.

### Usage
1.  **Upload Data**: Click "Upload TSV File" and select one or more `.tsv` files.
2.  **Select File**: If you uploaded multiple files, choose the one you want to analyze from the "Select File" dropdown.
3.  **Select Phenotype**: Choose the phenotype variable you wish to analyze from the dropdown menu.
4.  **Choose Test Method**: Select the desired statistical test. If you choose "T-test", you can select specific pairs for comparison.
5.  **Customize Plot**: Adjust options for significance markers, colors, and plot format (PNG/PDF).
6.  **Analyze**: Click "Analyze Data" to generate the results.
7.  **View & Download**: The results will be displayed in the tabs on the right. You can download the statistical tables and the violin plot using the download buttons.

### Data Format
The input file must be in TSV (Tab-Separated Values) format with the following columns:
- `Genotype`: A column with the genotype information.
- Phenotype Columns: At least one numeric column representing a phenotype (e.g., `Height`, `Weight`).

**Example `gwas_sample.tsv`:**
```tsv
Sample	Genotype	Individual	Height	Weight	LeafArea	RootLength	SeedCount	Yield
S01	A/A	Ind01	25.1	10.5	15.2	5.1	300	50.2
S02	A/G	Ind02	28.3	11.2	16.8	5.5	320	55.1
S03	G/G	Ind03	22.0	9.8	14.5	4.8	280	48.9
...
```

---

## 2. Gene Expression Visualization

### Functionality
This tool visualizes gene expression data, comparing expression levels between two groups. It can generate violin plots, box plots, or a combination of both.

### Usage
1.  **Upload Data**: Upload a `.csv` file containing gene expression data. If no file is uploaded, a default sample dataset is used.
2.  **Enter Genes**: Type the names of the genes you want to visualize, separated by commas or spaces.
3.  **Define Groups**: Specify the prefixes that identify the columns for each of the two groups (e.g., `HF2` and `WJ`).
4.  **Customize Plot**: Select the statistical test, Y-axis label, data transformation method, plot type, and colors.
5.  **Generate Plot**: Click "Generate Plot" to create the visualization.
6.  **Download**: Use the download button to save the plot as a PNG or PDF.

### Data Format
The input file must be a CSV (Comma-Separated Values) file. It must contain:
- `gene_symbol`: A column with gene names.
- Sample Columns: Numeric columns for each sample, with names that start with the prefixes you define for your groups.

**Example `expression_sample.csv`:**
```csv
gene_symbol,HF2_1,HF2_2,HF2_3,WJ_1,WJ_2,WJ_3
GeneA,10.1,10.5,10.3,12.2,12.5,12.4
GeneB,15.2,15.4,15.3,14.1,14.3,14.2
GeneC,8.5,8.7,8.6,10.1,10.3,10.2
...
```

---

## 3. Phenotype Data Visualization

### Functionality
This tool is for visualizing general phenotype data. It creates boxplots to compare a categorical variable against one or more numeric phenotype variables and supports various statistical tests.

### Usage
1.  **Upload Data**: Upload a `.csv` file with your phenotype data.
2.  **Select Variables**: Choose the categorical variable for the X-axis and one or more numeric phenotype variables for the Y-axis.
3.  **Choose Test**: Select a statistical test (t-test, ANOVA, or Kruskal-Wallis).
4.  **Customize Plot**: Customize the plot by showing/hiding points, adjusting colors, and setting font sizes.
5.  **Layout Settings**: If you are plotting multiple phenotypes, you can specify the number of rows and columns for the layout.
6.  **Download**: Save the final plot as a PNG or PDF.

### Data Format
The input file must be a CSV file containing at least one categorical column and one numeric column.

**Example `phenotype_sample.csv`:**
```csv
Treatment,PlantHeight,Yield,LeafWidth
Control,30.5,100.2,5.5
TreatmentA,25.2,120.5,5.1
TreatmentB,28.9,110.8,5.9
...
```

---

## 4. LocusZoom Plotting

### Functionality
This tool generates LocusZoom-style plots, which are used to visualize association data from GWAS results in the context of local gene information.

### Usage
1.  **Upload GWAS File**: Upload a `.csv` file containing GWAS results (SNP, chromosome, position, P-value, R2).
2.  **Upload Gene File**: Upload a `.csv` file with gene track information (chromosome, start, stop, gene name).
3.  **Enter Lead SNP**: Provide the ID of the lead SNP to be highlighted in the plot.
4.  **Generate Plot**: Click "Generate Plot" to create the LocusZoom plot.
5.  **Download**: Download the plot as a PNG or PDF.

### Data Formats

**GWAS File (`locuszoom_gwas_sample.csv`)**:
A CSV file with columns:
- `SNP`: SNP identifier.
- `chr`: Chromosome.
- `pos`: Position on the chromosome.
- `P`: P-value of association.
- `R2`: Linkage disequilibrium (R-squared) value with the lead SNP.

**Example:**
```csv
SNP,chr,pos,P,R2
rs1,1,1000,0.01,0.2
rs2,1,2000,0.005,0.85
rs3,1,3000,0.05,0.1
...
```

**Gene Track File (`locuszoom_genes_sample.csv`)**:
A CSV file with columns:
- `chr`: Chromosome.
- `start`: Gene start position.
- `stop`: Gene end position.
- `name`: Gene name.

**Example:**
```csv
chr,start,stop,name
1,1500,2500,GeneX
1,4000,5000,GeneY
...
```
