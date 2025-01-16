# gene-expression-analysis

## Download data (BAM file)
```bash
gdown --folder https://drive.google.com/drive/folders/1YN4FEJdWyKzPQziyFFCE7iyipXsVERwp?usp=drive_link
```

## Install requirements
```bash
pip install -r requirements.txt
conda install -q bioconda::subread -y
```

## Install R, VSCode and Jupyter
1. Prerequisites Installation:
   - Install R from the official R website (https://www.r-project.org)
   - Install VS Code if you haven't already
   - Install Python (required for Jupyter)

2. VS Code Extensions:
   - Install "Jupyter" extension from Microsoft
   - Install "R" extension for VS Code

3. Install Required R Packages:
Open R console or RStudio and run these commands:
```R
install.packages("IRkernel")
IRkernel::installspec()
install.packages("languageserver")
```

4. Create and Use R Jupyter Notebook:
   - Open VS Code
   - Click File > New File
   - Select "Jupyter Notebook"
   - When prompted to select a kernel, choose "R"

5. Basic Usage:
```R
# This is a code cell - you can run R code here
x <- c(1, 2, 3, 4, 5)
mean(x)

# Create a simple plot
plot(x, type = "l")

# Install and load packages
install.packages("ggplot2")
library(ggplot2)
```

6. Useful Keyboard Shortcuts:
   - Shift + Enter: Run current cell and move to next
   - Ctrl + Enter: Run current cell
   - B: Insert cell below
   - A: Insert cell above
   - DD: Delete current cell

7. Working with Data:
```R
# Read CSV file
data <- read.csv("your_file.csv")

# Basic data manipulation
head(data)
summary(data)

# Create visualizations
ggplot(data, aes(x = column1, y = column2)) +
  geom_point()
```

## Nextflow (nf-core)

**Install requirements**
```bash
conda install bioconda::nextflow -y
conda install nf-core -y
```


**Run example RNA fusion**
```bash
nextflow run nf-core/rnafusion \
   -profile docker,test \
   --outdir ./rnafusion\
   --build_references \
   -stub

nextflow run nf-core/rnafusion \
   -profile docker,test \
   --outdir ./rnafusion_analysis \
   -stub

   
nextflow run nf-core/rnafusion \
    --input nf-core/test-data/samplesheet_valid.csv \
    --outdir results \
    --genome GRCh38 \
    --tools arriba \
    --max_memory '10.GB' \
    --max_cpus 4 \
    -profile docker \
    -resume   

nextflow run nf-core/rnafusion \
    --input nf-core/test-data/samplesheet_valid.csv \
    --outdir results \
    --genome GRCh38 \
    --tools arriba \
    --max_memory '30.GB' \
    -profile singularity    
```   

Container-based
```bash

# Build container
docker build -t nextflow-container .

# Run container
docker run -it --rm \
    -v $(pwd):/workspace \
    nextflow-container bash

# Check Nextflow version
nextflow -version

# Run test pipeline
nextflow run hello

# Mount additional volumes if needed
docker run -it --rm \
    -v $(pwd):/workspace \
    -v /path/to/data:/data \
    nextflow-container bash

# Run pipeline
nextflow run nf-core/rnafusion \
    --input samples.csv \
    --outdir results
```    

## How to export mermaid graph
```bash
# Using mermaid-cli Docker image
docker run --rm -v $(pwd):/data minlag/mermaid-cli -i data/cancer-analysis.mermaid -o data/cancer-analysis.png

docker run --rm -v $(pwd):/data minlag/mermaid-cli \
    -i data/cancer-analysis.mermaid \
    -o data/cancer-analysis.svg


docker run --rm -v $(pwd):/data minlag/mermaid-cli \
    -i data/cancer-analysis.mermaid \
    -o data/cancer-analysis.png \
    -b transparent \
    -w 3840 \
    -H 2160 \
    -s 2 \
    -q 1.0    

docker run --rm -v $(pwd):/data minlag/mermaid-cli \
    -i data/cancer-analysis.mermaid \
    -o data/cancer-analysis.pdf \
    -w 4096
```