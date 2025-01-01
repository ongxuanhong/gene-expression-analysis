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