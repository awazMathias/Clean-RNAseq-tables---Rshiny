# CleanTableRNAseq

Allows table cleaning containing RNAseq data (number of counts per genes).
It takes as input files in txt, csv or xls format. The output format is either txt or csv.
It has two reactive outputs on the right: 
	- A table showing the merging of all the loaded files into a single table. The lines correspond to the genes
and the columns at the conditions (takes the name of the uploaded files) 
	- A boxplot allowing to observe the distribution of the number of genes according to the conditions (at log10)
Finally, on the left side, the upload, the download and all the options are present. 

## SOFTWARE AND PACKAGES

R version 4 and packages : 

    -shiny

    -shinydashboard

    -readxl

    -data.table

	-tools

	-ggplot2

	-reshape2

	-plotly
    

For install in R4: 

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("rtracklayer")

install.packages(c("shiny", "shinydashboard", "readxl", "data.table", "tools", "ggplot2", "reshape2", "plotly"))
```


## WIDGETS ##

1 - Uploader of tables. Tables must be uploaded one by one, the format may differ but must be
.txt, .csv or .xls. If the number of genes differs between the files, a warning message is displayed at the bottom
on the left of the screen and gives the difference. WARNING: Additional genes will be automatically deleted. 

2 - This widget is divided into two parts, the first is to give a value. This value corresponds
at minimum number of counts. In fact, if you leave 0, if a gene has a count of 0 under one of the conditions,
it will be deleted and will no longer appear in the table.
The deletion will be done after pressing the button "Removes lines < value".
After cleaning, a text is displayed under the button giving the number of lines / genes deleted.

3 - This widget is divided into two parts, the first is to give a threshold. This threshold corresponds to the minimum average
of the line. We can translate "average of the line" by the average of the counts under all the conditions for a given gene.
After pressing the "Clean by specific threshold" button, all lines with an average below the threshold will be deleted.
After cleaning, a text is displayed under the button giving the number of lines / genes deleted.

4 - Ce widget permet à l'utilisateur de définir le format de fichier de son (ses) tableau (s). Ces formats sont txt ou csv.

## PARAMETER - This part is divided into two sub-parts : 

5 - This first part allows you to define the character that will serve as a separator in your table (s)

6 - This second part allows to define if the names of your genes and the values ​​of the table are included, or not,
between single or double quotes

7 - This last widget corresponds to the downloader, two possibilities : 
	* Download the complete table which is displayed at the top right of your screen. So with all the counting
by condition, your genes
	* Download a file by condition. So with a file created by columns of the table displayed at the top right of your screen
and therefore, one file per condition. 
