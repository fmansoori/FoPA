FoPA v0.1
-------------------------------
Before running FoPA do the following:

- Install PRISM from http://www.prismmodelchecker.org/download.php
- Add the path of the lib and bin directory of PRISM (e.g. C:\Program Files\prism-4.2.1\lib and C:\Program Files\prism-4.2.1\bin) to environment variable Path
- Add prism.jar, all files of the lib directory and classes (e.g. C:\Program Files\prism-4.2.1\lib\prism.jar, C:\Program Files\prism-4.2.1\classes, C:\Program Files\prism-4.2.1\lib\*) to environment variable CP

- In the config file in FoPA folder:
      -Change the "FoPA_path" to the path of the FoPA in your system 
      -Change the "gene_path to the path of the experiment data
---------------------------------

How To Run "FoPA":
  Below is the general command:   
                                  ./FoPA  "DEGfile" "allfile" "outputfile address"
  
  - the first argument is the name of the differentially expressed genes file (the format of the file 
  is expressed below)
  - the second argument is the name of the all genes files 
     ( the address of these two files should be indicated in the config file as "gene_path")
   - the third argument is the name and the address of the output file, the file that contains the analysis' results)
   
  Example:  ./FoPA "GSE4107_DEG.txt" "GSE4107_all.txt" "./result.txt" 

-------------------------------------------------------------------
Format of the experiment data:

The "DEGfile" is a file containing the Entrez ID of genes found to be differentially expressed between disease and normal samples, 
the "allfile" is a tab-delemited file with the Entrez ID of all genes profiled on the microarray and the moderated t-score of each gene when comparing the desease and normal samples.
Some samples of these files could be found in ./data/input_genes folder.
