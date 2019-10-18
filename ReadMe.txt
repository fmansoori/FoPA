FoPA v0.1
-------------------------------
Before running FoPA do the following:

- Install PRISM from http://www.prismmodelchecker.org/download.php
- Add the path of the PRISM in your system (e.g. C:\Program Files\prism-4.5) to environment variable PRISM_DIR
- Add the "%PRISM_DIR%\lib;%PRISM_DIR%\bin" to the environment variable Path
- Add "%PRISM_DIR%\lib\prism.jar;%PRISM_DIR%\classes;%PRISM_DIR;%PRISM_DIR\lib\*; to the environment variable CP
- Copy the 'myprism' file to the [PRISM directory]\bin ('myprism' file is existed in the [FoPA directory]\PRISM folder)

- In the config file in FoPA folder change the paths to reflect the locations in your system:
      -Change the "FoPA_path" to the path of the FoPA in your system 
      -Change the "gene_path to the path of the experiment data
---------------------------------

How To Run "FoPA":
  Below is the general command:   

                                 > FoPA -[n/c] "DEGfile" "allfile" "outputfile"
                             or
                                 >> python.exe ".\FoPa.py" -[n/c] "DEGfile" "allfile"  "statefile"

  
  - the first argument is the name of the differentially expressed genes file (the format of the file 
  is expressed below)
  - the second argument is the name of the all genes files 
     ( the address of these two files should be indicated in the config file as "gene_path")
  - the third argument is the name and the address of the output file (This file contains the results)
  
options:  
   - -n is used when you want a new run.
   - -c is used when you want to start from the last previous state.
   
  Example:  >FoPA -n "GSE4107_DEG.txt" "GSE4107_all.txt" "./result.txt" 


-------------------------------------------------------------------
Format of the experiment data:

The "DEGfile" is a file containing the Entrez ID of genes found to be differentially expressed between disease and normal samples, 
the "allfile" is a tab-delemited file with the Entrez ID of all genes profiled on the microarray and the moderated t-score of each gene when comparing the desease and normal samples.
Some samples of these files could be found in ./data/input_genes folder.
