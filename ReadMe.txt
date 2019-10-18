FoPA v0.1
-------------------------------
Before running FoPA do the following:

- Install PRISM from http://www.prismmodelchecker.org/download.php
<<<<<<< HEAD
- Add the path of the PRISM in your system (e.g. C:\Program Files\prism-4.5) to environment variable PRISM_DIR
- Add the "%PRISM_DIR%\lib;%PRISM_DIR%\bin" to the environment variable Path
- Add "%PRISM_DIR%\lib\prism.jar;%PRISM_DIR%\classes;%PRISM_DIR;%PRISM_DIR\lib\*; to the environment variable CP
- Copy the 'myprism' file to the [PRISM directory]\bin ('myprism' file is existed in the [FoPA directory]\PRISM folder)

- In the config file in FoPA folder change the paths to reflect the locations in your system:
      -Change the "FoPA_path" to the path of the FoPA in your system 
      -Change the "gene_path to the path of the experiment data
	  
=======
- Add the path of the lib and bin directory of PRISM (e.g. C:\Program Files\prism-4.2.1\lib and C:\Program Files\prism-4.2.1\bin) to environment variable Path
- Add prism.jar, all files of the lib directory and classes (e.g. C:\Program Files\prism-4.2.1\lib\prism.jar, C:\Program Files\prism-4.2.1\classes, C:\Program Files\prism-4.2.1\lib\*) to environment variable CP

- In the config file in FoPA folder:
      -Change the "FoPA_path" to the path of the FoPA in your system 
      -Change the "gene_path to the path of the experiment data
>>>>>>> 8870f51a2d9f16029e04c9f98d679aeaec34197c
---------------------------------

How To Run "FoPA":
  Below is the general command:   
<<<<<<< HEAD
                                 > FoPA -[n/c] "DEGfile" "allfile" "outputfile"
                             or
                                 >> python.exe ".\FoPa.py" -[n/c] "DEGfile" "allfile"  "statefile"
=======
                                 >> ./build/FoPA/FoPA  "DEGfile" "allfile" "outputfile"
                             or
                                 >> python.exe ".\FoPa.py"  "DEGfile" "allfile"  "outputfile"
>>>>>>> 8870f51a2d9f16029e04c9f98d679aeaec34197c
  
  - the first argument is the name of the differentially expressed genes file (the format of the file 
  is expressed below)
  - the second argument is the name of the all genes files 
     ( the address of these two files should be indicated in the config file as "gene_path")
<<<<<<< HEAD
   - the third argument is the name and the address of the output file (This file contains the results)
  
options:  
   - -n is used when you want a new run.
   - -c is used when you want to start from the last previous state.
   
  Example:  >FoPA -n "GSE4107_DEG.txt" "GSE4107_all.txt" "./result.txt" 
=======
   - the third argument is the name and the address of the output file, the file that contains the analysis' results)
   
  Example:  ./FoPA "GSE4107_DEG.txt" "GSE4107_all.txt" "./result.txt" 
>>>>>>> 8870f51a2d9f16029e04c9f98d679aeaec34197c

-------------------------------------------------------------------
Format of the experiment data:

The "DEGfile" is a file containing the Entrez ID of genes found to be differentially expressed between disease and normal samples, 
the "allfile" is a tab-delemited file with the Entrez ID of all genes profiled on the microarray and the moderated t-score of each gene when comparing the desease and normal samples.
Some samples of these files could be found in ./data/input_genes folder.
