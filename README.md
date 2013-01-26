OmiX
====

This repo contains snippets to analyze lipidomics data.

How to use these snippets ?

To Analyze a set of lipidomes, follow these instructions :

First, place all you lipidomes in a folder. An example folder with name "Test" is given.

As you can see in the example folder there are a total of 10 files with the following names :

1. BY4741_24
2. BY4741_37
3. Elo1_24
4. Elo1_37
5. Elo2_24
6. Elo2_37
7. Elo3_24
8. Elo3_37
9. FA
10.UserSMILES

Files 1-8 are real lipidome data (Esjing et. al 2010 PNAS) included here for illustration. 9 and 10 are special files that have to included in all new folders you create.

Files 1-8 follow a tab delimited format. First column is the name of the lipid. Second column is the concentration of that lipid. If you do not have concentrations, have '0's in concentration column - do not leave that column empty.

Naming of each lipid follows a format. More than 20 classes of lipids are used in example lipidome datasets (files 1-8). Use those examples if you have to name a new lipid. Spacing and punctuation marks are important in writing lipid names.

Make sure there are no empty rows in any of the input files. Our programs do not check for empty rows/columns. In our experience, having empty rows/columns will lead to wrong results and we found it very diffucult to pin point source of wrong results when we had empty rows/columns.

File number 9, with name "FA" is used to draw fatty acid structures. In any given organism, including yeast, not all possible fatty acids are synthesized. User has to define what fatty acid possiblities have to used. Fatty Acid combinations we used for our analysis are given in the file. USer is free to add/remove some of them as per their requirement. But please note that lipidmaps structure drawing code has to modified accordingly. User may contact LipidMaps team (www.lipidmaps.org) or us to know how to edit structure drawing code. 

File number 10, with name "UserSMILES" contains User defined SMILES for lipid species. Our program uses Lipid Maps structure drawing tools to generate structures. Lipid Maps tools do not support all classes of lipids and it might happen that some of the lipids in your dataset are part of those not supported by Lipid Maps. In such cases, user can define their own SMILES. If you do not know which lipids are not supported, just run this program on your dataset(s) with default options. Lipids for which structures were not drawn are written to STDOUT. You can provide SMILES to those lipids by adding them to UserSMILES file and re-run the program.

make sure "bin" folder (provided with this package) is located in the same directory as this README file.

Last but not the least, there should be no other files in the folder. Our program treats all files in the folder (exluding files named "FA" and "UserSMILES") to be lipidome data that needs analysis.

run the program with following syntax

$ perl 130117_AnalyzeLipidomes.pl <folder_name>

e.g: $ perl 130117_AnalyzeLipidomes.pl Test

If there are any dependency issues, you will be notified through STDERR or STDOUT.

This program generates many result files, all of which will be located in the folder you have provided as input. In the example shown above, all results will be placed in the folder "Test". Taking example files and "Test" folder, we will descibe what to expect in the form of results. 

Our program generated the following (27) files with above example (excluding 10 input files)

01. ClassList
02. CentroidList.txt
03. SMILESlist.smi
04. Test.csv
05. Test.pca
06. Test_VariancePlot.pdf
07. Test_2D.svg
08. Test_3D.svg
09. BY4741_24_1.svg
10. BY4741_37_1.svg
11. Elo1_24_1.svg
12. Elo1_37_1.svg
13. Elo2_24_1.svg
14. Elo2_37_1.svg
15. Elo3_24_1.svg
16. Elo3_37_1.svg
17. WithConcMaxDist_from_PC1_PC2.csv
18. WithConcAvgDist_from_PC1_PC2.csv
19. WithOutConcMaxDist_from_PC1_PC2.csv
20. WithOutConcAvgDist_from_PC1_PC2.csv
21. WithConcMaxDist_from_SimilarityScore.csv
22. WithConcAvgDist_from_SimilarityScore.csv
23. WithOutConcMaxDist_from_SimilarityScore.csv
24. WithOutConcAvgDist_from_SimilarityScore.csv
25. Dendrograms.svg
26. Test.log
27. CalcStats

Files 1-8 contain results pertaining to all lipids combined. For the example, it means all lipid classes present in '8' lipidome datasets.

01. ClassList

This file containts list of lipid class (abbreviated name) present in input dataset. 

02. CentroidList.txt

This file contains three columns. First column is the lipid species (input). Second column is the representative isomer for that lipid. Third column is the number of isomers considered in choosing representative isomer of second column.

03. SMILESlist.smi

List of Non Canonical SMILES for the representative isomer for each lipid.

04. Test.csv

Pairwise distances between all lipid species using Levenshtein distance.

05. Test.pca

Coordinates of first three principal components for all lipid species.

06. Test_VariancePlot.pdf

Distribution of variance in the first 10 principal components.

07. Test.svg

Interactive PCA plot for first two principal compoents. Again, for all lipidomes combined.

08. Test.svg

Interactive PCA plot for first three principal compoents. Again, for all lipidomes combined.

09-16

Interactive PCA plot for first two principal compoents for each lipidomes seperately. 

17. WithConcMaxDist_from_PC1_PC2.csv

Pairwise maxiumum eucledian distance (also called Hausdroff distance in literature) for all pairs lipidomes as a comma seperated file; calculated with concentations.

18. WithConcAvgDist_from_PC1_PC2.csv

Pairwise average eucledian distance for all pairs lipidomes as a comma seperated file; calculated with concentations.

19. WithOutConcMaxDist_from_PC1_PC2.csv

Pairwise maxiumum eucledian distance for all pairs of lipidomes as a comma seperated file; calculated without concentations.

20. WithOutConcAvgDist_from_PC1_PC2.csv

Pairwise average eucledian distance for all pairs lipidomes as a comma seperated file; calculated without concentations.

21. WithConcMaxDist_from_SimilarityScore.csv

Pairwise maxiumum distance (not eucledian - distance here is (1-Similarty Score) - distance will be in aribitary units - but the range of distance value will be always 0 to 1) for all pairs lipidomes as a comma seperated file; calculated with concentations.

22. WithConcAvgDist_from_SimilarityScore.csv

Pairwise average distance (not eucledian) for all pairs lipidomes as a comma seperated file; calculated with concentations.

23. WithOutConcMaxDist_from_SimilarityScore.csv

Pairwise maxiumum distance (not eucledian - distance here is (1-Similarty Score) - distance will be in aribitary units - but the range of distance value will be always 0 to 1) for all pairs lipidomes as a comma seperated file; calculated without concentations.

24. WithOutConcAvgDist_from_SimilarityScore.csv

Pairwise average distance (not eucledian) for all pairs lipidomes as a comma seperated file; calculated without concentations.

25. Dendrograms.svg

Dendrograms generated using Pairwise Hausdroff Distance and Lipidome Overlap Scores.

26. Test.log

Log file for debugging.

27. CalcStats

Another log file for debugging purpose only.
