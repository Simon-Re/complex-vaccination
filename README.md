# complex-vaccination: Complex vaccination strategies prevent the emergence of vaccine resistance

This repository contains the relevant code used to simulate the results in https://pure.iiasa.ac.at/id/eprint/18732/1/WP-23-005.pdf.
10 Files can be found here. **pyrasaic_tools.py**, which contains most of the functionality and the main class. Other files, that carry "Figure" in their title, **Figure5.py**, **Figure6.py**, **Figure8.py**, **Figure9.py**, **SFigure1.py**, **SFigure2.py**, **SFigure3.py**, **Figure_Antibodies.py** and **Figure_Vacctypes.py**, contain utilities for recreating the main figures in the publication. They can also be adapted to perform calculations with different relevent epidemiological parameters, than those provided in the example.

## pyrasaic_tools.py

Contains several functions, which calculate the transmissibility, herd immunity threshold, probability of emergence and within host evolution for different pathogen variants. For details see https://pure.iiasa.ac.at/id/eprint/18732/1/WP-23-005.pdf.

## Figure5/6/8/9.py 

Reproduces Figures numbered 5, 6, 8 and 9 in the publication. Slight modifcations of these programs will allow reproducing Figure 4, 7, 10 and 11 aswell. 

## SFigure1/2/3.py

Reproduces Supplementary Figures. 

## Figure_Antibodies.py and Figure_Vacctypes.py

Can be used to create schematic figures of pathogen epitopes being detected by antibodies and schematics that show the immune status in a population after the vaccination campaign.
