# Initial Condition Generator for Relativistic Hydro

## Authors: José Hugo Elsas, Rafael Derradi 
## Colaborators: Takeshi Kodama, Tomoi Koide

Related publications: 
   - Souza, R. D. ; Kodama, T. ; Koide, T. ; Gaspar Elsas, J. H. C. ; Bratkovskaya, E.L. ; Cassing, W. .
Evaluating the collective properties of heavy-ion collisions from microscopic dynamics using the PHSD
model. In: The International Workshop on Collectivity in Relativistic Heavy Ion Collisions, 2015,
Kolymbari. PROCEEDINGS of The International Workshop on Collectivity in Relativistic Heavy Ion
Collisions. Barcelona: Edited by L. Bravina, V.K. Magas, A. Feijoo, 2014. v. 01. p. 111-110.* 

# Contents: 

- splitandfit .cpp and .h: Domain Splitting Algorithm implementation based on near equal integration value 
                           of each subdomain for the provided density function. Alternative algorithms will be Implemented in the future.

- trial-functions .cpp .h: Catalogue of Density and Velocity Functions to be fed to the splitandfit library. 
    Already Implemented:
      Density Profiles:
        -- Gaussian with Hotspot (winicon)
        -- Woodsaxon 
        -- Woodsaxon with Hotspot
        -- Gubser Entropy 
      Velocity Profiles:
        -- Null 
        -- Gubser 

- domain-spliting .ccp: Program that uses the splitandfit library

- sph-dens .cpp: Ploting Program, to check quality of fit from splitandfit library

- example .cpp: Example program of how to use the splitandfit library, 
                tutorial of how to use the progam included on the comments

- gaus2d.cpp: simple smooth profile

- drawIniConProfile: draw 2D plot of the initial
  condition profile (needs ROOT libraries)

Obs:

- In case makefile don't work, try creating a obj folder on the home directory. It's necessary to store the .o files temporarely
