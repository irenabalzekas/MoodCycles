# MoodCycles

*reference manuscript, dates, authors, and copywrites

Matlab code accompanying Manuscript: “A feasibility study of multiday mood cycles in three people with temporal lobe epilepsy undergoing ambulatory multimodal neuromonitoring and deep brain stimulation”
Author: Irena Balzekas, 2025. Please contact author with questions

Dependencies + functions
- BPWP package https://github.com/irenabalzekas/BPWP
- CVX matlab package http://cvxr.com/cvx/
- Matlab wavelet toolbox https://www.mathworks.com/products/wavelet.html
- Matlab circstat toolbox https://www.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics
- Matlab Watsons U2 https://www.mathworks.com/matlabcentral/fileexchange/43543-pierremegevand-watsons_u2

Main script and additional functions
 - MoodCyclesSeizureCyclesForPublication * start here - this works through simulated a case and generates figures and statistics as seen in the manuscript
 - simsignal() 
 - firls_nick()
 - MoodCycles_getEventPhase_Nobuiltinfilter()
 - MoodCycles_getHighLowEvents()
 - MoodSpikes_FilterFilip
 - MoodCycles_getHighLowEvents_withnoisfloor_NObuiltinfilter

Analyses were run on a computer with Intel® Core™ i5-9500 CPU @ 3.00GHz with 64 GB RAM, operating a 64-bit operating system, x64-based processor, version 22H2. 
