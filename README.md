This repository contains a fortran 90 code (miRNA_gillespie_single_history.f90), regarding the Gillespie simulation (Stochastic simulation) for time evolution of 'miRNA-mediated negative feedback loop' in the manuscript, "Effects of microRNA-mediated negative feedback on gene expression noise", Raunak Adhikary1, Arnab Roy1, Mohit Kumar Jolly2*, Dipjyoti Das1*

1 Department of Biological Sciences, Indian Institute of Science Education and Research Kolkata,Mohanpur, Nadia, West Bengal, 741246, India.
2 Centre for BioSystems Science and Engineering, Indian Institute of Science, Bengaluru 560012, India.

Email: ra19rs093@iiserkol.ac.in

*Corresponding Authors: dipjyoti.das@iiserkol.ac.in, mkjolly@iisc.ac.in

Instruction to get the code: Click on the 'Code' menu (green colour button) and from the drop down menu, choose the option 'Download Zip'. Then extract it in any directory of your desktop.

Instruction to run the code: Open the f90 code in editor. Model parameter values and initial molecule numbers are already given in the code. User can manipulate those values according to their need. Execute the code using fortran 90 compiler. After execution two datafiles will be created, 'parameter.dat' and  'output.dat'. Here, 'parameter.dat' will contain the model parameter values given as input and 'output.dat' will contain molecule numbers of mRNA, protein, miRNA and mRNA-miRNA complex at any time upto stopping time.
