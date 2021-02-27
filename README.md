# RMVB Read Me
Read Me. This readme document describes software codes which have been developed under support from the National Institutes of Health via grants EB021027, NS096761, MH114233, AT009263 and NSF CAREER CCF-1651825 awarded to Prof. Bin He. 
The source codes are provided as a service to the scientific community and may be used for any non-commercial purposes.  Users should use the codes or data at their own risks with no guarantees provided, as is. If anyone should use the codes provided here in its entirety or partially, we ask them to cite the following publication in any of their publications or presentations:
“S. A. Hossein Hosseini, A. Sohrabpour, M. Akçakaya and B. He, "Electromagnetic Brain Source Imaging by Means of a Robust Minimum Variance Beamformer," in IEEE Transactions on Biomedical Engineering, vol. 65, no. 10, pp. 2365-2374, Oct. 2018, doi: 10.1109/TBME.2018.2859204.”
------------------------------------------------------------------------------------------------------------------------------------------

This folder contains the codes used to analyze (simulated) EEG data within the robust minimum variance beamformer (RMVB) framework. The same codes provided here, were run to obtain results presented in the aforementioned paper and the results presented here are derived with this exact code, to ensure accuracy.
In preparation of the codes, the following (modified) packages have been borrowed:
- group_lasso (https://web.stanford.edu/~boyd/papers/admm/group_lasso/group_lasso.html)
- NewtonRaphson_0.5 (https://github.com/mikofski/NewtonRaphson)
- regu (http://www2.compute.dtu.dk/~pcha/Regutools/)
- CVX package (http://cvxr.com/cvx/)
-assignmentoptimal (https://www.mathworks.com/matlabcentral/fileexchange/6543-functions-for-the-rectangular-assignment-problem)

Except the CVX package, which should be installed on your computer, all other packages/codes are already provided in the project’s folders for your convenience.

How to Run the Program:
To run these codes, you need to have MATLAB installed on your computer (the current code has been tested on MATLAB 2017a and 2019b). The plots in the paper are also prepared in MATLAB. 
Solutions to the simulated imaging analysis, are saved in a subfolder of folder “Resource” named “DAT Files”. An empty “Resource\DAT Files” folder should be created before running the codes. It is also necessary to download folder “Resource” and place it next to the codes, since it contains important data (such as the lead-field matrices) that are necessary for the proper execution of the program. 

The file called “main.m” summarizes the pipeline for generating an example. This file calls 5 other main scripts as follows:
- Setup_Electrode.m (to set up electrodes configuration) 
- Setup_LFD.m (to generate uncertainty and inverse lead-field matrices)
- Setup_Source_vExample.m (to simulate source activity in the brain) 
- Main_Forward_vExample.m (to solve the forward problem)
- Main_Inverse_vExample.m (to solve the inverse problem)

Additionally, “main_setup.m” defines some of the high-level settings for the simulation. 

The figures generated by running the “Main_Inverse_vExample.m” plot the reconstructed activity at the first point of the time course by default. You can use arrow keys to change the time point at which the solution is depicted in the current figure. In addition, pressing ctrl + p will show the time point corresponding to the maximum activity. For other options, see function “setup_display.m”.

Finally, the provided codes implement a fast algorithm based on Lagrange multipliers to solve the second-oder cone programming (SOCP) of the RMVB by default. To change this algorithm to CVX, you need to uncomment/comment line 1387/1388 in the “solve_inverse_problem.m”.  
  
Seyed Amir Hossein Hosseini

February 24th, 2020
