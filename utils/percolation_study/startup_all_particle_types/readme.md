# Percolation simulation pipeline 

percolation_startup.sh

Input: T, phis and particle types to be initalized 
Output: directory structures ready for submission on vsc3. Initally the simluations run with energy_level = 0 for 1x10^10 steps and without cluster moves. This is necessary to properly intialize configurations.

Synopsis: Copies files into directories with naming scheme and sets parameters in copied para.ini and COMMITT.sh.

Note: 
For further parallel runs, we suggest to use the generated configurations as startup and sample more configs subsequently, after autocorr of positions and orientations has sufficiently decayed. 


 
