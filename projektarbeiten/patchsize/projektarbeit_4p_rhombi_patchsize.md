## Projektarbeit: Self assembly of patchy rhombi -- The role of patch size in feq-as2 systems. 

The aim of this project is to study the influence of the patch size on the self assembly of patchy rhombi.

#### Motivation

- **Molceular Level**
For some tetracarboxylic acids, the respective rhombi models overlap where the molecules would not overlap. Larger patches allow a larger edge to edge distance and therefore eleminate the problematic overlap.

![08583e3a.png](https://gitlab.com/catnat/rhombi/blob/master/projektarbeiten/patchsize/4p-overlap.png =250x)
**FIG.1** Overlap in rhombi model of a tetra-acid.
- **Colloidal Level**
In experiments with patchy colloids it is often easier to craft larger patches.
Too large patches however tend to yield a gel or other types of unordered structures instead of an ordered tiling.
In this case, computer simulations can be used to find the critical patch size after which the 
the system changes to an unordered structure.

Pleniminary results for **feq-as2** show 
that the critical patch size, depending on $\Delta$, is around $r_{p}=0.1$ for an edge length $L_x = 1.0$.

![88454ed0.png](https://gitlab.com/catnat/rhombi/blob/master/projektarbeiten/patchsize/delta_03_r01.png =250x)
**FIG.2** Onp tiling at $r_p = 0.1, \Delta=0.3, \epsilon=-4.2k_{B}T, \mu=0.25$

![538d2ea8.png](https://gitlab.com/catnat/rhombi/blob/master/projektarbeiten/patchsize/dma-as2-rp01.png =250x)
**FIG.3** Unordered assembly product at $r_p = 0.1, \Delta=0.2, \epsilon=-4.2 k_{B}T, \mu=0.25$

#### Relevant literature

[Design of patchy rhombi: from close-packed tilings to open lattices](https://arxiv.org/abs/1906.10938)

[Phys. Rev. Lett. 108, 035702 (2012) - Random and Ordered Phases of Off-Lattice Rhombus Tiles](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.108.035702)

[Broken symmetry and the variation of critical properties in the phase behaviour of supramolecular rhombus tilings \| Nature Chemistry](https://www.nature.com/articles/nchem.1199)


#### Task Description

##### Main Reserach Question
What is the critical patch size after which we observe an unordered structure?

##### Sub tasks
1. Run self-assembly simulations (grand canonical Monte-Carlo) on the Vienna Scientific Cluster (VSC-3) for patch radii (sizes) between 0.02 and 0.2 and interaction energies $\epsilon \in \{-3.8,-4.2,-4.8,-5.2\}k_{B}T$.Start a new simulation at $\mu = \mu_{eq} = 0.1$. After sufficient equilibration, change $\mu$ to $\mu* = 0.25$. 
2. Produce several runs per patch size/interaction strength to collect statistics on the self assembly product (at least 16 runs/per parameter state).
3. Automate the simulation runs (automatic generation of directories/simulation startpoints/checkpoints for each simulation run).
4. Regulary check the status of the assembly (the equilibration). Relevant quantities to check are: the ratio largest cluster/system size and packing fraction and total energy as function of MC Sweeps. It is also advisable to visually inspect some of the produced snapshots.
5. Collect statistics for checkpoints after equilibration of the system is reached. We are especially interested in the **degree of randomness $\Psi$**. The relevant output file is **psi_op.dat**.
6. Calculate and plot average randomness and randomness distribution for each patch size/ interaction strength.

##### Detailed Research Questions
- What is the critical patch size, after which the tilings become unordered?
- Is the transition to the unordered structure sharp or sequential? How would you quantify the sharpness of this transistion?
- Does the transition depend on the patch position $\Delta$? 
- What are the properties of the random network at big patch size? Is there a difference to the random tiling of feq-center at small patch size? If yes, how could this difference be quantified? 

#### Code Documentation
- **to clone** source from gitlab repository:
```bash
git clone git@gitlab.com:catnat/rhombi.git
```
###### Self assembly simulation (McPoly)
- **pre-requisits**: g++ -8, boost, gsl
- **name of binary**: McPoly
- **location**: in rhombi/source

- **to compile**:
```bash 
make all
```
- **required input to run**:
-- the compiled binary **McPoly**
-- the parameter file **para.ini** set with the correct paramters (automated with a python script over all paramter choices and runs)
- intial run:
```bash 
./McPoly -n <start time> <end time> 
``` 
- **run from specific checkpoint**
```bash 
./McPoly -f <checkpoint> <end time> <random var 1> <random var 2>
```
- **relevant output**:
-- positions_<checkpoint>.bin: binary file, contains center postions. 
-- orientations<checkpoint>.bin: binary file, contains rhombi orientations.
-- Box_<checkpoint>.bin: binary file, contains Box size and center position. 
-- RANDOM_STATE_R1_<checkpoint>.bin: binary file, contains state of the RNG 1. 
-- RANDOM_STATE_R2_<checkpoint>.bin: binary file, contains state of the RNG 2.
-- psi_op.dat: contains the randomness paramter for largest cluster, and whole system for all checkpoints.
-- NPT_OUT.txt: contains packing fraction for all checkpoints. 
-- All_Clusters_info.dat: sizes and particle ids for all clusters for all checkpoints.
-- All_Clusters_info.dat: sizes for all clusters  for all checkpoints. 
-- positions_.xyz: xyz-file for visual inspection
-- color_op.dat: orientation file for visual inspection
-- pos_c.xyz: xyz-file of center positions
###### Visual Inspection of Snapshots 
**pre-requisits**: vmd
**name of script**: vmd_patchy_rhombi.tcl 
- **location**: in rhombi/utils/visual_inspection

- **to run vmd to inspect rhombi**: vmd_patchy_rhombi.tcl and 
```bash 
vmd -e vmd_patchy_rhombi.tcl -args positions_.xyz 
``` 
( required extra input: color_op.dat)

- **to run vmd to inspect center positions**:
```bash
vmd pos_c.xyz
```
###### Running code on the VSC-3 cluster 
- **name of script**: submit.slrm
- **location**: in rhombi/submit_vsc3
- **to load modules** for compilation and running: 
```bash
module load <name of module>
```
- **to submit jobs**:
copy submit.slrm into simulation run directory. Then type:
```bash
sbatch submit.slrm
```
- **To check state** of simulation
```bash 
squeue -u <your username> 
```
- **To cancel** jobs
```bash 
scancel <job_id>
```

