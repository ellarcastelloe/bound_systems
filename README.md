# The Boundness Method - documentation and code

Most galaxy group finders identify "settled" groups that share a single dark matter halo and miss groups like our own Local Group that are gravitationally bound but do not yet share a common halo. The boundness method identifies bound multi-group systems consisting of multiple settled groups that are gravitationally bound. The boundness method allows us to identify "proto-groups"  in the early formation stages, as well as "Local Group analogues" that are similar to our own galaxy group. Identifying bound multi-group systems enables comparisons between groups in different formation stages and allows us to study  groups like our own are in the Universe.  

This repository contains code to test for gravitational boundness and calculate virialization state metrics for bound systems, as well as  documentation on the boundness method and virialization state metrics used.

# Using Catalogs of Bound Systems
<details>
  <summary>Details</summary>
  We test for gravitational boundess between neighboring settled groups (galaxy groups that share a common dark matter halo) to identify groups in early formation stages, including Local Group analogues. Below is a list of new attributes included in the catalog of bound systems in RESOLVE, ECO, and the mock catalogs. 
  
  * `boundFlag`: 1/0 flag for whether a galaxy is a member of a bound multi-group systems.
  * `lga`: 1/0 flag for whether a galaxy is a member of a Local Group analogue
  * `boundID`: Group ID for bound system. If `boundFlag = 0`, `boundID` and settled group ID `grp` match. If `boundFlag = 1`, `boundID` is a unique ID for the bound mutli-group system.
  * `grpR337`: Virial radius of settled group. 
  
    R<sub>337</sub> = (3 * 10<sup>logmh</sup>/4&pi; &Delta;<sub>mean</sub> &Omega;<sub>m</sub> &rho;<sub>crit</sub>)<sup>1/3</sup>
  
    Calculated using h=0.7, &Delta;<sub>mean</sub>=337, &rho;<sub>crit</sub> = 2.787e11 h<sup>2</sup> Msun/Mpc<sup>3</sup> and &Omega;<sub>m</sub> = 0.3075.
  
 For all attributes below, the quantity listed applies to the settled group if `boundFlag` = 0. If `boundFlag` = 1 the quantity applies to the bound multi-group system.
  * `boundN`: Number of galaxies in bound system.
  * `boundLog337`: Summed total of all halo masses in a bound system, using a halo mass convention of 337 times the background density.
  * `boundRproj`: Projected radius of bound system, calculated using method from Eckert+2017
  * `boundADalpha`: Alpha value obtained from Anderson-Darling test for bound systems with more than five members (`boundADalpha` = 0 if N<6). A higher `boundADalpha` means a system is more virialized. 
  * `boundTCross` (Gyr): System crossing time for bound systems with more than one member (`boundTCross` = 0 if N=1). We calculate crossing time following Firth+2006 as the average projected distance of group members from the group's coordinate center divided by the average velocity of group members. Systems with shorter crossing times are more virialized.
  * `boundLogG`: Log of group integrated gas mass. We use the `logmgas` column in RESOLVE and ECO for gas masses. 
  * `boundLogS`: Log of group integrated stellar mass
  * `boundURcolorgap`: Difference in u-r color between group central (galaxy with brightest r-magnitude) and brightest satellite, as in Eckert+2017             
  * `boundDSpval`: p-value from the Dressler & Shectman test for bound systems with more than 10 members (`boundDSpval` = 0 if N<11). A low p-value means that there is a high amount of subclustering, suggesting that the bound system is less virialized.
  * `boundR337overlap`: 1/0 flag for whether a bound multi-group system contains overlapping R337 of constituent settled groups. If two or more settled groups have overlapping R337 (using projection effect corrections to calculate the distance between group centers and the sum of R337 values for the two groups) then `boundR337overlap` = 1.

  For the mock catalogs, the above columns contain information for bound multi-group systems identified with projected data in the mocks using projection effect corrections (see section on testing for gravitational boundness between two settled groups below). There are also columns that have the same name  with `3d` on the end -- these columns contain information for "true" bound multi-group systems that were identified using 3D data in the mocks.  
  </details>
  
# Testing for Gravitational Boundness between two Settled Groups
  <details>
  We calculate whether two settled groups are gravitationally bound by comparing the relative velocity between the two groups to the escape velocity from one group at the location of the other group. If the relative velocity is smaller than the escape velocity, then the groups are bound. 
    
    
  
  **Calculating Escape Velocity**
    <details>
      We calculate the escape velocity from a chosen settled group at the distance of a neighboring settled group. Each settled group is treated as a point particle. 
      We use the equation
    v<sub>esc</sub> = (2GM/R<sub>grp-grp</sub>) <sup>1/2</sup>
    to calculate the escape velocity. M is the mass of the chosen group, calculated by summing the stellar and HI masses of each galaxy in the settled group with the halo mass of the group estimated using halo abundance matching (HAM). R<sub>grp-grp</sub> is the distance between the two groups. We calculate the projected distance between groups, R<sub>grp-grp (2D)</sub> using the Haversine formula. 
      
  We use the mock catalogs to correct for projection effects and approximate the 3D distance between group centers. For every pair of nearest neighbor groups in the mocks, we calculate the projected distance between groups using the Haversine formula. We also the true 3D distance between groups using the undistorted line-of-sight positions available in the mocks:
      
  R<sub>grp-grp (3D)</sub> = (R<sub>grp-grp(LOS)<sup>2</sup> + R<sub>grp-grp(2D)</sub><sup>2</sup>)<sup>1/2</sup>. 
      
  We create a distribution of R<sub>grp-grp (3D)</sub> / R<sub>grp-grp (2D)</sub> for every pair of nearest neighbor groups in the mocks. When testing for boundness in RESOLVE or ECO, we multiply the distribution of R<sub>grp-grp (3D)</sub> / R<sub>grp-grp (2D)</sub> from the mocks by the calculated R<sub>grp-grp (2D)</sub> for the pair of groups we're testing in RESOVLE or ECO, creating a distribution of possible R<sub>grp-grp (3D)</sub> values for that pair of groups. We use this distribution to calculate a probability of boundess, as described below. 
      
<img width="397" alt="forgithub_rgrpgrp" src="https://user-images.githubusercontent.com/46827591/124322644-0c40e680-db3d-11eb-9ee3-12b8f52c9f52.png">

  </details>
    
  **Calculating Relative Velocity between Groups**
    <details>
      We calculate the relative velocity between a chosen settled group and its neighbor group along the line-of-sight,
      v<sub>grp-grp (LOS)</sub> = |cz<sub>LOS, chosen group</sub>| - |cz<sub>LOS, neighbor group</sub>|. To approximate the 3D relative velocity between groups, we correct for projection effects using a similar method as described above for the escape velocity. For each pair of nearest neighbor groups in the mocks, we calculate v<sub>grp-grp (LOS)</sub>. We use the 3D velocity componenets for each galaxy in the mocks to calculate the 3D velocity of each settled group, and then calculate the relative 3D velocity between nearest neigbhbor groups. We create a distribution of v<sub>grp-grp (3D)</sub> / v<sub>grp-grp (LOS)</sub>. 
      
<img width="374" alt="forgithub_vgrpgrp" src="https://user-images.githubusercontent.com/46827591/124322719-31355980-db3d-11eb-9131-f1371ebfe041.png">

      
When testing whether a pair of groups in RESOLVE or ECO are bound, we multiply the calculated v<sub>grp-grp (LOS)</sub> for that pair of groups by the v<sub>grp-grp (3D)</sub> / v<sub>grp-grp (LOS)</sub> distribution and use the resulting distribution of possible v<sub>grp-grp (3D)</sub> for that pair of groups to calculate a probability of boundness, as described below. 
    </details>
    
  **Calculating Probability of Boundness**
    <details>
      Two settled groups are gravitationally bound if v<sub>esc (3D)</sub> > v<sub>grp-grp (3D)</sub>. We calculate the probability that a pair of groups is gravitationally bound using the distributions of v<sub>esc (3D)</sub> and v<sub>grp-grp (3D)</sub>.  We use Monte Carlo sampling with 10,000 samples each from the two distributions. We compare each pair of samples, and calculate the probability of boundness as the fraction of all samples where v<sub>esc (3D)</sub> > v<sub>grp-grp (3D)</sub>. In order for the groups to be considered part of the same bound multi-group system (`bound` = 1), the probability that the groups are bound must be greater than 90%. 
    </details>
    
    
  </details>
  

# Identifying Bound Multi-group Systems in a Catalog of Settled Groups
  <details>

  **Step 1: Identify nearest neighbor settled groups**
    
  We start with a catalog of "settled" groups that share a common dark matter halo, identified with FoF, the RESOLVE-G3 group finding algorithm (https://github.com/zhutchens1/g3groups) or another settled group finder. 
    
  We use a KD-Tree nearest neighbor search (https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.KDTree.html) to identify the 40 nearest neighbor settled groups to each settled group in the survey. 
  
  **Step 2: Test for gravitational boundness between neighboring settled groups**
    
  For each settled group in the group catalog:
    
  2a) Test whether "chosen" group is bound to nearest settled group
    
  If nearest neighbor groups are bound:
    
  2b) Assign groups the same `grp_bound` value and set `bound` = 1
    
  2c) Repeat steps 2a and 2b with incresasingly distant neighbors until chosen group is not bound to the neighboring group
    
  If nearest neighbor groups are not bound, the chosen group is not a member of a bound multi-group system, and the algorithm proceeds to the next settled group to test for boundness. 
    
  Settled groups are added to a bound multi-group system if they are bound to any other settled group already part of the bound multi-group system, so bound multi-group systems can continue to grow after they are first defined.
      
  **Step 3: Calculate properties of bound multi-group systems**
    
  Once boundness testing is finished, we calculate the properties of bound multi-group systems listed in the section "Using Catalogs of Bound Systems"
</details>
    
# Identifying Local Group Analogues in a Catalog of Bound Multi-group Systems
  <details>
    To identify galaxy group slike our own, we identify Local Group (LG) analogues that are a subset of bound multi-group systems. Each contains two giant galaxies, analogues for the Milky Way (MW) and Andromeda (M31), and their satellites identified by the settled group finder. To qualify as a LG analogue, a bound multi-group system must satisfy the following constraints:
    
  * Mass constraint (following Carlesi et al., 2019): the settled groups containing the MW and M31 analogues must each have a halo mass of at least 5x10<sup>11</sup>/h, and the two groups must have a combined halo mass of no more than 5x10<sup>12</sup>/h. The halo mass of the M31 analogue must be no more than 3 times greater than the halo mass of the MW analogue. 
  * The MW and M31 analogues must be separated by between 0.35-1.25 Mpc/h. We use our method for correcting for projection effects (see above) to estimate the 3D distance between galaxies, taking the distance between galaxies to be the median of the distribution of possible R<sub>grp-grp (3D)</sub> values. 
  * The settled groups containing the MW and M31 analogues must be gravitationally bound.
  * To ensure that the MW and M31 analogues are isolated from nearby large groups, the settled groups containing the MW and M31 analogues must not be bound to any other groups with halo mass above the gas-richness threshhold scale of 10<sup>11.5</sup> M<sub>sun</sub>.
  </details>

