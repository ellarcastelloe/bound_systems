# The Boundness Method - documentation, catalogs of bound systems, and code
This repository contains code to test for gravitational boundness and calculate virialization state metrics for bound systems, as well as catalogs of bound systems for RESOLVE, ECO, and mock catalogs and documentation on the boundness method and virialization state metrics used.

The boundness method can be used to identify systems of galaxies containing multiple "settled" groups that share a common dark matter halo. Bound systems include "Local Group analogues" that are similar to our own galaxy group, as well as low-mass proto-groups that are in early formation stages. 


# Using Catalogs of Bound Systems
<details>
  <summary>Details</summary>
  We test for gravitational boundess between neighboring settled groups (galaxy groups that share a common dark matter halo) to identify groups in early formation stages, including Local Group analogues. Below is a list of new attributes included in the catalog of bound systems in RESOLVE, ECO, and the mock catalogs. 
  
  * `bound`: 1/0 flag for whether a galaxy is a member of a bound multi-group systems.
  * `lga`: 1/0 flag for whether a galaxy is a member of a Local Group analogue
  * `grp_bound`: Group ID for bound system. If `bound = 0`, `grp_bound` and settled group ID `grp` match. If `bound = 1`, `grp_bound` is a unique ID for the bound mutli-group system.
  * `R337`: Virial radius of settled (FoF) group. 
  
    R<sub>337</sub> = (3 * 10<sup>logmh</sup>/4&pi; &Delta;<sub>mean</sub> &Omega;<sub>m</sub> &rho;<sub>crit</sub>)<sup>1/3</sup>
  
    Calculated using h=0.7, &Delta;<sub>mean</sub>=337, &rho;<sub>crit</sub> = 2.787e11 h<sup>2</sup> Msun/Mpc<sup>3</sup> and &Omega;<sub>m</sub> = 0.3075.
  
 For all attributes below, the quantity listed applies to the settled FoF group if `bound` = 0. If `bound` = 1 the quantity applies to the bound multi-group system.
  * `grpn_bound`: Number of galaxies in bound system.
  * `logmh337_bound`: Summed total of all halo masses in a bound system, using a halo mass convention of 337 times the background density.
  * `R337_bound`: R337 calculated using `logmh337_bound`. Note that this doesn't represent the radius of a bound multi-group system, but rather the radius of a halo containing all of the mass in the bound multi-group system. 
  * `Rproj_bound`: Projected radius of bound system, calculated using method from Eckert+2017
  * `ad_alpha`: Alpha value obtained from Anderson-Darling test for bound systems with more than five members (`ad_alpha` = 0 if N<6). A higher `ad_alpha` means the system is more virialized. 
  * `t_cross` (Gyr): System crossing time for bound systems with more than one member (`t_cross` = 0 if N=1). Systems with shorter crossing times are more virialized.
  * `grp_loggascontent`: Log of group integrated gas-to-stellar mass ratio.
  * `ur_colorgap`: Difference in u-r color between group central (galaxy with brightest r-magnitude) and brightest satellite, as in Eckert+2017             </details>
  
# Testing for Gravitational Boundness between two Settled Groups
  <details>
  We calculate whether two settled groups are gravitatinally bound by comparing the relative velocity between the two groups to the escape velocity from one group at the location of the other group. If the relative velocity is smaller than the escape velocity, then the groups are bound. 
    
    
  
  **Calculating Escape Velocity**
    <details>
      We calculate the escape velocity from a chosen settled group at the distance of a neighboring settled group. Each settled group is treated as a point particle. 
      We use the equation
    v<sub>esc</sub> = (2GM/R<sub>grp-grp</sub>) <sup>1/2</sup>
    to calculate the escape velocity. M is the mass of the chosen group, calculated by summing the stellar and HI masses of each galaxy in the settled group with the halo mass of the group estimated using halo abundance matching (HAM). R<sub>grp-grp</sub> is the distance between the two groups. We calculate the projected distance between groups, R<sub>grp-grp (2D)</sub> using the Haversine formula. 
      
  We use the mock catalogs to correct for projection effects and approximate the 3D distance between group centers. For every pair of nearest neighbor groups in the mocks, we calculate the projected distance between groups using the Haversine formula. We also the true 3D distance between groups using the undistorted line-of-sight positions available in the mocks:
      
  R<sub>grp-grp (3D)</sub> = (R<sub>grp-grp(LOS)<sup>2</sup> + R<sub>grp-grp(2D)</sub><sup>2</sup>)<sup>1/2</sup>. 
      
  We create a distribution of R<sub>grp-grp (3D)</sub> / R<sub>grp-grp (2D)</sub> for every pair of nearest neighbor groups in the mocks. When testing for boundness in RESOLVE or ECO, we multiply the distribution of R<sub>grp-grp (3D)</sub> / R<sub>grp-grp (2D)</sub> from the mocks by the calculated R<sub>grp-grp (2D)</sub> for the pair of groups we're testing in RESOVLE or ECO, creating a distribution of possible R<sub>grp-grp (3D)</sub> values for that pair of groups. We use this distribution to calculate a probability of boundess, as described below. 
      
![forgithub_rgrpgrpdist](https://user-images.githubusercontent.com/46827591/123892771-70cf2c00-d918-11eb-86b5-c5550bf0f77f.png)
    </details>
    
  **Calculating Relative Velocity between Groups**
    <details>
      We calculate the relative velocity between a chosen settled group and its neighbor group along the line-of-sight,
      v<sub>grp-grp (LOS)</sub> = |cz<sub>LOS, chosen group</sub>| - |cz<sub>LOS, neighbor group</sub>|. To approximate the 3D relative velocity between groups, we correct for projection effects using a similar method as described above for the escape velocity. For each pair of nearest neighbor groups in the mocks, we calculate v<sub>grp-grp (LOS)</sub>. We use the 3D velocity componenets for each galaxy in the mocks to calculate the 3D velocity of each settled group, and then calculate the relative 3D velocity between nearest neigbhbor groups. We create a distribution of v<sub>grp-grp (3D)</sub> / v<sub>grp-grp (LOS)</sub>. 
      
![forgithub_vgrpgrpdist](https://user-images.githubusercontent.com/46827591/123893005-d6bbb380-d918-11eb-88eb-ee9a9d294ffa.png)

      
When testing whether a pair of groups in RESOLVE or ECO are bound, we multiply the calculated v<sub>grp-grp (LOS)</sub> for that pair of groups by the v<sub>grp-grp (3D)</sub> / v<sub>grp-grp (LOS)</sub> distribution and use the resulting distribution of possible v<sub>grp-grp (3D)</sub> for that pair lof groups to calculate a probability of boundness, as described below. 
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
    Local Group (LG) analogues are a subset of bound multi-group systems. Each contains two giant galaxies, analogues for the Milky Way (MW) and Andromeda (M31), and their satellites identified by FoF. To qualify as a LG analogue, a bound multi-group system must satisfy the following constraints:
    
  * Mass constraint (following Carlesi et al., 2019): the FoF groups containing the MW and M31 analogues must each have a halo mass of at lease 5x10<sup>11</sup>/h, and the two groups must have a combined halo mass of no more than 5x10<sup>12</sup>/h. The halo mass of the M31 analogue must be no more than 3 times greater than the halo mass of the MW analogue. 
  * The MW and M31 analogues must be separated by between 0.35-1.25 Mpc/h. We use our method for correcting for projection effects (see above) to estimate the 3D distance between galaxies, taking the distance between galaxies to be the median of the distribution of possible R<sub>grp-grp (3D)</sub> values. 
  * The FoF groups containing the MW and M31 analogues must be gravitationally bound.
  * To ensure that the MW and M31 analogues are isolated from nearby large groups, the FoF groups containing the MW and M31 analogues must not be bound to any other groups with halo mass above the gas-richness threshhold scale of 10<sup>11.5</sup> M<sub>sun</sub>.
  </details>

