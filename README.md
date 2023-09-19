# Three-stage model of stochastic gene expression
This library provides data and codes for the manuscript [Exact solution of a three-stage model of stochastic gene expression including
cell-cycle dynamics](https://www.biorxiv.org/content/10.1101/2023.08.29.555255v2.full.pdf).

________________________________________________________________________________________________________
## Requirements
- Mathematica v13.2.1.0
- pandas v1.3.2
- numpy  v1.16.5～v1.23.0
- scipy v1.7.1
## File description
- "SSA_model _II _t10.csv" stores example data of Model II simulated by SSA.
- "exact_solution_Model_II.nb" computes the exact distributions of protein numbers of Model II in Mathematica.
- "population_SSA_IV.ipynb" generates SSA data for cell population in Model IV.
- "stationary_statistics_Model_III.nb"  computes the stationary gene
  product number statistics of Models I & III.

## Examples
__1. Exact solution of Model II.__  

The time-dependent exact distribution of protein counts of Model II is calculated from the generating function by the following piece code.
```
G = G0 + G1 /. param;
Bins = 80;
Gp = G /. w1 -> 0;
PP = ResourceFunction["NSeries"][Gp, {w2, -1, Bins}][[3]];
v = Table[{i - Bins - 1, Re[PP[[i]]]}, {i, Bins + 1, 2*Bins + 1}];
pG = ListPlot[v, PlotRange -> All]
```
`G` is the generating function. If one is interested in the probability distribution of mRNA, switch the positions of w1 and w2 in the codes above. For more details, please refer to the notes in [exact_solution_Model_II.nb]().

__2. SSA for a unsynchronized cell population of Model IV.__  

 `population_SSA` generates the gene product counts for a unsynchronized cell population of Model IV, which was adapted from the codes of Ref. [1]
```
data=population_SSA(m0,G0,G1,p0,t0,phase0,age0,Tmax,Ncycle,Tcycle,son,soff,rho,lam,dm)
```
`G0`,`G1`,`m0`,`p0`,`t0`,`phase0`,`age0` are the initial conditions.  
`Tmax` is the total simulation time.  
The parameter `k` for the exponential distribution is calculated by using `k`=`Ncycle`/`Tcycle`.  
`son`,`sof`,`rho`,`lam`,`dm` are the kinetic parameters.  
`data` is a matrix storing the information of
absolute time、 acitve gene、 inactive gene、 mRNA、 protein、 cell age and cell phase, each of which is represented by a rwo in `data`. For more details, please refer to the notes in [population_SSA_IV.ipynb]().  

__3. The exact statistics of Model III.__  

[stationary_statistics_Model_III.nb]() compares the stationary mRNA and protein statistics of Models I & III.


## Reference
[1] Beentjes, C. H., Perez-Carrasco, R., & Grima, R. (2020). [Exact solution of stochastic gene expression models with bursting, cell cycle and replication dynamics](https://journals.aps.org/pre/pdf/10.1103/PhysRevE.101.032403?casa_token=7ZbwkB0N77wAAAAA%3AFoM22TR7q45nanerhg1LpWNx4WfMr1Uk5Db0BV5er0s6i1kC0V_2m_claH3F7NQMR1pgUeJvEulG5Qmr). *Physical Review E*, *101*(3), 032403.

If you found this library useful in your research, please consider citing.

```
@article{wang2023exact,
  title={Exact solution of a three-stage model of stochastic gene expression model including cell-cycle dynamics},
  author={Wang, Yiling and Yu, Zhenhua and Cao, Zhixing and Grima, Ramon},
  journal={bioRxiv},
  pages={2023--08},
  year={2023},
  publisher={Cold Spring Harbor Laboratory}
}
```
