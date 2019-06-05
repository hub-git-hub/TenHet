Project
=====
TenHet, a system for performing tensor completion with  multi-view side information. The Alternating Direction Method of Multipliers (ADMM) is used to solve optimization problem. 

Code and Dataset
=====
`GenerateSyntheticData`: generate the synthtic tensor with multiple matrices.

`TenHet`: Tensor completion based on ADMM.

`solveE`: Compute the error tensor.

`SoverM`: compute the variable M.



`DrugID`, `DrugChem`, `DrugSide`: information about the drug in DrugBank.

`TargetID`, `TargetGO`, `TargetSW`: information about target protein in DrugBank.

`DiseaseID`, `DiseaseMin`, `DiseaseHPO`: information about disease in DrugBank.


Running
======
The implementation of TenHet has been tested on MATLAB2015. An example of running can be found in `Demo`. [Matlab Tensor Toolbox](http://www.sandia.gov/~tgkolda/TensorToolbox/index-2.6.html) is need and has been included in the package. 

Demo: 

(1) run `runSyntheticData` 

(2) run `runDrugBank`
