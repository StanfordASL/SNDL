# SNDL: Stabilizable Nonlinear Dynamics Learning

Code to accompany: [Learning Stabilizable Nonlinear Dynamics with Contraction-Based Regularization](https://arxiv.org/pdf/1907.13122.pdf)

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. All code is written in MATLAB.

### Prerequisites

There are five packages required: 
* [YALMIP](https://yalmip.github.io/) : parsing the learning problems (convex), 
* [Mosek](https://www.mosek.com/) : SDP solver for the learning problems, 
* [TOMLAB](https://tomopt.com/tomlab/) : parsing the trajectory optimization and CCM controller problems (non-convex), 
* [SNOPT, NPSOL](https://ccom.ucsd.edu/~optimizers/) : solvers for the non-convex problems,
* [cvx](http://cvxr.com/cvx/download/) : auxiliary functions. 

YALMIP and cvx are freely available. Mosek is freely available on an academic license. TOMLAB must be purchased (a 21-day trial license is available on their website). Evaluation copies of SNOPT and NPSOL are included with the trial TOMLAB download. In the near future, I will also release open-source alternatives for TOMLAB and SNOPT; additional details below. 

### Installing

Having installed the prerequisites (and adding them to the MATLAB path), download the repo, navigate to the repo in MATLAB, and execute the following in the command window:

```
run sndl_startup.m
```

This will add all necessary sub-directories of the repo to the MATLAB path.

## Workflow

This will take you through the various steps involved in reproducing the results in the paper. PVTOL simulation example in main branch; hardware version in pvtol_h branch. 

### Dataset generation (for Simulation)

Navigate to the data_gen folder, automated script: *generate_PVTOL_data.m*. Choose number of trajectories to generate and additional samples for constraint set as indicated within the script. Output files: *PVTOL_data_train.mat* and *PVTOL_data_val.mat*. 

Dataset files for hardware branch already included. 

### Learning

Navigate to main repo directory, automated script: *Main.m*. Adjustable parameters: 
* *Main.m*: *N_tr*: # points to use for training (i.e., included in regression loss). This must be less than *N_max* (total number of demonstration tuples).
* *Main.m*: *Nc*: # points to use for initial constraint set. This constraint set will dynamically change within the algorithm to ensure constraints hold at all points in X. 
* *load_PVTOL_params.m*: Parameter file for PVTOL. See definitions within.  

After this, the code enters the script *learn_loop.m* corresponding to Algorithm 1 in the paper. It will first compute and save two baseline models (unconstrained-unregularized, and unconstrained-l_2 regularized least squares). Following this, the CCM-regularized model will be computed iteratively. 

All saved functions are stored in *learned_functions/* and marked with identifiers 'unconu' and 'uncon' for the two baselines, and 'ccm' for the CCM-regularized model. The additional identifier is the number *N_tr*. Example solutions are provided.

### Evaluation (for Simulation)

Navigate to the test folder. In order:
* *gen_test_traj_PVTOL.m*: Generate test trajectories for the various models. Trajectory optimization uses the Pseudospectral Method with Lagrange interpolating polynomials, CGL nodes, and Clenshaw-Curtis quadrature; see also [Fahroo et. al.](https://arc-aiaa-org.stanford.idm.oclc.org/doi/pdfplus/10.2514/2.4862). Outputs: trajectory datasets *PVTOL_test_traj_* with identifiers: model type ('unconu', 'uncon', 'ccm') and training dataset size. 
* *test_dynamics_PVTOL_bulk.m*: Simulate and record data for each model on the generated test trajectories in bulk. Alternatively, load up the trajectory data from any of the generated trajectory datasets, and directly use the simulation functions *test_dynamics_PVTOL_LQR* or *test_dynamics_PVTOL_CCM*, to evaluate with TV-LQR or tracking MPC.
* *analyze_PVTOL_results.m*: Compare all models with shared (regression) training datasets. 

Notes: 
* Trajectory optimizer uses TOMLAB and SNOPT. To switch to alternative package, you can still use the same setup functions (constraints, gradients, Hessians, etc) and simply change the optimization call. An example using fmincon will be uploaded soon. 

* Trajectory generation file for hardware experiments: *gen_test_traj_PVTOL_hardware.m*

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

