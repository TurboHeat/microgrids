# Economic dispatch of a single micro gas turbine under CHP operation with uncertain demands

This contains the source code for the article in the title, found in https://doi.org/10.1016/j.apenergy.2021.118391

```tex
@article{Sharf_2022,
	title        = {Economic dispatch of a single micro gas turbine under {CHP} operation with uncertain demands},
	author       = {Miel Sharf and Iliya Romm and Michael Palman and Daniel Zelazo and Beni Cukurel},
	year         = 2022,
	month        = {mar},
	journal      = {Applied Energy},
	publisher    = {Elsevier {BV}},
	volume       = 309,
	pages        = 118391,
	doi          = {10.1016/j.apenergy.2021.118391},
	url          = {https://doi.org/10.1016%2Fj.apenergy.2021.118391}
}
```

**Contributors (in chronological order): Johannes F. Rist, Miguel de Freitas Dias, Iliya Romm, Anoop Jain, Miel Sharf.**

See also:

> Johannes F. Rist, et al. "_Economic dispatch of a single micro-gas turbine under CHP operation._" Applied Energy 200 (2017): 1-18. [DOI](https://doi.org/10.1016/j.apenergy.2017.05.064)

---

## Preparations

1) All data is saved in the folder called "Data", one level above this code package. This folder must:

    - ... contain a file with maps of a given turbine. This file must be called either `CHP.mat` or `sv_mappings.mat`.
    - ... have about 3.5GB of free space available for the saving of intermediate results.

1) `nakeinterp1.c` must be [compiled to MEX](https://www.mathworks.com/help/matlab/ref/mex.html).
1) Some version of the [`progressbar` library](https://github.com/fsaxen/ParforProgMon) must be on the MATLAB path.

---

## Description of Functions and Scripts

<span style="font-size:9pt;">Below are the functions and scripts needed to generate the transition table and case studies. Entries are ordered according to the sequential steps needed to go generate the econ report.</span>

---

**`createTransitionNetwork`** - generates the transition table. Does not need any external data files. The script is divided into several parts:

1. Build a small adjacency matrix, defining all possibile transitions from all states `(t,S,V)` for a fixed `t`, `S=1,...,s` and `V = 0,...,v-1`.
2. Replicate the adjacency matrix `T` times, where `T` is the number of time instances in the economic dispatch problem.
3. Add the source and terminal nodes.
4. Build the graph from the composite adjacency matrix.
Here, the nodes are indexed as `1,2,3,...`:

    - Node `1` is the source node. 
    - The next `(s*v+1)` nodes correspond to time `0`, followed by `(s*v+1)` nodes for time `1`, etc.
    - This repeats a total of `T` times, culminating in a single terminal node.  

    Inside each layer, the nodes are organized as follows:

    - The first node is the `'off'` state.
    - The rest are ordered in the following way:

       `(1,0),(2,0),...(s,0),(1,1),(2,1),...,(s,1),......,(1,v-1),(2,v-1),...,(s,v-1)`

- This script does not need to be run every time, only when changes to the structure of the network are required. 
- All the other scripts use the output of this one, called `graph_24h.mat`.

---

**`assignAllCosts`** - assigns the costs of the edges for all days and buildings combinations. Also generates the demand profiles and electricity tariffs and stores them
in a matrix. Smoothing of the demand data is performed here.

All the data needed to run the shortest path solver is assembled and generated here and the results are
saved in the file `graph_data_all_days.mat`. This script uses the function `assignCosts`, which is the building block of edge cost assignment.

Furthermore, the function `transformAdjacencyMatrix` is used to conduct auxiliary computations, allowing for a faster computation.

---

**`generate_all_figures`** - runs the shortest path solver, plots figures for all days, buildings and fuel costs, and generates all the economic metrics.

---

**`econAnalysisTables`** - with the data generated from `generate_all_figures`, performs analysis of economic metrics and outputs the results as an Excel spreadsheet.

---

**`createDemandProfileVector`** - Johannes' code, returns the electricity and heat demand profiles (note that the large hotel is the first building).

---

**`createElectricityTariffProfile`** - generates the electricity tariff profiles based on the pfd references, for all buildings.

---

**`generateDemandCharges`** - computes the demand charges for a given building and day. Used by `generate_all_figures`.

---

**`extractPath`** - retrieves the solution of the shortest path algorithm as power and heat production profiles, as well as the fuel consumption at each time step.

Uses the function `nakeinterp1` to interpolate the production when `w` the speed is ramped up so that the final length of the vector is correct.

---

**`transformAdjacencyMatrix`** - transforms the transitions table generated by `createTransitionNetwork` into a more memory-efficient representation.

---

**`runAll`** - runs `assignAllCosts`, `generate_all_figures` and `econAnalysisTables` in sequence.

---

### Notes

If there is a need to speed up the code, some computations can be moved to the GPU.
