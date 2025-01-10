[![arXiv][arxiv-shield]][arxiv-url]
[![DOI][doi-shield]][doi-url]
[![Documentation Status][docs-shield]][docs-url]
[![MIT License][license-shield]][license-url]

[doi-shield]: https://img.shields.io/badge/DOI-10.5281%20%2F%20zenodo.14617036-blue.svg?style=for-the-badge
[doi-url]: https://doi.org/10.5281/zenodo.14617036
[arxiv-shield]: https://img.shields.io/badge/arXiv-2501.05178-b31b1b.svg?style=for-the-badge
[arxiv-url]: https://doi.org/10.48550/arXiv.2501.05178
[license-shield]: https://img.shields.io/github/license/Jonas-Nicodemus/klap.svg?style=for-the-badge
[license-url]: https://github.com/Jonas-Nicodemus/klap/blob/main/LICENSE
[project-url]:https://github.com/Jonas-Nicodemus/klap/
[docs-shield]:https://img.shields.io/badge/docs-online-blue.svg?style=for-the-badge
[docs-url]:https://jonas-nicodemus.github.io/klap/

# KLAP: KYP Lemma based low-rank approximation for $\mathcal{H}_2$-optimal passivation
This repository contains the code for the paper [KLAP: KYP Lemma based low-rank approximation for H2-optimal passivation][arxiv-url].

<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary><h2 style="display: inline-block">Table of Contents</h2></summary>
  <ol>
    <li><a href="#citing">Citing</a></li>
    <li><a href="#installation">Installation</a></li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
  </ol>
</details>

## Citing
If you use this project for academic work, please consider citing our
[publication][arxiv-url]:

    J. Nicodemus, M. Voigt, S. Gugercin, and B. Unger
    KLAP: KYP Lemma based low-rank approximation for $\mathcal{H}_2$-optimal passivation
    ArXiv e-print 2501.05178, 2025.

## Installation
This code base is using the [Julia Language](https://julialang.org/) and
[DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> klap

To (locally) reproduce this project, do the following:

1. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
2. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

You may notice that most scripts start with the commands:
```julia
using DrWatson
@quickactivate "klap"
```
which auto-activate the project and enable local path handling from DrWatson.

## Usage
There are two executable scripts [`main.jl`](https://github.com/Jonas-Nicodemus/klap/tree/main/scripts/main.jl) and [`smartphone.jl`](https://github.com/Jonas-Nicodemus/klap/tree/main/scripts/smartphone.jl) located in the `scripts` directory:

- [`main.jl`](https://github.com/Jonas-Nicodemus/klap/tree/main/scripts/main.jl) applies and compares KLAP and LMI-based passivation methods on the ACC or CD player benchmark models.
- [`smartphone.jl`](https://github.com/Jonas-Nicodemus/klap/tree/main/scripts/smartphone.jl) applies KLAP to the smartphone benchmark model and compares the results with the model obtained by perturbation of Hamiltonian eigenvalues.

Furthermore, an [`example.ipynb`](https://github.com/Jonas-Nicodemus/klap/tree/main/notebooks/example.ipynb) notebook is provided in the `notebooks` directory, which generates the contour plots from the Examples 3.7, 3.8, 3.11 and 3.12.

## License
Distributed under the MIT License. See `LICENSE` for more information.

## Contact
Jonas Nicodemus - jonas.nicodemus@simtech.uni-stuttgart.de

Matthias Voigt - matthias.voigt@fernuni.ch\
Serkan Gugercin - gugercin@vt.edu\
Benjamin Unger - benjamin.unger@simtech.uni-stuttgart.de

Project Link: [https://github.com/Jonas-Nicodemus/klap][project-url]