

___StatusSwitchingQP.jl___


[![Build Status](https://github.com/PharosAbad/StatusSwitchingQP.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/PharosAbad/StatusSwitchingQP.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://github.com/PharosAbad/StatusSwitchingQP.jl/wiki)

<h1 align="center" margin=0px>
  Status Switching Method for Quadratic Programming and Linear Programming
</h1>

<p align="center">
  <a href="#features">Features</a> •
  <a href="#installation">Installation</a> •
  <a href="#license-">License</a> •
  <a href="https://github.com/PharosAbad/StatusSwitchingQP.jl/wiki">Documentation</a>
</p>

**StatusSwitchingQP.jl** solves the following convex quadratic programming problems (called `SSQP`):


$$
\begin{array}
[c]{cl}
\min & \frac{1}{2}\mathbf{x}^{\prime}\mathbf{Vx}+\mathbf{x}^{\prime}
\mathbf{q}\\
s.t. & \mathbf{Ax}=\mathbf{b}\in\mathbb{R}^{M}\\
& \mathbf{Gx}\leq\mathbf{g}\in\mathbb{R}^{J}\\
& \boldsymbol{d}\leq\mathbf{x}\leq\boldsymbol{u}\in\mathbb{R}^{N}
\end{array}
$$

with positive semi-definite symmetric matrix $\mathbf{V}\in\mathbb{R}^{N\times N}$.

## Features

* __Fast and Accurate__:  [Speed and Accuracy](https://github.com/PharosAbad/StatusSwitchingQP.jl/wiki/Speed-and-Accuracy)
* __Open Source__: Our code is available on [GitHub](https://github.com/PharosAbad/StatusSwitchingQP.jl) and distributed under the MIT License
* __Arbitrary Precision Arithmetic__: fully support for `BigFloat`


## Installation
__StatusSwitchingQP.jl__ can be added by

- `import Pkg; Pkg.add("StatusSwitchingQP")`
- `pkg> add StatusSwitchingQP`
- `pkg> dev StatusSwitchingQP` for testing nightly version. To use the registered version again `pkg> free StatusSwitchingQP`

## License 🔍
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

