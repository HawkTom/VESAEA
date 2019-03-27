# VESAEA

## Voronoi-based Efficient Surrogate-assisted Evolutionary Algorithm for Very Expensive Problems 

This repository is created for the source code of VESAEA. 

The VESAEA is an algorithm proposed for very computationally expensive problems, which the number of fitness evaluations is only 5D. (D is the dimensionality of the problems.)

The paper is accepted by CEC2019. If you are interested in this paper, you can contact me. 

**If you want use the code or VESAEA method in your paper, please cite my paper. The paper have not published in CEC2019 but you can cite the arXiv version now. The normal version will be updated after the paper is published by IEEE.**

The CALSAPSO is employed from [this](https://github.com/HandingWang/CALSAPSO). And the surrogate model used in the paper is from [this](https://sites.google.com/site/srgtstoolbox/). 

## Files 
```
surrogates -> surrogates model

CALSAPSO, EGO, GPEME, SSLAPSO -> four comparitive algorithms

TF, bounds -> test functions and problems' bounds

PSO, cmaes -> search algorithms

VESAEA -> proposed algorithm

VESAEA_woVLS -> VESAEA without Voronoi local search

compare -> compare VESAEA with other four algorithms

compare1 -> compare VESAEA with VESAEA_woVLS

```
## [*arXiv version*](https://arxiv.org/abs/1901.05755)

```
@article{tong2019voronoi,
  title={Voronoi-based Efficient Surrogate-assisted Evolutionary Algorithm for Very Expensive Problems},
  author={Tong, Hao and Huang, Changwu and Liu, Jialin and Yao, Xin},
  journal={arXiv preprint arXiv:1901.05755},
  year={2019}
}
```