
## 📌 Project Overview
It aims to provide a  structured workflow for heavy-ion collision simulations based on the following hybrid framework of:

1. Glauber model (GLISSANDO) that simulates  collision geometry and stored as `.dat` files.
2. That read by **VHLLE** to carried out Hydrodynamic evolution.
5. Followed by the  particle sampling and hadronic transport using **SMASH**.

The initial conditions used in this workflow are generated using a **Glauber collision model**.As a pre-generated `.dat` files describing the collision geometry and energy deposition.

---

##  Environment Setup:

This hybrid setup has been inspired by yukarpenko who develope the vhlle smash respository ,
Therefore,  the initial environment are performed by cloning the original repository:

```bash
git clone https://github.com/yukarpenko/run-vhlle-smash 
