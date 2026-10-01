# Simplified MITC4+ Shell Element MATLAB Implementation

This repository contains a full MATLAB module and benchmark suite implementing the **Simplified MITC4+ Shell Element** described in the research paper:

> **The simplified MITC4+ shell element and its performance in linear and nonlinear analysis**
> Hyung-Gyu Choi and Phill-Seung Lee
> *Computers & Structures*, Vol. 290, 2024, 107177.

---

## 📌 Repository Overview

```
.
├── main_demo.m                         # Main entry point script running all benchmarks
├── docs/
│   └── paper.pdf                       # Reference paper PDF
├── src/                                # Core FE & Element functions
│   ├── compute_director_vectors.m      # Calculates nodal director triads (Vn, V1, V2)
│   ├── compute_geometry_parameter.m    # Calculates skewness parameter mu (Eq. 55)
│   ├── shape_functions.m               # Bilinear quad shape functions
│   ├── element_stiffness_simplified_mitc4.m  # Simplified MITC4+ element stiffness matrix
│   ├── element_nonlinear_simplified_mitc4.m  # Geometric nonlinear Total Lagrangian formulation
│   ├── solve_linear_fe.m               # Linear FE system solver
│   └── solve_nonlinear_fe.m            # Newton-Raphson nonlinear FE solver
├── utils/                              # Mesh & Visualization utilities
│   ├── generate_quad_mesh.m            # Mesh generator (regular & distorted patterns)
│   ├── plot_mesh.m                     # Plots 3D finite element mesh
│   ├── plot_deformation.m              # Plots initial & deformed shell shapes
│   └── plot_convergence.m              # Plots log-log error convergence curves
└── benchmarks/                         # Benchmark test scripts
    ├── test_basic_tests.m              # Patch test, Zero-energy mode test, Isotropy test
    ├── run_linear_benchmarks.m         # L-shaped structure, Scordelis-Lo roof
    └── run_nonlinear_benchmarks.m      # Large deflection cantilever, Hemispherical shell
```

---

## 🚀 How to Run in MATLAB

1. Open MATLAB.
2. Set the working directory to this project root directory.
3. Run `main_demo`:
   ```matlab
   main_demo
   ```

---

## 🛠 Features Implemented

1. **Simplified Assumed Membrane Strains (Section 3.1):**
   - Implements Matrices $M(r,s)$ and $C$ (Eqs. 41 & 47) using 5 tying points (A–E).
   - Eliminates complex $n_i$ and $m_i$ coefficients from the improved MITC4+ element.
   - Converts strain in element-center covariant coordinate system using transformation $Q(r,s)$.

2. **Assumed Transverse Shear Strains (Section 2.2):**
   - MITC4 transverse shear strain interpolation at 4 tying points to alleviate transverse shear locking.

3. **Geometry-Dependent Gauss Integration (Section 3.2):**
   - Skew angle dependent parameter $\mu = \cos^2\theta$ for reduced sensitivity to mesh distortion.

4. **Linear & Geometric Nonlinear FE Analysis:**
   - Total Lagrangian formulation and incremental-iterative Newton-Raphson solver.

5. **Visualization Suite:**
   - 3D mesh plotting (`plot_mesh`).
   - Deformed geometry with color-mapped displacement magnitude (`plot_deformation`).
   - Log-log mesh refinement convergence curves (`plot_convergence`).

---

## 🧪 Benchmark Problems Included

- **Basic Tests:** Patch Test, Zero-Energy Mode Test, Isotropy Test (Section 4).
- **Linear Analysis:** L-shaped structure, Scordelis-Lo roof, Hyperboloid shell, Parabolic cylinder (Section 5.1–5.4).
- **Nonlinear Analysis:** Large deflection cantilever beam/plate, Slit annular plate, Hemispherical shell with hole (Section 5.5–5.7).
