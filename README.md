# Electro FDTD

*Electro FDTD* is a browser-based electromagnetic field simulation written in TypeScript, employing the *finite-difference in time-domain* (FDTD) method to solve Maxwell's equations in $TM_z$ mode with a $TE_{10}$ like source profile.

## Demo

https://github.com/user-attachments/assets/d88c6594-31b8-4bb5-9595-63d9288f104a

## What does this mean?

- The **Maxwell equations** are the governing equations for electromagnetic fields, that is the electric field $E$ and magnetic field $H$. The important equations for this project are:
  - $\nabla \times E(t) = -\mu \frac{\partial H(t)}{\partial t}$, which means intuitively: "the rotation of the $E$-field is equvialent to the negative change in time of the $H$-field"
  - $\nabla \times H(t) = \epsilon \frac{\partial E(t)}{\partial t}$, which means intuitively: "the rotation of the $H$-field is equivalent to the positive change in time of the $E$-field"
- Using the **FDTD** method, we solve these equations.
  - We approximate the derivatives in space with the central-difference method $f'(x) \approx \frac{f(x+h)-f(x-h)}{2h}$ which can be derived from the Taylor series. In practice this means that the simulation domain gets discretized into a grid where the cells have a distance of $h_x = h_y$ for the $x$- and $y$-axis respectively. The distances govern the resolution and accuracy of the simulation.
  - We discretize the change in time by timestepping, that is we alternatingly solve the $E$- and $H$-field in tiny steps where the $E$-field in step $k+1$ depends on the $H$-field in step $k+0.5$ and vice versa.
- **$TM_{z}$ mode** in this case means that the $H$-field propagates only on the $x$- and $y$-axis of the 2D-domain, while the $E$-field only changes on the $z$-axis. This value $E_z(t)$ is displayed by the simulator.

### Where can I learn more?

If I captured your interest and inspired you to build a FDTD sim for yourself, I highly recommend you watch this pretty amazing [lecture](https://www.youtube.com/watch?v=3s51YphMp04) :)

## What works
- ahead of time simulation of a fixed amount of epochs
- Yee-grid staggering
- TFSF/Huygens source boundary for x-directed propagation
- PEC objects
- relative permeability and permittivity
- UPML absorbing boundaries

## Goals
- C rewrite compiled to WebAssembly
- SIMD experiments
- live simulation
- a domain editor
- $S_{11}$ analysis

## Non-goals
- moving to 3D
- CPML

## Running locally
```bash
npm install
npm run dev
```

## Example usage
```typescript
const epochs = 2000;  // number of timesteps
const height = 200;   // number of rows of the grid
const width = 400;    // number of columns of the grid
const huygens_x = 70; // where on the x-axis the source should be injected
const sim = new Simulator(epochs, height, width, huygens_x);
sim.injectRectMaterial(150, 1, 40, 198, 1, 10); // x0, y0, width, height, mu_r, eps_r
sim.injectRectPEC(50, 0, 350, 1); // x0, y0, width, height
sim.buildUPML(30, 30, 0, 0); // left, right, upper, lower boundary width

sim.runSimulation();  // runs with simulation synchronously (blocks the main thread)
```
