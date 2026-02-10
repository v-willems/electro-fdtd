const C = 299792458;
const EPS0 = 8.8541878188e-12;
const MU0 = 4 * Math.PI * 1e-7;

export class Simulator {
  public rows: number;
  public cols: number;
  public epochs: number;

  private Ez: Float32Array;
  private Ezx: Float32Array;
  private Ezy: Float32Array;
  private Hx: Float32Array;
  private Hy: Float32Array;

  // Dampening coefficients for PML
  private cEx: Float32Array;
  private cEy: Float32Array;
  private bEx: Float32Array;
  private bEy: Float32Array;

  private cHx: Float32Array;
  private cHy: Float32Array;
  private bHx: Float32Array;
  private bHy: Float32Array;

  // derivative discretizations
  private dy: number;
  private dx: number;
  private dt: number;

  public mu_x: Float32Array;
  public mu_y: Float32Array;
  public eps: Float32Array;

  // PEC materials
  public pec: Uint8Array;

  // Huygens surface for TFSF
  huygens_j: number;
  // Center of the gaussian source
  public E0: number; // amplitude
  private f0: number; // example frequency
  private t0: number; // center time
  private tau: number; // pulse width

  private gaussSineSourcePositions: [number, number][] = [];
  private sineSourcePositions: [number, number][] = [];

  constructor(
    epochs: number,
    rows: number,
    cols: number,
    huygens_vertical: number,
  ) {
    this.epochs = epochs;
    this.rows = rows;
    this.cols = cols;

    const EzPlane = rows * cols;
    const HxPlane = rows * (cols + 2);
    const HyPlane = (rows + 2) * cols;

    this.Ez = new Float32Array(epochs * EzPlane);
    this.Ezx = new Float32Array(epochs * EzPlane);
    this.Ezy = new Float32Array(epochs * EzPlane);
    this.Hx = new Float32Array(epochs * HxPlane);
    this.Hy = new Float32Array(epochs * HyPlane);

    // PML
    this.cEx = new Float32Array(EzPlane);
    this.cEy = new Float32Array(EzPlane);
    this.bEx = new Float32Array(EzPlane);
    this.bEy = new Float32Array(EzPlane);

    this.cHx = new Float32Array(HxPlane);
    this.cHy = new Float32Array(HyPlane);
    this.bHx = new Float32Array(HxPlane);
    this.bHy = new Float32Array(HyPlane);

    this.dx = 1e-3;
    this.dy = 1e-3;
    this.dt =
      0.99 *
      (1 / (C * Math.sqrt(1 / (this.dx * this.dx) + 1 / (this.dy * this.dy))));

    // materials (prefer spatial-only)
    this.eps = new Float32Array(EzPlane);
    this.eps.fill(1.0);
    this.mu_x = new Float32Array(HxPlane);
    this.mu_x.fill(1.0);
    this.mu_y = new Float32Array(HyPlane);
    this.mu_y.fill(1.0);

    this.pec = new Uint8Array(EzPlane);
    this.pec.fill(1); // 1 means not perfect conductor, 0 means perfect conductor

    // TFSF boundary
    this.huygens_j = huygens_vertical;

    this.E0 = 1.0;
    this.f0 = 1e10;
    this.t0 = 10 * this.dt;
    this.tau = 30 * this.dt;
  }

  // Setup the PML coefficient and time update information
  // PML region is forced to be air - prevent touching other material
  public buildCoefficientArrays(
    left: number,
    right: number,
    up: number,
    down: number,
  ) {
    const m = 3;

    const R = 1e-8; // Target reflectance

    this.cEx = this.eps.map((v, _) => this.dt / (EPS0 * v));
    this.cEy = this.eps.map((v, _) => this.dt / (EPS0 * v));
    this.cHx = this.mu_x.map((v, _) => this.dt / (MU0 * v));
    this.cHy = this.mu_y.map((v, _) => this.dt / (MU0 * v));

    this.bEx.fill(1.0);
    this.bEy.fill(1.0);
    this.bHx.fill(1.0);
    this.bHy.fill(1.0);

    // For now use free space only
    // If that changes EPS0 has to be changed to a EPS0 * eps_r
    const sigma_max_left =
      (-((m + 1) * EPS0 * C) / (2 * this.dx * left)) * Math.log(R);
    const sigma_max_right =
      (-((m + 1) * EPS0 * C) / (2 * this.dx * right)) * Math.log(R);
    const sigma_max_up =
      (-((m + 1) * EPS0 * C) / (2 * this.dy * up)) * Math.log(R);
    const sigma_max_down =
      (-((m + 1) * EPS0 * C) / (2 * this.dy * down)) * Math.log(R);

    for (let j = 0; j < left; j++) {
      const s = Math.pow((left - (j + 1)) / left, m);
      for (let i = 0; i < this.rows; i++) {
        const sigmaX = s * sigma_max_left;
        const sigmaX_m = (MU0 / EPS0) * sigmaX; // magnetic field

        this.cEx[this.idx_g(i, j)] =
          this.dt / EPS0 / (1 + (sigmaX * this.dt) / (2 * EPS0));
        this.bEx[this.idx_g(i, j)] =
          (1 - (sigmaX * this.dt) / (2 * EPS0)) /
          (1 + (sigmaX * this.dt) / (2 * EPS0));
        this.bHx[this.idx_cHx(i, j)] =
          (1 - (sigmaX_m * this.dt) / (2 * MU0)) /
          (1 + (sigmaX_m * this.dt) / (2 * MU0));

        this.cHx[this.idx_cHx(i, j)] =
          this.dt / MU0 / (1 + (sigmaX_m * this.dt) / (2 * MU0));
      }
    }

    for (let j = this.cols + 2 - right; j < this.cols + 2; j++) {
      const s = Math.pow((j + 1 - (this.cols + 2 - right)) / right, m);

      for (let i = 0; i < this.rows; i++) {
        const sigmaX = s * sigma_max_right;
        const sigmaX_m = (MU0 / EPS0) * sigmaX; // magnetic field
        this.cEx[this.idx_g(i, j)] =
          this.dt / EPS0 / (1 + (sigmaX * this.dt) / (2 * EPS0));
        this.bEx[this.idx_g(i, j)] =
          (1 - (sigmaX * this.dt) / (2 * EPS0)) /
          (1 + (sigmaX * this.dt) / (2 * EPS0));
        this.bHx[this.idx_cHx(i, j)] =
          (1 - (sigmaX_m * this.dt) / (2 * MU0)) /
          (1 + (sigmaX_m * this.dt) / (2 * MU0));

        this.cHx[this.idx_cHx(i, j)] =
          this.dt / MU0 / (1 + (sigmaX_m * this.dt) / (2 * MU0));
      }
    }

    for (let i = 0; i < up; i++) {
      const s = Math.pow((up - (i + 1)) / up, m);

      for (let j = 0; j < this.cols; j++) {
        const sigmaY = s * sigma_max_up;
        const sigmaY_m = (MU0 / EPS0) * sigmaY; // magnetic field

        this.cEy[this.idx_g(i, j)] =
          this.dt / EPS0 / (1 + (sigmaY * this.dt) / (2 * EPS0));
        this.bEy[this.idx_g(i, j)] =
          (1 - (sigmaY * this.dt) / (2 * EPS0)) /
          (1 + (sigmaY * this.dt) / (2 * EPS0));

        this.bHy[this.idx_cHy(i, j)] =
          (1 - (sigmaY_m * this.dt) / (2 * MU0)) /
          (1 + (sigmaY_m * this.dt) / (2 * MU0));

        this.cHy[this.idx_cHy(i, j)] =
          this.dt / MU0 / (1 + (sigmaY_m * this.dt) / (2 * MU0));
      }
    }

    for (let i = this.rows + 2 - down; i < this.rows + 2; i++) {
      const s = Math.pow((i + 1 - (this.rows + 2 - down)) / down, m);

      for (let j = 0; j < this.cols; j++) {
        const sigmaY = s * sigma_max_down;
        const sigmaY_m = (MU0 / EPS0) * sigmaY; // magnetic field

        this.cEy[this.idx_g(i, j)] =
          this.dt / EPS0 / (1 + (sigmaY * this.dt) / (2 * EPS0));
        this.bEy[this.idx_g(i, j)] =
          (1 - (sigmaY * this.dt) / (2 * EPS0)) /
          (1 + (sigmaY * this.dt) / (2 * EPS0));

        this.bHy[this.idx_cHy(i, j)] =
          (1 - (sigmaY_m * this.dt) / (2 * MU0)) /
          (1 + (sigmaY_m * this.dt) / (2 * MU0));

        this.cHy[this.idx_cHy(i, j)] =
          this.dt / MU0 / (1 + (sigmaY_m * this.dt) / (2 * MU0));
      }
    }
  }

  private Hy_inc_at(k: number, i: number, j: number) {
    const H = -(1 / Math.sqrt(MU0 / EPS0)) * this.Ez_inc_at(k - 0.5, i, j);

    return H;
  }

  private Ez_inc_at(k: number, i: number, j: number) {
    const t = k * this.dt - (j * this.dx) / C;
    const m = 1;
    const y_sine = Math.sin((m * Math.PI * (i + 0.5)) / this.rows);
    const E = this.E0 * y_sine * this.gaussSine(t, this.f0, this.t0, this.tau);

    return E;
  }

  /**
   * next := k+1/2
   */
  private computeNextH_x(k: number, i: number, j: number) {
    const H_x_k = this.Hx[this.idx_Hx(k - 1, i, j)];

    const x_curl =
      (this.Ez[this.idx_Ez(k - 1, i + 1, j)] -
        this.Ez[this.idx_Ez(k - 1, i, j)]) /
      this.dy;

    const c = this.cHx[this.idx_cHx(i, j)];
    const b = this.bHx[this.idx_cHx(i, j)];

    return b * H_x_k - c * x_curl;
  }

  /**
   * next := k+1/2
   */
  private computeNextH_y(k: number, i: number, j: number) {
    const H_y_k = this.Hy[this.idx_Hy(k - 1, i, j)];

    // Correction term for TS/SF boundary
    let corr_E = 0;
    if (j === this.huygens_j - 1) {
      // Outside (scattered field)
      corr_E = -this.Ez_inc_at(k, i, j) / this.dx;
    }
    const y_curl =
      (this.Ez[this.idx_Ez(k - 1, i, j + 1)] -
        this.Ez[this.idx_Ez(k - 1, i, j)]) /
      this.dx;
    const c = this.cHy[this.idx_cHy(i, j)];
    const b = this.bHy[this.idx_cHy(i, j)];

    return b * H_y_k + c * (y_curl + corr_E);
  }

  /**
   * next := k+1
   */
  private computeNextE_zx(k: number, i: number, j: number) {
    const Ezx_k = this.Ezx[this.idx_Ez(k - 1, i, j)];

    // The grids might be staggered, but the indexing is not.
    // I.e. i+1/2 in E grid <=> i in H grid
    let corr_H = 0;
    if (j == this.huygens_j) {
      // Outside (scattered field)
      corr_H = -this.Hy_inc_at(k, i, j - 0.5) / this.dx;
    }
    const z_curl =
      (this.Hy[this.idx_Hy(k, i, j)] - this.Hy[this.idx_Hy(k, i, j - 1)]) /
      this.dx;

    // maybe idx_g is the wrong index function here
    const c = this.cEx[this.idx_g(i, j)];
    const b = this.bEx[this.idx_g(i, j)];

    return b * Ezx_k + c * (z_curl + corr_H);
  }

  private computeNextE_zy(k: number, i: number, j: number) {
    const Ezy_k = this.Ezy[this.idx_Ez(k - 1, i, j)];

    const z_curl =
      (this.Hx[this.idx_Hx(k, i, j)] - this.Hx[this.idx_Hx(k, i - 1, j)]) /
      this.dy;
    const c = this.cEy[this.idx_g(i, j)];
    const b = this.bEy[this.idx_g(i, j)];

    return b * Ezy_k - c * z_curl;
  }

  // Soft source
  gaussian(t: number, t0: number, tau: number) {
    const x = (t - t0) / tau;
    return Math.exp(-x * x);
  }

  gaussSine(t: number, f0: number, t0: number, tau: number) {
    return Math.sin(2 * Math.PI * f0 * t) * this.gaussian(t, t0, tau);
  }

  sine(t: number, f0: number) {
    return Math.sin(2 * Math.PI * f0 * t);
  }

  injectSoftEz(epoch: number) {
    const t = epoch * this.dt;
    const gs = this.E0 * this.gaussSine(t, this.f0, this.t0, this.tau);

    for (const [i, j] of this.gaussSineSourcePositions) {
      this.Ezx[this.idx_Ez(epoch, i, j)] += 0.5 * gs;
      this.Ezy[this.idx_Ez(epoch, i, j)] += 0.5 * gs;
    }

    const s = this.E0 * this.sine(t, this.f0);
    for (const [i, j] of this.sineSourcePositions) {
      this.Ezx[this.idx_Ez(epoch, i, j)] += 0.5 * s;
      this.Ezy[this.idx_Ez(epoch, i, j)] += 0.5 * s;
    }
  }

  public injectGaussSineSource(x: number, y: number) {
    this.gaussSineSourcePositions.push([x, y]);
  }
  public injectSineSource(x: number, y: number) {
    this.sineSourcePositions.push([x, y]);
  }

  public injectRectMaterial(
    x0: number,
    y0: number,
    w: number,
    h: number,
    mu_r: number,
    eps_r: number,
  ) {
    for (let i = y0; i < y0 + w; i++) {
      for (let j = x0; j < x0 + h; j++) {
        this.eps[this.idx_g(i, j)] = eps_r;
        this.mu_x[this.idx_cHx(i, j)] = mu_r;
        this.mu_y[this.idx_cHy(i, j)] = mu_r;
      }
    }
  }
  public injectRectPEC(x0: number, y0: number, w: number, h: number) {
    for (let i = y0; i < y0 + w; i++) {
      for (let j = x0; j < x0 + h; j++) {
        this.pec[this.idx_g(i, j)] = 0; // mark as perfect conductor
      }
    }
  }

  public runSimulation() {
    for (let epoch = 1; epoch < this.epochs; epoch++) {
      for (let i = 0; i < this.rows; i++) {
        for (let j = 0; j < this.cols; j++) {
          this.Hx[this.idx_Hx(epoch, i, j)] = this.computeNextH_x(epoch, i, j);
        }
      }
      for (let i = 0; i < this.rows; i++) {
        for (let j = 0; j < this.cols; j++) {
          this.Hy[this.idx_Hy(epoch, i, j)] = this.computeNextH_y(epoch, i, j);
        }
      }
      for (let i = 1; i < this.rows - 1; i++) {
        for (let j = 1; j < this.cols - 1; j++) {
          // In rewrite multiply to get branchless code
          if (this.pec[this.idx_g(i, j)] === 0) continue; // skip perfect conductor

          const Ezx = this.computeNextE_zx(epoch, i, j);
          this.Ezx[this.idx_Ez(epoch, i, j)] = Ezx;
          const Ezy = this.computeNextE_zy(epoch, i, j);
          this.Ezy[this.idx_Ez(epoch, i, j)] = Ezy;
          this.Ez[this.idx_Ez(epoch, i, j)] = Ezx + Ezy;
        }
      }
      this.injectSoftEz(epoch);
    }
  }

  public getE_z(epoch: number, i: number, j: number) {
    return this.Ez[this.idx_Ez(epoch, i, j)];
  }

  public idx_Ez(t: number, i: number, j: number): number {
    // i: 0..rows-1, j: 0..cols-1
    return t * (this.rows * this.cols) + i * this.cols + j;
  }

  private idx_Hx(t: number, i: number, j: number) {
    // Hx has column padding: j stored as (j+1)
    // i: 0..rows-1, j: -1..cols-1 is addressable
    const strideJ = this.cols + 2;
    return t * (this.rows * strideJ) + i * strideJ + (j + 1);
  }

  private idx_Hy(t: number, i: number, j: number) {
    // Hy has row padding: i stored as (i+1)
    // i: -1..rows-1 is addressable, j: 0..cols-1
    const strideJ = this.cols;
    return t * ((this.rows + 2) * strideJ) + (i + 1) * strideJ + j;
  }

  // Access cell information
  private idx_g(row: number, col: number) {
    return row * this.cols + col;
  }

  // Index coefficient arrays
  private idx_cHx(i: number, j: number) {
    return i * (this.cols + 2) + (j + 1);
  }
  private idx_cHy(i: number, j: number) {
    return (i + 1) * this.cols + j;
  }
  public idx_cEz(i: number, j: number) {
    return i * this.cols + j;
  }
}
