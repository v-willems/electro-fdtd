import { create, all } from "mathjs";
import type { Complex } from "mathjs";
const math = create(all);

const C = 299792458;
const EPS0 = 8.8541878188e-12;
const MU0 = 4 * Math.PI * 1e-7;

export class Simulator {
  public rows: number;
  public cols: number;
  public epochs: number;

  private Ez: Float32Array;
  private Hx: Float32Array;
  private Hy: Float32Array;

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

  // UPML
  public sigma_x: Float32Array;
  public sigma_y: Float32Array;

  public aEz: Float32Array;
  public bEz: Float32Array;
  public cEz: Float32Array;

  public aHx: Float32Array;
  public bHx: Float32Array;
  public cHx: Float32Array;

  public aHy: Float32Array;
  public bHy: Float32Array;
  public cHy: Float32Array;

  public Ez_acc: Float32Array;
  public Ez_acc_dx: Float32Array;
  public Ez_acc_dy: Float32Array;

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
    this.Hx = new Float32Array(epochs * HxPlane);
    this.Hy = new Float32Array(epochs * HyPlane);

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

    this.sigma_x = new Float32Array(EzPlane);
    this.sigma_x.fill(0);
    this.sigma_y = new Float32Array(EzPlane);
    this.sigma_y.fill(0);

    this.aEz = new Float32Array(EzPlane);
    this.bEz = new Float32Array(EzPlane);
    this.cEz = new Float32Array(EzPlane);
    this.aEz.fill(0);
    this.bEz.fill(0);
    this.cEz.fill(0);

    this.aHx = new Float32Array(HxPlane);
    this.bHx = new Float32Array(HxPlane);
    this.cHx = new Float32Array(HxPlane);
    this.aHx.fill(0);
    this.bHx.fill(0);
    this.cHx.fill(0);

    this.aHy = new Float32Array(HyPlane);
    this.bHy = new Float32Array(HyPlane);
    this.cHy = new Float32Array(HyPlane);
    this.aHy.fill(0);
    this.bHy.fill(0);
    this.cHy.fill(0);

    this.Ez_acc = new Float32Array(EzPlane);
    this.Ez_acc.fill(0);
    this.Ez_acc_dx = new Float32Array(HyPlane);
    this.Ez_acc_dx.fill(0);
    this.Ez_acc_dy = new Float32Array(HxPlane);
    this.Ez_acc_dy.fill(0);
  }

  public buildUPML(left: number, right: number, up: number, down: number) {
    const m = 3;

    for (let j = 0; j < left; j++) {
      const sigma = (EPS0 / (2 * this.dt)) * Math.pow((j + 1) / left, m);
      for (let i = 0; i < this.rows; i++) {
        this.sigma_x[this.idx_g(i, j)] = sigma;
      }
    }
    for (let j = 0; j < right; j++) {
      const sigma = (EPS0 / (2 * this.dt)) * Math.pow((right - j) / right, m);
      const idx_j = this.cols - right + j;
      for (let i = 0; i < this.rows; i++) {
        this.sigma_x[this.idx_g(i, idx_j)] = sigma;
      }
    }
    for (let i = 0; i < up; i++) {
      const sigma = (EPS0 / (2 * this.dt)) * Math.pow((i + 1) / up, m);
      for (let j = 0; j < this.cols; j++) {
        this.sigma_y[this.idx_g(i, j)] = sigma;
      }
    }
    for (let i = 0; i < down; i++) {
      const sigma = (EPS0 / (2 * this.dt)) * Math.pow((down - i) / down, m);
      const idx_i = this.rows - down + i;
      for (let j = 0; j < this.cols; j++) {
        this.sigma_y[this.idx_g(idx_i, j)] = sigma;
      }
    }
    // TODO when iterating over Hx, Hy dont forget to start at 1,
    // end at rows/cols to stay consistent with sigma_x,y based on Ez grid.

    // Fill in coeefficient arrays for Ez
    for (let i = 0; i < this.rows; i++) {
      for (let j = 0; j < this.cols; j++) {
        const idx = this.idx_g(i, j);
        const b =
          1 / this.dt +
          (this.sigma_x[idx] + this.sigma_y[idx]) / (2 * EPS0) +
          (this.sigma_x[idx] * this.sigma_y[idx] * this.dt) / (4 * EPS0 * EPS0);
        const a =
          (1 / this.dt -
            (this.sigma_x[idx] + this.sigma_y[idx]) / (2 * EPS0) -
            (this.sigma_x[idx] * this.sigma_y[idx] * this.dt) /
              (4 * EPS0 * EPS0)) /
          b;
        const c =
          ((1 / b) * (this.sigma_x[idx] * this.sigma_y[idx] * this.dt)) /
          (EPS0 * EPS0);

        const d = (1 / b) * (1 / (EPS0 * this.eps[idx]));
        this.aEz[idx] = a;
        this.bEz[idx] = c;
        this.cEz[idx] = d;
      }
    }

    // Fill in coefficient arrays for Hx
    for (let i = 0; i < this.rows - 1; i++) {
      for (let j = 0; j < this.cols - 1; j++) {
        const idx = this.idx_g(i, j);
        const f = 1 / this.dt + this.sigma_y[idx] / (2 * EPS0);
        const mu_r = MU0 * this.mu_x[this.idx_cHx(i, j)];

        const a = (1 / this.dt - this.sigma_y[idx] / (2 * EPS0)) / f;
        const b = 1 / mu_r / f;
        const c = (this.dt * this.sigma_x[idx]) / (mu_r * EPS0) / f;

        this.aHx[this.idx_cHx(i, j)] = a;
        this.bHx[this.idx_cHx(i, j)] = b;
        this.cHx[this.idx_cHx(i, j)] = c;
      }
    }

    // Fill in coefficient arrays for Hy
    for (let i = 0; i < this.rows - 1; i++) {
      for (let j = 0; j < this.cols - 1; j++) {
        const idx = this.idx_g(i, j);
        const f = 1 / this.dt + this.sigma_x[idx] / (2 * EPS0);
        const mu_r = MU0 * this.mu_y[this.idx_cHy(i, j)];

        const a = (1 / this.dt - this.sigma_x[idx] / (2 * EPS0)) / f;
        const b = 1 / mu_r / f;
        const c = (this.dt * this.sigma_y[idx]) / (mu_r * EPS0) / f;

        this.aHy[this.idx_cHy(i, j)] = a;
        this.bHy[this.idx_cHy(i, j)] = b;
        this.cHy[this.idx_cHy(i, j)] = c;
      }
    }
  }

  // UPML implementation
  private computeNextE_z(k: number, i: number, j: number): number {
    const Ez_k = this.Ez[this.idx_Ez(k - 1, i, j)];

    const idx = this.idx_g(i, j);
    const acc = this.Ez_acc[idx];

    // The grids might be staggered, but the indexing is not.
    // I.e. i+1/2 in E grid <=> i in H grid
    let corr_H = 0;
    if (j == this.huygens_j) {
      // Inside (total field)
      corr_H = -this.Hy_inc_at(k + 0.5, i, this.huygens_j - 1 + 0.5) / this.dx;
    }
    const z_curl_base =
      (this.Hy[this.idx_Hy(k, i, j)] - this.Hy[this.idx_Hy(k, i, j - 1)]) /
        this.dx -
      (this.Hx[this.idx_Hx(k, i, j)] - this.Hx[this.idx_Hx(k, i - 1, j)]) /
        this.dy;

    const z_curl = z_curl_base + corr_H;

    const a = this.aEz[idx];
    const b = this.bEz[idx];
    const c = this.cEz[idx];

    return a * Ez_k + b * acc + c * z_curl;
  }

  /**
   * next := k+1/2
   */
  private computeNextH_x(k: number, i: number, j: number) {
    const H_x_k = this.Hx[this.idx_Hx(k - 1, i, j)];

    const acc_dy = this.Ez_acc_dy[this.idx_cHx(i, j)];
    const x_curl =
      (this.Ez[this.idx_Ez(k - 1, i + 1, j)] -
        this.Ez[this.idx_Ez(k - 1, i, j)]) /
      this.dy;

    const a = this.aHx[this.idx_cHx(i, j)];
    const b = this.bHx[this.idx_cHx(i, j)];
    const c = this.cHx[this.idx_cHx(i, j)];

    return a * H_x_k - b * x_curl - c * acc_dy;
  }

  /**
   * next := k+1/2
   */
  private computeNextH_y(k: number, i: number, j: number) {
    const H_y_k = this.Hy[this.idx_Hy(k - 1, i, j)];

    const acc_dx = this.Ez_acc_dx[this.idx_cHy(i, j)];

    // Correction term for TS/SF boundary
    let corr_E = 0;
    if (j === this.huygens_j - 1) {
      // Outside (scattered field)
      corr_E = -this.Ez_inc_at(k, i, this.huygens_j) / this.dx;
    }
    const y_curl_base =
      (this.Ez[this.idx_Ez(k - 1, i, j + 1)] -
        this.Ez[this.idx_Ez(k - 1, i, j)]) /
      this.dx;

    const y_curl = y_curl_base + corr_E;

    const a = this.aHy[this.idx_cHy(i, j)];
    const b = this.bHy[this.idx_cHy(i, j)];
    const c = this.cHy[this.idx_cHy(i, j)];

    return a * H_y_k + b * y_curl - c * acc_dx;
  }

  private Hy_inc_at(k: number, i: number, j: number) {
    const H = -(1 / Math.sqrt(MU0 / EPS0)) * this.Ez_inc_at(k, i, j);

    return H;
  }

  private Ez_inc_at(k: number, i: number, j: number) {
    const t = k * this.dt - (j * this.dx) / C;
    const m = 1;
    const y_sine = Math.sin((m * Math.PI * (i + 0.5)) / (this.rows - 2));
    const E = this.E0 * y_sine * this.gaussSine(t, this.f0, this.t0, this.tau);

    return E;
  }

  public runSimulation() {
    console.log(`dt / EPS0: ${this.dt / EPS0}`);
    console.log(`dt / MU0: ${this.dt / MU0}`);
    for (let epoch = 1; epoch < this.epochs; epoch++) {
      // update Ez_acc for current epoch based on previous Ez
      for (let i = 0; i < this.rows; i++) {
        for (let j = 0; j < this.cols; j++) {
          this.Ez_acc[this.idx_g(i, j)] +=
            this.Ez[this.idx_Ez(epoch - 1, i, j)];
        }
      }
      for (let i = 0; i < this.rows; i++) {
        for (let j = 0; j < this.cols - 1; j++) {
          const CEy =
            (this.Ez[this.idx_Ez(epoch - 1, i, j + 1)] -
              this.Ez[this.idx_Ez(epoch - 1, i, j)]) /
            this.dx;
          this.Ez_acc_dx[this.idx_cHy(i, j)] += CEy;
        }
      }
      for (let i = 0; i < this.rows - 1; i++) {
        for (let j = 0; j < this.cols; j++) {
          const CEx =
            (this.Ez[this.idx_Ez(epoch - 1, i + 1, j)] -
              this.Ez[this.idx_Ez(epoch - 1, i, j)]) /
            this.dy;
          this.Ez_acc_dy[this.idx_cHx(i, j)] += CEx;
        }
      }

      // We leave out the first and last row as they have no Ez neighbours
      // and we ignore the outer boundary of Ez anyway
      for (let i = 0; i < this.rows - 1; i++) {
        for (let j = 0; j < this.cols; j++) {
          this.Hx[this.idx_Hx(epoch, i, j)] = this.computeNextH_x(epoch, i, j);
          if (!Number.isFinite(this.Hx[this.idx_Hx(epoch, i, j)]))
            throw new Error(`Hx NaN at epoch ${epoch}, i ${i}, j ${j}`);
        }
      }
      // We leave out the first and last column as they have no Ez neighbours
      // and we ignore the outer boundary of Ez anyway
      for (let i = 0; i < this.rows; i++) {
        for (let j = 0; j < this.cols - 1; j++) {
          this.Hy[this.idx_Hy(epoch, i, j)] = this.computeNextH_y(epoch, i, j);
          if (!Number.isFinite(this.Hy[this.idx_Hy(epoch, i, j)]))
            throw new Error(`Hy NaN at epoch ${epoch}, i ${i}, j ${j}`);
        }
      }
      for (let i = 1; i < this.rows - 1; i++) {
        for (let j = 1; j < this.cols - 1; j++) {
          // In rewrite multiply to get branchless code
          if (this.pec[this.idx_g(i, j)] === 0) continue; // skip perfect conductor

          this.Ez[this.idx_Ez(epoch, i, j)] = this.computeNextE_z(epoch, i, j);
          if (!Number.isFinite(this.Ez[this.idx_Ez(epoch, i, j)]))
            throw new Error(`Ez NaN at epoch ${epoch}, i ${i}, j ${j}`);
        }
      }
      // console.log(`Ez_acc ${this.Ez_acc[this.idx_g(50, 20)]}`);
      // console.log(
      //   `Ez ${this.Ez[this.idx_Ez(epoch, 50, 10)]} at epoch ${epoch}`,
      // );
    }
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

  public injectRectMaterial(
    x0: number,
    y0: number,
    w: number,
    h: number,
    mu_r: number,
    eps_r: number,
  ) {
    for (let i = y0; i < y0 + h; i++) {
      for (let j = x0; j < x0 + w; j++) {
        this.eps[this.idx_g(i, j)] = eps_r;
        this.mu_x[this.idx_cHx(i, j)] = mu_r;
        this.mu_y[this.idx_cHy(i, j)] = mu_r;
      }
    }
  }
  public injectRectPEC(x0: number, y0: number, w: number, h: number) {
    for (let i = y0; i < y0 + h; i++) {
      for (let j = x0; j < x0 + w; j++) {
        this.pec[this.idx_g(i, j)] = 0; // mark as perfect conductor
      }
    }
  }

  public getEzStatsFrequencies(i: number, j: number) {
    const Ez_time = new Float32Array(this.epochs);
    for (let n = 0; n < this.epochs; n++) {
      Ez_time[n] = this.Ez[this.idx_Ez(n, i, j)];
    }

    const X = math.fft(Array.from(Ez_time)) as unknown as Complex[];

    const N = this.epochs;
    const half = Math.floor(N / 2) + 1;
    const fs = 1 / this.dt;

    const freqs = new Float32Array(half);
    const mags = new Float32Array(half);

    for (let k = 0; k < half; k++) {
      freqs[k] = (k * fs) / N;

      // magnitude (non-negative)
      let m = Math.hypot(X[k].re, X[k].im);

      // normalize and convert to single-sided amplitude spectrum
      m /= N;
      if (k !== 0 && !(N % 2 === 0 && k === half - 1)) m *= 2;

      mags[k] = m;
    }

    return { freqs, mags };
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
  public idx_g(row: number, col: number) {
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
