// Visualizer.ts
type Stop = { t: number; r: number; g: number; b: number };

type VisualizerOpts = {
  x0: number;
  y0: number;
  gridSpacingV: number;
  gridSpacingH: number;
  pixelSize?: number;
  lutSize?: number;
  stops?: Stop[];
  clearStyle?: string;

  // material overlay colors
  pecColor?: string; // default "black"
  dielectricColor?: string; // default "rgb(210 210 210)"
  dielectricThreshold?: number; // eps_r != 1 test tolerance, default 1e-6
};

export class Visualizer {
  private ctx: CanvasRenderingContext2D;

  private x0: number;
  private y0: number;
  private gridSpacingV: number;
  private gridSpacingH: number;
  private pixelSize: number;

  private lut: Uint8ClampedArray; // RGBRGB...
  private lutN: number;

  private clearStyle: string;

  private pecColor: string;
  private dielectricColor: string;
  private dielectricThreshold: number;

  static readonly blueWhiteRed: Stop[] = [
    { t: 0.0, r: 20, g: 60, b: 220 },
    { t: 0.5, r: 255, g: 255, b: 255 },
    { t: 1.0, r: 220, g: 40, b: 40 },
  ];

  constructor(ctx: CanvasRenderingContext2D, opts: VisualizerOpts) {
    this.ctx = ctx;

    this.x0 = opts.x0;
    this.y0 = opts.y0;
    this.gridSpacingV = opts.gridSpacingV;
    this.gridSpacingH = opts.gridSpacingH;

    this.pixelSize = opts.pixelSize ?? 1;

    const stops = opts.stops ?? Visualizer.blueWhiteRed;
    this.lutN = opts.lutSize ?? 256;
    this.lut = Visualizer.buildLUT(stops, this.lutN);

    this.clearStyle = opts.clearStyle ?? "lightgray";

    this.pecColor = opts.pecColor ?? "black";
    this.dielectricColor = opts.dielectricColor ?? "rgb(210 210 210)";
    this.dielectricThreshold = opts.dielectricThreshold ?? 1e-6;
  }

  getColorRangeFunc(maxVal: number, minVal: number): (E: number) => number {
    const span = maxVal - minVal;
    if (span === 0) return () => 0.5;
    return (E: number) => Visualizer.clamp((E - minVal) / span, 0, 1);
  }

  clear() {
    const ctx = this.ctx;
    const c = ctx.canvas;
    ctx.fillStyle = this.clearStyle;
    ctx.fillRect(0, 0, c.width, c.height);
  }

  // Draw material mask as background overlay.
  // - PEC (pec[idx]==0) drawn black
  // - eps_r != 1 drawn light gray (unless PEC)
  drawMaterials(sim: {
    rows: number;
    cols: number;
    eps: Float32Array;
    pec: Uint8Array;
  }) {
    const ctx = this.ctx;
    const ps = this.pixelSize;

    const rows = sim.rows;
    const cols = sim.cols;

    // eps is EzPlane (rows*cols), pec is EzPlane
    for (let i = 0; i < rows; i++) {
      const x = this.x0 + i * this.gridSpacingV;
      const base = i * cols;
      for (let j = 0; j < cols; j++) {
        const idx = base + j;
        const y = this.y0 + j * this.gridSpacingH;

        // pec: 1 means not PEC, 0 means PEC
        if (sim.pec[idx] === 0) {
          ctx.fillStyle = this.pecColor;
          ctx.fillRect(x, y, ps, ps);
          continue;
        }

        const er = sim.eps[idx];
        if (Math.abs(er - 1.0) > this.dielectricThreshold) {
          ctx.fillStyle = this.dielectricColor;
          ctx.fillRect(x, y, ps, ps);
        }
      }
    }
  }

  drawEField(
    epoch: number,
    sim: {
      rows: number;
      cols: number;
      getE_z(t: number, i: number, j: number): number;
      idx_Ez(t: number, i: number, j: number): number;
      idx_cEz(i: number, j: number): number;
      // optional material overlays (if present, PEC overrides field color)
      eps?: Float32Array;
      pec?: Uint8Array;
    },
    colorRange: (E: number) => number,
  ) {
    const ctx = this.ctx;
    const ps = this.pixelSize;

    const hasMask = !!sim.pec && !!sim.eps;

    for (let i = 0; i < sim.rows; i++) {
      const y = this.y0 + i * this.gridSpacingH;
      for (let j = 0; j < sim.cols; j++) {
        const x = this.x0 + j * this.gridSpacingV;

        const idx = sim.idx_cEz(i, j);

        if (hasMask) {
          // PEC overrides everything
          if (sim.pec![idx] === 0) {
            ctx.fillStyle = this.pecColor;
            ctx.fillRect(x, y, ps, ps);
            continue;
          }

          // eps_r overlay (optional: keep it as background by drawing it first instead)
          const er = sim.eps![idx];
          if (Math.abs(er - 1.0) > this.dielectricThreshold) {
            ctx.fillStyle = this.dielectricColor;
            ctx.fillRect(x, y, ps, ps);
            // fall through to draw field on top if you prefer; comment next line to keep only gray
          }
        }

        const e = sim.getE_z(epoch, i, j);
        const t = colorRange(e);
        const lutIdx = this.lutIndex(t);

        const r = this.lut[lutIdx + 0];
        const g = this.lut[lutIdx + 1];
        const b = this.lut[lutIdx + 2];

        ctx.fillStyle = `rgb(${r} ${g} ${b} / 50%)`;
        ctx.fillRect(x, y, ps, ps);
      }
    }
  }

  animateEField(
    sim: {
      epochs: number;
      rows: number;
      cols: number;
      getE_z(t: number, i: number, j: number): number;
      idx_Ez(t: number, i: number, j: number): number;
      idx_cEz(i: number, j: number): number;
      eps?: Float32Array;
      pec?: Uint8Array;
    },
    colorRange: (E: number) => number,
    opts?: {
      onFrame?: (epoch: number) => void;
      loop?: boolean;
      startEpoch?: number;
      maxFps?: number;
      drawMaterialsEachFrame?: boolean; // default false; if true, draws masks every frame
    },
  ): () => void {
    let running = true;
    let epoch = opts?.startEpoch ?? 0;
    const loop = opts?.loop ?? false;

    const maxFps = opts?.maxFps ?? 0;
    const minDt = maxFps > 0 ? 1000 / maxFps : 0;
    let lastT = 0;

    const tick = (t: number) => {
      if (!running) return;

      if (minDt > 0 && t - lastT < minDt) {
        requestAnimationFrame(tick);
        return;
      }
      lastT = t;

      this.clear();

      // If materials are static, draw once outside the animation loop instead (faster).

      this.drawEField(epoch, sim, colorRange);
      opts?.onFrame?.(epoch);

      epoch++;
      if (epoch >= sim.epochs) {
        if (loop) epoch = 0;
        else {
          running = false;
          return;
        }
      }

      requestAnimationFrame(tick);
    };

    requestAnimationFrame(tick);

    return () => {
      running = false;
    };
  }

  private lutIndex(t01: number) {
    const t = Visualizer.clamp(t01, 0, 1);
    const i = Math.min(this.lutN - 1, (t * (this.lutN - 1)) | 0);
    return i * 3;
  }

  private static clamp(val: number, lower: number, upper: number) {
    return Math.max(lower, Math.min(val, upper));
  }

  private static lerp(a: number, b: number, t: number) {
    return a + (b - a) * t;
  }

  private static clamp01(x: number) {
    return Math.max(0, Math.min(1, x));
  }

  private static colorFromStops(stops: Stop[], t: number) {
    t = Visualizer.clamp01(t);

    let k = 0;
    while (k < stops.length - 2 && t > stops[k + 1].t) k++;

    const a = stops[k];
    const b = stops[k + 1];
    const span = b.t - a.t || 1;
    const u = (t - a.t) / span;

    return {
      r: Math.round(Visualizer.lerp(a.r, b.r, u)),
      g: Math.round(Visualizer.lerp(a.g, b.g, u)),
      b: Math.round(Visualizer.lerp(a.b, b.b, u)),
    };
  }

  private static buildLUT(stops: Stop[], n = 256) {
    const lut = new Uint8ClampedArray(n * 3);
    for (let i = 0; i < n; i++) {
      const t = i / (n - 1);
      const c = Visualizer.colorFromStops(stops, t);
      lut[3 * i + 0] = c.r;
      lut[3 * i + 1] = c.g;
      lut[3 * i + 2] = c.b;
    }
    return lut;
  }
}
