import { Visualizer } from "./Visualizer.ts";
import { Simulator } from "./Simulator.ts";
import { Chart } from "chart.js/auto";

const canvas = document.getElementById("grid") as HTMLCanvasElement;
const frameCounter = document.getElementById("framecounter") as HTMLDivElement;
const ctx = canvas.getContext("2d");
if (!ctx) throw new Error("2D context not available");

const height = canvas.height;
const width = canvas.width;

ctx.fillStyle = "lightgray";
ctx.fillRect(0, 0, width, height);

const startBtn = document.getElementById(
  "startBtn",
) as HTMLButtonElement | null;
if (!startBtn) throw new Error("startBtn not found");

const epochs = 200;
const simulator = new Simulator(epochs, 200, 400, 50);
// simulator.injectGaussSineSource(100, 100);
simulator.injectRectMaterial(150, 80, 40, 40, 2.8, 2.8);
simulator.injectRectPEC(50, 1, 1, 350);
simulator.injectRectPEC(50, 198, 1, 350);
simulator.buildCoefficientArrays(30, 30, 0, 0);
simulator.runSimulation();

const viz = new Visualizer(ctx, {
  x0: 0,
  y0: 0,
  gridSpacingV: 1,
  gridSpacingH: 1,
  pixelSize: 1,
  clearStyle: "white",
});

let minE = +simulator.E0,
  maxE = -simulator.E0;
// for (let t = 0; t < epochs; t++) {
//   for (let i = 0; i < simulator.rows; i++) {
//     for (let j = 0; j < simulator.cols; j++) {
//       const e = simulator.getE_z(t, i, j);
//       if (e < minE) minE = e;
//       if (e > maxE) maxE = e;
//     }
//   }
// }

const colorRange = viz.getColorRangeFunc(maxE, minE);

let stopAnim: (() => void) | null = null;

startBtn.addEventListener("click", () => {
  if (stopAnim) return; // already running

  stopAnim = viz.animateEField(simulator, colorRange, {
    onFrame: (e) => (frameCounter.textContent = String(e)),
    loop: true, // set true to repeat
    maxFps: 60, // optional throttle
    drawMaterialsEachFrame: true,
  });
});
