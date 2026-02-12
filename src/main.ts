import { Visualizer } from "./Visualizer.ts";
import { Simulator } from "./Simulator.ts";
import {
  Chart,
  type ChartConfiguration,
  type ScatterDataPoint,
  LineController,
  LineElement,
  PointElement,
  LinearScale,
  Title,
  Tooltip,
  Legend,
} from "chart.js";

Chart.register(
  LineController,
  LineElement,
  PointElement,
  LinearScale,
  Title,
  Tooltip,
  Legend,
);

const sim_canvas = document.getElementById("grid") as HTMLCanvasElement;
const frameCounter = document.getElementById("framecounter") as HTMLDivElement;
const sim_ctx = sim_canvas.getContext("2d");
if (!sim_ctx) throw new Error("2D context not available");

const chart_canvas = document.getElementById("chart") as HTMLCanvasElement;

const height = sim_canvas.height;
const width = sim_canvas.width;

sim_ctx.fillStyle = "lightgray";
sim_ctx.fillRect(0, 0, width, height);

const startBtn = document.getElementById(
  "startBtn",
) as HTMLButtonElement | null;
if (!startBtn) throw new Error("startBtn not found");

const epochs = 1200;
const sim = new Simulator(epochs, 200, 400, 70);
// simulator.injectRectMaterial(150, 1, 40, 198, 1, 10);
// simulator.injectRectPEC(50, 0, 350, 1);
// simulator.injectRectPEC(50, 198, 350, 1);
// simulator.injectRectPEC(150, 1, 40, 198);
sim.buildUPML(30, 30, 0, 0);
// console.log(`aEz ${sim.aEz[sim.idx_g(50, 50)]}`);
// console.log(`bEz ${sim.bEz[sim.idx_g(50, 50)]}`);
// console.log(`cEz ${sim.cEz[sim.idx_g(50, 50)]}`);
// console.log(`aHy ${sim.aHy[sim.idx_g(50, 50)]}`);
// console.log(`bHy ${sim.bHy[sim.idx_g(50, 50)]}`);
// console.log(`cHy ${sim.cHy[sim.idx_g(50, 50)]}`);
console.log(`aHx ${sim.aHx[sim.idx_g(50, 50)]}`);
console.log(`bHx ${sim.bHx[sim.idx_g(50, 50)]}`);
console.log(`cHx ${sim.cHx[sim.idx_g(50, 50)]}`);
sim.runSimulation();

const { freqs, mags } = sim.getEzStatsFrequencies(40, 50);

const points: ScatterDataPoint[] = Array.from(freqs).map((f, idx) => ({
  x: f,
  y: mags[idx],
}));

const config: ChartConfiguration<"line", ScatterDataPoint[], number> = {
  type: "line",
  data: {
    datasets: [
      {
        label: "FFT magnitude",
        data: points,
        parsing: false as const, // <- key change
        pointRadius: 0,
        borderWidth: 1,
      },
    ],
  },
  options: {
    animation: false,
    scales: {
      x: {
        type: "linear",
        title: { display: true, text: "Frequency (Hz)" },
        ticks: {
          callback: (value) => (Number(value) / 1e9).toString(),
        },
        max: 30e9,
        min: 0,
      },
      y: {
        title: { display: true, text: "Magnitude" },
      },
    },
  },
};

const chart = new Chart(chart_canvas, config);

const viz = new Visualizer(sim_ctx, {
  x0: 0,
  y0: 0,
  gridSpacingV: 1,
  gridSpacingH: 1,
  pixelSize: 1,
  clearStyle: "white",
});

let minE = +sim.E0,
  maxE = -sim.E0;
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

  stopAnim = viz.animateEField(sim, colorRange, {
    onFrame: (e) => (frameCounter.textContent = String(e)),
    loop: true, // set true to repeat
    maxFps: 60, // optional throttle
    drawMaterialsEachFrame: true,
  });
});
