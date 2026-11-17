/*
  The molecules of one phase in a thin layer of water, for panel (b).

  The number of molecules is the real one: the number density of the phase times the volume of the layer.
  The motion is schematic. In a fluid the molecules are disks that move in the plane and repel one another when
  they touch, at speeds that follow the temperature. In ice the molecules sit on a honeycomb lattice, in which
  each molecule has three neighbors at 0.28 nm, and vibrate about their sites.

  When the layer shown would hold less than FEW molecules, a larger piece of the layer is followed and the
  picture shows one corner of it, so that a molecule crosses the picture as often as it should.
*/
export const THICK = 0.33;     // nm, thickness of the layer, about one molecule
export const RADIUS = 0.115;   // nm, radius of the disk drawn for a molecule
const SIGMA = 0.26;            // nm, two molecules closer than this repel
const STIFF = 300;             // 1/ps^2, strength of the repulsion divided by the mass of a molecule
const DT = 0.004;              // ps, time step
const FEW = 8;
const KT_M = 4.615e-4;         // nm^2/ps^2 per K, k_B T / m for a water molecule

const hash = (i, j, b, k) => { const s = Math.sin(i * 127.1 + j * 311.7 + b * 74.7 + k * 19.3) * 43758.5453; return s - Math.floor(s); };

export class Sample {
  constructor(random) {
    this.random = random;
    this.kind = null; this.fluid = false; this.time = 0; this.mean = 0;
    this.x = []; this.y = []; this.vx = []; this.vy = []; this.ax = []; this.ay = [];
    this.lx = 1; this.ly = 1; this.w = 1; this.h = 1; this.T = 300;
    this.sites = []; this.a = 0;
  }

  // kind: 'solid', 'liquid', 'gas', ...; number: molecules per nm3; T in K; w, h: size of the picture in nm
  set(kind, number, T, w, h) {
    const areal = number * THICK;                     // molecules per nm2 of the layer
    this.mean = areal * w * h;
    this.T = T;
    if (kind === 'solid') {
      const a = Math.sqrt(4 / (3 * Math.sqrt(3) * areal));          // distance between neighbors on the honeycomb lattice
      if (this.fluid || this.kind !== 'solid' || Math.abs(a - this.a) > 1e-7 || w !== this.w || h !== this.h) this.lattice(a, w, h);
      this.fluid = false;
    } else {
      const scale = this.mean >= FEW ? 1 : Math.sqrt(FEW / this.mean);
      const count = this.mean >= FEW ? Math.round(this.mean) : FEW;
      const lx = w * scale, ly = h * scale;
      if (!this.fluid) { this.fluid = true; this.fill(count, lx, ly); }
      else if (count !== this.x.length || lx !== this.lx || ly !== this.ly) this.adjust(count, lx, ly);
      this.fluid = true;
    }
    this.kind = kind; this.w = w; this.h = h;
  }

  lattice(a, w, h) {
    this.a = a;
    const sites = [], cx = w / 2, cy = h / 2, r3 = Math.sqrt(3);
    const nj = Math.ceil((h / 2 + 1) / (1.5 * a)) + 1, ni = Math.ceil((w / 2 + 1) / (r3 * a)) + nj;
    for (let j = -nj; j <= nj; j++) for (let i = -ni; i <= ni; i++) for (let b = 0; b < 2; b++) {
      const x = cx + (i + j / 2) * r3 * a, y = cy + j * 1.5 * a + b * a - a / 2;
      if (x < -RADIUS || x > w + RADIUS || y < -RADIUS || y > h + RADIUS) continue;
      sites.push({ x, y, sx: 0.6 + 0.8 * hash(i, j, b, 1), sy: 0.6 + 0.8 * hash(i, j, b, 2),
                   fx: 2 + 2 * hash(i, j, b, 3), fy: 2 + 2 * hash(i, j, b, 4), px: 6.2832 * hash(i, j, b, 5), py: 6.2832 * hash(i, j, b, 6) });
    }
    this.sites = sites;
  }

  normal() { return Math.sqrt(-2 * Math.log(1 - this.random())) * Math.cos(2 * Math.PI * this.random()); }

  add(x, y) {
    const s = Math.sqrt(KT_M * this.T);
    this.x.push(x); this.y.push(y); this.vx.push(s * this.normal()); this.vy.push(s * this.normal()); this.ax.push(0); this.ay.push(0);
  }

  fill(count, lx, ly) {                               // new fluid: the molecules start on a grid with some disorder
    this.x = []; this.y = []; this.vx = []; this.vy = []; this.ax = []; this.ay = [];
    this.lx = lx; this.ly = ly;
    const nx = Math.max(1, Math.ceil(Math.sqrt(count * lx / ly))), ny = Math.ceil(count / nx);
    const cells = [];
    for (let k = 0; k < nx * ny; k++) cells.push(k);
    for (let k = cells.length - 1; k > 0; k--) { const m = Math.floor(this.random() * (k + 1)); [cells[k], cells[m]] = [cells[m], cells[k]]; }
    for (let k = 0; k < count; k++) {
      const i = cells[k] % nx, j = Math.floor(cells[k] / nx);
      this.add((i + 0.2 + 0.6 * this.random()) / nx * lx, (j + 0.2 + 0.6 * this.random()) / ny * ly);
    }
    this.relax(80);
    for (let k = 0; k < 60; k++) this.step(0.064);    // 4 ps of motion, so that the first picture shows no trace of the grid
  }

  adjust(count, lx, ly) {                             // the density or the size of the picture changed
    const fx = lx / this.lx, fy = ly / this.ly;
    for (let i = 0; i < this.x.length; i++) { this.x[i] *= fx; this.y[i] *= fy; }
    this.lx = lx; this.ly = ly;
    while (this.x.length > count) {
      const i = Math.floor(this.random() * this.x.length);
      for (const list of [this.x, this.y, this.vx, this.vy, this.ax, this.ay]) { list[i] = list[list.length - 1]; list.pop(); }
    }
    const added = this.x.length < count;
    while (this.x.length < count) this.add(this.random() * lx, this.random() * ly);
    if (added || fx < 1 || fy < 1) this.relax(25);
  }

  forces() {
    const { x, y, ax, ay, lx, ly } = this, n = x.length, hx = lx / 2, hy = ly / 2, s2 = SIGMA * SIGMA;
    for (let i = 0; i < n; i++) { ax[i] = 0; ay[i] = 0; }
    for (let i = 0; i < n; i++) {
      const xi = x[i], yi = y[i];
      for (let j = i + 1; j < n; j++) {
        let dx = xi - x[j];
        if (dx > hx) dx -= lx; else if (dx < -hx) dx += lx;
        if (dx > SIGMA || dx < -SIGMA) continue;
        let dy = yi - y[j];
        if (dy > hy) dy -= ly; else if (dy < -hy) dy += ly;
        const r2 = dx * dx + dy * dy;
        if (r2 >= s2 || r2 < 1e-12) continue;
        const r = Math.sqrt(r2), f = STIFF * (SIGMA - r) / r;
        ax[i] += f * dx; ay[i] += f * dy; ax[j] -= f * dx; ay[j] -= f * dy;
      }
    }
  }

  wrap() {
    const { x, y, lx, ly } = this;
    for (let i = 0; i < x.length; i++) { x[i] = ((x[i] % lx) + lx) % lx; y[i] = ((y[i] % ly) + ly) % ly; }
  }

  relax(steps) {                                      // push overlapping molecules apart, without motion in time
    for (let s = 0; s < steps; s++) {
      this.forces();
      for (let i = 0; i < this.x.length; i++) { this.x[i] += 0.25 * this.ax[i] / STIFF; this.y[i] += 0.25 * this.ay[i] / STIFF; }
      this.wrap();
    }
  }

  step(ps) {                                          // advance by a time in ps
    this.time += ps;
    if (!this.fluid) return;
    const { x, y, vx, vy, ax, ay } = this, n = x.length;
    const steps = Math.max(1, Math.min(16, Math.round(ps / DT)));
    for (let s = 0; s < steps; s++) {
      this.forces();
      for (let i = 0; i < n; i++) { vx[i] += ax[i] * DT; vy[i] += ay[i] * DT; x[i] += vx[i] * DT; y[i] += vy[i] * DT; }
      this.wrap();
    }
    let sum = 0;                                      // keep the speeds at the temperature
    for (let i = 0; i < n; i++) sum += vx[i] * vx[i] + vy[i] * vy[i];
    if (sum > 0) {
      const ratio = Math.min(4, Math.max(0.25, 2 * n * KT_M * this.T / sum));
      const factor = Math.sqrt(1 + 0.2 * (ratio - 1));
      for (let i = 0; i < n; i++) { vx[i] *= factor; vy[i] *= factor; }
    }
  }

  // Calls draw(x, y) for every molecule that touches the picture, 0 <= x <= w and 0 <= y <= h (nm)
  each(draw) {
    const { w, h } = this;
    if (!this.fluid) {
      const amp = 0.02 * Math.sqrt(this.T / 273.16), t = this.time;
      for (const s of this.sites) draw(s.x + amp * s.sx * Math.sin(6.2832 * s.fx * t + s.px), s.y + amp * s.sy * Math.sin(6.2832 * s.fy * t + s.py));
      return;
    }
    const { x, y, lx, ly } = this;
    for (let i = 0; i < x.length; i++) {
      for (const X of [x[i], x[i] - lx, x[i] + lx]) for (const Y of [y[i], y[i] - ly, y[i] + ly]) {
        if (X > -RADIUS && X < w + RADIUS && Y > -RADIUS && Y < h + RADIUS) draw(X, Y);
      }
    }
  }
}
