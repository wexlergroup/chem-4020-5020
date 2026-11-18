import React, { useEffect, useRef, useState } from 'react';
import { LineChart, Line, XAxis, YAxis, ResponsiveContainer, ReferenceLine, ReferenceDot } from 'recharts';
import { useBook, Studio, Panel, Slider, V, axis, seeded, LINE, DASH, STATIC, onPaper } from '../../studio-kit/kit.jsx';
import { reaction } from './reaction.js';

/*
  Chemical equilibrium (Section 7.1). Panel (a) is G(xi) - G(0) for 2 NO2 <=> N2O4 with the tangent at the extent that
  the student chooses; its slope is the reaction Gibbs energy. Panel (b) shows the molecules of the mixture at that
  extent (40 NO2 at the start, schematic) in a flask whose brown color follows the concentration of NO2.
  Panel (c) lists K, Q_p, the slope, and the composition (src/reaction.js).
*/
const START = { T: 298.15, logP: 0, xi: 0.5 };
const UNITS = 40;                                     // molecules of NO2 at xi = 0, for 2 mol
const SIZE = 300;                                     // px, height of the flask
const FONT = 'Inter, ui-sans-serif, system-ui, sans-serif';
const BROWN = [150, 75, 20];                          // color of NO2 gas
const ABSORB = 8.85;                                  // L/mol: pure NO2 at 1 bar and 298 K gives 30 % of the darkest shade

function ticksOf(lo, hi, n) {
  const rough = (hi - lo) / n, power = 10 ** Math.floor(Math.log10(rough));
  const step = [1, 2, 5, 10].map((m) => m * power).find((s) => s >= rough - 1e-12);
  const ticks = [];
  for (let v = Math.ceil(lo / step - 1e-9) * step; v <= hi + 1e-9; v += step) ticks.push(Number((Math.round(v / step) * step).toPrecision(10)));
  return { ticks, step };
}
const minus = (s) => String(s).replace('-', '−');
function num(x, digits = 3) {                         // 3 significant figures, as m × 10^n outside 0.01 to 10000
  if (x === 0) return '0';
  const a = Math.abs(x);
  if (a >= 0.01 && a < 1e4) return minus(Number(x.toPrecision(digits)).toString());
  const p = Math.floor(Math.log10(a)), m = x / 10 ** p;
  return <>{minus(m.toFixed(digits - 1))} × 10<sup>{minus(p)}</sup></>;
}

function Flask({ xi, c1, playing, color }) {
  const canvas = useRef(null), now = useRef({}), units = useRef(null);
  if (!units.current) {                               // each NO2 unit has a home and a partner for dimerization
    const random = seeded(7010), list = [], cells = [];  // start on a 7 x 6 grid with some disorder, so that no two overlap
    for (let k = 0; k < 42; k++) cells.push(k);
    for (let k = cells.length - 1; k > 0; k--) { const m = Math.floor(random() * (k + 1)); [cells[k], cells[m]] = [cells[m], cells[k]]; }
    for (let i = 0; i < UNITS; i++) {
      const cx = cells[i] % 7, cy = Math.floor(cells[i] / 7);
      list.push({ x: (cx + 0.25 + 0.5 * random()) / 7, y: (cy + 0.25 + 0.5 * random()) / 6, vx: random() - 0.5, vy: random() - 0.5, a: 6.28 * random(), w: random() - 0.5 });
    }
    units.current = list;
  }
  now.current = { xi, c1, playing, color };

  useEffect(() => {
    let frame = 0, last = performance.now();
    const draw = (time) => {
      const { xi, c1, playing, color } = now.current;
      const element = canvas.current;
      if (!element) return;
      const width = element.clientWidth, ratio = window.devicePixelRatio || 1;
      if (element.width !== Math.round(width * ratio)) { element.width = Math.round(width * ratio); element.height = Math.round(SIZE * ratio); }
      const seconds = Math.min((time - last) / 1000, 0.05);
      last = time;
      const ctx = element.getContext('2d');
      ctx.setTransform(ratio, 0, 0, ratio, 0, 0);
      ctx.clearRect(0, 0, width, SIZE);
      const side = Math.min(width - 4, SIZE - 30), left = (width - side) / 2, top = 2;
      // The gas, brown in proportion to the light that NO2 absorbs
      const alpha = 0.55 * (1 - Math.exp(-ABSORB * c1));
      ctx.fillStyle = color.fill; ctx.fillRect(left, top, side, side);
      ctx.fillStyle = `rgba(${BROWN.join(',')},${alpha.toFixed(3)})`; ctx.fillRect(left, top, side, side);
      ctx.strokeStyle = color.axis; ctx.lineWidth = 1.5; ctx.strokeRect(left, top, side, side);
      // Molecules: pairs (2k, 2k+1) form N2O4 first, so that 2 xi of every 2 units are bound
      const dimers = Math.round(xi * UNITS / 2), list = units.current, r = 5.5, pad = 10;
      ctx.save(); ctx.beginPath(); ctx.rect(left, top, side, side); ctx.clip();
      ctx.strokeStyle = color.ground; ctx.lineWidth = 1.5;                       // a thin rim keeps the molecules clear of the brown
      for (let i = 0; i < UNITS; i++) {
        const u = list[i];
        if (playing) {
          u.x += 0.05 * u.vx * seconds; u.y += 0.05 * u.vy * seconds; u.a += 1.5 * u.w * seconds;
          if (u.x < 0 || u.x > 1) { u.vx = -u.vx; u.x = Math.min(1, Math.max(0, u.x)); }
          if (u.y < 0 || u.y > 1) { u.vy = -u.vy; u.y = Math.min(1, Math.max(0, u.y)); }
        }
        const pair = Math.floor(i / 2), bound = pair < dimers;
        if (bound && i % 2 === 1) continue;           // the partner draws the dimer
        const x = left + pad + u.x * (side - 2 * pad), y = top + pad + u.y * (side - 2 * pad);
        if (bound) {
          const dx = r * Math.cos(u.a), dy = r * Math.sin(u.a);
          ctx.fillStyle = color.blue;
          ctx.beginPath(); ctx.arc(x - dx, y - dy, r, 0, 6.2832); ctx.stroke(); ctx.beginPath(); ctx.arc(x + dx, y + dy, r, 0, 6.2832); ctx.stroke();
          ctx.beginPath(); ctx.arc(x - dx, y - dy, r, 0, 6.2832); ctx.arc(x + dx, y + dy, r, 0, 6.2832); ctx.fill();
        } else {
          ctx.fillStyle = color.vermilion;
          ctx.beginPath(); ctx.arc(x, y, r, 0, 6.2832); ctx.stroke(); ctx.fill();
        }
      }
      ctx.restore();
      // Key
      const y = top + side + 16;
      ctx.font = `12.5px ${FONT}`; ctx.textBaseline = 'middle'; ctx.textAlign = 'left';
      ctx.fillStyle = color.vermilion; ctx.beginPath(); ctx.arc(left + 6, y, r, 0, 6.2832); ctx.fill();
      ctx.fillStyle = color.text; ctx.fillText('NO₂', left + 16, y);
      ctx.fillStyle = color.blue; ctx.beginPath(); ctx.arc(left + 72, y, r, 0, 6.2832); ctx.arc(left + 82, y, r, 0, 6.2832); ctx.fill();
      ctx.fillStyle = color.text; ctx.fillText('N₂O₄', left + 94, y);
      frame = requestAnimationFrame(draw);
    };
    frame = requestAnimationFrame(draw);
    const off = onPaper(() => { cancelAnimationFrame(frame); draw(performance.now()); });
    return () => { off(); cancelAnimationFrame(frame); };
  }, []);
  const dimers = Math.round(xi * UNITS / 2);
  return <canvas ref={canvas} style={{ height: SIZE }} role="img"
                 aria-label={`A flask with ${UNITS - 2 * dimers} molecules of NO2 and ${dimers} molecules of N2O4. The gas is brown in proportion to the concentration of NO2.`} />;
}

export default function Equilibrium() {
  const { color } = useBook();
  const [T, setT] = useState(START.T);
  const [logP, setLogP] = useState(START.logP);
  const [xi, setXi] = useState(START.xi);
  const [playing, setPlaying] = useState(() => {
    if (STATIC) return false;
    try { return !window.matchMedia('(prefers-reduced-motion: reduce)').matches; } catch (e) { return true; }
  });
  const P = 10 ** logP;
  const r = reaction(T, P);
  const s = r.at(xi);

  const curve = [];
  for (let i = 0; i <= 200; i++) { const x = i / 200; curve.push({ xi: x, G: r.g(x) }); }
  const gs = curve.map((p) => p.G), lo = Math.min(...gs), hi = Math.max(...gs), pad = 0.08 * (hi - lo);
  const { ticks: yTicks, step } = ticksOf(lo - pad, hi + pad, 5);
  const domain = [Math.min(yTicks[0], lo - pad), Math.max(yTicks[yTicks.length - 1], hi + pad)];
  const half = 0.16;
  const tangent = [
    { xi: Math.max(0, xi - half), T: s.G - Math.min(half, xi) * s.slope },
    { xi: Math.min(1, xi + half), T: s.G + Math.min(half, 1 - xi) * s.slope },
  ];

  const dr = (sym) => <>Δ<sub>r</sub><V>{sym}</V></>;
  const rows = [
    [<>Standard reaction Gibbs energy</>, <>{dr('G')}<sup>∘</sup> = {dr('H')}<sup>∘</sup> − <V>T</V>{dr('S')}<sup>∘</sup></>, <>{minus(r.dG0.toFixed(2))} kJ mol<sup>−1</sup></>],
    [<>Equilibrium constant</>, <><V>K</V> = exp(−{dr('G')}<sup>∘</sup>/<V>RT</V>)</>, num(r.K)],
    [<>Reaction quotient at ξ</>, <><V sub="p">Q</V> = <V sub={<>N<sub>2</sub>O<sub>4</sub></>}>y</V><V>P</V><sup>∘</sup>/(<V sub={<>NO<sub>2</sub></>}>y</V><sup>2</sup><V>P</V>)</>, num(s.Q)],
    [<>Slope of <V>G</V>(ξ) at ξ</>, <>{dr('G')} = <V>RT</V> ln(<V sub="p">Q</V>/<V>K</V>)</>, <>{minus(s.slope.toFixed(2))} kJ mol<sup>−1</sup></>],
    [<>Mole fractions at ξ</>, <><V sub={<>NO<sub>2</sub></>}>y</V> and <V sub={<>N<sub>2</sub>O<sub>4</sub></>}>y</V></>, <>{s.y1.toFixed(3)} and {s.y2.toFixed(3)}</>],
    [<>Extent at the minimum</>, <>ξ<sub>eq</sub>, from <V>K</V></>, <>{r.xiEq.toFixed(3)} mol</>],
  ];

  return (
    <Studio title="Chemical equilibrium"
            lede="Follow the Gibbs free energy of a mixture of NO₂ and N₂O₄ as the reaction 2 NO₂ ⇌ N₂O₄ proceeds, and change the temperature and the total pressure.">
      <p className="key">
        <span style={{ color: color.blue }}><i className="line" /><span><V>G</V>(ξ) − <V>G</V>(0)</span></span>
        <span style={{ color: color.vermilion }}><i className="line" /><span>tangent at ξ</span></span>
        <span style={{ color: color.gray }}><i className="dash" /><span>ξ<sub>eq</sub></span></span>
      </p>
      <div className="panels two left-wide" style={{ marginTop: 8 }}>
        <Panel letter="a" title={<>Gibbs free energy, <V>G</V>(ξ) − <V>G</V>(0) (kJ)</>}>
          <div className="chart" style={{ height: SIZE - 24 }}>
            <ResponsiveContainer width="100%" height="100%">
              <LineChart margin={{ top: 12, right: 14, bottom: -6, left: 0 }}>
                <XAxis dataKey="xi" type="number" domain={[0, 1]} ticks={[0, 0.2, 0.4, 0.6, 0.8, 1]} {...axis(color)} />
                <YAxis type="number" domain={domain} ticks={yTicks} width={44} allowDataOverflow
                       tickFormatter={(v) => minus(v.toFixed(step < 1 ? 1 : 0))} {...axis(color)} />
                <ReferenceLine x={r.xiEq} stroke={color.gray} strokeWidth={1.2} strokeDasharray={DASH} />
                <Line data={curve} dataKey="G" stroke={color.blue} strokeWidth={LINE} dot={false} isAnimationActive={false} />
                <Line data={tangent} dataKey="T" stroke={color.vermilion} strokeWidth={2.2} dot={false} isAnimationActive={false} />
                <ReferenceDot x={r.xiEq} y={r.g(r.xiEq)} r={4} fill={color.text} stroke="none" />
                <ReferenceDot x={xi} y={s.G} r={6} fill={color.vermilion} stroke={color.ground} strokeWidth={2} />
              </LineChart>
            </ResponsiveContainer>
          </div>
          <p className="xlabel">Extent of reaction, ξ (mol)</p>
        </Panel>
        <Panel letter="b" title={<>The mixture at ξ</>}>
          <Flask xi={xi} c1={s.c1} playing={playing} color={color} />
        </Panel>
      </div>
      <div className="panels" style={{ marginTop: 14 }}>
        <Panel letter="c" title={<>2 NO<sub>2</sub>(g) ⇌ N<sub>2</sub>O<sub>4</sub>(g), from 2 mol of NO<sub>2</sub>, at ξ = {xi.toFixed(2)} mol</>}>
          <table className="values">
            <tbody>
              {rows.map(([name, formula, value], i) => (
                <tr key={i}><td>{name}<span className="formula">{formula}</span></td><td>{value}</td></tr>
              ))}
            </tbody>
          </table>
        </Panel>
      </div>

      <div className="controls">
        <Slider label={<>Extent of reaction, ξ</>} show={<>{xi.toFixed(2)} mol</>} value={xi} min={0.01} max={0.99} step={0.01} onChange={setXi} />
        <Slider label={<>Temperature, <V>T</V></>} show={<>{T.toFixed(T % 1 ? 2 : 0)} K</>} value={T} min={250} max={400} step={1} onChange={setT} />
        <Slider label={<>Total pressure, <V>P</V></>} show={<>{num(P)} bar</>} value={logP} min={-1} max={1} step={0.02}
                onChange={(v) => setLogP(Math.abs(v) < 1e-9 ? 0 : v)} />
        <div className="row" style={{ gridColumn: 'auto', alignSelf: 'end' }}>
          <button type="button" onClick={() => setXi(Number(r.xiEq.toFixed(2)))}>Go to ξ<sub>eq</sub></button>
          <button type="button" onClick={() => setPlaying(!playing)}>{playing ? 'Pause' : 'Play'}</button>
          <button type="button" onClick={() => { setT(START.T); setLogP(START.logP); setXi(START.xi); }}>Reset</button>
        </div>
      </div>
      <p className="hint">
        The curve is Eq. (7.1.32) of the notes with the total pressure <V>P</V>: <V>G</V>(ξ) − <V>G</V>(0) = ξ({dr('G')}<sup>∘</sup> − <V>RT</V> ln(<V>P</V>/<V>P</V><sup>∘</sup>))
        + <V>RT</V>[(2 − 2ξ) ln <V sub={<>NO<sub>2</sub></>}>y</V> + ξ ln <V sub={<>N<sub>2</sub>O<sub>4</sub></>}>y</V>]. {dr('H')}<sup>∘</sup> = −57.11 kJ mol<sup>−1</sup> and
        {' '}{dr('S')}<sup>∘</sup> = −175.69 J mol<sup>−1</sup> K<sup>−1</sup> are the values at 298.15 K from the NIST-JANAF tables, taken as constant.
        Panel (b) shows 40 units of NO<sub>2</sub> for the 2 mol, so that each N<sub>2</sub>O<sub>4</sub> is two of them bound together. NO<sub>2</sub> is brown and
        N<sub>2</sub>O<sub>4</sub> is colorless, and the shade of the flask follows the concentration of NO<sub>2</sub>, <V sub={<>NO<sub>2</sub></>}>y</V><V>P</V>/<V>RT</V>.
      </p>
    </Studio>
  );
}
