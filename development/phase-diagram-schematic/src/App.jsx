import React, { useEffect, useRef, useState } from 'react';
import { useBook, Studio, Panel, Slider, Sci, V, seeded, LINE, DASH, STATIC } from '../../studio-kit/kit.jsx';
import { TRIPLE, CRITICAL, saturationPressure, sublimationPressure, meltingPressure,
         saturationTemperature, sublimationTemperature, meltingTemperature, phases } from './water.js';
import { Sample, THICK, RADIUS } from './layer.js';

/*
  Phase diagram of water (Section 6.1). Panel (a) is the P-T diagram with a state that the student moves.
  Panel (b) shows the molecules of each phase present at that state, in a thin layer and in their real number,
  and panel (c) lists the density, the molar volume, and the mean distance between molecules (src/water.js).

  Near a coexistence curve the state moves onto the curve (`on` is 'lg', 'sg', 'sl', 'triple', or 'critical'),
  and panels (b) and (c) then show the phases that coexist. The motion in panel (b) is schematic (src/layer.js);
  one second shows 1 ps.
*/
const T_MIN = 200, T_MAX = 800;                       // K
const P_MIN = 1, P_MAX = 1e8;                         // Pa: from 0.00001 bar to 1000 bar
const START = { T: 600, P: saturationPressure(600), on: 'lg' };
const EDGE = 5, TALL = 6;                             // nm, width and height of the layer in panel (b)
const GAP = 0.16;                                     // nm, space between the pictures of two phases
const PS_PER_SECOND = 1;
const HEIGHT = 392;                                   // px, height of the diagram
const PICTURE = 300;                                  // px, largest width of the picture of the layer
const FONT = 'Inter, ui-sans-serif, system-ui, sans-serif';

// The three coexistence curves as lists of [T, P]
const T_LOW = sublimationTemperature(P_MIN);
const SG = Array.from({ length: 81 }, (_, i) => { const T = T_LOW + (TRIPLE.T - T_LOW) * i / 80; return [T, sublimationPressure(T)]; });
const SL = Array.from({ length: 81 }, (_, i) => { const P = TRIPLE.P * (P_MAX / TRIPLE.P) ** (i / 80); return [meltingTemperature(P), P]; });
const LG = Array.from({ length: 161 }, (_, i) => { const T = TRIPLE.T + (CRITICAL.T - TRIPLE.T) * i / 160; return [T, saturationPressure(T)]; });
const T_TOP = SL[80][0];                              // where the solid-liquid curve leaves the diagram

const NAME = { solid: 'solid', liquid: 'liquid', gas: 'gas', supercritical: 'supercritical fluid', critical: 'fluid at the critical point' };
const hues = (color) => ({ solid: color.blue, liquid: color.green, gas: color.vermilion, supercritical: color.purple, critical: color.purple });
const clamp = (x, lo, hi) => Math.min(hi, Math.max(lo, x));
const SUPER = { '-': '⁻', 0: '⁰', 1: '¹', 2: '²', 3: '³', 4: '⁴', 5: '⁵', 6: '⁶', 7: '⁷', 8: '⁸', 9: '⁹' };
function plain(value) {                               // a small number as text for the canvas: 0.19, or 5 × 10⁻⁴
  if (value >= 0.01) return value.toPrecision(2);
  const power = Math.floor(Math.log10(value)), mantissa = Math.round(value / 10 ** power);
  return mantissa === 10 ? `10${[...String(power + 1)].map((c) => SUPER[c]).join('')}`
                         : `${mantissa} × 10${[...String(power)].map((c) => SUPER[c]).join('')}`;
}

// A number with `digits` significant figures, as m x 10^n when it is very small or very large
function Num({ value, digits = 4 }) {
  if (value >= 10 ** digits || value < 0.01) return <Sci value={value} digits={Math.min(digits, 3)} />;
  return <>{value.toPrecision(digits)}</>;
}

function Diagram({ state, onMove, color, dark }) {
  const box = useRef(null), dragging = useRef(false);
  const [width, setWidth] = useState(430);
  useEffect(() => {
    const element = box.current, measure = () => setWidth(element.clientWidth);
    const sizes = new ResizeObserver(measure);
    sizes.observe(element); measure();
    return () => sizes.disconnect();
  }, []);

  const m = { l: 58, r: 12, t: 10, b: 24 };
  const pw = Math.max(120, width - m.l - m.r), ph = HEIGHT - m.t - m.b;
  const right = m.l + pw, bottom = m.t + ph;
  const sx = (T) => m.l + (T - T_MIN) / (T_MAX - T_MIN) * pw;
  const sy = (P) => m.t + (1 - Math.log10(P / P_MIN) / Math.log10(P_MAX / P_MIN)) * ph;
  const xy = ([T, P]) => sx(T).toFixed(1) + ' ' + sy(P).toFixed(1);
  const line = (points) => points.map((p, i) => (i ? 'L' : 'M') + xy(p)).join(' ');
  const more = (points) => points.map((p) => 'L' + xy(p)).join(' ');
  const back = (points) => more([...points].reverse());
  const hue = hues(color), narrow = pw < 330;

  // The regions in which one phase is stable
  const regions = [
    [hue.solid, `M${xy([T_MIN, P_MAX])} L${xy([T_TOP, P_MAX])} ${back(SL)} ${back(SG)} L${xy([T_MIN, P_MIN])} Z`],
    [hue.liquid, `M${xy([T_TOP, P_MAX])} L${xy([CRITICAL.T, P_MAX])} L${xy([CRITICAL.T, CRITICAL.P])} ${back(LG)} ${more(SL)} Z`],
    [hue.supercritical, `M${xy([CRITICAL.T, P_MAX])} L${xy([T_MAX, P_MAX])} L${xy([T_MAX, CRITICAL.P])} L${xy([CRITICAL.T, CRITICAL.P])} Z`],
    [hue.gas, `${line(SG)} ${more(LG)} L${xy([T_MAX, CRITICAL.P])} L${xy([T_MAX, P_MIN])} Z`],
  ];

  const place = (clientX, clientY) => {
    const frame = box.current.getBoundingClientRect();
    const x = clamp(clientX - frame.left, m.l, right), y = clamp(clientY - frame.top, m.t, bottom);
    const T = T_MIN + (x - m.l) / pw * (T_MAX - T_MIN), P = P_MIN * (P_MAX / P_MIN) ** (1 - (y - m.t) / ph);
    for (const [on, point] of [['triple', TRIPLE], ['critical', CRITICAL]]) {
      if (Math.hypot(x - sx(point.T), y - sy(point.P)) < 10) return onMove({ T: point.T, P: point.P, on });
    }
    let best = null;                                  // the nearest coexistence curve within 8 px
    for (const [on, points] of [['lg', LG], ['sg', SG], ['sl', SL]]) for (const [Tc, Pc] of points) {
      const d = Math.hypot(x - sx(Tc), y - sy(Pc));
      if (d < 8 && (!best || d < best.d)) best = { d, on };
    }
    if (!best) return onMove({ T, P, on: null });
    if (best.on === 'sl') { const p = clamp(P, TRIPLE.P, P_MAX); return onMove({ T: meltingTemperature(p), P: p, on: 'sl' }); }
    if (best.on === 'lg') { const t = clamp(T, TRIPLE.T, CRITICAL.T); return onMove({ T: t, P: saturationPressure(t), on: 'lg' }); }
    const t = clamp(T, T_LOW, TRIPLE.T);
    return onMove({ T: t, P: sublimationPressure(t), on: 'sg' });
  };

  const small = { fill: color.text, fontSize: 12 }, label = { fontSize: narrow ? 11.5 : 12.5, fontWeight: 600, textAnchor: 'middle' };
  const cx = sx(state.T), cy = sy(state.P);
  return (
    <div ref={box} className="chart">
      <svg width={width} height={HEIGHT} role="img" style={{ display: 'block', cursor: 'crosshair', userSelect: 'none' }}
           aria-label={`Phase diagram of water, pressure on a logarithmic scale against temperature. The state is at ${state.T.toFixed(1)} kelvin and ${(state.P / 1e5).toPrecision(3)} bar.`}
           onPointerDown={(e) => {
             if (e.pointerType !== 'mouse' && !e.target.dataset.handle) return;        // a finger drags only the point
             dragging.current = true; e.currentTarget.setPointerCapture(e.pointerId); place(e.clientX, e.clientY);
           }}
           onPointerMove={(e) => { if (dragging.current) place(e.clientX, e.clientY); }}
           onPointerUp={() => { dragging.current = false; }} onPointerCancel={() => { dragging.current = false; }}
           onClick={(e) => place(e.clientX, e.clientY)}>
        {regions.map(([fill, d], i) => <path key={i} d={d} fill={fill} fillOpacity={dark ? 0.16 : 0.1} />)}

        {/* 1 bar, and the conventional limits of the supercritical fluid */}
        <path d={`M${m.l} ${sy(1e5)} H${right}`} stroke={color.gray} strokeWidth={1.2} strokeDasharray={DASH} fill="none" />
        <text x={right - 4} y={sy(1e5) - 5} textAnchor="end" {...small} fill={color.muted}>1 bar</text>
        <path d={`M${xy([CRITICAL.T, P_MAX])} L${xy([CRITICAL.T, CRITICAL.P])} L${xy([T_MAX, CRITICAL.P])}`} stroke={color.gray} strokeWidth={1.2} strokeDasharray="1.5 3" fill="none" />

        {[SG, SL, LG].map((points, i) => <path key={i} d={line(points)} stroke={color.text} strokeWidth={LINE} strokeLinecap="round" strokeLinejoin="round" fill="none" />)}

        <text x={sx(236)} y={sy(2e6)} fill={hue.solid} {...label}>solid</text>
        <text x={sx(400)} y={sy(6e6)} fill={hue.liquid} {...label}>liquid</text>
        <text x={sx(560)} y={sy(1.5e3)} fill={hue.gas} {...label}>gas</text>
        <text x={sx(0.5 * (CRITICAL.T + T_MAX))} y={sy(P_MAX) + 12.5} fill={hue.supercritical} {...label} fontSize={11.5}>supercritical</text>
        <text x={sx(0.5 * (CRITICAL.T + T_MAX))} y={sy(P_MAX) + 25.5} fill={hue.supercritical} {...label} fontSize={11.5}>fluid</text>

        <circle cx={sx(TRIPLE.T)} cy={sy(TRIPLE.P)} r={3.6} fill={color.text} />
        <text x={sx(TRIPLE.T) + 9} y={sy(TRIPLE.P) + 15} {...small}>triple point</text>
        <circle cx={sx(CRITICAL.T)} cy={sy(CRITICAL.P)} r={3.6} fill={color.text} />
        <text x={sx(CRITICAL.T) + 9} y={sy(CRITICAL.P) + 16} {...small}>critical point</text>

        {/* Axes */}
        <path d={`M${m.l} ${m.t} V${bottom} H${right}`} stroke={color.axis} fill="none" />
        {[200, 300, 400, 500, 600, 700, 800].filter((T) => !narrow || T % 200 === 0).map((T) => (
          <g key={T}>
            <path d={`M${sx(T)} ${bottom} v4`} stroke={color.axis} />
            <text x={sx(T)} y={bottom + 18} textAnchor="middle" {...small}>{T}</text>
          </g>
        ))}
        {[-5, -4, -3, -2, -1, 0, 1, 2, 3].map((power) => {
          const y = sy(1e5 * 10 ** power);
          return (
            <g key={power}>
              <path d={`M${m.l} ${y} h-4`} stroke={color.axis} />
              <text x={m.l - 8} y={y + 4} textAnchor="end" {...small}>
                {power === 0 ? '1' : power === 1 ? '10' : <>10<tspan dy={-5} fontSize={9}>{String(power).replace('-', '−')}</tspan></>}
              </text>
            </g>
          );
        })}
        <text x={13} y={m.t + ph / 2} textAnchor="middle" transform={`rotate(-90 13 ${m.t + ph / 2})`} fill={color.text} fontSize={12.5}>
          Pressure, <tspan fontStyle="italic">P</tspan> (bar), logarithmic scale
        </text>

        {/* The state */}
        <circle cx={cx} cy={cy} r={6.5} fill={color.ground} stroke={color.text} strokeWidth={2.4} />
        <circle cx={cx} cy={cy} r={2} fill={color.text} />
        <circle data-handle="1" cx={cx} cy={cy} r={20} fill="transparent" style={{ touchAction: 'none', cursor: 'grab' }} />
      </svg>
      <p className="xlabel" style={{ paddingLeft: m.l }}>Temperature, <V>T</V> (K)</p>
    </div>
  );
}

function Layer({ list, T, playing, color }) {
  const canvas = useRef(null), samples = useRef(null), now = useRef({});
  if (!samples.current) { const random = seeded(6010); samples.current = [0, 1, 2].map(() => new Sample(random)); }
  now.current = { list, T, playing, color };

  useEffect(() => {
    let frame = 0, last = performance.now();
    const draw = (time) => {
      const { list, T, playing, color } = now.current;
      const element = canvas.current;
      if (!element) return;
      const width = element.clientWidth, ratio = window.devicePixelRatio || 1;
      const side = Math.min(width, PICTURE), tall = side * TALL / EDGE, height = tall + 50;
      if (element.width !== Math.round(width * ratio) || element.height !== Math.round(height * ratio)) {
        element.width = Math.round(width * ratio); element.height = Math.round(height * ratio); element.style.height = height + 'px';
      }
      const seconds = Math.min((time - last) / 1000, 0.05);
      last = time;
      const ctx = element.getContext('2d');
      ctx.setTransform(ratio, 0, 0, ratio, 0, 0);
      ctx.clearRect(0, 0, width, height);
      const scale = side / EDGE, left = (width - side) / 2, top = 22, hue = hues(color);
      const n = list.length, w = (EDGE - GAP * (n - 1)) / n;                     // nm, width of the picture of one phase
      list.forEach((phase, k) => {
        const sample = samples.current[k], x0 = left + k * (w + GAP) * scale, wide = w * scale;
        sample.set(phase.name, phase.number, T, w, TALL);
        if (playing) sample.step(seconds * PS_PER_SECOND);
        ctx.fillStyle = color.fill; ctx.fillRect(x0, top, wide, tall);
        ctx.save();
        ctx.beginPath(); ctx.rect(x0, top, wide, tall); ctx.clip();
        ctx.fillStyle = hue[phase.name];
        sample.each((x, y) => { ctx.beginPath(); ctx.arc(x0 + x * scale, top + y * scale, RADIUS * scale, 0, 2 * Math.PI); ctx.fill(); });
        ctx.restore();
        if (sample.mean < 1) {                        // fewer than one molecule on average: say how many
          ctx.fillStyle = color.muted; ctx.font = `12px ${FONT}`; ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
          const mid = x0 + wide / 2, lines = wide > 150 ? [`${plain(sample.mean)} molecules`, 'on average'] : [plain(sample.mean), 'molecules', 'on average'];
          lines.forEach((text, i) => ctx.fillText(text, mid, top + tall / 2 + (i - (lines.length - 1) / 2) * 16));
        }
        ctx.strokeStyle = color.axis; ctx.lineWidth = 1.5; ctx.strokeRect(x0, top, wide, tall);
        ctx.fillStyle = hue[phase.name]; ctx.font = `600 12.5px ${FONT}`; ctx.textAlign = 'center'; ctx.textBaseline = 'alphabetic';
        ctx.fillText(NAME[phase.name], x0 + wide / 2, top - 7);
      });
      for (let k = n; k < 3; k++) { samples.current[k].kind = null; samples.current[k].fluid = false; }
      // Scale bar of 1 nm under the left corner
      const y = top + tall + 16;
      ctx.strokeStyle = color.text; ctx.lineWidth = 1.5;
      ctx.beginPath(); ctx.moveTo(left, y); ctx.lineTo(left + scale, y); ctx.moveTo(left, y - 4); ctx.lineTo(left, y + 4);
      ctx.moveTo(left + scale, y - 4); ctx.lineTo(left + scale, y + 4); ctx.stroke();
      ctx.fillStyle = color.text; ctx.font = `12px ${FONT}`; ctx.textAlign = 'left'; ctx.textBaseline = 'middle';
      ctx.fillText('1 nm', left + scale + 8, y);
      frame = requestAnimationFrame(draw);
    };
    frame = requestAnimationFrame(draw);
    return () => cancelAnimationFrame(frame);
  }, []);

  return <canvas ref={canvas} style={{ height: PICTURE * TALL / EDGE + 50 }} role="img"
                 aria-label={`A thin layer of water as ${list.map((p) => NAME[p.name]).join(', ')}, with the molecules drawn as disks in their real number.`} />;
}

export default function PhaseDiagram() {
  const { color, dark } = useBook();
  const [state, setState] = useState(START);
  const [playing, setPlaying] = useState(() => {
    if (STATIC) return false;
    try { return !window.matchMedia('(prefers-reduced-motion: reduce)').matches; } catch (e) { return true; }
  });
  const { T, P, on } = state;
  const list = phases(T, P, on);
  const hue = hues(color);
  const width = (EDGE - GAP * (list.length - 1)) / list.length;                // nm, as in panel (b)

  // A slider moves the state along its coexistence curve as far as the curve goes, and off it beyond the end.
  const setT = (t) => {
    if (on === 'lg' && t >= TRIPLE.T && t <= CRITICAL.T) return setState({ T: t, P: saturationPressure(t), on });
    if (on === 'sg' && t >= T_LOW && t <= TRIPLE.T) return setState({ T: t, P: sublimationPressure(t), on });
    if (on === 'sl' && t >= T_TOP && t <= TRIPLE.T) return setState({ T: t, P: meltingPressure(t), on });
    return setState({ T: t, P, on: null });
  };
  const setP = (p) => {
    if (on === 'lg' && p >= TRIPLE.P && p <= CRITICAL.P) return setState({ T: saturationTemperature(p), P: p, on });
    if (on === 'sg' && p >= P_MIN && p <= TRIPLE.P) return setState({ T: sublimationTemperature(p), P: p, on });
    if (on === 'sl' && p >= TRIPLE.P && p <= P_MAX) return setState({ T: meltingTemperature(p), P: p, on });
    return setState({ T, P: p, on: null });
  };

  const where = {
    lg: 'On the liquid–gas curve, liquid and gas coexist.',
    sg: 'On the solid–gas curve, solid and gas coexist.',
    sl: 'On the solid–liquid curve, solid and liquid coexist.',
    triple: 'At the triple point, solid, liquid, and gas coexist.',
    critical: 'At the critical point, liquid and gas are one phase.',
  }[on] || `One phase is stable at this state: ${list[0].name === 'solid' ? 'solid (ice Ih)' : NAME[list[0].name]}.`;
  const vbar = <V bar>v</V>;
  const rows = [
    [<>Density, <V>ρ</V> (g cm<sup>−3</sup>)</>, (p) => <Num value={p.rho / 1000} />],
    [<>Molar volume, {vbar} (cm<sup>3</sup> mol<sup>−1</sup>)</>, (p) => <Num value={p.v * 1e6} />],
    [<>Mean distance between molecules, ({vbar}/<V sub="A">N</V>)<sup>1/3</sup> (nm)</>, (p) => <Num value={p.spacing} digits={3} />],
    ['Molecules in the picture in panel (b), on average', (p) => <Num value={p.number * THICK * width * TALL} digits={3} />],
  ];

  return (
    <Studio title="Phase diagram of water"
            lede="Move a state through the phase diagram of water, and look at the molecules of the phases that are stable there.">
      <div className="panels two left-wide">
        <Panel letter="a" title="Phase diagram of water">
          <Diagram state={state} onMove={setState} color={color} dark={dark} />
        </Panel>
        <Panel letter="b" title="Molecules in a layer 0.33 nm thick">
          <Layer list={list} T={T} playing={playing} color={color} />
        </Panel>
      </div>
      <div className="panels" style={{ marginTop: 14 }}>
        <Panel letter="c" title={where}>
          <table className="values">
            <thead>
              <tr><th></th>{list.map((p) => <th key={p.name} style={{ color: hue[p.name] }}>{p.name === 'critical' ? 'fluid' : NAME[p.name]}</th>)}</tr>
            </thead>
            <tbody>
              {rows.map(([name, value], i) => (
                <tr key={i}><td>{name}</td>{list.map((p) => <td key={p.name} className="num">{value(p)}</td>)}</tr>
              ))}
            </tbody>
          </table>
        </Panel>
      </div>

      <div className="controls">
        <Slider label={<>Temperature, <V>T</V></>} show={<>{T.toFixed(on === 'triple' || on === 'critical' ? 2 : 1)} K</>} value={T} min={T_MIN} max={T_MAX} step={1} onChange={setT} />
        <Slider label={<>Pressure, <V>P</V></>} show={<><Num value={P / 1e5} /> bar</>} value={Math.log10(P / 1e5)} min={-5} max={3} step={0.01}
                onChange={(x) => setP(clamp(1e5 * 10 ** x, P_MIN, P_MAX))} />
        <div className="row">
          <button type="button" onClick={() => setState({ T: TRIPLE.T, P: TRIPLE.P, on: 'triple' })}>Triple point</button>
          <button type="button" onClick={() => setState({ T: CRITICAL.T, P: CRITICAL.P, on: 'critical' })}>Critical point</button>
          <button type="button" onClick={() => setPlaying(!playing)}>{playing ? 'Pause' : 'Play'}</button>
          <button type="button" onClick={() => setState(START)}>Reset</button>
        </div>
      </div>
      <p className="hint">
        Drag the point in panel (a), or use the sliders. Near a coexistence curve the point moves onto the curve, and the sliders then move it
        along the curve. Panel (b) shows a layer of water 5 nm wide, 6 nm high, and 0.33 nm thick, which is about one molecule. The number of molecules in it
        is the real one, from the density of the phase. Their motion is schematic, and one second shows 1 ps. In ice Ih each molecule has four
        neighbors, and the picture shows three of them in a plane. The dotted lines at the critical temperature and the critical pressure are a
        convention and not coexistence curves. The coexistence curves and the densities follow the equations of the International Association
        for the Properties of Water and Steam.
      </p>
    </Studio>
  );
}
