import React, { useEffect, useMemo, useRef, useState } from 'react';
import { BarChart, Bar, LineChart, Line, XAxis, YAxis, Tooltip, ResponsiveContainer, ReferenceLine } from 'recharts';
import { useBook, Studio, Panel, Slider, Sci, V, axis, niceTicks, DASH, STATIC } from '../../studio-kit/kit.jsx';
import { MOLECULES, rotor } from './rotor.js';

/*
  Linear rigid rotor (Section 2.7). A diatomic molecule of the masses m_A and m_B and the bond length r turns
  about its center of mass. Panel (a) shows the molecule to scale, panel (b) the quantities of the notes
  (src/rotor.js), panel (c) the probabilities of the levels, and panel (d) their energies beside the temperature.

  The picture turns with the classical root mean square angular velocity, (2 k_B T / I)^(1/2), and one second
  shows 0.25 ps.
*/
const CO = MOLECULES.find((m) => m.name === 'CO');
const START = { mA: CO.mA, mB: CO.mB, r: CO.r, T: 300, same: false };
const PS_PER_SECOND = 0.25;
const SIZE = 340;                                     // px, height of the picture
const SCALE = 80;                                     // px per angstrom

function Molecule({ mA, mB, r, omega, same, playing, color }) {
  const canvas = useRef(null), now = useRef({}), angle = useRef(0.5);
  now.current = { mA, mB, r, omega, same, playing, color };

  useEffect(() => {
    let frame = 0, last = performance.now();
    const draw = (time) => {
      const { mA, mB, r, omega, same, playing, color } = now.current;
      const element = canvas.current;
      if (!element) return;
      const width = element.clientWidth, ratio = window.devicePixelRatio || 1;
      if (element.width !== Math.round(width * ratio)) { element.width = Math.round(width * ratio); element.height = Math.round(SIZE * ratio); }
      const seconds = Math.min((time - last) / 1000, 0.05);
      last = time;
      if (playing) angle.current += omega * PS_PER_SECOND * 1e-12 * seconds;
      const ctx = element.getContext('2d');
      ctx.setTransform(ratio, 0, 0, ratio, 0, 0);
      ctx.clearRect(0, 0, width, SIZE);
      const cx = width / 2, cy = (SIZE - 26) / 2;
      const a = -r * SCALE * mB / (mA + mB), b = r * SCALE * mA / (mA + mB);    // distances of A and B from the center of mass
      const cos = Math.cos(angle.current), sin = Math.sin(angle.current);
      const radius = (m) => 6 + 2.4 * Math.cbrt(m);
      ctx.strokeStyle = color.axis; ctx.lineWidth = 4; ctx.lineCap = 'round';
      ctx.beginPath(); ctx.moveTo(cx + a * cos, cy + a * sin); ctx.lineTo(cx + b * cos, cy + b * sin); ctx.stroke();
      ctx.font = '600 12px Inter, ui-sans-serif, system-ui, sans-serif'; ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
      for (const [d, m, fill, name] of [[a, mA, color.blue, 'A'], [b, mB, same ? color.blue : color.vermilion, 'B']]) {
        ctx.fillStyle = fill; ctx.beginPath(); ctx.arc(cx + d * cos, cy + d * sin, radius(m), 0, 2 * Math.PI); ctx.fill();
        ctx.fillStyle = color.ground; ctx.fillText(name, cx + d * cos, cy + d * sin + 0.5);
      }
      // The center of mass, which stays in place
      ctx.strokeStyle = color.text; ctx.lineWidth = 1.5; ctx.lineCap = 'butt';
      ctx.beginPath(); ctx.moveTo(cx - 5, cy); ctx.lineTo(cx + 5, cy); ctx.moveTo(cx, cy - 5); ctx.lineTo(cx, cy + 5); ctx.stroke();
      // Scale bar of 1 angstrom
      const y = SIZE - 12;
      ctx.beginPath(); ctx.moveTo(2, y); ctx.lineTo(2 + SCALE, y); ctx.moveTo(2, y - 4); ctx.lineTo(2, y + 4);
      ctx.moveTo(2 + SCALE, y - 4); ctx.lineTo(2 + SCALE, y + 4); ctx.stroke();
      ctx.fillStyle = color.text; ctx.font = '12px Inter, ui-sans-serif, system-ui, sans-serif'; ctx.textAlign = 'left';
      ctx.fillText('1 Å', 2 + SCALE + 8, y);
      frame = requestAnimationFrame(draw);
    };
    frame = requestAnimationFrame(draw);
    return () => cancelAnimationFrame(frame);
  }, []);

  return <canvas ref={canvas} style={{ height: SIZE }} role="img"
                 aria-label="A molecule of two atoms, A and B, joined by a bond, turns about its center of mass, which is marked by a cross." />;
}

function LevelTip({ active, payload, name, digits, unit }) {
  if (!active || !payload || !payload.length) return null;
  const point = payload[0].payload;
  return <div className="tip"><div><V>J</V> = {point.J}</div><div>{name} = {point[payload[0].dataKey].toFixed(digits)}{unit}</div></div>;
}

export default function RigidRotor() {
  const { color } = useBook();
  const [mA, setMA] = useState(START.mA);               // u
  const [mBfree, setMB] = useState(START.mB);           // u
  const [r, setR] = useState(START.r);                  // angstrom
  const [T, setT] = useState(START.T);                  // K
  const [same, setSame] = useState(START.same);         // two identical atoms: sigma = 2
  const [playing, setPlaying] = useState(() => {
    if (STATIC) return false;
    try { return !window.matchMedia('(prefers-reduced-motion: reduce)').matches; } catch (e) { return true; }
  });
  const mB = same ? mA : mBfree, sigma = same ? 2 : 1;

  const g = useMemo(() => rotor(mA, mB, r, T, sigma), [mA, mB, r, T, sigma]);
  // The levels shown: up to the level at which 99.95 % of the probability is reached, and at least J = 0 to 6
  const shown = useMemo(() => {
    let total = 0, last = g.levels.length - 1;
    for (let i = 0; i < g.levels.length; i++) { total += g.levels[i].p; if (total > 0.9995) { last = i; break; } }
    return g.levels.slice(0, Math.max(last, 6) + 1);
  }, [g]);
  const top = shown.length - 1, most = shown.reduce((best, level) => (level.p > best.p ? level : best), shown[0]);
  const jTicks = niceTicks(top, 6).filter((t) => Number.isInteger(t)), jStep = jTicks.length > 1 ? jTicks[1] - jTicks[0] : 1;
  const pTicks = niceTicks(most.p * 1.02, 4), pStep = pTicks[1] - pTicks[0], pTop = Math.ceil(most.p / pStep - 1e-9) * pStep;
  const eMax = Math.max(shown[top].e, T * 1.15), eTicks = niceTicks(eMax, 4), eStep = eTicks[1] - eTicks[0], eTop = Math.ceil(eMax / eStep - 1e-9) * eStep;
  const margin = { top: 14, right: 12, bottom: -6, left: 0 };
  const xTitle = <p className="xlabel">Level, <V>J</V></p>;
  const kB = <V sub="B">k</V>, Th = <>Θ<sub>rot</sub></>, q = <V sub="rot">q</V>, EJ = <V sub={<i className="v">J</i>}>E</V>, gJ = <V sub={<i className="v">J</i>}>g</V>;
  const digits4 = (x) => Number(x.toPrecision(4)).toString();
  const big = (x) => (x >= 1000 ? x.toFixed(0) : x >= 100 ? x.toFixed(1) : x.toFixed(2));
  // The division by sigma = 2 holds only well above the rotational temperature; below 2 Theta_rot the values are not shown
  const low = same && T < 2 * g.theta;
  const rows = [
    ['Reduced mass', <><V>μ</V> = <V sub="A">m</V><V sub="B">m</V>/(<V sub="A">m</V> + <V sub="B">m</V>)</>, <>{digits4(g.mu)} u</>],
    ['Moment of inertia', <><V>I</V> = <V>μ</V><V>r</V><sup>2</sup></>, <><Sci value={g.I} /> kg m<sup>2</sup></>],
    ['Rotational constant', <><V>B̃</V> = <V>h</V>/8π<sup>2</sup><V>c</V><V>I</V></>, <>{digits4(g.B)} cm<sup>−1</sup></>],
    ['Rotational temperature', <>{Th} = ℏ<sup>2</sup>/2{kB}<V>I</V></>, <>{digits4(g.theta)} K</>],
    ['Partition function, sum over J', <>{q} = (1/<V>σ</V>) Σ {gJ} e<sup>−{EJ}/{kB}<V>T</V></sup></>, low ? '—' : big(g.q)],
    [<>Partition function for <V>T</V> ≫ {Th}</>, <>{q} ≈ <V>T</V>/<V>σ</V>{Th}</>, low ? '—' : big(g.high)],
    ['Molar entropy of rotation', low ? <span style={{ color: color.vermilion }}>needs <V>T</V> ≫ {Th} when <V>σ</V> = 2</span>
      : <><V>S</V> = <V>R</V>[ln {q} + <V sub="rot">U</V>/{kB}<V>T</V>]</>, low ? '—' : <>{g.S.toFixed(2)} J K<sup>−1</sup> mol<sup>−1</sup></>],
  ];
  const pick = (m) => { setMA(m.mA); setMB(m.mB); setR(m.r); setSame(m.same); };
  const picked = (m) => Math.abs(m.mA - mA) < 1e-9 && Math.abs(m.mB - mB) < 1e-9 && Math.abs(m.r - r) < 1e-9 && m.same === same;
  const reset = () => { pick(CO); setT(START.T); };

  return (
    <Studio title="Linear rigid rotor"
            lede="A molecule of two atoms turns about its center of mass. Choose the masses, the bond length, and the temperature, and follow the rotational temperature, the partition function, and the probabilities of the rotational levels.">
      <div className="panels two left-narrow">
        <Panel letter="a" title="The molecule, to scale">
          <Molecule mA={mA} mB={mB} r={r} omega={g.omega} same={same} playing={playing} color={color} />
        </Panel>
        <Panel letter="b" title={<>Quantities of the rotor (<V>σ</V> = {sigma})</>}>
          <table className="values">
            <tbody>
              {rows.map(([name, formula, value], i) => (
                <tr key={i}><td>{name}<span className="formula">{formula}</span></td><td>{value}</td></tr>
              ))}
            </tbody>
          </table>
        </Panel>
      </div>

      <div className="panels two" style={{ marginTop: 14 }}>
        <Panel letter="c" title={<>Probabilities of the levels, <V sub={<i className="v">J</i>}>p</V></>}>
          <div className="chart" style={{ height: 190 }}>
            <ResponsiveContainer width="100%" height="100%">
              <BarChart data={shown} margin={margin} barCategoryGap={shown.length > 40 ? 0 : '18%'}>
                <XAxis dataKey="J" interval={jStep - 1} {...axis(color)} />
                <YAxis domain={[0, pTop]} ticks={niceTicks(pTop, Math.round(pTop / pStep))} width={44} {...axis(color)} />
                <Tooltip content={<LevelTip name={<V sub={<i className="v">J</i>}>p</V>} digits={4} unit="" />} cursor={{ fill: color.fill }} isAnimationActive={false} />
                <Bar dataKey="p" fill={color.blue} isAnimationActive={false} />
              </BarChart>
            </ResponsiveContainer>
          </div>
          {xTitle}
          <p className="note"><span>Most probable level: <V>J</V> = {most.J}</span></p>
        </Panel>
        <Panel letter="d" title={<>Energies of the levels, {EJ}/{kB} (K)</>}>
          <div className="chart" style={{ height: 190 }}>
            <ResponsiveContainer width="100%" height="100%">
              <LineChart data={shown} margin={margin}>
                <XAxis dataKey="J" type="number" domain={[0, top]} ticks={jTicks} {...axis(color)} />
                <YAxis domain={[0, eTop]} ticks={niceTicks(eTop, Math.round(eTop / eStep))} width={44} {...axis(color)} />
                <Tooltip content={<LevelTip name={<>{EJ}/{kB}</>} digits={1} unit=" K" />} isAnimationActive={false} />
                <ReferenceLine y={T} stroke={color.gray} strokeDasharray={DASH} strokeWidth={1.4} />
                <Line dataKey="e" stroke="none" dot={{ r: shown.length > 60 ? 1.5 : 2.6, fill: color.text, stroke: 'none' }} activeDot={{ r: 4 }} isAnimationActive={false} />
              </LineChart>
            </ResponsiveContainer>
          </div>
          {xTitle}
          <p className="key"><span style={{ color: color.gray }}><i className="dash" /><span>the temperature, <V>T</V> = {T.toFixed(0)} K</span></span></p>
        </Panel>
      </div>

      <div className="controls">
        <div style={{ gridColumn: '1 / -1' }}>
          <div className="row" style={{ justifyContent: 'flex-start', flexWrap: 'wrap' }}>
            {MOLECULES.map((m) => <button type="button" key={m.name} aria-pressed={picked(m)} onClick={() => pick(m)}>{m.name}</button>)}
          </div>
        </div>
        <Slider label={<>Mass of atom A, <V sub="A">m</V></>} show={<>{mA.toFixed(2)} u</>} value={mA} min={1} max={130} step={0.01} onChange={setMA} />
        <Slider label={<>Mass of atom B, <V sub="B">m</V></>} show={<>{mB.toFixed(2)} u</>} value={mB} min={1} max={130} step={0.01} onChange={setMB} disabled={same} />
        <Slider label={<>Bond length, <V>r</V></>} show={<>{r.toFixed(3)} Å</>} value={r} min={0.5} max={3} step={0.001} onChange={setR} />
        <Slider label={<>Temperature, <V>T</V></>} show={<>{T.toFixed(0)} K</>} value={T} min={10} max={1000} step={10} onChange={setT} />
        <div className="checks" style={{ marginTop: 0 }}>
          <label><input type="checkbox" checked={same} onChange={(e) => setSame(e.target.checked)} /><span>Two identical atoms (<V>σ</V> = 2)</span></label>
        </div>
        <div className="row" style={{ gridColumn: 'auto' }}>
          <button type="button" onClick={() => setPlaying(!playing)}>{playing ? 'Pause' : 'Play'}</button>
          <button type="button" onClick={reset}>Reset</button>
        </div>
      </div>
      <p className="hint">
        The levels are {EJ} = {kB}{Th}<V>J</V>(<V>J</V> + 1) with the degeneracy {gJ} = 2<V>J</V> + 1. The buttons set the masses and the bond length to those of
        the molecules in the table of Section 2.7. In the picture, the molecule turns with the classical root mean square angular velocity, and one second shows 0.25 ps.
        For two identical atoms, the spins of the nuclei change the weights of the levels with even and odd <V>J</V>. The symmetry number <V>σ</V> = 2 accounts
        for this at high temperature only. The probabilities in panel (c) do not include it, and below 2{Th} the studio does not show {q} and <V>S</V> for two identical atoms.
      </p>
    </Studio>
  );
}
