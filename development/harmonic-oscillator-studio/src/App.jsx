import React, { useMemo, useState } from 'react';
import { BarChart, Bar, Cell, LineChart, Line, XAxis, YAxis, Tooltip, ResponsiveContainer, ReferenceLine, ReferenceDot } from 'recharts';
import { useBook, Studio, Panel, Slider, V, axis, niceTicks, LINE, DASH } from '../../studio-kit/kit.jsx';

/*
  Harmonic oscillator (Section 2.6). One oscillator with the levels E_n = hbar omega (n + 1/2) and the
  vibrational temperature Theta_vib = hbar omega / k_B. With x = Theta_vib / T:

      q   = exp(-x/2) / (1 - exp(-x))               partition function
      p_n = exp(-n x) (1 - exp(-x))                 probability of level n
      U   = hbar omega [1/2 + 1/(exp(x) - 1)]       internal energy
      C_V = k_B x^2 exp(x) / (exp(x) - 1)^2         heat capacity

  The molecules and their vibrational temperatures are those of the table in Section 2.6
  (harmonic wavenumbers from the NIST Chemistry WebBook, Theta_vib = 1.43878 cm K times the wavenumber).
*/
const MOLECULES = [
  ['H₂', 6332], ['HCl', 4303], ['N₂', 3393], ['CO', 3122], ['O₂', 2274], ['Cl₂', 805], ['Br₂', 468], ['I₂', 309],
];
const START = { theta: 805, T: 298, range: 1000 };
const LEVELS = 10;                                    // bars for n = 0 to 9, and one bar for all higher levels

const energy = (x) => 0.5 + 1 / Math.expm1(x);                       // U / hbar omega
const heat = (x) => (x > 500 ? 0 : x * x * Math.exp(x) / Math.expm1(x) ** 2);   // C_V / k_B

function Tip({ active, payload, label, name, digits }) {
  if (!active || !payload || !payload.length) return null;
  return (
    <div className="tip">
      <div><V>T</V> = {Number(label).toFixed(0)} K</div>
      <div>{name} = {payload[0].payload[payload[0].dataKey].toFixed(digits)}</div>
    </div>
  );
}

function LevelTip({ active, payload }) {
  if (!active || !payload || !payload.length) return null;
  const { n, p } = payload[0].payload;
  return <div className="tip"><div><V>n</V> = {n}</div><div>probability = {p.toFixed(4)}</div></div>;
}

export default function HarmonicOscillator() {
  const { color } = useBook();
  const [theta, setTheta] = useState(START.theta);      // K
  const [T, setT] = useState(START.T);                  // K, the temperature of panel (a) and of the markers
  const [range, setRange] = useState(START.range);      // K, largest temperature on the axes

  const x = theta / T;
  const levels = useMemo(() => {
    const bars = [];
    for (let n = 0; n < LEVELS; n++) bars.push({ n: String(n), p: Math.exp(-n * x) * -Math.expm1(-x) });
    bars.push({ n: `${LEVELS}+`, p: Math.exp(-LEVELS * x) });         // sum over all levels from n = 10
    return bars;
  }, [x]);
  const curves = useMemo(() => {
    const points = [];
    for (let i = 1; i <= 240; i++) {
      const t = range * i / 240;
      points.push({ t, u: energy(theta / t), kT: t / theta, zero: 0.5, c: heat(theta / t), one: 1 });
    }
    return points;
  }, [theta, range]);

  const uTicks = niceTicks(Math.max(1.5, energy(theta / range) * 1.02), 4);
  const uStep = uTicks[1] - uTicks[0], uTop = uTicks[uTicks.length - 1] + (uTicks[uTicks.length - 1] < energy(theta / range) ? uStep : 0);
  const xAxis = { dataKey: 't', type: 'number', domain: [0, range], ticks: niceTicks(range, 5), ...axis(color) };
  const margin = { top: 14, right: 12, bottom: -6, left: 0 };
  const xTitle = <p className="xlabel">Temperature, <V>T</V> (K)</p>;
  const hw = <>ℏ<V>ω</V></>, kB = <V sub="B">k</V>, Cv = <V sub={<i className="v">V</i>}>C</V>, Th = <>Θ<sub>vib</sub></>;
  const molecule = MOLECULES.find(([, value]) => value === theta);
  const changeRange = (value) => { setRange(value); if (T > value) setT(value); };
  const pick = (value) => setTheta(value);
  // What each plot shows, for a reader who cannot see it
  const said = {
    a: `Bar chart of the probabilities of the levels n = 0 to ${LEVELS - 1}, with one more bar for all higher levels together, at `
       + `${T.toFixed(0)} kelvin for a vibrational temperature of ${theta.toFixed(0)} kelvin. The probability of level 0 is `
       + `${levels[0].p.toFixed(3)} and the probability of level 1 is ${levels[1].p.toFixed(3)}.`,
    b: `Internal energy in units of h-bar omega against temperature from 0 to ${range.toFixed(0)} kelvin. The curve starts at the `
       + `zero-point energy, one half, and approaches the line of kBT at high temperature. At ${T.toFixed(0)} kelvin the internal `
       + `energy is ${energy(x).toFixed(3)}.`,
    c: `Heat capacity in units of the Boltzmann constant against temperature from 0 to ${range.toFixed(0)} kelvin. The curve rises `
       + `from zero toward 1, its limit at high temperature. At ${T.toFixed(0)} kelvin the heat capacity is ${heat(x).toFixed(3)}.`,
  };

  return (
    <Studio title="Harmonic oscillator"
            lede="One oscillator, such as the vibration of a diatomic molecule, in contact with a heat bath. Choose its vibrational temperature and a temperature, and follow the probabilities of its levels, its internal energy, and its heat capacity.">
      <div className="panels">
        <Panel letter="a" title={<>Probabilities of the levels at <V>T</V> = {T.toFixed(0)} K</>}>
          <div className="chart" style={{ height: 190 }}>
            <ResponsiveContainer width="100%" height="100%">
              <BarChart data={levels} margin={margin} barCategoryGap="22%" aria-label={said.a}>
                <XAxis dataKey="n" {...axis(color)} />
                <YAxis domain={[0, 1]} ticks={[0, 0.25, 0.5, 0.75, 1]} width={44} {...axis(color)} />
                <Tooltip content={<LevelTip />} cursor={{ fill: color.fill }} isAnimationActive={false} />
                <Bar dataKey="p" isAnimationActive={false}>
                  {levels.map((bar, i) => <Cell key={bar.n} fill={i < LEVELS ? color.blue : color.gray} />)}
                </Bar>
              </BarChart>
            </ResponsiveContainer>
          </div>
          <p className="xlabel">Level, <V>n</V> (the last bar is the sum over all levels from <V>n</V> = {LEVELS})</p>
          <p className="note">
            <span>{Th}/<V>T</V> = {x.toFixed(2)}</span>
            <span><V sub="0">p</V> = {levels[0].p.toFixed(3)}</span>
            <span><V sub="1">p</V> = {levels[1].p.toFixed(3)}</span>
          </p>
        </Panel>
      </div>

      <div className="panels two" style={{ marginTop: 14 }}>
        <Panel letter="b" title={<>Internal energy, <V>U</V>/{hw}</>}>
          <p className="key">
            <span style={{ color: color.gray }}><i className="dash" /><span>zero-point energy and {kB}<V>T</V></span></span>
          </p>
          <div className="chart" style={{ height: 190 }}>
            <ResponsiveContainer width="100%" height="100%">
              <LineChart data={curves} margin={margin} aria-label={said.b}>
                <XAxis {...xAxis} />
                <YAxis domain={[0, uTop]} ticks={niceTicks(uTop, Math.round(uTop / uStep))} width={44} allowDataOverflow {...axis(color)} />
                <Tooltip content={<Tip name={<><V>U</V>/{hw}</>} digits={3} />} isAnimationActive={false} />
                <Line dataKey="zero" stroke={color.gray} strokeWidth={1.4} strokeDasharray={DASH} dot={false} activeDot={false} isAnimationActive={false} />
                <Line dataKey="kT" stroke={color.gray} strokeWidth={1.4} strokeDasharray={DASH} dot={false} activeDot={false} isAnimationActive={false} />
                <Line dataKey="u" stroke={color.text} strokeWidth={LINE} dot={false} isAnimationActive={false} />
                <ReferenceDot x={T} y={energy(x)} r={4.5} fill={color.vermilion} stroke={color.ground} strokeWidth={1.5} />
              </LineChart>
            </ResponsiveContainer>
          </div>
          {xTitle}
          <p className="note"><span><span style={{ color: color.vermilion }}>●</span> at {T.toFixed(0)} K: <V>U</V> = {energy(x).toFixed(3)} {hw}</span></p>
        </Panel>
        <Panel letter="c" title={<>Heat capacity, {Cv}/{kB}</>}>
          <p className="key">
            <span style={{ color: color.gray }}><i className="dash" /><span>limit at high temperature</span></span>
          </p>
          <div className="chart" style={{ height: 190 }}>
            <ResponsiveContainer width="100%" height="100%">
              <LineChart data={curves} margin={margin} aria-label={said.c}>
                <XAxis {...xAxis} />
                <YAxis domain={[0, 1.2]} ticks={[0, 0.25, 0.5, 0.75, 1]} width={44} {...axis(color)} />
                <Tooltip content={<Tip name={<>{Cv}/{kB}</>} digits={3} />} isAnimationActive={false} />
                <Line dataKey="one" stroke={color.gray} strokeWidth={1.4} strokeDasharray={DASH} dot={false} activeDot={false} isAnimationActive={false} />
                <Line dataKey="c" stroke={color.text} strokeWidth={LINE} dot={false} isAnimationActive={false} />
                <ReferenceDot x={T} y={heat(x)} r={4.5} fill={color.vermilion} stroke={color.ground} strokeWidth={1.5} />
              </LineChart>
            </ResponsiveContainer>
          </div>
          {xTitle}
          <p className="note"><span><span style={{ color: color.vermilion }}>●</span> at {T.toFixed(0)} K: {Cv} = {heat(x).toFixed(3)} {kB}</span></p>
        </Panel>
      </div>

      <div className="controls">
        <div style={{ gridColumn: '1 / -1' }}>
          <Slider label={<>Vibrational temperature, {Th}</>}
                  show={<>{theta.toFixed(0)} K{molecule ? `, as for ${molecule[0]}` : ''} (wavenumber {(theta / 1.43878).toFixed(0)} cm<sup>−1</sup>)</>}
                  value={theta} min={100} max={6500} step={1} onChange={setTheta} />
          <div className="row" style={{ justifyContent: 'flex-start', marginTop: 6, flexWrap: 'wrap' }}>
            {MOLECULES.map(([name, value]) => (
              <button type="button" key={name} aria-pressed={theta === value} onClick={() => pick(value)}>{name}</button>
            ))}
          </div>
        </div>
        <Slider label={<>Temperature, <V>T</V></>} show={<>{T.toFixed(0)} K</>} value={T} min={10} max={range} step={1} onChange={setT} />
        <Slider label="Largest temperature shown" show={<>{range.toFixed(0)} K</>} value={range} min={200} max={6000} step={100} onChange={changeRange} />
        <div className="row">
          <button type="button" onClick={() => { setTheta(START.theta); setRange(START.range); setT(START.T); }}>Reset</button>
        </div>
      </div>
      <p className="hint">
        The levels are <V sub={<i className="v">n</i>}>E</V> = {hw}(<V>n</V> + 1/2), and {Th} = {hw}/{kB}. With <V>x</V> = {Th}/<V>T</V>, the
        probability of level <V>n</V> is <V sub={<i className="v">n</i>}>p</V> = e<sup>−<V>n</V><V>x</V></sup>(1 − e<sup>−<V>x</V></sup>), the internal energy
        is <V>U</V> = {hw}[1/2 + 1/(e<sup><V>x</V></sup> − 1)], and the heat capacity is {Cv} = {kB}<V>x</V><sup>2</sup>e<sup><V>x</V></sup>/(e<sup><V>x</V></sup> − 1)<sup>2</sup>.
        The buttons set {Th} to the values of the table in Section 2.6.
      </p>
    </Studio>
  );
}
