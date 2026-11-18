import React, { useMemo, useState } from 'react';
import { LineChart, Line, XAxis, YAxis, Tooltip, ResponsiveContainer, ReferenceLine } from 'recharts';
import { useBook, Studio, Panel, Slider, V, axis, niceTicks, LINE, DASH } from '../../studio-kit/kit.jsx';

/*
  Two-state system (Section 2.2). One particle has two states: state 1 with the energy 0 and state 2 with the
  energy dE. With x = dE / (k_B T):

      Q   = 1 + exp(-x)                         partition function
      p_1 = 1 / Q,  p_2 = exp(-x) / Q           probabilities of the two states
      C_V = k_B x^2 exp(-x) / (1 + exp(-x))^2   heat capacity

  Energies and k_B T are in units of epsilon, an arbitrary unit of energy.
*/

// C_V is largest where (x/2) tanh(x/2) = 1, that is, at x = 2.39936
const PEAK_T = 1 / 2.39936;        // k_B T / dE at the maximum of C_V
const PEAK_C = 0.43923;            // C_V / k_B at the maximum
const START = { gap: 2.0, range: 5.0 };

function Tip({ active, payload, label, rows }) {
  if (!active || !payload || !payload.length) return null;
  const point = payload[0].payload;
  return (
    <div className="tip">
      <div><V sub="B">k</V><V>T</V> = {Number(label).toFixed(2)} <V>ε</V></div>
      {rows.map(([name, key, digits]) => <div key={key}>{name} = {point[key].toFixed(digits)}</div>)}
    </div>
  );
}

export default function TwoStateSystem() {
  const { color } = useBook();
  const [gap, setGap] = useState(START.gap);          // energy of state 2 above state 1
  const [range, setRange] = useState(START.range);    // largest k_B T on the axes

  const data = useMemo(() => {
    const points = [];
    for (let i = 1; i <= 240; i++) {
      const t = range * i / 240, e = Math.exp(-gap / t), q = 1 + e;
      points.push({ t, q, p1: 1 / q, p2: e / q, c: (gap / t) ** 2 * e / (q * q) });
    }
    return points;
  }, [gap, range]);

  const peak = PEAK_T * gap;                          // k_B T at the maximum of C_V
  const x = { dataKey: 't', type: 'number', domain: [0, range], ticks: niceTicks(range, 5), ...axis(color) };
  const margin = { top: 14, right: 12, bottom: -6, left: 0 };
  const xTitle = <p className="xlabel">Temperature, <V sub="B">k</V><V>T</V>/<V>ε</V></p>;
  const p1 = <><V sub="1">p</V></>, p2 = <><V sub="2">p</V></>;
  // What each plot shows, for a reader who cannot see it
  const end = data[data.length - 1];
  const axes = `against temperature, with kBT from 0 to ${range.toFixed(1)} epsilon, for an energy difference of ${gap.toFixed(1)} epsilon`;
  const said = {
    a: `Probabilities of the two states ${axes}. The probability of state 1 falls from 1 and the probability of state 2 rises `
       + `from 0, and both approach one half. At the largest temperature shown they are ${end.p1.toFixed(2)} and ${end.p2.toFixed(2)}.`,
    b: `Heat capacity in units of the Boltzmann constant ${axes}. `
       + (peak <= range
         ? `It rises from zero to its largest value, ${PEAK_C.toFixed(3)}, at kBT = ${peak.toFixed(2)} epsilon, and then falls.`
         : `It rises from zero to ${end.c.toFixed(3)} at the largest temperature shown. Its largest value, ${PEAK_C.toFixed(3)}, lies at kBT = ${peak.toFixed(2)} epsilon, beyond the temperatures shown.`),
    c: `Partition function ${axes}. It rises from 1 toward 2, and it is ${end.q.toFixed(2)} at the largest temperature shown.`,
  };

  return (
    <Studio title="Two-state system"
            lede="A particle has two states, state 1 and state 2. Change the energy difference between them, and follow the probabilities of the two states, the heat capacity, and the partition function as the temperature rises.">
      <div className="panels">
        <Panel letter="a" title="Probabilities of the two states">
          <p className="key">
            <span style={{ color: color.blue }}><i className="line" /><span>state 1, {p1}</span></span>
            <span style={{ color: color.vermilion }}><i className="line" /><span>state 2, {p2}</span></span>
          </p>
          <div className="chart" style={{ height: 220 }}>
            <ResponsiveContainer width="100%" height="100%">
              <LineChart data={data} margin={margin} aria-label={said.a}>
                <XAxis {...x} />
                <YAxis domain={[0, 1]} ticks={[0, 0.25, 0.5, 0.75, 1]} width={44} {...axis(color)} />
                <Tooltip content={<Tip rows={[[p1, 'p1', 3], [p2, 'p2', 3]]} />} isAnimationActive={false} />
                <ReferenceLine y={0.5} stroke={color.gray} strokeDasharray={DASH} strokeWidth={1.2} />
                <Line dataKey="p1" stroke={color.blue} strokeWidth={LINE} dot={false} isAnimationActive={false} />
                <Line dataKey="p2" stroke={color.vermilion} strokeWidth={LINE} dot={false} isAnimationActive={false} />
              </LineChart>
            </ResponsiveContainer>
          </div>
          {xTitle}
        </Panel>
      </div>

      <div className="panels two" style={{ marginTop: 14 }}>
        <Panel letter="b" title={<>Heat capacity, <V sub={<i className="v">V</i>}>C</V>/<V sub="B">k</V></>}>
          <div className="chart" style={{ height: 190 }}>
            <ResponsiveContainer width="100%" height="100%">
              <LineChart data={data} margin={margin} aria-label={said.b}>
                <XAxis {...x} />
                <YAxis domain={[0, 0.5]} ticks={[0, 0.1, 0.2, 0.3, 0.4, 0.5]} width={44} {...axis(color)} />
                <Tooltip content={<Tip rows={[[<><V sub={<i className="v">V</i>}>C</V>/<V sub="B">k</V></>, 'c', 3]]} />} isAnimationActive={false} />
                {peak <= range ? <ReferenceLine x={peak} stroke={color.gray} strokeDasharray={DASH} strokeWidth={1.2} /> : null}
                <Line dataKey="c" stroke={color.text} strokeWidth={LINE} dot={false} isAnimationActive={false} />
              </LineChart>
            </ResponsiveContainer>
          </div>
          {xTitle}
          <p className="note">
            <span>Largest at <V sub="B">k</V><V>T</V> = {PEAK_T.toFixed(3)} Δ<V>E</V> = {peak.toFixed(2)} <V>ε</V> (dashed line), where <V sub={<i className="v">V</i>}>C</V> = {PEAK_C.toFixed(3)} <V sub="B">k</V></span>
          </p>
        </Panel>
        <Panel letter="c" title={<>Partition function, <V>Q</V></>}>
          <div className="chart" style={{ height: 190 }}>
            <ResponsiveContainer width="100%" height="100%">
              <LineChart data={data} margin={margin} aria-label={said.c}>
                <XAxis {...x} />
                <YAxis domain={[1, 2]} ticks={[1, 1.25, 1.5, 1.75, 2]} width={44} {...axis(color)} />
                <Tooltip content={<Tip rows={[[<V>Q</V>, 'q', 3]]} />} isAnimationActive={false} />
                <Line dataKey="q" stroke={color.text} strokeWidth={LINE} dot={false} isAnimationActive={false} />
              </LineChart>
            </ResponsiveContainer>
          </div>
          {xTitle}
        </Panel>
      </div>

      <div className="controls">
        <Slider label={<>Energy of state 2 above state 1, Δ<V>E</V></>} show={<>{gap.toFixed(1)} <V>ε</V></>}
                value={gap} min={0.5} max={10} step={0.1} onChange={setGap} />
        <Slider label={<>Largest temperature shown, <V sub="B">k</V><V>T</V></>} show={<>{range.toFixed(1)} <V>ε</V></>}
                value={range} min={1} max={20} step={0.5} onChange={setRange} />
        <div className="row"><button type="button" onClick={() => { setGap(START.gap); setRange(START.range); }}>Reset</button></div>
      </div>
      <p className="hint">
        Energies are in units of an arbitrary energy <V>ε</V>. State 1 has the energy 0, and state 2 has the energy Δ<V>E</V>.
        With <V>x</V> = Δ<V>E</V>/<V sub="B">k</V><V>T</V>, the partition function is <V>Q</V> = 1 + e<sup>−<V>x</V></sup>, the
        probabilities are <V sub="1">p</V> = 1/<V>Q</V> and <V sub="2">p</V> = e<sup>−<V>x</V></sup>/<V>Q</V>, and the heat capacity
        of one particle is <V sub={<i className="v">V</i>}>C</V> = <V sub="B">k</V> <V>x</V><sup>2</sup> e<sup>−<V>x</V></sup>/(1 + e<sup>−<V>x</V></sup>)<sup>2</sup>.
        The maximum of this heat capacity is known as the Schottky anomaly.
      </p>
    </Studio>
  );
}
