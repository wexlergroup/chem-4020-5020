import React, { useMemo, useState } from 'react';
import { LineChart, Line, XAxis, YAxis, Tooltip, ResponsiveContainer } from 'recharts';
import { useBook, Studio, Panel, Slider, Select, Checks, V, axis, niceTicks, logTicks, LINE, DASH } from '../../studio-kit/kit.jsx';
import { R, GASES, vanDerWaals, pengRobinson, coexistence, stableZ } from './eos.js';

/*
  Real gas (Section 1.4). One mole of a gas at a fixed temperature, described by the ideal-gas law, the
  van der Waals equation, and the Peng-Robinson equation (src/eos.js). Panel (a) shows the compressibility
  factor Z = P V_m / (R T) against the pressure, and panel (b) shows the isotherm, P against V_m.

  The temperature and the largest pressure are kept as ratios to the critical values of the gas, so that a
  change of the gas keeps the reduced conditions.
*/
const START = { gas: 'co2', tr: 1.1, pr: 8, shown: { ideal: true, vdw: true, pr: true } };
const NAMES = { ideal: 'ideal gas', vdw: 'van der Waals', pr: 'Peng–Robinson' };
const POINTS = 320;

const digits3 = (x) => Number(x.toPrecision(3)).toString();              // three significant figures
const kelvin = (T) => (T < 100 ? T.toFixed(1) : T.toFixed(0));
const digits4 = (x) => Number(x.toPrecision(4)).toString();              // critical constants, as tabulated
const bar = (P) => (P < 100 ? digits3(P) : P.toFixed(0));
const Vm = <V sub="m">V</V>;
const Tc = <V sub="c">T</V>;
const Pc = <V sub="c">P</V>;

function Tip({ active, payload, label, head, unit, value, after, shown }) {
  if (!active || !payload || !payload.length) return null;
  const point = payload[0].payload;
  return (
    <div className="tip">
      <div>{head} = {digits3(Number(label))} {unit}</div>
      {['ideal', 'vdw', 'pr'].filter((key) => shown[key] && point[key] != null).map((key) => (
        <div key={key}>{NAMES[key]}: {value} = {digits3(point[key])}{after}</div>
      ))}
    </div>
  );
}

export default function RealGas() {
  const { color } = useBook();
  const [gasKey, setGasKey] = useState(START.gas);
  const [tr, setTr] = useState(START.tr);               // T / T_c
  const [pr, setPr] = useState(START.pr);               // largest pressure shown / P_c
  const [shown, setShown] = useState(START.shown);

  const gas = GASES.find((g) => g.key === gasKey);
  const T = tr * gas.Tc, top = pr * gas.Pc;
  const unit = R * gas.Tc / gas.Pc;                     // L/mol; the volume axis runs from 0.1 to 10 times this volume
  const volumes = [0.1 * unit, 10 * unit];

  const { rows, isotherm, sat, constants } = useMemo(() => {
    const models = { vdw: vanDerWaals(gas, T), pr: pengRobinson(gas, T) };
    const sat = { vdw: coexistence(models.vdw), pr: coexistence(models.pr) };

    // (a) Z against P. The pressures of coexistence are added to the grid, so that the jumps of Z are vertical.
    const pressures = [];
    for (let i = 1; i <= POINTS; i++) pressures.push(top * i / POINTS);
    for (const key of ['vdw', 'pr']) if (sat[key] && sat[key].P < top) pressures.push(sat[key].P * (1 - 1e-7), sat[key].P * (1 + 1e-7));
    pressures.sort((x, y) => x - y);
    const rows = pressures.map((P) => ({ P, ideal: 1, vdw: stableZ(models.vdw, P, sat.vdw), pr: stableZ(models.pr, P, sat.pr) }));

    // (b) P against V_m. Between the volumes of the liquid and the vapor, the stable isotherm is the tie line.
    const grid = [];
    for (let i = 0; i <= POINTS; i++) grid.push(volumes[0] * 100 ** (i / POINTS));
    for (const key of ['vdw', 'pr']) if (sat[key]) grid.push(sat[key].Vl, sat[key].Vg);
    grid.sort((x, y) => x - y);
    const isotherm = grid.map((volume) => {
      const row = { V: volume, ideal: R * T / volume };
      for (const key of ['vdw', 'pr']) {
        const model = models[key], s = sat[key];
        if (volume <= model.b * 1.001) { row[key] = null; row[key + 'Loop'] = null; continue; }
        const inside = s && volume >= s.Vl && volume <= s.Vg;
        row[key] = inside ? s.P : model.pressure(volume);
        row[key + 'Loop'] = inside ? model.pressure(volume) : null;
      }
      return row;
    });
    return { rows, isotherm, sat, constants: { a: models.vdw.a, b: models.vdw.b } };
  }, [gasKey, tr, pr]);

  // Top of the Z axis: above the largest value of the lines shown, and at least 1.2
  const largest = Math.max(1.2, ...rows.map((row) => Math.max(shown.vdw ? row.vdw : 0, shown.pr ? row.pr : 0)));
  const zTicks = niceTicks(largest * 1.0001, 5);
  const zStep = zTicks[1] - zTicks[0], zTop = Math.ceil(largest / zStep - 1e-9) * zStep;
  const zAll = niceTicks(zTop, Math.round(zTop / zStep));
  const pTicks = niceTicks(top, 5);

  const margin = { top: 14, right: 14, bottom: -6, left: 0 };
  const lines = (suffix) => [
    shown.ideal && !suffix ? <Line key="ideal" dataKey="ideal" stroke={color.gray} strokeWidth={1.6} strokeDasharray={DASH} dot={false} isAnimationActive={false} /> : null,
    shown.vdw ? <Line key={'vdw' + suffix} dataKey={'vdw' + suffix} stroke={color.blue} strokeWidth={suffix ? 1.1 : LINE} dot={false} activeDot={!suffix} isAnimationActive={false} /> : null,
    shown.pr ? <Line key={'pr' + suffix} dataKey={'pr' + suffix} stroke={color.vermilion} strokeWidth={suffix ? 1.1 : LINE} dot={false} activeDot={!suffix} isAnimationActive={false} /> : null,
  ];
  const below = sat.vdw || sat.pr;
  const coexist = ['vdw', 'pr'].filter((key) => shown[key] && sat[key]);
  // What each plot shows, for a reader who cannot see it
  const real = ['vdw', 'pr'].filter((key) => shown[key]);
  const names = [shown.ideal ? 'the ideal-gas law' : null, ...real.map((key) => `the ${NAMES[key]} equation`)].filter(Boolean);
  const list = names.length > 1 ? `${names.slice(0, -1).join(', ')} and ${names[names.length - 1]}` : names[0] || 'no equation';
  const lowest = (key) => rows.reduce((best, row) => (row[key] < best[key] ? row : best), rows[0]);
  const said = {
    a: `Compressibility factor Z of ${gas.name.toLowerCase()} at ${kelvin(T)} kelvin against pressure from 0 to ${bar(top)} bar, for ${list}.`
       + (shown.ideal ? ' For the ideal gas, Z = 1 at every pressure.' : '')
       + real.map((key) => ` With the ${NAMES[key]} equation, Z is lowest, ${lowest(key)[key].toFixed(2)}, at ${bar(lowest(key).P)} bar, and it is `
                           + `${rows[rows.length - 1][key].toFixed(2)} at the largest pressure shown.`).join('')
       + (coexist.length ? ' Z jumps at the pressure at which liquid and vapor coexist.' : ''),
    b: `Isotherm of ${gas.name.toLowerCase()} at ${kelvin(T)} kelvin: pressure from 0 to ${bar(top)} bar against molar volume on a logarithmic `
       + `scale from ${digits3(volumes[0])} to ${digits3(volumes[1])} liters per mole, for ${list}. The pressure falls as the molar volume grows.`
       + (real.length ? ' At large molar volume the curves of the real-gas equations approach the ideal-gas isotherm, and at small molar volume they rise steeply.' : '')
       + (coexist.length ? ` The isotherm is horizontal where liquid and vapor coexist, at ${coexist.map((key) => `${bar(sat[key].P)} bar (${NAMES[key]})`).join(' and ')}.` : ''),
  };

  return (
    <Studio title="Real gas"
            lede="Choose a gas and a temperature, and compare the ideal-gas law with two equations of state for real gases, the van der Waals equation and the Peng–Robinson equation.">
      <p className="key">
        {shown.ideal ? <span style={{ color: color.gray }}><i className="dash" /><span>ideal gas</span></span> : null}
        {shown.vdw ? <span style={{ color: color.blue }}><i className="line" /><span>van der Waals</span></span> : null}
        {shown.pr ? <span style={{ color: color.vermilion }}><i className="line" /><span>Peng–Robinson</span></span> : null}
      </p>
      <div className="panels two" style={{ marginTop: 8 }}>
        <Panel letter="a" title={<>Compressibility factor, <V>Z</V></>}>
          <div className="chart" style={{ height: 250 }}>
            <ResponsiveContainer width="100%" height="100%">
              <LineChart data={rows} margin={margin} aria-label={said.a}>
                <XAxis dataKey="P" type="number" domain={[0, top]} ticks={pTicks} {...axis(color)} />
                <YAxis domain={[0, zTop]} ticks={zAll} width={44} allowDataOverflow {...axis(color)} />
                <Tooltip content={<Tip head={<V>P</V>} unit="bar" value={<V>Z</V>} shown={shown} />} isAnimationActive={false} />
                {lines('')}
              </LineChart>
            </ResponsiveContainer>
          </div>
          <p className="xlabel">Pressure, <V>P</V> (bar)</p>
        </Panel>
        <Panel letter="b" title={<>Isotherm: pressure, <V>P</V> (bar)</>}>
          <div className="chart" style={{ height: 250 }}>
            <ResponsiveContainer width="100%" height="100%">
              <LineChart data={isotherm} margin={margin} aria-label={said.b}>
                <XAxis dataKey="V" type="number" scale="log" domain={volumes} ticks={logTicks(volumes[0], volumes[1])}
                       allowDataOverflow {...axis(color)} />
                <YAxis domain={[0, top]} ticks={pTicks} width={44} allowDataOverflow {...axis(color)} />
                <Tooltip content={<Tip head={Vm} unit="L/mol" value={<V>P</V>} after=" bar" shown={shown} />} isAnimationActive={false} />
                {lines('Loop')}
                {lines('')}
              </LineChart>
            </ResponsiveContainer>
          </div>
          <p className="xlabel">Molar volume, {Vm} (L/mol), logarithmic scale</p>
        </Panel>
      </div>
      <p className="note">
        <span>{gas.name}: {Tc} = {digits4(gas.Tc)} K, {Pc} = {digits4(gas.Pc)} bar</span>
        <span>van der Waals constants: <V sub="m">a</V> = {digits3(constants.a)} L<sup>2</sup> bar mol<sup>−2</sup>, <V sub="m">b</V> = {digits3(constants.b)} L mol<sup>−1</sup></span>
      </p>
      <p className="note lines3 live">
        {below ? (
          <span>
            Below {Tc}, liquid and vapor coexist at one pressure
            {coexist.length ? <>: {coexist.map((key) => `${bar(sat[key].P)} bar (${NAMES[key]})`).join(' and ')}</> : null}.
            There <V>Z</V> jumps in (a) and the isotherm is horizontal in (b). The thin lines in (b) are the loops of the two equations.
          </span>
        ) : (
          <span>At this temperature, each equation gives one molar volume at every pressure. Below {Tc}, liquid and vapor coexist at one pressure.</span>
        )}
      </p>

      <div className="controls">
        <Select label="Gas" value={gasKey} options={GASES.map((g) => [g.key, g.name])} onChange={setGasKey} />
        <Checks label="Equations shown" onChange={(key, on) => setShown({ ...shown, [key]: on })}
                items={[['ideal', 'ideal gas', shown.ideal], ['vdw', 'van der Waals', shown.vdw], ['pr', 'Peng–Robinson', shown.pr]]} />
        <Slider label={<>Temperature, <V>T</V></>} show={<>{kelvin(T)} K = {tr.toFixed(2)} {Tc}</>}
                value={tr} min={0.7} max={4} step={0.01} onChange={setTr} />
        <Slider label={<>Largest pressure shown</>} show={<>{bar(top)} bar = {pr.toFixed(1)} {Pc}</>}
                value={pr} min={1} max={12} step={0.5} onChange={setPr} />
        <div className="row">
          <button type="button" onClick={() => { setGasKey(START.gas); setTr(START.tr); setPr(START.pr); setShown(START.shown); }}>Reset</button>
        </div>
      </div>
      <p className="hint">
        The compressibility factor is <V>Z</V> = <V>P</V>{Vm}/<V>R</V><V>T</V>, which is 1 for an ideal gas. The van der Waals
        equation is (<V>P</V> + <V sub="m">a</V>/{Vm}<sup>2</sup>)({Vm} − <V sub="m">b</V>) = <V>R</V><V>T</V>, with <V sub="m">a</V> = 27<V>R</V><sup>2</sup>{Tc}<sup>2</sup>/64{Pc} and <V sub="m">b</V> = <V>R</V>{Tc}/8{Pc}.
        The Peng–Robinson equation, <V>P</V> = <V>R</V><V>T</V>/({Vm} − <V>b</V>) − <V>a</V>(<V>T</V>)/({Vm}<sup>2</sup> + 2<V>b</V>{Vm} − <V>b</V><sup>2</sup>), is
        a later equation of the same kind. Its constants also follow from {Tc} and {Pc}, and <V>a</V>(<V>T</V>) also depends on the acentric factor of the gas.
        The temperature and the largest pressure are set as multiples of {Tc} and {Pc}, and they keep these multiples when you change the gas.
      </p>
    </Studio>
  );
}
