import React, { useState } from 'react';
import { ComposedChart, Scatter, Line, XAxis, YAxis, Tooltip, ResponsiveContainer, ReferenceLine } from 'recharts';
import { useBook, Studio, Panel, Slider, Select, V, axis, LINE, DASH } from '../../studio-kit/kit.jsx';
import { LIQUIDS } from './data.js';
import { fit, referenceAt } from './fit.js';

/*
  Clausius–Clapeyron (Section 6.2). The vapor pressure of a liquid (a), the Clausius–Clapeyron plot of ln P_sat
  against 1/T with a straight line fitted to the points that the student chooses (b), the residuals of the fit (c),
  and the enthalpy of vaporization from the slope beside the reference value (d). The data come from the reference
  equation of state of each liquid (src/data.js); for water they are the table of the mini-lab of the section.
*/
const START = { liquid: 'water' };
const HEIGHT = 230;
const margin = { top: 12, right: 14, bottom: -6, left: 0 };

function ticksBetween(lo, hi, n) {                    // round ticks that cover lo to hi
  const rough = (hi - lo) / n, power = 10 ** Math.floor(Math.log10(rough));
  const step = [1, 2, 5, 10].map((m) => m * power).find((s) => s >= rough - 1e-12);
  const ticks = [];
  for (let v = Math.floor(lo / step) * step; v <= hi + step * 1e-6; v += step) ticks.push(Math.round(v / step) * step);
  if (ticks[ticks.length - 1] < hi - 1e-9) ticks.push(ticks[ticks.length - 1] + step);
  return ticks.map((v) => Number(v.toPrecision(10)));
}
const fixed = (step) => Math.max(0, -Math.floor(Math.log10(step) + 1e-9));

function Dot({ cx, cy, payload, color }) {            // filled when the point is in the fit, open when it is not
  if (cx == null || cy == null) return null;
  return payload.used
    ? <circle cx={cx} cy={cy} r={4.6} fill={color.blue} />
    : <circle cx={cx} cy={cy} r={4.2} fill={color.ground} stroke={color.blue} strokeWidth={1.6} />;
}

function Tip({ active, payload }) {
  if (!active || !payload || !payload.length) return null;
  const p = payload[0].payload;
  if (p.t == null) return null;
  return (
    <div className="tip">
      <div>{String(p.t).replace('-', '−')} °C, <V>T</V> = {p.T.toFixed(2)} K</div>
      <div><V sub="sat">P</V> = {p.P.toPrecision(5)} kPa</div>
      <div>ln(<V sub="sat">P</V>/1 bar) = {p.y.toFixed(3)}</div>
      <div>Residual = {p.residual >= 0 ? '+' : '−'}{Math.abs(p.residual).toFixed(4)}</div>
    </div>
  );
}

export default function ClausiusClapeyron() {
  const { color } = useBook();
  const [key, setKey] = useState(START.liquid);
  const liquid = LIQUIDS.find((l) => l.key === key);
  const last = liquid.points.length - 1;
  const [range, setRange] = useState([0, last]);
  const [from, to] = range;
  const f = fit(liquid.points, from, to);
  const points = f.all.map((p, i) => ({ ...p, used: i >= from && i <= to, inv: 1000 * p.x }));
  const reference = referenceAt(liquid.reference, f.Tmean);

  const choose = (k) => { setKey(k); setRange([0, LIQUIDS.find((l) => l.key === k).points.length - 1]); };
  const setFrom = (i) => setRange([Math.min(i, to - 1), to]);
  const setTo = (i) => setRange([from, Math.max(i, from + 1)]);

  // (a) P_sat against t, with the fitted line turned back into a curve over the range of the fit
  const tTicks = ticksBetween(points[0].t, points[last].t, 5);
  const pTicks = ticksBetween(0, Math.max(...points.map((p) => p.P)) * 1.03, 6);
  const curve = [];
  for (let i = 0; i <= 60; i++) {
    const t = points[from].t + (points[to].t - points[from].t) * i / 60;
    curve.push({ t, fitP: 100 * Math.exp(f.slope / (t + 273.15) + f.intercept) });
  }
  // (b) ln P against 1000/T; (c) the residuals
  const xs = points.map((p) => p.inv), ys = points.map((p) => p.y);
  const xTicks = ticksBetween(Math.min(...xs), Math.max(...xs), 5);
  const yTicks = ticksBetween(Math.min(...ys), Math.max(...ys), 5);
  const line = [points[to], points[from]].map((p) => ({ inv: p.inv, fitY: f.slope * p.x + f.intercept }));
  const big = Math.max(...points.map((p) => Math.abs(p.residual)), 1e-4);
  const rTicks = ticksBetween(-big, big, 4);
  const xStep = xTicks[1] - xTicks[0], yStep = yTicks[1] - yTicks[0], rStep = rTicks[1] - rTicks[0];

  const x = <><V>T</V></>, dh = <>Δ<V bar>h</V><sub>vap</sub></>;
  const pt = (i) => `${String(points[i].t).replace('-', '−')} °C`;
  const rows = [
    ['Points in the fit', <>{pt(from)} to {pt(to)}</>, f.n],
    ['Slope of the line', <>−{dh}/<V>R</V></>, <>{Math.round(f.slope).toLocaleString('en-US').replace('-', '−')} K</>],
    [<>Coefficient of determination, <V>R</V><sup>2</sup></>, null, f.r2.toFixed(6)],
    [<>Enthalpy of vaporization from the slope</>, <>{dh} = −<V>R</V> × slope</>, <>{f.dh.toFixed(2)} kJ mol<sup>−1</sup></>],
    [<>Reference value of {dh}</>,
     <>at {f.Tmean.toFixed(1)} K ({(f.Tmean - 273.15).toFixed(1)} °C), where 1/{x} is the mean of 1/{x} over the points in the fit</>,
     <>{reference.toFixed(2)} kJ mol<sup>−1</sup></>],
  ];

  const scatter = (dataKey) => (
    <Scatter data={points} dataKey={dataKey} shape={(p) => <Dot {...p} color={color} />} isAnimationActive={false} />
  );

  return (
    <Studio title="Clausius–Clapeyron"
            lede="Fit a straight line to the Clausius–Clapeyron plot of a liquid over a range of temperatures that you choose, and compare the enthalpy of vaporization from its slope with the reference value.">
      <p className="key">
        <span style={{ color: color.blue }}><svg width="12" height="12"><circle cx="6" cy="6" r="4.6" fill={color.blue} /></svg><span>point in the fit</span></span>
        <span style={{ color: color.blue }}><svg width="12" height="12"><circle cx="6" cy="6" r="4.2" fill="none" stroke={color.blue} strokeWidth="1.6" /></svg><span>point not in the fit</span></span>
        <span style={{ color: color.vermilion }}><i className="line" /><span>fitted line</span></span>
      </p>
      <div className="panels two" style={{ marginTop: 8 }}>
        <Panel letter="a" title={<>Vapor pressure, <V sub="sat">P</V> (kPa)</>}>
          <div className="chart" style={{ height: HEIGHT }}>
            <ResponsiveContainer width="100%" height="100%">
              <ComposedChart margin={margin}>
                <XAxis dataKey="t" type="number" domain={[tTicks[0], tTicks[tTicks.length - 1]]} ticks={tTicks} allowDataOverflow
                       tickFormatter={(v) => String(v).replace('-', '−')} {...axis(color)} />
                <YAxis type="number" domain={[0, pTicks[pTicks.length - 1]]} ticks={pTicks} width={48} allowDataOverflow {...axis(color)} />
                <Tooltip content={<Tip />} isAnimationActive={false} cursor={false} />
                <Line data={curve} dataKey="fitP" stroke={color.vermilion} strokeWidth={LINE} dot={false} activeDot={false} isAnimationActive={false} />
                {scatter('P')}
              </ComposedChart>
            </ResponsiveContainer>
          </div>
          <p className="xlabel">Temperature, <V>T</V> (°C)</p>
        </Panel>
        <Panel letter="b" title={<>Clausius–Clapeyron plot, ln(<V sub="sat">P</V>/1 bar)</>}>
          <div className="chart" style={{ height: HEIGHT }}>
            <ResponsiveContainer width="100%" height="100%">
              <ComposedChart margin={margin}>
                <XAxis dataKey="inv" type="number" domain={[xTicks[0], xTicks[xTicks.length - 1]]} ticks={xTicks} allowDataOverflow
                       tickFormatter={(v) => v.toFixed(fixed(xStep))} {...axis(color)} />
                <YAxis type="number" domain={[yTicks[0], yTicks[yTicks.length - 1]]} ticks={yTicks} width={48} allowDataOverflow
                       tickFormatter={(v) => v.toFixed(fixed(yStep)).replace('-', '−')} {...axis(color)} />
                <Tooltip content={<Tip />} isAnimationActive={false} cursor={false} />
                <Line data={line} dataKey="fitY" stroke={color.vermilion} strokeWidth={LINE} dot={false} activeDot={false} isAnimationActive={false} />
                {scatter('y')}
              </ComposedChart>
            </ResponsiveContainer>
          </div>
          <p className="xlabel">1/<V>T</V> (10<sup>−3</sup> K<sup>−1</sup>)</p>
        </Panel>
      </div>
      <div className="panels two" style={{ marginTop: 14 }}>
        <Panel letter="c" title={<>Residuals, data − fit, in ln(<V sub="sat">P</V>/1 bar)</>}>
          <div className="chart" style={{ height: HEIGHT }}>
            <ResponsiveContainer width="100%" height="100%">
              <ComposedChart margin={margin}>
                <XAxis dataKey="inv" type="number" domain={[xTicks[0], xTicks[xTicks.length - 1]]} ticks={xTicks} allowDataOverflow
                       tickFormatter={(v) => v.toFixed(fixed(xStep))} {...axis(color)} />
                <YAxis type="number" domain={[rTicks[0], rTicks[rTicks.length - 1]]} ticks={rTicks} width={48} allowDataOverflow
                       tickFormatter={(v) => (Math.abs(v) < 1e-12 ? '0' : v.toFixed(fixed(rStep)).replace('-', '−'))} {...axis(color)} />
                <Tooltip content={<Tip />} isAnimationActive={false} cursor={false} />
                <ReferenceLine y={0} stroke={color.gray} strokeWidth={1.2} strokeDasharray={DASH} />
                {scatter('residual')}
              </ComposedChart>
            </ResponsiveContainer>
          </div>
          <p className="xlabel">1/<V>T</V> (10<sup>−3</sup> K<sup>−1</sup>)</p>
        </Panel>
        <Panel letter="d" title={`${liquid.name}: the fit`}>
          <table className="values">
            <tbody>
              {rows.map(([name, formula, value], i) => (
                <tr key={i}><td>{name}{formula ? <span className="formula">{formula}</span> : null}</td><td>{value}</td></tr>
              ))}
            </tbody>
          </table>
        </Panel>
      </div>

      <div className="controls">
        <Slider label="Lowest temperature in the fit" show={pt(from)} value={from} min={0} max={last} step={1} onChange={setFrom} />
        <Slider label="Highest temperature in the fit" show={pt(to)} value={to} min={0} max={last} step={1} onChange={setTo} />
        <Select label="Liquid" value={key} options={LIQUIDS.map((l) => [l.key, l.name])} onChange={choose} />
        <div className="row" style={{ gridColumn: 'auto', alignSelf: 'end' }}>
          <button type="button" onClick={() => setRange([0, last])}>All points</button>
        </div>
      </div>
      <p className="hint">
        The fit is a straight line through the filled points by least squares, and the enthalpy of vaporization comes from its slope,
        as in Part A of the mini-lab. The data are vapor pressures from the reference equation of state of each liquid, evaluated with
        CoolProp: for water the IAPWS-95 formulation, which gives the table of the mini-lab, for ethanol the equation of Schroeder and
        co-workers (2014), and for benzene that of Thol and co-workers (2012). The data are smooth, so the residuals show the curvature of
        the plot and no scatter of measurement.
      </p>
    </Studio>
  );
}
