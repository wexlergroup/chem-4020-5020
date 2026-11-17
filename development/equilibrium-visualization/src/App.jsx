import React, { useState } from 'react';
import { LineChart, Line, XAxis, YAxis, CartesianGrid, ReferenceLine, ResponsiveContainer, ReferenceDot } from 'recharts';

const EquilibriumVis = () => {
  // State for user controls
  const [temperature, setTemperature] = useState(298.15); // Kelvin
  const [pressure, setPressure] = useState(1); // Bar
  const [currentXi, setCurrentXi] = useState(0.5); // Extent of reaction (0 to 1)

  // 2 NO2(g) <=> N2O4(g), starting from 2 mol of NO2.
  // Standard reaction enthalpy and entropy at 298.15 K (NIST-JANAF tables),
  // taken as constant over the temperature range of the slider.
  const dH_r = -57.112; // kJ/mol (gives -4.73 kJ/mol for the standard Gibbs energy at 298.15 K)
  const dS_r = -0.17569; // kJ/(mol K)
  const R = 0.008314462618; // kJ/(mol K)

  const RT = R * temperature;
  const dG_r_not = dH_r - temperature * dS_r;
  const K = Math.exp(-dG_r_not / RT);

  // G(xi) - G(0) in kJ, where G(0) is pure NO2 at the same T and P:
  // xi * (dG_r_not - RT ln P) + RT * sum_i n_i ln y_i   (P in bar, standard pressure 1 bar)
  const gRel = (xi) => {
    const n_NO2 = 2 - 2 * xi;
    const n_N2O4 = xi;
    const n_total = 2 - xi;
    let mix = 0;
    if (n_NO2 > 0) mix += n_NO2 * Math.log(n_NO2 / n_total);
    if (n_N2O4 > 0) mix += n_N2O4 * Math.log(n_N2O4 / n_total);
    return xi * (dG_r_not - RT * Math.log(pressure)) + RT * mix;
  };

  // The G vs xi curve
  const data = [];
  for (let i = 0; i <= 200; i += 1) {
    const xi = i / 200;
    data.push({ xi, G: gRel(xi) });
  }

  // Equilibrium extent from K P = xi (2 - xi) / (2 - 2 xi)^2
  const xiEq = 1 - 1 / Math.sqrt(1 + 4 * K * pressure);
  const gEq = gRel(xiEq);

  // Vertical range of the plot, fixed by the curve so that the tangent line does not rescale it
  const gMin = Math.min(...data.map((p) => p.G));
  const gMax = Math.max(...data.map((p) => p.G));
  const pad = 0.08 * (gMax - gMin);
  const yDomain = [gMin - pad, gMax + pad];
  // Round tick values (1, 2, or 5 times a power of ten apart)
  const rough = (yDomain[1] - yDomain[0]) / 5;
  const mag = Math.pow(10, Math.floor(Math.log10(rough)));
  const step = [1, 2, 5, 10].map((m) => m * mag).find((v) => v >= rough);
  const yTicks = [];
  for (let v = Math.ceil(yDomain[0] / step) * step; v <= yDomain[1] + 1e-9; v += step) {
    yTicks.push(Math.abs(v) < 1e-9 ? 0 : v);
  }

  // Current state
  const y_NO2 = (2 - 2 * currentXi) / (2 - currentXi);
  const y_N2O4 = currentXi / (2 - currentXi);
  const Q = y_N2O4 / (y_NO2 * y_NO2 * pressure);
  const delta_r_G = RT * Math.log(Q / K); // slope of G(xi), kJ/mol
  const gNow = gRel(currentXi);

  let slopeText = 'At the minimum';
  let slopeColor = 'text-blue-600';
  if (currentXi < xiEq - 0.005) {
    slopeText = 'Forward is spontaneous';
    slopeColor = 'text-green-700';
  } else if (currentXi > xiEq + 0.005) {
    slopeText = 'Reverse is spontaneous';
    slopeColor = 'text-red-600';
  }

  // Tangent line at the current point: its slope is the reaction Gibbs energy
  const half = 0.12;
  const tangent = [
    { xi: Math.max(0, currentXi - half), G: gNow - Math.min(half, currentXi) * delta_r_G },
    { xi: Math.min(1, currentXi + half), G: gNow + Math.min(half, 1 - currentXi) * delta_r_G },
  ];

  const fmt = (x) => (x >= 1000 || x < 0.01 ? x.toExponential(1) : x >= 100 ? x.toFixed(0) : x.toFixed(2));

  return (
    <div className="w-full max-w-4xl mx-auto p-2 bg-white font-sans text-slate-800">
      <div className="flex flex-wrap items-baseline gap-x-3 border-b border-slate-200 pb-1 mb-2">
        <h2 className="text-base font-bold">Minimizing Gibbs Free Energy</h2>
        <span className="text-xs text-slate-500">2 NO₂(g) ⇌ N₂O₄(g), starting from 2 mol of NO₂</span>
      </div>

      <div className="grid grid-cols-1 sm:grid-cols-5 gap-3">
        {/* Left column: plot and extent slider */}
        <div className="sm:col-span-3">
          <div className="bg-slate-50 rounded-lg border border-slate-200 relative" style={{ height: 284 }}>
            <ResponsiveContainer width="100%" height="100%">
              <LineChart margin={{ top: 12, right: 16, left: 6, bottom: 18 }}>
                <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0" />
                <XAxis
                  dataKey="xi"
                  type="number"
                  domain={[0, 1]}
                  ticks={[0, 0.2, 0.4, 0.6, 0.8, 1]}
                  tick={{ fontSize: 11 }}
                  label={{ value: 'Extent of reaction ξ (mol)', position: 'bottom', offset: 2, fontSize: 12 }}
                />
                <YAxis
                  domain={yDomain}
                  allowDataOverflow={true}
                  tick={{ fontSize: 11 }}
                  ticks={yTicks}
                  tickFormatter={(v) => v.toFixed(step < 1 ? 1 : 0)}
                  width={46}
                  label={{ value: 'G(ξ) − G(0) (kJ)', angle: -90, position: 'insideLeft', offset: 4, fontSize: 12, style: { textAnchor: 'middle' } }}
                />
                <ReferenceLine x={xiEq} stroke="#059669" strokeDasharray="4 3" />
                <Line data={data} type="linear" dataKey="G" stroke="#2563eb" strokeWidth={2.5} dot={false} isAnimationActive={false} />
                <Line data={tangent} type="linear" dataKey="G" stroke="#ef4444" strokeWidth={1.5} dot={false} isAnimationActive={false} />
                <ReferenceDot x={xiEq} y={gEq} r={4} fill="#059669" stroke="#fff" strokeWidth={1} />
                <ReferenceDot x={currentXi} y={gNow} r={6} fill="#ef4444" stroke="#fff" strokeWidth={2} />
              </LineChart>
            </ResponsiveContainer>
            <div className="absolute top-2 left-1/2 -translate-x-1/2 bg-white/90 px-2 py-0.5 rounded text-xs text-emerald-700 border border-emerald-200">
              minimum at ξ = {xiEq.toFixed(2)}
            </div>
          </div>

          <div className="bg-slate-100 px-3 py-2 rounded-lg mt-2">
            <div className="flex justify-between items-center">
              <label htmlFor="xi" className="text-sm font-semibold">Extent ξ <span className="font-normal text-slate-500">(moves the red point)</span></label>
              <span className="font-mono text-sm">{currentXi.toFixed(2)} mol</span>
            </div>
            <input
              id="xi"
              type="range"
              min="0.01"
              max="0.99"
              step="0.01"
              value={currentXi}
              onChange={(e) => setCurrentXi(parseFloat(e.target.value))}
              className="w-full h-2 bg-slate-300 rounded-lg appearance-none cursor-pointer accent-red-500"
            />
            <div className="flex justify-between text-xs text-slate-500">
              <span>pure NO₂</span>
              <span>pure N₂O₄</span>
            </div>
          </div>
        </div>

        {/* Right column: conditions and readouts */}
        <div className="sm:col-span-2 space-y-2">
          <div className="border border-slate-200 rounded-lg px-3 py-2">
            <div className="flex justify-between">
              <label htmlFor="T" className="text-sm">Temperature</label>
              <span className="text-sm font-mono font-bold">{temperature} K</span>
            </div>
            <input
              id="T"
              type="range"
              min="250"
              max="400"
              step="1"
              value={temperature}
              onChange={(e) => setTemperature(parseFloat(e.target.value))}
              className="w-full h-2 bg-slate-200 rounded-lg appearance-none cursor-pointer accent-blue-600"
            />
            <div className="flex justify-between mt-2">
              <label htmlFor="P" className="text-sm">Total pressure</label>
              <span className="text-sm font-mono font-bold">{pressure.toFixed(1)} bar</span>
            </div>
            <input
              id="P"
              type="range"
              min="0.1"
              max="10"
              step="0.1"
              value={pressure}
              onChange={(e) => setPressure(parseFloat(e.target.value))}
              className="w-full h-2 bg-slate-200 rounded-lg appearance-none cursor-pointer accent-purple-600"
            />
          </div>

          <div className="bg-slate-50 border border-slate-200 rounded-lg px-3 py-2 text-sm">
            <div className="flex justify-between py-0.5">
              <span>Equilibrium constant K</span>
              <span className="font-mono font-bold text-emerald-700">{fmt(K)}</span>
            </div>
            <div className="flex justify-between py-0.5">
              <span>Reaction quotient Q<sub>p</sub></span>
              <span className="font-mono font-bold">{fmt(Q)}</span>
            </div>
            <div className="flex justify-between py-0.5 border-t border-slate-200 mt-1 pt-1">
              <span>Slope Δ<sub>r</sub>G</span>
              <span className="font-mono font-bold">{delta_r_G.toFixed(2)} kJ/mol</span>
            </div>
            <div className="text-xs text-slate-500">Δ<sub>r</sub>G = RT ln(Q<sub>p</sub>/K), the slope of the red line</div>
            <div className={`text-center font-semibold mt-1 ${slopeColor}`}>{slopeText}</div>

            <div className="relative h-3 bg-slate-200 rounded-full w-full overflow-hidden mt-2">
              <div className="absolute top-0 bottom-0 left-1/2 w-px bg-slate-500" />
              <div
                className="absolute top-0 bottom-0 w-1.5 bg-red-500"
                style={{ left: `clamp(2%, ${50 + Math.log10(Q / K) * 15}%, 97%)` }}
              />
            </div>
            <div className="flex justify-between text-xs text-slate-500 font-mono">
              <span>Q &lt; K</span>
              <span>Q = K</span>
              <span>Q &gt; K</span>
            </div>
          </div>
        </div>
      </div>
    </div>
  );
};

export default EquilibriumVis;
