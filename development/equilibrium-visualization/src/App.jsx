import React, { useState, useMemo, useEffect } from 'react';
import { LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, ReferenceLine, ResponsiveContainer, ReferenceDot } from 'recharts';
import { Info, ArrowRight, ArrowLeft, Minus } from 'lucide-react';

const EquilibriumVis = () => {
  // State for user controls
  const [temperature, setTemperature] = useState(298); // Kelvin
  const [pressure, setPressure] = useState(1); // Bar
  const [currentXi, setCurrentXi] = useState(0.5); // Extent of reaction (0 to 1)

  // Constants for NO2 -> N2O4 dimerization
  // 2 NO2 <=> N2O4
  // Enthalpies and Entropies of formation (approximate values)
  const Hf_NO2 = 33.18; // kJ/mol
  const S_NO2 = 0.240; // kJ/(mol K)
  
  const Hf_N2O4 = 9.16; // kJ/mol
  const S_N2O4 = 0.304; // kJ/(mol K)

  const R = 0.008314; // kJ/(mol K)

  // -- Calculation Helper Functions --

  // Calculate standard Chemical Potential (mu_not) at T
  // mu^o = H^o - T*S^o (assuming H and S are const over small T range)
  const calcMuNot = (H, S, T) => H - T * S;

  // Generate the G vs Xi curve data
  const data = useMemo(() => {
    const points = [];
    const muNot_NO2 = calcMuNot(Hf_NO2, S_NO2, temperature);
    const muNot_N2O4 = calcMuNot(Hf_N2O4, S_N2O4, temperature);
    
    // We scan xi from 0.01 to 0.99 to avoid log(0) errors
    for (let xi = 0.01; xi <= 0.99; xi += 0.01) {
      // Moles
      const n_NO2 = 2 - 2 * xi;
      const n_N2O4 = xi;
      const n_total = n_NO2 + n_N2O4;

      // Mole fractions
      const y_NO2 = n_NO2 / n_total;
      const y_N2O4 = n_N2O4 / n_total;

      // Partial Pressures (P_i = y_i * P_total)
      const P_NO2 = y_NO2 * pressure;
      const P_N2O4 = y_N2O4 * pressure;

      // Chemical Potentials (mu_i = mu_i^o + RT ln(P_i/P^o))
      // P^o is 1 bar
      const mu_NO2 = muNot_NO2 + R * temperature * Math.log(P_NO2);
      const mu_N2O4 = muNot_N2O4 + R * temperature * Math.log(P_N2O4);

      // Total Gibbs Energy of the System
      // G_sys = sum(n_i * mu_i)
      const G_sys = n_NO2 * mu_NO2 + n_N2O4 * mu_N2O4;

      points.push({
        xi: parseFloat(xi.toFixed(2)),
        G: G_sys,
      });
    }
    return points;
  }, [temperature, pressure]);

  // Find the theoretical minimum (Equilibrium) from the dataset
  const equilibriumPoint = useMemo(() => {
    return data.reduce((min, p) => p.G < min.G ? p : min, data[0]);
  }, [data]);

  // Calculate current state properties based on user slider 'currentXi'
  const currentState = useMemo(() => {
    const xi = currentXi;
    const n_NO2 = 2 - 2 * xi;
    const n_N2O4 = xi;
    const n_total = n_NO2 + n_N2O4;
    
    // Reaction Quotient Qp
    // Qp = P_N2O4 / (P_NO2)^2
    //    = (y_N2O4 * P) / (y_NO2 * P)^2
    //    = (y_N2O4) / (y_NO2^2 * P)
    const y_NO2 = n_NO2 / n_total;
    const y_N2O4 = n_N2O4 / n_total;
    const Q = y_N2O4 / (Math.pow(y_NO2, 2) * pressure);

    // Equilibrium Constant Kp
    // delta G_r_not = delta H_r - T * delta S_r
    const dH_r = Hf_N2O4 - 2 * Hf_NO2;
    const dS_r = S_N2O4 - 2 * S_NO2;
    const dG_r_not = dH_r - temperature * dS_r;
    const K = Math.exp(-dG_r_not / (R * temperature));

    // Slope (Delta_r G) = RT ln(Q/K)
    // Also = dG/dXi
    const delta_r_G = R * temperature * Math.log(Q / K);
    
    // For the tangent line visualization
    // y = mx + b -> b = y - mx
    // We pick the G value from our pre-calculated data to ensure visual consistency
    const dataPoint = data.find(p => Math.abs(p.xi - xi) < 0.005);
    const G_val = dataPoint ? dataPoint.G : 0;
    
    return {
      Q,
      K,
      delta_r_G,
      G_val,
      dH_r,
      dS_r
    };
  }, [currentXi, temperature, pressure, data]);


  // Determine Slope Direction Text
  let slopeText = "";
  let directionIcon = <Minus size={20} />;
  let slopeColor = "text-gray-500";
  
  if (currentState.delta_r_G < -0.5) {
    slopeText = "Forward Spontaneous";
    directionIcon = <ArrowRight size={20} />;
    slopeColor = "text-green-600";
  } else if (currentState.delta_r_G > 0.5) {
    slopeText = "Reverse Spontaneous";
    directionIcon = <ArrowLeft size={20} />;
    slopeColor = "text-red-500";
  } else {
    slopeText = "Equilibrium";
    directionIcon = <div className="w-5 h-5 rounded-full border-2 border-blue-500 flex items-center justify-center text-xs font-bold text-blue-500">E</div>;
    slopeColor = "text-blue-600";
  }

  // Calculate tangent line points for visualization
  const tangentPoints = [
    { xi: Math.max(0, currentXi - 0.15), G: currentState.G_val - 0.15 * currentState.delta_r_G },
    { xi: Math.min(1, currentXi + 0.15), G: currentState.G_val + 0.15 * currentState.delta_r_G }
  ];

  return (
    <div className="w-full max-w-4xl mx-auto p-4 bg-white rounded-xl shadow-lg font-sans">
      <div className="mb-6 border-b pb-4">
        <h2 className="text-2xl font-bold text-slate-800">Minimizing Gibbs Free Energy</h2>
        <p className="text-slate-500 text-sm mt-1">
          Visualizing equilibrium for <span className="font-mono font-medium text-slate-700">2 NO₂(g) ⇌ N₂O₄(g)</span>
        </p>
      </div>

      <div className="grid grid-cols-1 lg:grid-cols-3 gap-8">
        {/* Left Column: Visualization */}
        <div className="lg:col-span-2 space-y-4">
          
          {/* Main Chart */}
          <div className="h-80 bg-slate-50 rounded-lg border border-slate-200 p-2 relative">
            <ResponsiveContainer width="100%" height="100%">
              <LineChart margin={{ top: 20, right: 30, left: 20, bottom: 20 }}>
                <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0" />
                <XAxis 
                  dataKey="xi" 
                  type="number" 
                  domain={[0, 1]} 
                  label={{ value: 'Extent of Reaction (ξ)', position: 'bottom', offset: 0 }} 
                  allowDataOverflow={true}
                />
                <YAxis 
                  domain={['auto', 'auto']} 
                  label={{ value: 'G (kJ)', angle: -90, position: 'insideLeft' }} 
                  hide={false}
                />
                <Tooltip 
                  formatter={(value) => value.toFixed(2)} 
                  labelFormatter={(label) => `ξ = ${label}`}
                  contentStyle={{ backgroundColor: 'rgba(255,255,255,0.9)', borderRadius: '8px', border: '1px solid #cbd5e1' }}
                />
                
                {/* The Main G Curve */}
                <Line 
                  data={data} 
                  type="monotone" 
                  dataKey="G" 
                  stroke="#3b82f6" 
                  strokeWidth={3} 
                  dot={false} 
                  animationDuration={300}
                />

                {/* Theoretical Equilibrium Point (Minimum) */}
                <ReferenceDot 
                  x={equilibriumPoint.xi} 
                  y={equilibriumPoint.G} 
                  r={6} 
                  fill="transparent" 
                  stroke="#10b981" 
                  strokeWidth={2}
                  strokeDasharray="3 3" 
                />

                {/* Current State Point (User Controlled) */}
                <ReferenceDot 
                  x={currentXi} 
                  y={currentState.G_val} 
                  r={8} 
                  fill="#ef4444" 
                  stroke="#fff" 
                  strokeWidth={2} 
                />

                {/* Tangent Line (Visualizing Slope) */}
                {/* We simulate a tangent by drawing a reference line between two calc'd points */}
                {/* Note: Recharts doesn't support arbitrary line segments easily, so we overlay an SVG or use ReferenceLine if axis matches. 
                    However, strictly plotting a segment in Recharts is tricky. 
                    Simplification: We won't draw the tangent line inside Recharts to avoid artifacting, 
                    instead we rely on the slope readout below. */}
              </LineChart>
            </ResponsiveContainer>

            {/* Overlay Text for Min */}
            <div className="absolute top-4 right-4 bg-white/80 p-2 rounded text-xs text-green-700 border border-green-200">
              Min at ξ ≈ {equilibriumPoint.xi}
            </div>
          </div>

          {/* Slider for Extent (Moving the Red Dot) */}
          <div className="bg-slate-100 p-4 rounded-lg">
            <div className="flex justify-between items-center mb-2">
              <label className="font-semibold text-slate-700 flex items-center gap-2">
                Current Extent (ξ)
                <span className="text-xs font-normal text-slate-500 bg-white px-2 py-0.5 rounded border">
                  Interactive
                </span>
              </label>
              <span className="font-mono text-slate-800">{currentXi.toFixed(2)}</span>
            </div>
            <input
              type="range"
              min="0.01"
              max="0.99"
              step="0.01"
              value={currentXi}
              onChange={(e) => setCurrentXi(parseFloat(e.target.value))}
              className="w-full h-2 bg-slate-300 rounded-lg appearance-none cursor-pointer accent-red-500"
            />
            <div className="flex justify-between text-xs text-slate-500 mt-1">
              <span>Pure Reactants</span>
              <span>Pure Products</span>
            </div>
          </div>

          {/* Slope Visualization Panel */}
          <div className={`flex items-center justify-between p-4 rounded-lg border-l-4 ${
              Math.abs(currentState.delta_r_G) < 0.5 ? 'bg-blue-50 border-blue-500' : 'bg-slate-50 border-slate-300'
            }`}>
            <div>
              <div className="text-xs uppercase tracking-wide text-slate-500 font-bold mb-1">Current Slope</div>
              <div className="font-mono text-lg font-bold text-slate-800">
                ΔᵣG = {currentState.delta_r_G.toFixed(2)} <span className="text-sm font-normal text-slate-500">kJ/mol</span>
              </div>
            </div>
            <div className={`flex items-center gap-2 px-4 py-2 rounded-full bg-white shadow-sm font-bold ${slopeColor}`}>
              {directionIcon}
              {slopeText}
            </div>
          </div>

        </div>

        {/* Right Column: Controls & Data */}
        <div className="space-y-6">
          
          {/* Environmental Controls */}
          <div className="bg-white border border-slate-200 rounded-lg p-5 shadow-sm">
            <h3 className="font-bold text-slate-800 mb-4 flex items-center gap-2">
              <Info size={16} /> Conditions
            </h3>
            
            <div className="space-y-6">
              {/* Temperature Control */}
              <div>
                <div className="flex justify-between mb-1">
                  <label className="text-sm text-slate-600">Temperature</label>
                  <span className="text-sm font-mono font-bold text-slate-800">{temperature} K</span>
                </div>
                <input
                  type="range"
                  min="250"
                  max="400"
                  step="1"
                  value={temperature}
                  onChange={(e) => setTemperature(parseInt(e.target.value))}
                  className="w-full h-2 bg-slate-200 rounded-lg appearance-none cursor-pointer accent-blue-600"
                />
                <p className="text-xs text-slate-400 mt-1">Higher T favors reactants (exothermic).</p>
              </div>

              {/* Pressure Control */}
              <div>
                <div className="flex justify-between mb-1">
                  <label className="text-sm text-slate-600">Pressure</label>
                  <span className="text-sm font-mono font-bold text-slate-800">{pressure} bar</span>
                </div>
                <input
                  type="range"
                  min="0.1"
                  max="10"
                  step="0.1"
                  value={pressure}
                  onChange={(e) => setPressure(parseFloat(e.target.value))}
                  className="w-full h-2 bg-slate-200 rounded-lg appearance-none cursor-pointer accent-purple-600"
                />
                <p className="text-xs text-slate-400 mt-1">Higher P favors products (lower volume).</p>
              </div>
            </div>
          </div>

          {/* Q vs K Comparison */}
          <div className="bg-slate-50 border border-slate-200 rounded-lg p-5">
            <h3 className="font-bold text-slate-800 mb-4">Ratio Analysis</h3>
            
            <div className="space-y-3">
              <div className="flex justify-between items-center pb-2 border-b border-slate-200">
                <span className="text-sm text-slate-600">Reaction Quotient (Q)</span>
                <span className="font-mono font-bold text-slate-800">{currentState.Q.toFixed(2)}</span>
              </div>
              <div className="flex justify-between items-center pb-2 border-b border-slate-200">
                <span className="text-sm text-slate-600">Equilibrium Const (K)</span>
                <span className="font-mono font-bold text-green-700">{currentState.K.toFixed(2)}</span>
              </div>
              
              {/* Visual Scale */}
              <div className="pt-2">
                <div className="text-xs text-center mb-1 text-slate-500">
                  Driving Force: RT ln(Q/K)
                </div>
                <div className="relative h-4 bg-gray-200 rounded-full w-full overflow-hidden flex">
                  <div className="w-1/2 border-r border-white/50 bg-gray-300"></div>
                  <div className="w-1/2 bg-gray-300"></div>
                  
                  {/* Marker */}
                  <div 
                    className="absolute top-0 bottom-0 w-1 bg-blue-600 transition-all duration-300"
                    style={{ 
                      left: `clamp(5%, ${50 + (Math.log10(currentState.Q/currentState.K) * 20)}%, 95%)` 
                    }}
                  />
                </div>
                <div className="flex justify-between text-xs text-slate-400 mt-1 font-mono">
                  <span>Q &lt; K</span>
                  <span>Q = K</span>
                  <span>Q &gt; K</span>
                </div>
              </div>
            </div>
          </div>

        </div>
      </div>
    </div>
  );
};

export default EquilibriumVis;