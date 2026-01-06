import React, { useState, useMemo } from 'react';
import { LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, Legend, ResponsiveContainer, ReferenceLine } from 'recharts';
import { Settings, Thermometer, Activity, Zap, Calculator } from 'lucide-react';

const TwoLevelSystem = () => {
  // --- State ---
  const [deltaE, setDeltaE] = useState(2.0); // Energy gap
  const [maxTemp, setMaxTemp] = useState(5.0); // Max temperature for simulation
  const [points, setPoints] = useState(100); // Resolution

  // --- Physics Constants ---
  const kB = 1.0; // Boltzmann constant (in arbitrary units for simplicity)

  // --- Calculation Logic ---
  const data = useMemo(() => {
    const dataPoints = [];
    const step = maxTemp / points;

    // Start from a small epsilon to avoid T=0 division by zero
    for (let i = 1; i <= points; i++) {
      const T = i * step;
      const beta = 1 / (kB * T);
      const x = beta * deltaE; // Dimensionless parameter

      // 1. Partition Function Z = 1 + exp(-beta * deltaE)
      // (Assuming ground state E1 = 0, excited state E2 = deltaE)
      const expTerm = Math.exp(-x);
      const Z = 1 + expTerm;

      // 2. Populations
      const P_ground = 1 / Z;
      const P_excited = expTerm / Z;

      // 3. Heat Capacity Cv
      const Cv = kB * (x * x) * expTerm / Math.pow(1 + expTerm, 2);

      dataPoints.push({
        T: parseFloat(T.toFixed(2)),
        Z: parseFloat(Z.toFixed(3)),
        P_ground: parseFloat(P_ground.toFixed(3)),
        P_excited: parseFloat(P_excited.toFixed(3)),
        Cv: parseFloat(Cv.toFixed(4)),
      });
    }
    return dataPoints;
  }, [deltaE, maxTemp, points]);

  // Find Peak Cv (Schottky Anomaly)
  const peakCvData = useMemo(() => {
    if (data.length === 0) return { T: 0, Cv: 0 };
    return data.reduce((max, current) => (current.Cv > max.Cv ? current : max), data[0]);
  }, [data]);

  // Theoretical Peak Location: kB * T_peak / deltaE ≈ 0.417
  const theoreticalPeakT = 0.417 * deltaE / kB;

  return (
    <div className="min-h-screen bg-slate-100 p-4 md:p-8 font-sans">
      <header className="max-w-7xl mx-auto mb-8">
        <div className="flex items-center gap-3 mb-2">
          <Activity className="w-8 h-8 text-indigo-600" />
          <h1 className="text-3xl font-bold text-slate-900">Two-Level System Studio</h1>
        </div>
        <p className="text-slate-600 max-w-2xl">
          An interactive computational physics studio exploring the statistical mechanics of a system with two discrete energy levels.
        </p>
      </header>

      <main className="max-w-7xl mx-auto grid grid-cols-1 lg:grid-cols-12 gap-6">

        {/* Left Column: Controls */}
        <div className="lg:col-span-3">
          <div className="bg-slate-50 p-6 rounded-xl border border-slate-200 shadow-sm h-full">
            <div className="flex items-center gap-2 mb-4 text-slate-800">
              <Settings className="w-5 h-5" />
              <h2 className="font-semibold text-lg">System Parameters</h2>
            </div>

            <div className="space-y-6">
              <div>
                <div className="flex justify-between mb-2">
                  <label className="text-sm font-medium text-slate-600 flex items-center gap-2">
                    <Zap className="w-4 h-4 text-amber-500" />
                    Energy Gap (ΔE)
                  </label>
                  <span className="text-sm font-bold text-slate-900">{deltaE.toFixed(1)} ε</span>
                </div>
                <input
                  type="range"
                  min="0.5"
                  max="10"
                  step="0.1"
                  value={deltaE}
                  onChange={(e) => setDeltaE(parseFloat(e.target.value))}
                  className="w-full h-2 bg-slate-200 rounded-lg appearance-none cursor-pointer accent-indigo-600"
                />
                <p className="text-xs text-slate-500 mt-1">Difference between ground and excited state energies.</p>
              </div>

              <div>
                <div className="flex justify-between mb-2">
                  <label className="text-sm font-medium text-slate-600 flex items-center gap-2">
                    <Thermometer className="w-4 h-4 text-rose-500" />
                    Max Temperature (T)
                  </label>
                  <span className="text-sm font-bold text-slate-900">{maxTemp.toFixed(1)} K</span>
                </div>
                <input
                  type="range"
                  min="1"
                  max="20"
                  step="0.5"
                  value={maxTemp}
                  onChange={(e) => setMaxTemp(parseFloat(e.target.value))}
                  className="w-full h-2 bg-slate-200 rounded-lg appearance-none cursor-pointer accent-rose-600"
                />
                <p className="text-xs text-slate-500 mt-1">Range of the thermal simulation.</p>
              </div>

              <div className="pt-4 border-t border-slate-200">
                 <div className="bg-indigo-50 p-4 rounded-lg border border-indigo-100">
                    <h3 className="text-xs font-bold text-indigo-900 uppercase tracking-wider mb-2 flex items-center gap-2">
                      <Calculator className="w-3 h-3" />
                      Schottky Analysis
                    </h3>
                    <div className="flex justify-between items-center mb-1">
                      <span className="text-xs text-indigo-700">Peak T (Simulated):</span>
                      <span className="text-sm font-mono font-bold text-indigo-900">{peakCvData?.T} K</span>
                    </div>
                    <div className="flex justify-between items-center mb-1">
                       <span className="text-xs text-indigo-700">Peak T (Theory):</span>
                       <span className="text-sm font-mono font-bold text-indigo-900">~{theoreticalPeakT.toFixed(2)} K</span>
                    </div>
                    <div className="flex justify-between items-center">
                      <span className="text-xs text-indigo-700">Max C<sub>v</sub>/k<sub>B</sub>:</span>
                      <span className="text-sm font-mono font-bold text-indigo-900">{peakCvData?.Cv}</span>
                    </div>
                 </div>
              </div>
            </div>
          </div>
        </div>

        {/* Right Column: Graphs */}
        <div className="lg:col-span-9 space-y-6">

          {/* Graph 1: Populations */}
          <div className="bg-white p-6 rounded-xl border border-slate-200 shadow-sm">
            <div className="flex items-center justify-between mb-4">
              <h3 className="font-semibold text-slate-800 flex items-center gap-2">
                <div className="w-2 h-2 rounded-full bg-emerald-500"></div>
                Level Populations (Probability)
              </h3>
              <span className="text-xs bg-slate-100 px-2 py-1 rounded text-slate-500">P vs Temperature</span>
            </div>
            <div className="h-64 w-full">
              <ResponsiveContainer width="100%" height="100%">
                <LineChart data={data} margin={{ top: 10, right: 30, bottom: 20, left: 0 }}>
                  <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0" />
                  <XAxis 
                    dataKey="T" 
                    label={{ value: 'Temperature (kBT)', position: 'insideBottom', offset: -10 }} 
                    stroke="#64748b"
                    fontSize={12}
                    tickMargin={5}
                  />
                  <YAxis 
                    domain={[0, 1.1]}
                    stroke="#64748b"
                    fontSize={12}
                  />
                  <Tooltip 
                    contentStyle={{ borderRadius: '8px', border: 'none', boxShadow: '0 4px 6px -1px rgb(0 0 0 / 0.1)' }}
                  />
                  <Legend 
                    layout="vertical"
                    verticalAlign="top"
                    align="right"
                    wrapperStyle={{
                        top: 10,
                        right: 10,
                        backgroundColor: 'rgba(255, 255, 255, 0.9)',
                        border: '1px solid #e2e8f0',
                        borderRadius: '8px',
                        padding: '10px',
                        boxShadow: '0 2px 4px rgba(0,0,0,0.05)'
                    }}
                  />
                  <Line 
                    type="monotone" 
                    dataKey="P_ground" 
                    name="Ground State (P0)" 
                    stroke="#3b82f6" 
                    strokeWidth={2}
                    dot={false}
                  />
                  <Line 
                    type="monotone" 
                    dataKey="P_excited" 
                    name="Excited State (P1)" 
                    stroke="#ef4444" 
                    strokeWidth={2}
                    dot={false}
                  />
                  <ReferenceLine y={0.5} stroke="#cbd5e1" strokeDasharray="3 3" />
                </LineChart>
              </ResponsiveContainer>
            </div>
            <p className="text-sm text-slate-500 mt-2">
              At T=0, the system is frozen in the ground state (P<sub>0</sub>=1). As T→∞, both populations approach 0.5.
            </p>
          </div>

          <div className="grid grid-cols-1 md:grid-cols-2 gap-6">

            {/* Graph 2: Heat Capacity */}
            <div className="bg-white p-6 rounded-xl border border-slate-200 shadow-sm">
              <div className="flex items-center justify-between mb-4">
                <h3 className="font-semibold text-slate-800 flex items-center gap-2">
                  <div className="w-2 h-2 rounded-full bg-indigo-500"></div>
                  Heat Capacity (C<sub>v</sub>)
                </h3>
              </div>
              <div className="h-48 w-full">
                <ResponsiveContainer width="100%" height="100%">
                  <LineChart data={data} margin={{ top: 5, right: 10, bottom: 5, left: 0 }}>
                    <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0" />
                    <XAxis dataKey="T" stroke="#64748b" fontSize={12} />
                    <YAxis stroke="#64748b" fontSize={12} />
                    <Tooltip contentStyle={{ borderRadius: '8px', border: 'none', boxShadow: '0 4px 6px -1px rgb(0 0 0 / 0.1)' }} />
                    <Line 
                      type="monotone" 
                      dataKey="Cv" 
                      name="Cv / kB" 
                      stroke="#6366f1" 
                      strokeWidth={2}
                      dot={false}
                    />
                    <ReferenceLine x={peakCvData?.T} stroke="#fbbf24" strokeDasharray="3 3" label={{ value: 'Peak', position: 'top', fill:'#d97706', fontSize: 10 }} />
                  </LineChart>
                </ResponsiveContainer>
              </div>
              <div className="mt-2 bg-indigo-50 p-3 rounded text-xs text-indigo-800 border border-indigo-100">
                <strong>Schottky Anomaly:</strong> The peak represents the temperature where the system most effectively absorbs energy to populate the excited state.
              </div>
            </div>

            {/* Graph 3: Partition Function */}
            <div className="bg-white p-6 rounded-xl border border-slate-200 shadow-sm">
              <div className="flex items-center justify-between mb-4">
                <h3 className="font-semibold text-slate-800 flex items-center gap-2">
                  <div className="w-2 h-2 rounded-full bg-slate-500"></div>
                  Partition Function (Z)
                </h3>
              </div>
              <div className="h-48 w-full">
                <ResponsiveContainer width="100%" height="100%">
                  <LineChart data={data} margin={{ top: 5, right: 10, bottom: 5, left: 0 }}>
                    <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0" />
                    <XAxis dataKey="T" stroke="#64748b" fontSize={12} />
                    <YAxis domain={['auto', 'auto']} stroke="#64748b" fontSize={12} />
                    <Tooltip contentStyle={{ borderRadius: '8px', border: 'none', boxShadow: '0 4px 6px -1px rgb(0 0 0 / 0.1)' }} />
                    <Line 
                      type="monotone" 
                      dataKey="Z" 
                      name="Z" 
                      stroke="#64748b" 
                      strokeWidth={2}
                      dot={false}
                    />
                  </LineChart>
                </ResponsiveContainer>
              </div>
               <p className="text-sm text-slate-500 mt-2">
                Z counts the thermally accessible states. <br/>Z → 1 as T → 0. Z → 2 as T → ∞.
              </p>
            </div>
          </div>

        </div>
      </main>
    </div>
  );
};

export default TwoLevelSystem;
