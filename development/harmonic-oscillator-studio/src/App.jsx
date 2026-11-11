import React, { useState, useMemo } from 'react';
import { LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer, AreaChart, Area } from 'recharts';
import { Settings, Thermometer, Activity, Info, Zap } from 'lucide-react';

const HarmonicOscillatorStudio = () => {
  // State for the simulation
  const [theta, setTheta] = useState(100); // Characteristic Temperature (Theta_E) in Kelvin
  const [maxTemp, setMaxTemp] = useState(500); // Max T for the graph
  const [showLimits, setShowLimits] = useState(true);
  const [activeTab, setActiveTab] = useState('cv'); // 'cv' or 'u'

  // Constants
  // kB is used conceptually but plots are normalized
  
  // Generate data points
  const data = useMemo(() => {
    const points = [];
    const step = maxTemp / 100;
    
    // Avoid T=0 division by zero, start small
    for (let T = 1; T <= maxTemp; T += step) {
      const x = theta / T; // x = hbar*omega / (kB*T)
      
      // Internal Energy U (normalized by kB*Theta)
      // U = kB*Theta * (1/2 + 1/(e^x - 1))
      // Plotting U / (kB*Theta)
      const u_reduced = 0.5 + 1 / (Math.exp(x) - 1);
      
      // High T limit for U: U approx kB*T -> U/(kB*Theta) approx T/Theta
      const u_high_T = T / theta;
      
      // Low T limit for U: U approx 0.5 * kB*Theta -> normalized 0.5
      const u_low_T = 0.5;

      // Heat Capacity Cv (normalized by kB)
      // Cv = kB * x^2 * e^x / (e^x - 1)^2
      let cv_reduced = 0;
      if (x < 100) { // Avoid overflow for very low T
        const numerator = Math.pow(x, 2) * Math.exp(x);
        const denominator = Math.pow(Math.exp(x) - 1, 2);
        cv_reduced = numerator / denominator;
      }
      
      // High T limit for Cv: 1 (Dulong-Petit)
      const cv_high_T = 1;

      // Low T limit for Cv: x^2 * e^-x
      const cv_low_T = Math.pow(x, 2) * Math.exp(-x);

      points.push({
        T: Math.round(T),
        x: parseFloat(x.toFixed(2)),
        u: parseFloat(u_reduced.toFixed(4)),
        u_high: parseFloat(u_high_T.toFixed(4)),
        u_low: parseFloat(u_low_T.toFixed(4)),
        cv: parseFloat(cv_reduced.toFixed(4)),
        cv_high: 1, // Constant 1 kB
        cv_low: parseFloat(cv_low_T.toFixed(4)),
      });
    }
    return points;
  }, [theta, maxTemp]);

  // Calculations for current state display
  const currentX = theta / (maxTemp / 2); // Sample point in middle
  const isQuantum = currentX > 2; // Arbitrary threshold for UI feedback

  return (
    <div className="flex flex-col h-screen bg-slate-50 text-slate-800 font-sans overflow-hidden">
      {/* Header */}
      <header className="bg-slate-900 text-white p-4 shadow-lg flex justify-between items-center z-10">
        <div className="flex items-center gap-2">
          <Activity className="text-cyan-400" />
          <h1 className="text-xl font-bold tracking-tight">Harmonic Oscillator <span className="font-light text-cyan-400">Computational Studio</span></h1>
        </div>
        <div className="text-xs text-slate-400">
          Canonical Ensemble • 1D Oscillator
        </div>
      </header>

      <div className="flex flex-1 overflow-hidden">
        {/* Sidebar Controls */}
        <aside className="w-80 bg-white border-r border-slate-200 p-6 overflow-y-auto flex flex-col gap-8 shadow-sm z-10">
          
          {/* Parameter Section */}
          <section>
            <h2 className="flex items-center gap-2 text-sm font-bold uppercase text-slate-500 mb-4">
              <Settings size={16} /> System Parameters
            </h2>
            
            <div className="mb-6">
              <label className="block text-sm font-medium text-slate-700 mb-2">
                Characteristic Temp (Θ_E)
              </label>
              <div className="flex items-center gap-3">
                <input
                  type="range"
                  min="50"
                  max="1000"
                  step="10"
                  value={theta}
                  onChange={(e) => setTheta(Number(e.target.value))}
                  className="w-full h-2 bg-slate-200 rounded-lg appearance-none cursor-pointer accent-cyan-600"
                />
                <span className="text-cyan-700 font-mono font-bold w-12">{theta}K</span>
              </div>
              <p className="text-xs text-slate-500 mt-2">
                Proportional to frequency (ħω/k_B). Higher Θ_E means "more quantum" at room temp.
              </p>
            </div>

            <div className="mb-6">
              <label className="block text-sm font-medium text-slate-700 mb-2">
                Temperature Range (Max T)
              </label>
              <div className="flex items-center gap-3">
                <input
                  type="range"
                  min="100"
                  max="2000"
                  step="50"
                  value={maxTemp}
                  onChange={(e) => setMaxTemp(Number(e.target.value))}
                  className="w-full h-2 bg-slate-200 rounded-lg appearance-none cursor-pointer accent-slate-600"
                />
                <span className="text-slate-700 font-mono font-bold w-12">{maxTemp}K</span>
              </div>
            </div>

            <div className="flex items-center justify-between p-3 bg-slate-100 rounded-lg">
              <span className="text-sm font-medium text-slate-700">Show Limit Lines</span>
              <button 
                onClick={() => setShowLimits(!showLimits)}
                className={`w-12 h-6 rounded-full transition-colors relative ${showLimits ? 'bg-cyan-600' : 'bg-slate-300'}`}
              >
                <div className={`absolute top-1 w-4 h-4 rounded-full bg-white transition-transform ${showLimits ? 'left-7' : 'left-1'}`} />
              </button>
            </div>
          </section>

          {/* Theory Section */}
          <section className="bg-slate-50 rounded-xl p-4 border border-slate-100">
             <h2 className="flex items-center gap-2 text-sm font-bold uppercase text-slate-500 mb-4">
              <Info size={16} /> Physics Core
            </h2>
            <div className="space-y-4 text-sm text-slate-600">
              <div>
                <strong className="block text-slate-800 mb-1">Partition Function Variable</strong>
                <div className="font-mono text-xs bg-white p-2 rounded border border-slate-200">
                  x = <span className="text-purple-600">Θ_E</span> / T
                </div>
              </div>
              <div>
                <strong className="block text-slate-800 mb-1">Heat Capacity (C_V)</strong>
                <div className="font-mono text-xs bg-white p-2 rounded border border-slate-200">
                   k_B · x² · e^x / (e^x - 1)²
                </div>
              </div>
              <div>
                <strong className="block text-slate-800 mb-1">Energy (U)</strong>
                <div className="font-mono text-xs bg-white p-2 rounded border border-slate-200">
                  k_B·Θ_E · (1/2 + 1/(e^x - 1))
                </div>
              </div>
            </div>
          </section>

          {/* Regime Indicator */}
          <div className={`mt-auto p-4 rounded-lg border ${isQuantum ? 'bg-purple-50 border-purple-200 text-purple-800' : 'bg-orange-50 border-orange-200 text-orange-800'}`}>
            <div className="flex items-center gap-2 font-bold mb-1">
              <Zap size={16} />
              {isQuantum ? "Quantum Regime" : "Classical Regime"}
            </div>
            <p className="text-xs opacity-90">
              At mid-graph (T ≈ {maxTemp/2}K), thermal energy is {isQuantum ? "less than" : "greater than"} the oscillator spacing.
            </p>
          </div>
        </aside>

        {/* Main Content Area */}
        <main className="flex-1 flex flex-col min-w-0 bg-slate-50">
          
          {/* Tabs */}
          <div className="flex border-b border-slate-200 bg-white px-6 pt-4">
            <button
              onClick={() => setActiveTab('cv')}
              className={`pb-4 px-6 font-medium text-sm transition-colors relative ${activeTab === 'cv' ? 'text-cyan-600' : 'text-slate-500 hover:text-slate-700'}`}
            >
              Heat Capacity C_V(T)
              {activeTab === 'cv' && <div className="absolute bottom-0 left-0 right-0 h-0.5 bg-cyan-600" />}
            </button>
            <button
              onClick={() => setActiveTab('u')}
              className={`pb-4 px-6 font-medium text-sm transition-colors relative ${activeTab === 'u' ? 'text-cyan-600' : 'text-slate-500 hover:text-slate-700'}`}
            >
              Internal Energy U(T)
              {activeTab === 'u' && <div className="absolute bottom-0 left-0 right-0 h-0.5 bg-cyan-600" />}
            </button>
          </div>

          {/* Chart Container */}
          <div className="flex-1 p-6 relative">
            <div className="bg-white rounded-2xl shadow-sm border border-slate-200 h-full p-4 flex flex-col">
              
              <div className="mb-4 flex justify-between items-end">
                <div>
                   <h3 className="text-lg font-bold text-slate-800">
                    {activeTab === 'cv' ? 'Heat Capacity Analysis' : 'Internal Energy Analysis'}
                  </h3>
                  <p className="text-sm text-slate-500">
                    {activeTab === 'cv' 
                      ? 'Visualizing the "freezing out" of degrees of freedom at low temperatures.' 
                      : 'Visualizing the transition from zero-point energy to equipartition.'}
                  </p>
                </div>
                {/* Legend/Key */}
                <div className="flex gap-4 text-xs font-medium">
                  <div className="flex items-center gap-1.5">
                    <div className="w-3 h-3 rounded-full bg-cyan-500"></div>
                    <span>Quantum Model</span>
                  </div>
                  {showLimits && (
                    <>
                      <div className="flex items-center gap-1.5">
                        <div className="w-3 h-3 rounded-full bg-orange-400 opacity-50"></div>
                        <span>Classical Limit</span>
                      </div>
                      <div className="flex items-center gap-1.5">
                        <div className="w-3 h-3 rounded-full bg-indigo-400 opacity-50"></div>
                        <span>Low-T Approx</span>
                      </div>
                    </>
                  )}
                </div>
              </div>

              <div className="flex-1 w-full min-h-0">
                <ResponsiveContainer width="100%" height="100%">
                  {activeTab === 'cv' ? (
                    <AreaChart data={data} margin={{ top: 10, right: 30, left: 0, bottom: 0 }}>
                      <defs>
                        <linearGradient id="colorCv" x1="0" y1="0" x2="0" y2="1">
                          <stop offset="5%" stopColor="#0891b2" stopOpacity={0.1}/>
                          <stop offset="95%" stopColor="#0891b2" stopOpacity={0}/>
                        </linearGradient>
                      </defs>
                      <CartesianGrid strokeDasharray="3 3" vertical={false} stroke="#e2e8f0" />
                      <XAxis 
                        dataKey="T" 
                        label={{ value: 'Temperature (K)', position: 'insideBottomRight', offset: -5 }} 
                        stroke="#64748b"
                        tick={{fill: '#64748b', fontSize: 12}}
                      />
                      <YAxis 
                        label={{ value: 'Cv / kB', angle: -90, position: 'insideLeft' }} 
                        stroke="#64748b"
                        tick={{fill: '#64748b', fontSize: 12}}
                        domain={[0, 1.2]}
                      />
                      <Tooltip 
                        contentStyle={{ borderRadius: '8px', border: 'none', boxShadow: '0 4px 6px -1px rgb(0 0 0 / 0.1)' }}
                        itemStyle={{ fontSize: '12px' }}
                        labelStyle={{ color: '#64748b', marginBottom: '0.25rem' }}
                      />
                      
                      {showLimits && (
                        <Line 
                          type="monotone" 
                          dataKey="cv_high" 
                          stroke="#fb923c" 
                          strokeWidth={2} 
                          strokeDasharray="5 5" 
                          dot={false}
                          name="Classical Limit (Dulong-Petit)"
                        />
                      )}

                      {showLimits && (
                        <Line 
                          type="monotone" 
                          dataKey="cv_low" 
                          stroke="#818cf8" 
                          strokeWidth={2} 
                          strokeDasharray="5 5" 
                          dot={false}
                          name="Low-T Approximation"
                        />
                      )}

                      <Area 
                        type="monotone" 
                        dataKey="cv" 
                        stroke="#06b6d4" 
                        strokeWidth={3}
                        fillOpacity={1} 
                        fill="url(#colorCv)" 
                        name="Exact Quantum Model"
                      />
                    </AreaChart>
                  ) : (
                     <LineChart data={data} margin={{ top: 10, right: 30, left: 0, bottom: 0 }}>
                      <CartesianGrid strokeDasharray="3 3" vertical={false} stroke="#e2e8f0" />
                      <XAxis 
                        dataKey="T" 
                        label={{ value: 'Temperature (K)', position: 'insideBottomRight', offset: -5 }} 
                        stroke="#64748b"
                        tick={{fill: '#64748b', fontSize: 12}}
                      />
                      <YAxis 
                        label={{ value: 'U / (kB·Θ)', angle: -90, position: 'insideLeft' }} 
                        stroke="#64748b"
                        tick={{fill: '#64748b', fontSize: 12}}
                      />
                      <Tooltip 
                        contentStyle={{ borderRadius: '8px', border: 'none', boxShadow: '0 4px 6px -1px rgb(0 0 0 / 0.1)' }}
                        itemStyle={{ fontSize: '12px' }}
                        labelStyle={{ color: '#64748b', marginBottom: '0.25rem' }}
                      />
                      
                      {showLimits && (
                        <Line 
                          type="monotone" 
                          dataKey="u_high" 
                          stroke="#fb923c" 
                          strokeWidth={2} 
                          strokeDasharray="5 5" 
                          dot={false}
                          name="Classical Equipartition (kB T)"
                        />
                      )}
                       
                       {showLimits && (
                        <Line 
                          type="monotone" 
                          dataKey="u_low" 
                          stroke="#818cf8" 
                          strokeWidth={2} 
                          strokeDasharray="5 5" 
                          dot={false}
                          name="Zero Point Energy"
                        />
                      )}

                      <Line 
                        type="monotone" 
                        dataKey="u" 
                        stroke="#06b6d4" 
                        strokeWidth={3} 
                        dot={false}
                        name="Exact Quantum Model"
                      />
                    </LineChart>
                  )}
                </ResponsiveContainer>
              </div>

              {/* Annotation Area */}
              <div className="mt-6 grid grid-cols-2 gap-6 p-4 bg-slate-50 rounded-xl border border-slate-100">
                  <div>
                    <h4 className="font-bold text-sm text-slate-800 mb-2 flex items-center gap-2">
                       <Thermometer size={14} className="text-orange-500" /> High-Temperature Limit (T ≫ Θ_E)
                    </h4>
                    <p className="text-xs text-slate-600 leading-relaxed">
                      As T → ∞, x → 0. The exponential term e^x ≈ 1+x. 
                      The heat capacity approaches the constant value k_B (Law of Dulong-Petit), 
                      showing the system behaves classically.
                    </p>
                  </div>
                  <div>
                    <h4 className="font-bold text-sm text-slate-800 mb-2 flex items-center gap-2">
                      <Thermometer size={14} className="text-blue-500" /> Low-Temperature Limit (T ≪ Θ_E)
                    </h4>
                    <p className="text-xs text-slate-600 leading-relaxed">
                       As T → 0, x → ∞. The heat capacity drops to zero exponentially (e^(-Θ_E/T)). 
                       Degrees of freedom "freeze out" because thermal energy k_B T is insufficient 
                       to bridge the gap ħω.
                    </p>
                  </div>
              </div>
            </div>
          </div>
        </main>
      </div>
    </div>
  );
};

export default HarmonicOscillatorStudio;