import React, { useState, useEffect, useMemo, useRef } from 'react';
import { Play, Pause, RotateCcw, Info, Atom, Thermometer, Ruler, Scale } from 'lucide-react';
import { BarChart, Bar, XAxis, YAxis, CartesianGrid, Tooltip as RechartsTooltip, ResponsiveContainer, Cell, LineChart, Line } from 'recharts';

/**
 * ROTATIONAL STUDIO
 * A React application to simulate rigid rotor mechanics and thermodynamics.
 */

// --- Constants ---
const CONSTANTS = {
  h: 6.62607015e-34, // Planck constant (J⋅s)
  kb: 1.380649e-23,  // Boltzmann constant (J/K)
  c: 2.99792458e10,  // Speed of light (cm/s) for wavenumber conversion
  amu: 1.66053906660e-27, // Atomic mass unit (kg)
  na: 6.02214076e23  // Avogadro's number
};

// --- Helper Functions ---

// Calculate moment of inertia I = mu * r^2
const calculateInertia = (m1_amu, m2_amu, r_angstrom) => {
  const m1 = m1_amu * CONSTANTS.amu;
  const m2 = m2_amu * CONSTANTS.amu;
  const mu = (m1 * m2) / (m1 + m2);
  const r = r_angstrom * 1e-10; // convert to meters
  return mu * r * r;
};

// Calculate Rotational Constant B in Joules
// B = h^2 / (8 * pi^2 * I)
const calculateB_Joules = (I) => {
  return (CONSTANTS.h * CONSTANTS.h) / (8 * Math.PI * Math.PI * I);
};

// Convert B (Joules) to wavenumbers (cm^-1)
// B_cm = B_J / (h * c * 100)  -- c in m/s here usually, but let's be careful with units
// easier: B_J / (h * c_cm_s)
const joulesToWavenumbers = (E_joules) => {
  return E_joules / (CONSTANTS.h * CONSTANTS.c);
};

// Calculate Rotational Temperature Theta_rot = B / kb
const calculateThetaRot = (B_joules) => {
  return B_joules / CONSTANTS.kb;
};

const RotorSimulation = () => {
  // --- State ---
  const [mass1, setMass1] = useState(12); // amu (e.g., Carbon)
  const [mass2, setMass2] = useState(16); // amu (e.g., Oxygen)
  const [bondLength, setBondLength] = useState(1.13); // Angstroms (CO is ~1.128)
  const [temperature, setTemperature] = useState(300); // Kelvin
  const [isHomonuclear, setIsHomonuclear] = useState(false);
  const [isPlaying, setIsPlaying] = useState(true);
  const [showExplanation, setShowExplanation] = useState(false);

  // Animation Ref
  const canvasRef = useRef(null);
  const animationRef = useRef(null);
  const rotationAngle = useRef(0);

  // --- Physics Calculations ---
  
  // Enforce homonuclear constraint
  useEffect(() => {
    if (isHomonuclear) {
      setMass2(mass1);
    }
  }, [isHomonuclear, mass1]);

  const physicsData = useMemo(() => {
    // 1. Basic Molecular Properties
    const I = calculateInertia(mass1, mass2, bondLength);
    const B_Joules = calculateB_Joules(I);
    const B_cm = joulesToWavenumbers(B_Joules);
    const thetaRot = calculateThetaRot(B_Joules);
    const sigma = isHomonuclear ? 2 : 1;

    // 2. Statistical Mechanics (Partition Function & Populations)
    // We'll calculate up to a J where population is negligible
    const maxJ = 50; 
    let Z_sum = 0;
    const levels = [];

    // Nuclear Spin Statistics Approximation:
    // For this educational model, we will visualize the 'sigma' effect 
    // by modifying the partition function divisor. 
    // Note: In reality, homonuclear diatomics have spin isomers (ortho/para).
    // Here we show the thermodynamic symmetry number effect where Z -> Z/sigma.
    
    for (let J = 0; J <= maxJ; J++) {
      const degeneracy = 2 * J + 1;
      const energy_J = B_Joules * J * (J + 1);
      const boltzmannFactor = Math.exp(-energy_J / (CONSTANTS.kb * temperature));
      
      // For Homonuclear, half the states are forbidden (roughly). 
      // We simulate this by weighing the sum by 1/sigma effectively in the total Z,
      // but for individual populations in the plot, we'll keep the standard Boltzmann curve
      // and show the dropped Z value.
      
      Z_sum += degeneracy * boltzmannFactor;
      
      levels.push({
        J,
        energy_J,
        energy_cm: joulesToWavenumbers(energy_J),
        degeneracy,
        boltzmannFactor,
        weightedFactor: degeneracy * boltzmannFactor
      });
    }

    // Apply Symmetry Number to Z
    const Z_rot = Z_sum / sigma;

    // Calculate normalized populations
    const dataPoints = levels.map(level => ({
      ...level,
      population: (level.weightedFactor / (Z_rot * sigma)) // P_J = (g_J * exp) / (sigma * Z_rot_corrected)? 
      // Actually P(J) is probability of finding molecule in state J. 
      // Sum of P(J) must be 1.
      // If sigma=2, we effectively say "half the rotational states don't exist".
      // So the sum of allowed states would be Z_sum/2.
      // So P(J) for an *allowed* state is (2J+1)exp / (Z_sum/2).
      // For this vis, we will plot the envelope scaled to sum to 1.
    }));

    // Entropy (Rotational) S_rot = Nk [ ln(Z) + T(dlnZ/dT) ] -> High T approx: Nk [ ln(T/sigma*theta) + 1 ]
    // Using the exact sum for Z: S = k (ln Z + <E>/kT)
    // Average Energy <E>
    let totalEnergyWeighted = 0;
    levels.forEach(l => {
      totalEnergyWeighted += l.energy_J * l.weightedFactor;
    });
    const avgEnergy = totalEnergyWeighted / Z_sum;
    const S_rot_joules = CONSTANTS.kb * (Math.log(Z_rot) + avgEnergy / (CONSTANTS.kb * temperature));
    
    // Convert S to J/(mol K)
    const S_molar = S_rot_joules * CONSTANTS.na;

    return {
      I,
      B_cm,
      thetaRot,
      sigma,
      Z_rot,
      S_molar,
      dataPoints: dataPoints.slice(0, 30), // Only show first 30 for chart clarity
      maxPop: Math.max(...dataPoints.map(d => d.population))
    };
  }, [mass1, mass2, bondLength, temperature, isHomonuclear]);


  // --- Animation Loop ---
  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    const ctx = canvas.getContext('2d');
    
    const render = () => {
      if (!ctx) return;
      
      // Clear canvas
      ctx.clearRect(0, 0, canvas.width, canvas.height);
      const centerX = canvas.width / 2;
      const centerY = canvas.height / 2;

      // Physics-based rotation speed (heuristic visualization)
      // Classical frequency ~ sqrt(T / I)
      const baseSpeed = 0.05;
      const speedFactor = Math.sqrt(temperature) * 1e-22 / Math.sqrt(physicsData.I); 
      // The 1e-22 is a scaling factor to make the visual speed reasonable on screen
      
      if (isPlaying) {
        rotationAngle.current += baseSpeed + (speedFactor * 0.5); 
      }

      // Draw Molecule
      ctx.save();
      ctx.translate(centerX, centerY);
      ctx.rotate(rotationAngle.current);

      // Draw Bond
      const visualBondLength = bondLength * 60; // Scale for pixels
      ctx.beginPath();
      ctx.moveTo(-visualBondLength / 2, 0);
      ctx.lineTo(visualBondLength / 2, 0);
      ctx.lineWidth = 6;
      ctx.strokeStyle = '#94a3b8'; // Slate 400
      ctx.stroke();

      // Draw Atom 1
      const r1 = 15 + Math.sqrt(mass1) * 2; // Size based on mass
      ctx.beginPath();
      ctx.arc(-visualBondLength / 2, 0, r1, 0, 2 * Math.PI);
      ctx.fillStyle = isHomonuclear ? '#3b82f6' : '#ef4444'; // Blue if homo, Red if hetero
      ctx.fill();
      ctx.strokeStyle = 'white';
      ctx.lineWidth = 2;
      ctx.stroke();
      
      // Text Atom 1
      ctx.fillStyle = 'white';
      ctx.font = '12px Arial';
      ctx.textAlign = 'center';
      ctx.textBaseline = 'middle';
      ctx.save();
      ctx.translate(-visualBondLength / 2, 0);
      ctx.rotate(-rotationAngle.current); // Keep text upright-ish
      ctx.fillText(mass1, 0, 0);
      ctx.restore();

      // Draw Atom 2
      const r2 = 15 + Math.sqrt(mass2) * 2;
      ctx.beginPath();
      ctx.arc(visualBondLength / 2, 0, r2, 0, 2 * Math.PI);
      ctx.fillStyle = '#3b82f6'; // Blue
      ctx.fill();
      ctx.stroke();

      // Text Atom 2
      ctx.fillStyle = 'white';
      ctx.save();
      ctx.translate(visualBondLength / 2, 0);
      ctx.rotate(-rotationAngle.current);
      ctx.fillText(mass2, 0, 0);
      ctx.restore();

      ctx.restore();
      
      animationRef.current = requestAnimationFrame(render);
    };

    render();

    return () => cancelAnimationFrame(animationRef.current);
  }, [physicsData, bondLength, temperature, isPlaying, mass1, mass2, isHomonuclear]);


  return (
    <div className="flex flex-col h-screen bg-slate-50 text-slate-800 font-sans overflow-hidden">
      {/* HEADER */}
      <header className="flex items-center justify-between px-6 py-4 bg-white border-b border-slate-200 shadow-sm z-10">
        <div className="flex items-center gap-3">
          <div className="bg-indigo-600 p-2 rounded-lg text-white">
            <Atom size={24} />
          </div>
          <div>
            <h1 className="text-xl font-bold text-slate-900">Rigid Rotor Computational Studio</h1>
            <p className="text-xs text-slate-500 uppercase tracking-wider">Symmetry & Partition Functions</p>
          </div>
        </div>
        
        <div className="flex items-center gap-4">
          <div className="flex items-center gap-2 bg-indigo-50 px-3 py-1 rounded-full text-indigo-700 text-sm font-medium border border-indigo-100">
             <span className="text-xs uppercase text-indigo-400">Current Symmetry Number</span>
             <span className="text-lg">σ = {physicsData.sigma}</span>
          </div>
          <button 
            onClick={() => setShowExplanation(!showExplanation)}
            className="p-2 text-slate-400 hover:text-indigo-600 transition-colors"
          >
            <Info size={20} />
          </button>
        </div>
      </header>

      {/* MAIN CONTENT GRID */}
      <main className="flex-1 grid grid-cols-1 lg:grid-cols-12 overflow-hidden">
        
        {/* LEFT PANEL: VISUALIZATION & CONTROLS */}
        <section className="lg:col-span-4 bg-white border-r border-slate-200 flex flex-col overflow-y-auto">
          
          {/* 3D CANVAS */}
          <div className="relative h-64 bg-slate-900 flex items-center justify-center overflow-hidden shrink-0">
            <div className="absolute top-4 left-4 text-white/50 text-xs font-mono">
              ROTATION SIMULATION<br/>
              SPEED ∝ √T/I
            </div>
            <canvas ref={canvasRef} width={400} height={300} className="w-full h-full object-contain" />
            
            <div className="absolute bottom-4 right-4 flex gap-2">
              <button 
                onClick={() => setIsPlaying(!isPlaying)}
                className="p-2 bg-white/10 hover:bg-white/20 text-white rounded-full backdrop-blur-sm transition"
              >
                {isPlaying ? <Pause size={16} /> : <Play size={16} />}
              </button>
            </div>
          </div>

          {/* CONTROLS */}
          <div className="p-6 space-y-8">
            
            {/* Type Toggle */}
            <div className="bg-slate-50 p-4 rounded-xl border border-slate-200">
              <label className="text-xs font-bold text-slate-400 uppercase mb-3 block">Molecule Type</label>
              <div className="flex bg-white rounded-lg p-1 border border-slate-200 shadow-sm">
                <button
                  onClick={() => setIsHomonuclear(false)}
                  className={`flex-1 py-2 text-sm font-medium rounded-md transition-all ${!isHomonuclear ? 'bg-indigo-600 text-white shadow-md' : 'text-slate-500 hover:bg-slate-50'}`}
                >
                  Heteronuclear
                </button>
                <button
                  onClick={() => { setIsHomonuclear(true); setMass2(mass1); }}
                  className={`flex-1 py-2 text-sm font-medium rounded-md transition-all ${isHomonuclear ? 'bg-indigo-600 text-white shadow-md' : 'text-slate-500 hover:bg-slate-50'}`}
                >
                  Homonuclear
                </button>
              </div>
              <p className="text-xs text-slate-500 mt-2 leading-relaxed">
                {isHomonuclear 
                  ? "Symmetry number σ = 2. Rotation by 180° results in an indistinguishable configuration."
                  : "Symmetry number σ = 1. All orientations are unique."}
              </p>
            </div>

            {/* Sliders */}
            <div className="space-y-6">
              <div>
                <div className="flex justify-between mb-2">
                  <label className="flex items-center gap-2 text-sm font-semibold text-slate-700">
                    <Scale size={16} className="text-indigo-500"/> Mass 1 (amu)
                  </label>
                  <span className="text-sm font-mono bg-slate-100 px-2 rounded text-slate-600">{mass1}</span>
                </div>
                <input 
                  type="range" min="1" max="100" step="1" 
                  value={mass1} 
                  onChange={(e) => setMass1(Number(e.target.value))}
                  className="w-full h-2 bg-slate-200 rounded-lg appearance-none cursor-pointer accent-indigo-600"
                />
              </div>

              <div className={isHomonuclear ? "opacity-50 pointer-events-none" : ""}>
                <div className="flex justify-between mb-2">
                  <label className="flex items-center gap-2 text-sm font-semibold text-slate-700">
                    <Scale size={16} className="text-cyan-500"/> Mass 2 (amu)
                  </label>
                  <span className="text-sm font-mono bg-slate-100 px-2 rounded text-slate-600">{mass2}</span>
                </div>
                <input 
                  type="range" min="1" max="100" step="1" 
                  value={mass2} 
                  onChange={(e) => setMass2(Number(e.target.value))}
                  className="w-full h-2 bg-slate-200 rounded-lg appearance-none cursor-pointer accent-cyan-600"
                />
              </div>

              <div>
                <div className="flex justify-between mb-2">
                  <label className="flex items-center gap-2 text-sm font-semibold text-slate-700">
                    <Ruler size={16} className="text-emerald-500"/> Bond Length (Å)
                  </label>
                  <span className="text-sm font-mono bg-slate-100 px-2 rounded text-slate-600">{bondLength}</span>
                </div>
                <input 
                  type="range" min="0.5" max="3.0" step="0.01" 
                  value={bondLength} 
                  onChange={(e) => setBondLength(Number(e.target.value))}
                  className="w-full h-2 bg-slate-200 rounded-lg appearance-none cursor-pointer accent-emerald-600"
                />
              </div>

              <div>
                <div className="flex justify-between mb-2">
                  <label className="flex items-center gap-2 text-sm font-semibold text-slate-700">
                    <Thermometer size={16} className="text-rose-500"/> Temperature (K)
                  </label>
                  <span className="text-sm font-mono bg-slate-100 px-2 rounded text-slate-600">{temperature}</span>
                </div>
                <input 
                  type="range" min="10" max="1000" step="10" 
                  value={temperature} 
                  onChange={(e) => setTemperature(Number(e.target.value))}
                  className="w-full h-2 bg-slate-200 rounded-lg appearance-none cursor-pointer accent-rose-600"
                />
              </div>
            </div>

            {/* Reset */}
            <button 
              onClick={() => {
                setMass1(12);
                setMass2(16);
                setBondLength(1.13);
                setTemperature(300);
                setIsHomonuclear(false);
              }}
              className="flex items-center gap-2 text-sm text-slate-400 hover:text-slate-600 transition-colors mx-auto"
            >
              <RotateCcw size={14} /> Reset Defaults
            </button>

          </div>
        </section>

        {/* RIGHT PANEL: ANALYTICS */}
        <section className="lg:col-span-8 bg-slate-50 flex flex-col h-full overflow-hidden">
          
          {/* DASHBOARD STATS */}
          <div className="grid grid-cols-2 md:grid-cols-4 gap-4 p-6 shrink-0">
            <div className="bg-white p-4 rounded-xl shadow-sm border border-slate-100">
              <div className="text-slate-400 text-xs font-bold uppercase tracking-wider mb-1">Rotational Constant (B)</div>
              <div className="text-2xl font-bold text-slate-800">{physicsData.B_cm.toFixed(2)} <span className="text-sm text-slate-400 font-normal">cm⁻¹</span></div>
            </div>
            <div className="bg-white p-4 rounded-xl shadow-sm border border-slate-100">
              <div className="text-slate-400 text-xs font-bold uppercase tracking-wider mb-1">Char. Temp (Θrot)</div>
              <div className="text-2xl font-bold text-slate-800">{physicsData.thetaRot.toFixed(2)} <span className="text-sm text-slate-400 font-normal">K</span></div>
            </div>
            <div className="bg-white p-4 rounded-xl shadow-sm border border-slate-100 relative overflow-hidden group">
              <div className="absolute right-0 top-0 w-16 h-16 bg-indigo-50 rounded-bl-full -mr-8 -mt-8 transition-transform group-hover:scale-110"></div>
              <div className="text-indigo-500 text-xs font-bold uppercase tracking-wider mb-1 relative z-10">Partition Function (q)</div>
              <div className="text-2xl font-bold text-indigo-700 relative z-10">{physicsData.Z_rot.toFixed(2)}</div>
            </div>
            <div className="bg-white p-4 rounded-xl shadow-sm border border-slate-100">
              <div className="text-slate-400 text-xs font-bold uppercase tracking-wider mb-1">Molar Entropy (S)</div>
              <div className="text-2xl font-bold text-slate-800">{physicsData.S_molar.toFixed(2)} <span className="text-sm text-slate-400 font-normal">J/mol·K</span></div>
            </div>
          </div>

          {/* CHARTS */}
          <div className="flex-1 grid grid-rows-2 gap-6 px-6 pb-6 min-h-0">
            
            {/* Chart 1: Population Distribution */}
            <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-4 flex flex-col min-h-0">
              <div className="flex justify-between items-center mb-4">
                <h3 className="text-sm font-bold text-slate-700">Rotational Population Distribution P(J)</h3>
                <div className="text-xs text-slate-400">Higher T spreads population to higher J</div>
              </div>
              <div className="flex-1 min-h-0">
                <ResponsiveContainer width="100%" height="100%">
                  <BarChart data={physicsData.dataPoints} margin={{ top: 5, right: 30, left: 20, bottom: 5 }}>
                    <CartesianGrid strokeDasharray="3 3" vertical={false} stroke="#e2e8f0" />
                    <XAxis 
                      dataKey="J" 
                      label={{ value: 'Quantum Number J', position: 'insideBottomRight', offset: -5, fontSize: 12 }} 
                      tick={{fontSize: 10}}
                    />
                    <YAxis 
                      hide={false} 
                      label={{ value: 'Population Probability', angle: -90, position: 'insideLeft', fontSize: 12 }} 
                      tick={{fontSize: 10}}
                    />
                    <RechartsTooltip 
                      cursor={{fill: '#f1f5f9'}}
                      contentStyle={{borderRadius: '8px', border: 'none', boxShadow: '0 4px 6px -1px rgb(0 0 0 / 0.1)'}}
                    />
                    <Bar dataKey="population" name="Population" radius={[4, 4, 0, 0]}>
                      {physicsData.dataPoints.map((entry, index) => (
                        <Cell key={`cell-${index}`} fill={entry.population === physicsData.maxPop ? '#6366f1' : '#94a3b8'} />
                      ))}
                    </Bar>
                  </BarChart>
                </ResponsiveContainer>
              </div>
            </div>

            {/* Chart 2: Energy Levels */}
            <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-4 flex flex-col min-h-0">
               <div className="flex justify-between items-center mb-4">
                <h3 className="text-sm font-bold text-slate-700">Energy Levels (E_J)</h3>
                <div className="text-xs text-slate-400">Spacing increases as 2J + 1</div>
              </div>
              <div className="flex-1 min-h-0">
                <ResponsiveContainer width="100%" height="100%">
                  <LineChart data={physicsData.dataPoints} margin={{ top: 5, right: 30, left: 20, bottom: 5 }}>
                    <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0" />
                    <XAxis dataKey="J" label={{ value: 'Quantum Number J', position: 'insideBottomRight', offset: -5, fontSize: 12 }} tick={{fontSize: 10}} />
                    <YAxis label={{ value: 'Energy (cm⁻¹)', angle: -90, position: 'insideLeft', fontSize: 12 }} tick={{fontSize: 10}} />
                    <RechartsTooltip 
                      contentStyle={{borderRadius: '8px', border: 'none', boxShadow: '0 4px 6px -1px rgb(0 0 0 / 0.1)'}}
                    />
                    <Line type="monotone" dataKey="energy_cm" stroke="#0ea5e9" strokeWidth={2} dot={{r: 2}} activeDot={{r: 6}} />
                  </LineChart>
                </ResponsiveContainer>
              </div>
            </div>

          </div>
        </section>
      </main>
      
      {/* EXPLANATION MODAL */}
      {showExplanation && (
        <div className="absolute inset-0 bg-black/50 z-50 flex items-center justify-center p-4">
          <div className="bg-white rounded-2xl max-w-2xl w-full max-h-[90vh] overflow-y-auto shadow-2xl">
            <div className="p-6 border-b border-slate-100 flex justify-between items-center">
              <h2 className="text-xl font-bold text-slate-800">Symmetry & The Rotor</h2>
              <button onClick={() => setShowExplanation(false)} className="p-2 hover:bg-slate-100 rounded-full"><RotateCcw className="rotate-45" size={20}/></button>
            </div>
            <div className="p-6 space-y-4 text-slate-600 leading-relaxed">
              <h3 className="font-bold text-indigo-600">The Symmetry Number (σ)</h3>
              <p>
                The symmetry number, <strong>σ</strong>, corrects for the overcounting of indistinguishable molecular orientations in the partition function.
              </p>
              <ul className="list-disc pl-5 space-y-2">
                <li><strong>Heteronuclear (σ=1):</strong> Example: CO. Rotation by 360° is required to return to the original configuration. All rotational states are unique.</li>
                <li><strong>Homonuclear (σ=2):</strong> Example: N₂. Rotation by 180° brings the molecule to a configuration indistinguishable from the start. We must divide the partition function by 2 to avoid double counting states.</li>
              </ul>
              
              <h3 className="font-bold text-indigo-600 mt-4">Thermodynamic Impact</h3>
              <p>
                Since <strong>q<sub>rot</sub> ∝ 1/σ</strong>, a homonuclear molecule has exactly half the rotational partition function of a hypothetical heteronuclear molecule with identical mass/geometry.
              </p>
              <p>
                This leads to a reduction in Entropy: <strong>S = k ln(q) + ...</strong>. <br/>
                The difference is <strong>ΔS = -R ln(2) ≈ -5.76 J/(mol K)</strong>.
              </p>
              
              <div className="bg-slate-100 p-4 rounded-lg mt-4 text-sm">
                <strong>Note on Nuclear Spin:</strong> In a full quantum mechanical treatment, homonuclear molecules (like H₂) have "ortho" and "para" states based on nuclear spin, restricting J to only even or odd numbers. This simulation uses the high-temperature approximation (σ factor) which is standard for general thermodynamic calculations.
              </div>
            </div>
            <div className="p-4 bg-slate-50 border-t border-slate-100 text-right">
              <button onClick={() => setShowExplanation(false)} className="px-6 py-2 bg-indigo-600 text-white rounded-lg hover:bg-indigo-700">Got it</button>
            </div>
          </div>
        </div>
      )}
    </div>
  );
};

export default RotorSimulation;