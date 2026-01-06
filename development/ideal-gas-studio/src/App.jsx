import React, { useState, useEffect, useRef, useMemo } from 'react';
import { Play, Pause, RefreshCw, Thermometer, Box, Hash, Info, ChevronRight, ChevronDown } from 'lucide-react';

const IdealGasStudio = () => {
  // --- State ---
  const [isPlaying, setIsPlaying] = useState(true);
  const [particleCount, setParticleCount] = useState(50); // N
  const [volume, setVolume] = useState(200); // V (represented as box width in pixels)
  const [temperature, setTemperature] = useState(20); // T
  const [showMath, setShowMath] = useState(true);

  // Constants for "Simulation Units"
  const kB = 1.0; 
  const mass = 1.0;
  const h = 10.0; // Planck constant equivalent for scaling
  
  // --- Physics Engine Refs ---
  const canvasRef = useRef(null);
  const requestRef = useRef(null);
  const particlesRef = useRef([]);
  const boxSizeRef = useRef(volume);

  // --- Initialization ---
  useEffect(() => {
    // Initialize particles with random positions and velocities based on T
    const newParticles = [];
    for (let i = 0; i < 200; i++) { // Max pool size
      newParticles.push({
        x: Math.random() * volume,
        y: Math.random() * 300, // Fixed height for 2D visual simplicity (V ~ Width)
        vx: (Math.random() - 0.5) * Math.sqrt(temperature),
        vy: (Math.random() - 0.5) * Math.sqrt(temperature),
        color: `hsl(${Math.random() * 60 + 200}, 70%, 50%)`
      });
    }
    particlesRef.current = newParticles;
  }, []); // Run once on mount

  // --- Update Simulation Parameters ---
  useEffect(() => {
    // Update active particle velocities when T changes
    // We scale existing velocities to match new RMS velocity
    const currentParticles = particlesRef.current;
    const currentSpeedFactor = Math.sqrt(temperature);
    
    // Simple heuristic to adjust speed without resetting positions
    currentParticles.forEach(p => {
       const oldSpeed = Math.sqrt(p.vx*p.vx + p.vy*p.vy) || 1;
       const dirX = p.vx / oldSpeed;
       const dirY = p.vy / oldSpeed;
       // Add some randomness so they don't all move in perfect sync
       const speed = currentSpeedFactor * (0.8 + Math.random() * 0.4); 
       p.vx = dirX * speed;
       p.vy = dirY * speed;
    });
    
    boxSizeRef.current = volume;
  }, [temperature, volume]);

  // --- Animation Loop ---
  const animate = () => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    const ctx = canvas.getContext('2d');
    
    // Clear
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    
    // Draw Box
    const boxWidth = boxSizeRef.current;
    const boxHeight = 300; // Fixed height
    
    // Fill background of box
    ctx.fillStyle = '#f8fafc';
    ctx.fillRect(0, 0, boxWidth, boxHeight);
    
    // Draw Walls
    ctx.strokeStyle = '#334155';
    ctx.lineWidth = 4;
    ctx.strokeRect(0, 0, boxWidth, boxHeight);

    // Update & Draw Particles
    const activeParticles = particlesRef.current.slice(0, particleCount);
    
    activeParticles.forEach(p => {
      // Move
      p.x += p.vx;
      p.y += p.vy;

      // Wall Collisions (Bounce)
      if (p.x <= 2) { p.x = 2; p.vx *= -1; }
      if (p.x >= boxWidth - 2) { p.x = boxWidth - 2; p.vx *= -1; }
      if (p.y <= 2) { p.y = 2; p.vy *= -1; }
      if (p.y >= boxHeight - 2) { p.y = boxHeight - 2; p.vy *= -1; }

      // Draw
      ctx.beginPath();
      ctx.arc(p.x, p.y, 4, 0, Math.PI * 2);
      ctx.fillStyle = p.color;
      ctx.fill();
    });

    if (isPlaying) {
      requestRef.current = requestAnimationFrame(animate);
    }
  };

  useEffect(() => {
    requestRef.current = requestAnimationFrame(animate);
    return () => cancelAnimationFrame(requestRef.current);
  }, [isPlaying, particleCount, volume, temperature]); // Re-bind loop when state changes

  // --- Computational Physics (The "Derivation") ---

  // 1. Calculate Thermal De Broglie Wavelength (Lambda)
  // Lambda = h / sqrt(2 * pi * m * kB * T)
  const lambda = useMemo(() => {
    return h / Math.sqrt(2 * Math.PI * mass * kB * temperature);
  }, [temperature]);

  // 2. Partition Function (Z) for N non-interacting particles
  // Z_1 = V / Lambda^2 (Using 2D area/volume for this sim)
  // Z_N = (Z_1)^N / N!
  // We work with ln(Z) to avoid overflow
  const lnZ = useMemo(() => {
    const area = volume * 300; // V
    const Z1 = area / (lambda * lambda); // Single particle partition function (2D)
    
    // Sterling's approximation for ln(N!) approx N ln N - N
    const lnNFactorial = particleCount * Math.log(particleCount) - particleCount;
    
    return (particleCount * Math.log(Z1)) - lnNFactorial;
  }, [volume, lambda, particleCount]);

  // 3. Helmholtz Free Energy (F = -kB T ln Z)
  const F = -kB * temperature * lnZ;

  // 4. Internal Energy (U)
  // Equipartition theorem: U = N * (f/2) * kB * T. For 2D monatomic gas, f=2. So U = N * kB * T.
  // We can also derive it from - d(lnZ)/d(beta), but let's use the explicit result for the display
  // to show the "Exact" value vs the simulation scaling.
  const U_derived = particleCount * kB * temperature; // 2D specific (2/2 * NkT)

  // 5. Pressure (P)
  // P = - (dF / dV)_T = N * kB * T / V (Ideal Gas Law)
  // 2D "Pressure" is Force/Length, effectively P * Area = NkT
  // We will display P_2D = NkT / Area
  const area = volume * 300;
  const P_derived = (particleCount * kB * temperature) / area;

  // 6. Entropy (S)
  // S = - (dF / dT)_V
  // S = (U - F) / T
  const S_derived = (U_derived - F) / temperature;

  // Format helper
  const fmt = (num) => num.toFixed(2);

  return (
    <div className="flex flex-col h-screen bg-slate-50 font-sans text-slate-800 overflow-hidden">
      {/* Header */}
      <div className="bg-white border-b border-slate-200 p-4 shadow-sm flex items-center justify-between">
        <div className="flex items-center space-x-2">
          <div className="bg-indigo-600 p-2 rounded-lg">
            <Hash className="w-5 h-5 text-white" />
          </div>
          <div>
            <h1 className="text-lg font-bold text-slate-900">Ideal Gas Computational Studio</h1>
            <p className="text-xs text-slate-500">Microstate Simulation & Statistical Mechanics Derivation</p>
          </div>
        </div>
        <div className="flex items-center space-x-4 text-sm text-slate-600 bg-slate-100 px-3 py-1 rounded-full">
           <span>k<sub>B</sub> = 1.0</span>
           <span>m = 1.0</span>
        </div>
      </div>

      <div className="flex flex-1 overflow-hidden">
        
        {/* LEFT COLUMN: VISUAL SIMULATION */}
        <div className="flex-1 p-6 flex flex-col space-y-6 overflow-y-auto">
          
          {/* Simulation Container */}
          <div className="bg-white rounded-xl shadow-lg border border-slate-200 p-1 flex flex-col items-center justify-center relative min-h-[400px]">
            <canvas 
              ref={canvasRef} 
              width={800} 
              height={320} 
              className="w-full h-full rounded-lg bg-slate-50 cursor-crosshair"
            />
            
            {/* Overlay Stats */}
            <div className="absolute top-4 right-4 bg-white/90 backdrop-blur border border-slate-200 p-3 rounded-lg shadow-sm text-xs space-y-1">
              <div className="flex justify-between w-32">
                <span className="font-semibold text-indigo-600">Kinetic E:</span>
                <span>{fmt(U_derived)}</span>
              </div>
              <div className="flex justify-between w-32">
                <span className="font-semibold text-emerald-600">Particles:</span>
                <span>{particleCount}</span>
              </div>
            </div>
          </div>

          {/* Controls */}
          <div className="grid grid-cols-1 md:grid-cols-3 gap-6">
            
            {/* N Control */}
            <div className="bg-white p-5 rounded-xl border border-slate-200 shadow-sm transition-all hover:shadow-md">
              <div className="flex items-center space-x-2 mb-4">
                <div className="p-2 bg-emerald-100 rounded-lg">
                   <Hash className="w-5 h-5 text-emerald-600" />
                </div>
                <span className="font-semibold text-slate-700">Particle Count (N)</span>
              </div>
              <input 
                type="range" min="1" max="200" value={particleCount} 
                onChange={(e) => setParticleCount(parseInt(e.target.value))}
                className="w-full h-2 bg-emerald-100 rounded-lg appearance-none cursor-pointer accent-emerald-600"
              />
              <div className="flex justify-between mt-2 text-sm text-slate-500">
                <span>1</span>
                <span className="font-bold text-emerald-700">{particleCount}</span>
                <span>200</span>
              </div>
              <p className="mt-2 text-xs text-slate-400">Scales U and P linearly.</p>
            </div>

            {/* V Control */}
            <div className="bg-white p-5 rounded-xl border border-slate-200 shadow-sm transition-all hover:shadow-md">
              <div className="flex items-center space-x-2 mb-4">
                <div className="p-2 bg-blue-100 rounded-lg">
                   <Box className="w-5 h-5 text-blue-600" />
                </div>
                <span className="font-semibold text-slate-700">Volume (V)</span>
              </div>
              <input 
                type="range" min="100" max="800" value={volume} 
                onChange={(e) => setVolume(parseInt(e.target.value))}
                className="w-full h-2 bg-blue-100 rounded-lg appearance-none cursor-pointer accent-blue-600"
              />
              <div className="flex justify-between mt-2 text-sm text-slate-500">
                <span>Small</span>
                <span className="font-bold text-blue-700">{volume} units</span>
                <span>Large</span>
              </div>
              <p className="mt-2 text-xs text-slate-400">Inversely proportional to P.</p>
            </div>

            {/* T Control */}
            <div className="bg-white p-5 rounded-xl border border-slate-200 shadow-sm transition-all hover:shadow-md">
              <div className="flex items-center space-x-2 mb-4">
                <div className="p-2 bg-rose-100 rounded-lg">
                   <Thermometer className="w-5 h-5 text-rose-600" />
                </div>
                <span className="font-semibold text-slate-700">Temperature (T)</span>
              </div>
              <input 
                type="range" min="1" max="100" value={temperature} 
                onChange={(e) => setTemperature(parseInt(e.target.value))}
                className="w-full h-2 bg-rose-100 rounded-lg appearance-none cursor-pointer accent-rose-600"
              />
              <div className="flex justify-between mt-2 text-sm text-slate-500">
                <span>Cold</span>
                <span className="font-bold text-rose-700">{temperature} K</span>
                <span>Hot</span>
              </div>
              <p className="mt-2 text-xs text-slate-400">Increases Entropy & Pressure.</p>
            </div>

          </div>
        </div>

        {/* RIGHT COLUMN: COMPUTATIONAL DERIVATION */}
        <div className="w-1/3 bg-slate-900 text-slate-100 p-6 flex flex-col space-y-6 shadow-2xl overflow-y-auto border-l border-slate-700">
          
          <div className="flex items-center justify-between pb-4 border-b border-slate-700">
            <h2 className="text-xl font-bold text-indigo-400">Thermodynamic State</h2>
            <button 
              onClick={() => setShowMath(!showMath)}
              className="text-xs bg-slate-800 hover:bg-slate-700 px-3 py-1 rounded-full transition-colors"
            >
              {showMath ? 'Hide Math' : 'Show Math'}
            </button>
          </div>

          {/* Core Partition Function Box */}
          <div className="bg-slate-800/50 p-4 rounded-lg border border-slate-700 space-y-2">
            <div className="flex items-center space-x-2 text-yellow-400 mb-1">
              <Info className="w-4 h-4" />
              <span className="text-sm font-semibold uppercase tracking-wider">The Partition Function (Z)</span>
            </div>
            {showMath && (
              <div className="text-sm font-mono text-slate-400 bg-slate-950 p-3 rounded mb-3">
                 Z = V<sup>N</sup> / (N! &lambda;<sup>2N</sup>)
              </div>
            )}
            <div className="flex justify-between items-end">
              <span className="text-slate-400 text-sm">ln(Z) Value:</span>
              <span className="text-2xl font-mono text-yellow-400">{fmt(lnZ)}</span>
            </div>
          </div>

          {/* Derived Internal Energy */}
          <div className="space-y-3">
            <div className="flex items-center space-x-2">
               <div className="w-2 h-2 rounded-full bg-indigo-500"></div>
               <h3 className="font-bold">Internal Energy (U)</h3>
            </div>
            {showMath && (
               <div className="text-xs text-slate-400 font-mono pl-4 border-l-2 border-indigo-500/30">
                 U = - &part; ln(Z) / &part; &beta; <br/>
                 &nbsp; = N k<sub>B</sub> T
               </div>
            )}
            <div className="bg-indigo-900/20 p-4 rounded-lg border border-indigo-500/30 flex justify-between items-center">
              <span className="text-sm text-indigo-300">Energy Units</span>
              <span className="text-2xl font-bold text-indigo-400">{fmt(U_derived)}</span>
            </div>
          </div>

          {/* Derived Pressure */}
          <div className="space-y-3">
            <div className="flex items-center space-x-2">
               <div className="w-2 h-2 rounded-full bg-blue-500"></div>
               <h3 className="font-bold">Pressure (P)</h3>
            </div>
            {showMath && (
               <div className="text-xs text-slate-400 font-mono pl-4 border-l-2 border-blue-500/30">
                 P = k<sub>B</sub>T ( &part; ln(Z) / &part; V )<br/>
                 &nbsp; = N k<sub>B</sub> T / V
               </div>
            )}
            <div className="bg-blue-900/20 p-4 rounded-lg border border-blue-500/30 flex justify-between items-center">
              <span className="text-sm text-blue-300">Force / Area</span>
              <span className="text-2xl font-bold text-blue-400">{P_derived.toFixed(4)}</span>
            </div>
          </div>

          {/* Derived Entropy */}
          <div className="space-y-3">
            <div className="flex items-center space-x-2">
               <div className="w-2 h-2 rounded-full bg-rose-500"></div>
               <h3 className="font-bold">Entropy (S)</h3>
            </div>
            {showMath && (
               <div className="text-xs text-slate-400 font-mono pl-4 border-l-2 border-rose-500/30">
                 S = k<sub>B</sub> (ln Z + &beta;U)<br/>
                 <span className="opacity-70 italic">(Sackur-Tetrode Eq)</span>
               </div>
            )}
            <div className="bg-rose-900/20 p-4 rounded-lg border border-rose-500/30 flex justify-between items-center">
              <span className="text-sm text-rose-300">J / K</span>
              <span className="text-2xl font-bold text-rose-400">{fmt(S_derived)}</span>
            </div>
          </div>

          {/* Scaling Analysis - Dynamic Feedback */}
          <div className="mt-auto pt-6 border-t border-slate-700">
             <h4 className="text-sm font-semibold text-slate-400 mb-2">Scaling Analysis</h4>
             <div className="grid grid-cols-2 gap-2 text-xs">
                <div className="bg-slate-800 p-2 rounded">
                   <span className="text-slate-500 block">P vs V</span>
                   <span className={volume > 400 ? "text-red-400" : "text-green-400"}>
                     {volume > 400 ? "Low Pressure" : "High Pressure"}
                   </span>
                </div>
                <div className="bg-slate-800 p-2 rounded">
                   <span className="text-slate-500 block">S vs T</span>
                   <span className={temperature > 50 ? "text-rose-400" : "text-blue-400"}>
                     {temperature > 50 ? "High Disorder" : "Low Disorder"}
                   </span>
                </div>
             </div>
          </div>

        </div>
      </div>
    </div>
  );
};

export default IdealGasStudio;