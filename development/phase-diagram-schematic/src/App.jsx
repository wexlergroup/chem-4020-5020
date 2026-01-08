import React, { useState, useEffect, useRef, useCallback } from 'react';
import * as THREE from 'three';

// --- CONSTANTS ---
const MARGIN = { top: 40, right: 40, bottom: 60, left: 80 };
const MIN_T = 200, MAX_T = 700; // Kelvin
const MIN_P = 1, MAX_P = 100000000; // Pascals (Log scale)
const TRIPLE_POINT = { t: 273.16, p: 611.66 };
const CRITICAL_POINT = { t: 647.1, p: 22064000 };
const STD_PRESSURE = 101325; // 1 atm

// --- HELPER FUNCTIONS ---
const mapT = (t, width) => {
  const plotWidth = width - MARGIN.left - MARGIN.right;
  return MARGIN.left + ((t - MIN_T) / (MAX_T - MIN_T)) * plotWidth;
};

const mapP = (p, height) => {
  const plotHeight = height - MARGIN.top - MARGIN.bottom;
  const logMin = Math.log10(MIN_P);
  const logMax = Math.log10(MAX_P);
  const logP = Math.log10(p);
  return (height - MARGIN.bottom) - ((logP - logMin) / (logMax - logMin) * plotHeight);
};

const invMapT = (x, width) => {
  const plotWidth = width - MARGIN.left - MARGIN.right;
  return ((x - MARGIN.left) / plotWidth) * (MAX_T - MIN_T) + MIN_T;
};

const invMapP = (y, height) => {
  const plotHeight = height - MARGIN.top - MARGIN.bottom;
  const logMin = Math.log10(MIN_P);
  const logMax = Math.log10(MAX_P);
  const normY = ((height - MARGIN.bottom) - y) / plotHeight;
  return Math.pow(10, normY * (logMax - logMin) + logMin);
};

const getPhase = (t, p) => {
  if (p < TRIPLE_POINT.p) {
    if (t > TRIPLE_POINT.t) return 'vapor';
    // Sublimation curve approx
    const boundaryP = 611 * Math.exp(-20 * (1 - t / 273));
    return p < boundaryP ? 'vapor' : 'solid';
  } else {
    if (t < 273) return 'solid';
    if (t > 647 && p > 22064000) return 'supercritical';
    
    // Saturation curve approx
    const pSat = 611 * Math.exp((17.27 * (t - 273.16)) / (t - 35.86));
    if (p > pSat) return 'liquid';
    return 'vapor';
  }
};

// --- COMPONENTS ---

const PhaseDiagram = ({ onPhaseSelect, currentPhase }) => {
  const canvasRef = useRef(null);
  const containerRef = useRef(null);
  const [hoverInfo, setHoverInfo] = useState({ x: 0, y: 0, visible: false, text: '', coords: '--' });

  const draw = useCallback((width, height) => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    const ctx = canvas.getContext('2d');
    
    // Handle high-DPI displays
    const dpr = window.devicePixelRatio || 1;
    canvas.width = width * dpr;
    canvas.height = height * dpr;
    ctx.scale(dpr, dpr);
    canvas.style.width = `${width}px`;
    canvas.style.height = `${height}px`;

    ctx.clearRect(0, 0, width, height);

    // --- GRID ---
    ctx.strokeStyle = '#334155';
    ctx.lineWidth = 0.5;
    ctx.font = '10px sans-serif';
    ctx.fillStyle = '#94a3b8';
    ctx.textAlign = 'center';

    const plotBottom = height - MARGIN.bottom;
    const plotTop = MARGIN.top;
    const plotLeft = MARGIN.left;
    const plotRight = width - MARGIN.right;

    // X-Axis (Temp)
    for (let t = 200; t <= 700; t += 50) {
      let x = mapT(t, width);
      ctx.beginPath();
      ctx.moveTo(x, plotTop); ctx.lineTo(x, plotBottom + 5);
      ctx.stroke();
      ctx.fillText(t, x, plotBottom + 18);
    }
    ctx.fillText("Temperature (K)", (plotLeft + plotRight) / 2, plotBottom + 35);

    // Y-Axis (Pressure)
    ctx.textAlign = 'right';
    ctx.textBaseline = 'middle';
    for (let e = 0; e <= 8; e++) {
      let p = Math.pow(10, e);
      let y = mapP(p, height);
      ctx.beginPath();
      ctx.moveTo(plotLeft - 5, y); ctx.lineTo(plotRight, y);
      ctx.stroke();
      
      let label = `10^${e}`;
      if (e === 0) label = "1 Pa";
      if (e === 5) label = "1 bar";
      if (e === 8) label = "100 MPa";
      ctx.fillText(label, plotLeft - 10, y);
    }

    ctx.save();
    ctx.translate(20, (plotTop + plotBottom) / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.textAlign = 'center';
    ctx.fillText("Pressure (Log Scale)", 0, 0);
    ctx.restore();

    // Axes Box
    ctx.lineWidth = 2;
    ctx.strokeStyle = '#cbd5e1';
    ctx.strokeRect(plotLeft, plotTop, plotRight - plotLeft, plotBottom - plotTop);

    // --- CURVES ---
    ctx.lineWidth = 3;
    const tpX = mapT(TRIPLE_POINT.t, width);
    const tpY = mapP(TRIPLE_POINT.p, height);
    const cpX = mapT(CRITICAL_POINT.t, width);
    const cpY = mapP(CRITICAL_POINT.p, height);

    // Sublimation
    ctx.beginPath();
    ctx.strokeStyle = '#94a3b8';
    ctx.moveTo(mapT(MIN_T, width), mapP(1, height));
    ctx.quadraticCurveTo(mapT(250, width), mapP(100, height), tpX, tpY);
    ctx.stroke();

    // Melting (Negative Slope)
    ctx.beginPath();
    ctx.strokeStyle = '#60a5fa';
    ctx.moveTo(tpX, tpY);
    ctx.lineTo(mapT(250, width), mapP(MAX_P, height));
    ctx.stroke();

    // Vaporization
    ctx.beginPath();
    ctx.strokeStyle = '#f87171';
    ctx.moveTo(tpX, tpY);
    ctx.quadraticCurveTo(mapT(450, width), mapP(100000, height), cpX, cpY);
    ctx.stroke();

    // --- MARKERS ---
    // 1 atm
    ctx.setLineDash([4, 4]);
    ctx.lineWidth = 1;
    ctx.strokeStyle = 'rgba(255,255,255,0.3)';
    const y1atm = mapP(STD_PRESSURE, height);
    ctx.beginPath();
    ctx.moveTo(plotLeft, y1atm);
    ctx.lineTo(plotRight, y1atm);
    ctx.stroke();
    
    ctx.fillStyle = '#fbbf24';
    ctx.textAlign = 'left';
    ctx.fillText("1 atm", plotRight - 35, y1atm - 5);

    ctx.setLineDash([]);
    
    // Points
    const drawPoint = (x, y, color, label, offset = 0) => {
      ctx.fillStyle = color;
      ctx.beginPath();
      ctx.arc(x, y, 4, 0, Math.PI * 2);
      ctx.fill();
      ctx.fillText(label, x + 8, y + offset);
    };

    drawPoint(cpX, cpY, '#f87171', 'Critical Point');
    drawPoint(tpX, tpY, '#60a5fa', 'Triple Point', 15);

    // Annotations
    ctx.font = 'bold 16px sans-serif';
    ctx.fillStyle = 'rgba(255,255,255,0.8)';
    ctx.textAlign = 'center';
    ctx.fillText('SOLID (ICE)', mapT(215, width), mapP(100000, height));
    ctx.fillText('LIQUID', mapT(350, width), mapP(10000000, height));
    ctx.fillText('VAPOR', mapT(500, width), mapP(1000, height));
    
    // Slope Note
    ctx.font = 'italic 11px sans-serif';
    ctx.fillStyle = '#60a5fa';
    ctx.textAlign = 'right';
    ctx.fillText("Negative Slope", mapT(250, width) - 10, mapP(10000000, height));
    ctx.fillText("(Ice less dense)", mapT(250, width) - 10, mapP(10000000, height) + 12);

  }, []);

  useEffect(() => {
    const observer = new ResizeObserver((entries) => {
      window.requestAnimationFrame(() => {
        if (!entries.length) return;
        const { width, height } = entries[0].contentRect;
        draw(width, height);
      });
    });
    if (containerRef.current) observer.observe(containerRef.current);
    return () => observer.disconnect();
  }, [draw]);

  const handleMouseMove = (e) => {
    const rect = canvasRef.current.getBoundingClientRect();
    const x = e.clientX - rect.left;
    const y = e.clientY - rect.top;
    
    // Use container dimensions for mapping, not canvas internal resolution
    const width = rect.width;
    const height = rect.height;

    const t = invMapT(x, width);
    const p = invMapP(y, height);
    
    const phase = getPhase(t, p);
    
    setHoverInfo({
      x: x, 
      y: y, 
      visible: true, 
      text: phase.toUpperCase(), 
      coords: `${Math.round(t)} K / ${p.toExponential(1)} Pa`
    });
  };

  const handleMouseLeave = () => {
    setHoverInfo(prev => ({ ...prev, visible: false, coords: '--' }));
  };

  const handleClick = () => {
    if (hoverInfo.text) {
      onPhaseSelect(hoverInfo.text.toLowerCase());
    }
  };

  return (
    <div className="w-full md:w-1/2 h-1/2 md:h-full border-r border-slate-700 relative bg-slate-900 flex flex-col">
      <div className="p-4 border-b border-slate-700 bg-slate-800">
        <h1 className="text-xl font-bold text-blue-400">H₂O Phase Diagram</h1>
        <p className="text-sm text-slate-400 mt-1">Click a region to view molecular structure.</p>
      </div>
      
      <div ref={containerRef} className="relative flex-grow cursor-pointer">
        <canvas 
          ref={canvasRef} 
          className="w-full h-full block"
          onMouseMove={handleMouseMove}
          onMouseLeave={handleMouseLeave}
          onClick={handleClick}
        />
        {hoverInfo.visible && (
          <div 
            className="absolute bg-white text-slate-900 text-xs px-2 py-1 rounded shadow pointer-events-none transform -translate-x-1/2 -translate-y-full mt-[-10px]"
            style={{ left: hoverInfo.x, top: hoverInfo.y }}
          >
            {hoverInfo.text}
          </div>
        )}
      </div>

      <div className="p-4 bg-slate-800 text-xs text-slate-400 border-t border-slate-700">
        <div className="flex justify-between items-center mb-2">
            <span className="text-slate-500">Current Phase:</span>
            <span className={`font-bold uppercase text-sm ${
                currentPhase === 'solid' ? 'text-blue-300' : 
                currentPhase === 'liquid' ? 'text-blue-500' : 'text-red-400'
            }`}>
                {currentPhase}
            </span>
        </div>
        <div className="flex justify-between border-t border-slate-700 pt-2">
            <span>Temp/Press:</span>
            <span>{hoverInfo.coords}</span>
        </div>
      </div>
    </div>
  );
};

const MolecularSimulation = ({ phase, isPaused, togglePause }) => {
  const mountRef = useRef(null);
  const sceneRef = useRef(null);
  const moleculesRef = useRef([]);
  
  // These refs allow the animation loop to access the latest React state
  // without needing to be recreated on every render
  const phaseRef = useRef(phase);
  const isPausedRef = useRef(isPaused);

  useEffect(() => { phaseRef.current = phase; }, [phase]);
  useEffect(() => { isPausedRef.current = isPaused; }, [isPaused]);

  useEffect(() => {
    // --- THREE.JS SETUP ---
    const scene = new THREE.Scene();
    scene.background = new THREE.Color(0x000000);
    scene.fog = new THREE.FogExp2(0x000000, 0.02);
    
    const camera = new THREE.PerspectiveCamera(50, mountRef.current.clientWidth / mountRef.current.clientHeight, 0.1, 100);
    camera.position.z = 20;

    const renderer = new THREE.WebGLRenderer({ antialias: true });
    renderer.setSize(mountRef.current.clientWidth, mountRef.current.clientHeight);
    mountRef.current.appendChild(renderer.domElement);

    const ambientLight = new THREE.AmbientLight(0x404040, 2);
    scene.add(ambientLight);
    const pointLight = new THREE.PointLight(0xffffff, 1, 100);
    pointLight.position.set(10, 10, 10);
    scene.add(pointLight);
    const pointLight2 = new THREE.PointLight(0x5555ff, 0.5, 100);
    pointLight2.position.set(-10, -10, 5);
    scene.add(pointLight2);

    sceneRef.current = { scene, camera, renderer };

    // Resize Handler
    const handleResize = () => {
      if (!mountRef.current) return;
      camera.aspect = mountRef.current.clientWidth / mountRef.current.clientHeight;
      camera.updateProjectionMatrix();
      renderer.setSize(mountRef.current.clientWidth, mountRef.current.clientHeight);
    };
    window.addEventListener('resize', handleResize);

    // Animation Loop
    let animationId;
    const animate = () => {
      animationId = requestAnimationFrame(animate);
      
      const time = Date.now() * 0.001;
      const currentSimType = phaseRef.current;
      const paused = isPausedRef.current;

      if (!paused) {
        if (currentSimType === 'solid') {
          moleculesRef.current.forEach(m => {
            m.mesh.position.x = m.basePos.x + Math.sin(time * 10 + m.phaseOffset) * 0.05;
            m.mesh.position.y = m.basePos.y + Math.cos(time * 12 + m.phaseOffset) * 0.05;
            m.mesh.position.z = m.basePos.z + Math.sin(time * 8 + m.phaseOffset) * 0.05;
          });
          scene.rotation.y = time * 0.1;
        } else if (currentSimType === 'liquid') {
          moleculesRef.current.forEach(m => {
            m.mesh.position.add(m.velocity);
            m.mesh.rotation.x += m.rotVel.x;
            m.mesh.rotation.y += m.rotVel.y;
            const dist = m.mesh.position.length();
            if (dist > 8) m.velocity.add(m.mesh.position.clone().multiplyScalar(-0.002));
            if (Math.random() > 0.95) {
              m.velocity.add(new THREE.Vector3((Math.random()-0.5)*0.02, (Math.random()-0.5)*0.02, (Math.random()-0.5)*0.02));
            }
          });
          scene.rotation.y = time * 0.05;
        } else if (currentSimType === 'vapor' || currentSimType === 'supercritical') {
          const bounds = 12;
          moleculesRef.current.forEach(m => {
            m.mesh.position.add(m.velocity);
            m.mesh.rotation.x += m.rotVel.x;
            m.mesh.rotation.y += m.rotVel.y;
            if (m.mesh.position.x > bounds || m.mesh.position.x < -bounds) m.velocity.x *= -1;
            if (m.mesh.position.y > bounds || m.mesh.position.y < -bounds) m.velocity.y *= -1;
            if (m.mesh.position.z > bounds || m.mesh.position.z < -bounds) m.velocity.z *= -1;
          });
          scene.rotation.y += 0.001;
        }
      }
      
      renderer.render(scene, camera);
    };
    animate();

    return () => {
      window.removeEventListener('resize', handleResize);
      cancelAnimationFrame(animationId);
      if (mountRef.current && renderer.domElement) {
        mountRef.current.removeChild(renderer.domElement);
      }
      renderer.dispose();
    };
  }, []); // Run once on mount

  // Re-build molecules when phase changes
  useEffect(() => {
    if (!sceneRef.current) return;
    const { scene, camera } = sceneRef.current;

    // Cleanup old molecules
    moleculesRef.current.forEach(m => scene.remove(m.mesh));
    moleculesRef.current = [];

    const createMolecule = () => {
      const group = new THREE.Group();
      const oxy = new THREE.Mesh(
        new THREE.SphereGeometry(0.4, 16, 16),
        new THREE.MeshPhongMaterial({ color: 0xff0000, shininess: 100 })
      );
      group.add(oxy);
      const hydMat = new THREE.MeshPhongMaterial({ color: 0xffffff, shininess: 100 });
      const hydGeo = new THREE.SphereGeometry(0.25, 16, 16);
      const h1 = new THREE.Mesh(hydGeo, hydMat); h1.position.set(0.45, 0.3, 0); group.add(h1);
      const h2 = new THREE.Mesh(hydGeo, hydMat); h2.position.set(-0.45, 0.3, 0); group.add(h2);
      return group;
    };

    if (phase === 'solid') {
      camera.position.z = 12;
      const spacingX = 1.8, spacingY = 2.0, spacingZ = 1.8;
      for(let x = -2; x <= 2; x++) {
        for(let y = -2; y <= 2; y++) {
          for(let z = -1; z <= 1; z++) {
            const mol = createMolecule();
            const offsetX = (y % 2) * 0.5 * spacingX;
            const offsetZ = (y % 2) * 0.5 * spacingZ;
            mol.position.set(x * spacingX + offsetX, y * spacingY, z * spacingZ + offsetZ);
            mol.rotation.x = Math.random() * 0.5;
            mol.rotation.y = Math.random() * Math.PI * 2;
            scene.add(mol);
            moleculesRef.current.push({ mesh: mol, basePos: mol.position.clone(), phaseOffset: Math.random() * 100 });
          }
        }
      }
    } else if (phase === 'liquid') {
      camera.position.z = 15;
      for(let i=0; i<80; i++) {
        const mol = createMolecule();
        const r = 6 * Math.cbrt(Math.random());
        const theta = Math.random() * Math.PI * 2;
        const phi = Math.acos(2 * Math.random() - 1);
        mol.position.set(r * Math.sin(phi) * Math.cos(theta), r * Math.sin(phi) * Math.sin(theta), r * Math.cos(phi));
        mol.rotation.set(Math.random()*Math.PI, Math.random()*Math.PI, 0);
        scene.add(mol);
        moleculesRef.current.push({
          mesh: mol,
          velocity: new THREE.Vector3((Math.random()-0.5)*0.05, (Math.random()-0.5)*0.05, (Math.random()-0.5)*0.05),
          rotVel: new THREE.Vector3((Math.random()-0.5)*0.1, (Math.random()-0.5)*0.1, 0)
        });
      }
    } else { // Vapor or Supercritical
      camera.position.z = 20;
      const count = phase === 'supercritical' ? 100 : 30;
      const speedMult = phase === 'supercritical' ? 0.8 : 0.5;
      const boxSize = phase === 'supercritical' ? 18 : 20;
      
      for(let i=0; i<count; i++) {
        const mol = createMolecule();
        mol.position.set((Math.random()-0.5)*boxSize, (Math.random()-0.5)*boxSize, (Math.random()-0.5)*boxSize);
        scene.add(mol);
        moleculesRef.current.push({
          mesh: mol,
          velocity: new THREE.Vector3((Math.random()-0.5)*speedMult, (Math.random()-0.5)*speedMult, (Math.random()-0.5)*speedMult),
          rotVel: new THREE.Vector3((Math.random()-0.5)*0.2, (Math.random()-0.5)*0.2, (Math.random()-0.5)*0.2)
        });
      }
    }

  }, [phase]);

  const getSimInfo = () => {
    switch(phase) {
      case 'solid': return { title: "Solid Phase (Ice Ih)", desc: "Molecules locked in a lattice. Vibrating but not translating." };
      case 'liquid': return { title: "Liquid Phase", desc: "Molecules close together, high density, flowing freely." };
      case 'vapor': return { title: "Vapor Phase", desc: "Molecules far apart, high speeds, negligible forces." };
      case 'supercritical': return { title: "Supercritical Fluid", desc: "Hybrid state: Gas-like kinetic energy, Liquid-like density." };
      default: return { title: "Waiting...", desc: "Select a phase." };
    }
  };

  const info = getSimInfo();

  return (
    <div className="w-full md:w-1/2 h-1/2 md:h-full relative bg-black">
      <div ref={mountRef} className="w-full h-full" />
      <div id="overlay" className="absolute top-5 left-5 pointer-events-none z-10">
        <div className="bg-slate-900/80 p-4 rounded-lg border border-slate-700 max-w-xs">
          <h2 className="text-lg font-bold text-blue-300 mb-1">{info.title}</h2>
          <p className="text-sm text-gray-300 mb-3">{info.desc}</p>
          <button 
            onClick={togglePause}
            className={`pointer-events-auto px-3 py-1 rounded text-sm font-bold transition-colors ${
              isPaused ? "bg-blue-600 hover:bg-blue-500 text-white" : "bg-slate-600 hover:bg-slate-500 text-white"
            }`}
          >
            {isPaused ? "▶ Play Simulation" : "⏸ Pause Simulation"}
          </button>
        </div>
      </div>
    </div>
  );
};

const App = () => {
  const [phase, setPhase] = useState('liquid');
  const [isPaused, setIsPaused] = useState(true);

  return (
    <div className="flex flex-col h-screen md:flex-row font-sans bg-slate-900 text-white">
      <PhaseDiagram onPhaseSelect={setPhase} currentPhase={phase} />
      <MolecularSimulation phase={phase} isPaused={isPaused} togglePause={() => setIsPaused(!isPaused)} />
    </div>
  );
};

export default App;