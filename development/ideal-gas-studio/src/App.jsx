import React, { useEffect, useRef, useState } from 'react';
import { useBook, Studio, Panel, Slider, Select, Sci, V, seeded, STATIC } from '../../studio-kit/kit.jsx';
import { GASES, gas } from './gas.js';

/*
  Ideal gas (Section 2.5). N atoms of a noble gas move in a cube. Panel (a) shows the atoms, seen from one side
  of the cube. Panel (b) lists the quantities that follow from the partition function Q = q_trans^N / N!
  (src/gas.js), in SI units.

  The picture: the two components of each velocity in the plane of the picture come from the Maxwell-Boltzmann
  distribution. A change of the temperature or of the gas rescales all velocities, and a change of the volume
  rescales all positions. One second of the animation shows 5 ps.
*/
const START = { gas: 'he', N: 50, edge: 10, T: 300 };
const POOL = 200;                                     // largest number of atoms
const LARGEST = 15;                                   // nm, largest edge of the cube
const PS_PER_SECOND = 5;                              // time scale of the animation
const SIZE = 372;                                     // px, height of the picture

function makeAtoms() {                                // positions as fractions of the edge, velocities in units of the rms speed
  const random = seeded(4020), atoms = [];
  const normal = () => Math.sqrt(-2 * Math.log(1 - random())) * Math.cos(2 * Math.PI * random());
  for (let i = 0; i < POOL; i++) atoms.push({ x: random(), y: random(), vx: normal(), vy: normal() });
  return atoms;
}

function Atoms({ N, edge, speed, playing, color }) {
  const canvas = useRef(null), atoms = useRef(null), now = useRef({});
  if (!atoms.current) atoms.current = makeAtoms();
  now.current = { N, edge, speed, playing, color };

  useEffect(() => {
    let frame = 0, last = performance.now();
    const draw = (time) => {
      const { N, edge, speed, playing, color } = now.current;
      const element = canvas.current;
      if (!element) return;
      const width = element.clientWidth, ratio = window.devicePixelRatio || 1;
      if (element.width !== Math.round(width * ratio)) { element.width = Math.round(width * ratio); element.height = Math.round(SIZE * ratio); }
      const seconds = Math.min((time - last) / 1000, 0.05);
      last = time;
      // speed is in m/s = 0.001 nm/ps; positions are fractions of the edge
      const step = playing ? speed * 1e-3 * PS_PER_SECOND * seconds / edge : 0;
      const ctx = element.getContext('2d');
      ctx.setTransform(ratio, 0, 0, ratio, 0, 0);
      ctx.clearRect(0, 0, width, SIZE);
      const scale = Math.min(width - 4, SIZE - 34) / LARGEST;           // px per nm
      const side = edge * scale, left = (width - side) / 2, top = (SIZE - 34 - side) / 2 + 2;
      const R = 3;                                                      // px, radius of the dot of an atom
      ctx.fillStyle = color.fill; ctx.fillRect(left, top, side, side);
      ctx.strokeStyle = color.axis; ctx.lineWidth = 1.5; ctx.strokeRect(left, top, side, side);
      ctx.fillStyle = color.blue;
      for (let i = 0; i < N; i++) {
        const a = atoms.current[i];
        a.x += a.vx * step; a.y += a.vy * step;
        if (a.x < 0) { a.x = -a.x; a.vx = -a.vx; } else if (a.x > 1) { a.x = 2 - a.x; a.vx = -a.vx; }
        if (a.y < 0) { a.y = -a.y; a.vy = -a.vy; } else if (a.y > 1) { a.y = 2 - a.y; a.vy = -a.vy; }
        a.x = Math.min(1, Math.max(0, a.x)); a.y = Math.min(1, Math.max(0, a.y));
        ctx.beginPath(); ctx.arc(left + R + a.x * (side - 2 * R), top + R + a.y * (side - 2 * R), R, 0, 2 * Math.PI); ctx.fill();
      }
      // Scale bar of 5 nm under the left corner of the cube
      const y = top + side + 18, x0 = left;
      ctx.strokeStyle = color.text; ctx.lineWidth = 1.5;
      ctx.beginPath(); ctx.moveTo(x0, y); ctx.lineTo(x0 + 5 * scale, y); ctx.moveTo(x0, y - 4); ctx.lineTo(x0, y + 4);
      ctx.moveTo(x0 + 5 * scale, y - 4); ctx.lineTo(x0 + 5 * scale, y + 4); ctx.stroke();
      ctx.fillStyle = color.text; ctx.font = '12px Inter, ui-sans-serif, system-ui, sans-serif'; ctx.textBaseline = 'middle';
      ctx.fillText('5 nm', x0 + 5 * scale + 8, y);
      frame = requestAnimationFrame(draw);
    };
    frame = requestAnimationFrame(draw);
    return () => cancelAnimationFrame(frame);
  }, []);

  return <canvas ref={canvas} style={{ height: SIZE }} role="img"
                 aria-label={`${N} atoms move in straight lines inside a square box and bounce off its walls.`} />;
}

export default function IdealGas() {
  const { color } = useBook();
  const [gasKey, setGasKey] = useState(START.gas);
  const [N, setN] = useState(START.N);
  const [edge, setEdge] = useState(START.edge);       // nm
  const [T, setT] = useState(START.T);                // K
  const [playing, setPlaying] = useState(() => {
    if (STATIC) return false;
    try { return !window.matchMedia('(prefers-reduced-motion: reduce)').matches; } catch (e) { return true; }
  });

  const atom = GASES.find((g) => g.key === gasKey);
  const g = gas(N, edge, T, atom.mass);
  const kB = <V sub="B">k</V>, q = <V sub="trans">q</V>;
  const rows = [
    ['Thermal de Broglie wavelength', <>Λ = (<V>h</V><sup>2</sup>/2π<V>m</V>{kB}<V>T</V>)<sup>1/2</sup></>, <>{(g.lambda * 1e12).toFixed(1)} pm</>],
    ['Mean distance between atoms', <>(<V>V</V>/<V>N</V>)<sup>1/3</sup></>, <>{(g.spacing * 1e9).toFixed(2)} nm</>],
    ['Partition function of one atom', <>{q} = <V>V</V>/Λ<sup>3</sup></>, <Sci value={g.q} />],
    ['Partition function of the gas', <>ln <V>Q</V> = <V>N</V> ln {q} − ln <V>N</V>!</>, g.lnQ.toFixed(1)],
    ['Internal energy', <><V>U</V> = (3/2) <V>N</V>{kB}<V>T</V></>, <><Sci value={g.U} /> J</>],
    ['Pressure', <><V>P</V> = <V>N</V>{kB}<V>T</V>/<V>V</V></>, <>{g.P < 0.1 ? g.P.toFixed(4) : g.P < 10 ? g.P.toFixed(2) : g.P.toFixed(1)} bar</>],
    ['Entropy', <><V>S</V> = <V>U</V>/<V>T</V> + {kB} ln <V>Q</V></>, <><Sci value={g.S} /> J K<sup>−1</sup></>],
  ];
  const reset = () => { setGasKey(START.gas); setN(START.N); setEdge(START.edge); setT(START.T); };

  return (
    <Studio title="Ideal gas"
            lede="Atoms of a noble gas move in a cube. Change the number of atoms, the volume, and the temperature, and follow the quantities that the partition function of the gas gives.">
      <div className="panels two left-narrow">
        <Panel letter="a" title="Atoms in the cube, seen from one side">
          <Atoms N={N} edge={edge} speed={g.speed} playing={playing} color={color} />
        </Panel>
        <Panel letter="b" title="From the partition function">
          <table className="values">
            <tbody>
              {rows.map(([name, formula, value]) => (
                <tr key={name}><td>{name}<span className="formula">{formula}</span></td><td>{value}</td></tr>
              ))}
            </tbody>
          </table>
        </Panel>
      </div>

      <div className="controls">
        <Select label="Gas" value={gasKey} options={GASES.map((x) => [x.key, x.name])} onChange={setGasKey} />
        <Slider label={<>Number of atoms, <V>N</V></>} show={N} value={N} min={1} max={POOL} step={1} onChange={setN} />
        <Slider label={<>Volume, <V>V</V></>} show={<>{(edge ** 3).toFixed(0)} nm<sup>3</sup> (edge {edge.toFixed(1)} nm)</>}
                value={edge} min={5} max={LARGEST} step={0.5} onChange={setEdge} />
        <Slider label={<>Temperature, <V>T</V></>} show={<>{T.toFixed(0)} K</>} value={T} min={50} max={1000} step={10} onChange={setT} />
        <div className="row">
          <button type="button" onClick={() => setPlaying(!playing)}>{playing ? 'Pause' : 'Play'}</button>
          <button type="button" onClick={reset}>Reset</button>
        </div>
      </div>
      <p className="hint">
        The gas is ideal: the atoms are points that do not act on one another. In the picture, one second shows 5 ps of motion.
        The partition function of the gas is <V>Q</V> = {q}<sup><V>N</V></sup>/<V>N</V>!. Because the number of atoms is small, ln <V>N</V>! is
        the exact sum here and not the Stirling approximation. The notes derive <V>U</V> and <V>P</V> from ln <V>Q</V> in Section 2.5
        and the entropy in Section 4.3.
      </p>
    </Studio>
  );
}
