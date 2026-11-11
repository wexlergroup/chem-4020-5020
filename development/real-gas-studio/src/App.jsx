import React, { useState, useEffect, useRef, useMemo } from 'react';
import { FlaskConical } from 'lucide-react';

// --- Constants & Data ---
const R = 8.314; // J/(mol·K)

const GAS_PRESETS = [
    { name: "Ideal Gas (Reference)", Tc: 0, Pc: 1, omega: 0, type: "ideal" },
    { name: "Carbon Dioxide (CO2)", Tc: 304.13, Pc: 73.77e5, omega: 0.224, type: "real" },
    { name: "Nitrogen (N2)", Tc: 126.2, Pc: 33.9e5, omega: 0.037, type: "real" },
    { name: "Water Vapor (H2O)", Tc: 647.1, Pc: 220.6e5, omega: 0.344, type: "real" },
    { name: "Methane (CH4)", Tc: 190.56, Pc: 45.99e5, omega: 0.011, type: "real" },
    { name: "Helium (He)", Tc: 5.19, Pc: 2.27e5, omega: -0.385, type: "real" },
];

// --- Math Engine ---

// Cubic Solver Helper: z^3 + a2*z^2 + a1*z + a0 = 0
// Returns the largest real root (Gas phase Z)
function solveCubic(a2, a1, a0) {
    const Q = (3 * a1 - a2 * a2) / 9;
    const R_val = (9 * a2 * a1 - 27 * a0 - 2 * a2 * a2 * a2) / 54;
    const D = Q * Q * Q + R_val * R_val; // Discriminant

    if (D > 0) { // One real root
        const S = Math.cbrt(R_val + Math.sqrt(D));
        const T = Math.cbrt(R_val - Math.sqrt(D));
        return S + T - a2 / 3;
    } else { // Three real roots (take the largest for Vapor phase Z)
        const theta = Math.acos(R_val / Math.sqrt(-Q * Q * Q));
        const z1 = 2 * Math.sqrt(-Q) * Math.cos(theta / 3) - a2 / 3;
        const z2 = 2 * Math.sqrt(-Q) * Math.cos((theta + 2 * Math.PI) / 3) - a2 / 3;
        const z3 = 2 * Math.sqrt(-Q) * Math.cos((theta + 4 * Math.PI) / 3) - a2 / 3;
        return Math.max(z1, z2, z3);
    }
}

// --- EOS Implementations ---

const Models = {
    ideal: {
        getZ: () => 1.0,
        color: '#94a3b8', // Slate 400
        label: 'Ideal Gas'
    },
    vdw: {
        // Van der Waals
        // P = RT/(V-b) - a/V^2
        getZ: (P, T, gas) => {
            const a = (27 * (R ** 2) * (gas.Tc ** 2)) / (64 * gas.Pc);
            const b = (R * gas.Tc) / (8 * gas.Pc);
            
            const A = (a * P) / ((R * T) ** 2);
            const B = (b * P) / (R * T);
            
            const a2 = -(1 + B);
            const a1 = A;
            const a0 = -A * B;
            
            return solveCubic(a2, a1, a0);
        },
        getParams: (gas) => {
            const a = (27 * (R ** 2) * (gas.Tc ** 2)) / (64 * gas.Pc);
            const b = (R * gas.Tc) / (8 * gas.Pc);
            return { a, b };
        },
        color: '#3b82f6', // Blue 500
        label: 'Van der Waals'
    },
    pr: {
        // Peng-Robinson
        getZ: (P, T, gas) => {
            const Tr = T / gas.Tc;
            const kappa = 0.37464 + 1.54226 * gas.omega - 0.26992 * (gas.omega ** 2);
            const alpha = (1 + kappa * (1 - Math.sqrt(Tr))) ** 2;

            const a = (0.45724 * (R ** 2) * (gas.Tc ** 2)) / gas.Pc;
            const b = (0.07780 * R * gas.Tc) / gas.Pc;
            const a_alpha = a * alpha;

            const A = (a_alpha * P) / ((R * T) ** 2);
            const B = (b * P) / (R * T);

            const a2 = -(1 - B);
            const a1 = A - 3 * (B ** 2) - 2 * B;
            const a0 = -(A * B - (B ** 2) - (B ** 3));

            return solveCubic(a2, a1, a0);
        },
        getParams: (gas, T) => {
            const Tr = T / gas.Tc;
            const kappa = 0.37464 + 1.54226 * gas.omega - 0.26992 * (gas.omega ** 2);
            const alpha = (1 + kappa * (1 - Math.sqrt(Tr))) ** 2;
            const a = (0.45724 * (R ** 2) * (gas.Tc ** 2)) / gas.Pc;
            const b = (0.07780 * R * gas.Tc) / gas.Pc;
            return { a: a * alpha, b };
        },
        color: '#10b981', // Emerald 500
        label: 'Peng-Robinson'
    }
};

// --- Custom Hooks ---

// Hook to load scripts dynamically (for Plotly CDN in this standalone component)
const useExternalScript = (url) => {
    const [loaded, setLoaded] = useState(false);

    useEffect(() => {
        if (window.Plotly) {
            setLoaded(true);
            return;
        }
        const script = document.createElement('script');
        script.src = url;
        script.async = true;
        script.onload = () => setLoaded(true);
        document.body.appendChild(script);

        return () => {
            // Optional cleanup
        };
    }, [url]);

    return loaded;
};

// --- Components ---

const ControlPanel = ({ gas, setGas, temp, setTemp, pMax, setPMax, models, toggleModel }) => {
    const currentParamsVDW = gas.type === 'real' ? Models.vdw.getParams(gas) : {a:0, b:0};
    const currentParamsPR = gas.type === 'real' ? Models.pr.getParams(gas, temp) : {a:0, b:0};

    return (
        <div className="bg-slate-900 border-r border-slate-800 w-full md:w-80 p-6 flex flex-col gap-8 overflow-y-auto h-full flex-shrink-0">
            <div>
                <h1 className="text-xl font-bold text-blue-400 mb-2 flex items-center gap-2">
                    <FlaskConical className="w-6 h-6" /> Real Gas Studio
                </h1>
                <p className="text-xs text-slate-400">Computational EOS Comparison</p>
            </div>

            {/* Gas Selection */}
            <div className="space-y-3">
                <label className="text-sm font-semibold text-slate-300 uppercase tracking-wider">Substance</label>
                <select 
                    value={GAS_PRESETS.indexOf(gas)}
                    onChange={(e) => setGas(GAS_PRESETS[e.target.value])}
                    className="w-full bg-slate-950 border border-slate-700 rounded p-2 text-sm focus:border-blue-500 outline-none transition-colors text-slate-200"
                >
                    {GAS_PRESETS.map((g, idx) => (
                        <option key={g.name} value={idx}>{g.name}</option>
                    ))}
                </select>
                {gas.type === 'real' && (
                    <div className="text-xs text-slate-500 grid grid-cols-2 gap-2 mt-2 bg-slate-950 p-2 rounded">
                        <div>Tc: {gas.Tc} K</div>
                        <div>Pc: {(gas.Pc/1e5).toFixed(1)} bar</div>
                        <div>ω: {gas.omega}</div>
                    </div>
                )}
            </div>

            {/* Conditions */}
            <div className="space-y-4">
                <label className="text-sm font-semibold text-slate-300 uppercase tracking-wider">Conditions</label>
                
                <div>
                    <div className="flex justify-between text-xs text-slate-400 mb-1">
                        <span>Temperature (T)</span>
                        <span className="text-blue-400 font-mono">{temp} K</span>
                    </div>
                    <input 
                        type="range" min="100" max="1000" step="5" 
                        value={temp} onChange={(e) => setTemp(Number(e.target.value))}
                        className="w-full h-1 bg-slate-700 rounded-lg appearance-none cursor-pointer slider-thumb"
                    />
                </div>

                <div>
                    <div className="flex justify-between text-xs text-slate-400 mb-1">
                        <span>Max Pressure (P)</span>
                        <span className="text-blue-400 font-mono">{pMax} bar</span>
                    </div>
                    <input 
                        type="range" min="50" max="500" step="10" 
                        value={pMax} onChange={(e) => setPMax(Number(e.target.value))}
                        className="w-full h-1 bg-slate-700 rounded-lg appearance-none cursor-pointer slider-thumb"
                    />
                </div>
            </div>

            {/* Model Toggles */}
            <div className="space-y-3">
                <label className="text-sm font-semibold text-slate-300 uppercase tracking-wider">Models</label>
                <div className="flex flex-col gap-2">
                    {Object.keys(Models).map(key => (
                        <label key={key} className="flex items-center justify-between p-2 rounded hover:bg-slate-800 cursor-pointer transition-colors group">
                            <div className="flex items-center gap-2">
                                <div className={`w-3 h-3 rounded-full`} style={{background: Models[key].color}}></div>
                                <span className="text-sm text-slate-300 group-hover:text-white">{Models[key].label}</span>
                            </div>
                            <input 
                                type="checkbox" 
                                checked={models[key]} 
                                onChange={() => toggleModel(key)}
                                className="accent-blue-500"
                            />
                        </label>
                    ))}
                </div>
            </div>

            {/* Calculated Parameters Display */}
            {gas.type === 'real' && (
                <div className="mt-auto pt-6 border-t border-slate-800">
                     <label className="text-xs font-semibold text-slate-500 uppercase tracking-wider mb-2 block">Computed Parameters</label>
                     <div className="space-y-2 text-xs font-mono">
                        <div className="bg-slate-950 p-2 rounded border-l-2 border-blue-500">
                            <div className="text-slate-400 font-bold mb-1">Van der Waals</div>
                            <div className="flex justify-between"><span>a:</span> <span>{currentParamsVDW.a.toFixed(3)}</span></div>
                            <div className="flex justify-between"><span>b:</span> <span>{(currentParamsVDW.b * 1e5).toFixed(2)} e-5</span></div>
                        </div>
                        <div className="bg-slate-950 p-2 rounded border-l-2 border-emerald-500">
                            <div className="text-slate-400 font-bold mb-1">Peng-Robinson</div>
                            <div className="flex justify-between"><span>a(T):</span> <span>{currentParamsPR.a.toFixed(3)}</span></div>
                            <div className="flex justify-between"><span>b:</span> <span>{(currentParamsPR.b * 1e5).toFixed(2)} e-5</span></div>
                        </div>
                     </div>
                </div>
            )}
        </div>
    );
};

const Charts = ({ gas, temp, pMax, models }) => {
    const zPlotRef = useRef(null);
    const pvPlotRef = useRef(null);
    const plotlyLoaded = useExternalScript("https://cdn.plot.ly/plotly-2.24.1.min.js");

    // Generate Data
    const data = useMemo(() => {
        const steps = 100;
        const pStep = (pMax * 1e5) / steps; 
        
        // Initialize arrays
        const datasets = {
            ideal: { xP: [], yZ: [], xV: [], yP: [] },
            vdw: { xP: [], yZ: [], xV: [], yP: [] },
            pr: { xP: [], yZ: [], xV: [], yP: [] }
        };

        // Avoid P=0 singularity
        for(let i = 1; i <= steps; i++) {
            const P = i * pStep;
            const P_bar = P / 1e5;

            // Ideal
            if(models.ideal) {
                datasets.ideal.xP.push(P_bar);
                datasets.ideal.yZ.push(1);
                datasets.ideal.yP.push(P_bar);
                datasets.ideal.xV.push((R * temp) / P * 1000); // L/mol
            }

            // VdW
            if(models.vdw) {
                if(gas.type === 'ideal') {
                     datasets.vdw.xP.push(P_bar); datasets.vdw.yZ.push(1);
                } else {
                    const z = Models.vdw.getZ(P, temp, gas);
                    const v = (z * R * temp) / P;
                    datasets.vdw.xP.push(P_bar);
                    datasets.vdw.yZ.push(z);
                    datasets.vdw.yP.push(P_bar);
                    datasets.vdw.xV.push(v * 1000); // L/mol
                }
            }

            // PR
            if(models.pr) {
                 if(gas.type === 'ideal') {
                     datasets.pr.xP.push(P_bar); datasets.pr.yZ.push(1);
                } else {
                    const z = Models.pr.getZ(P, temp, gas);
                    const v = (z * R * temp) / P;
                    datasets.pr.xP.push(P_bar);
                    datasets.pr.yZ.push(z);
                    datasets.pr.yP.push(P_bar);
                    datasets.pr.xV.push(v * 1000); // L/mol
                }
            }
        }
        return datasets;
    }, [gas, temp, pMax, models]);

    // Render Plots
    useEffect(() => {
        if (!plotlyLoaded || !window.Plotly) return;

        const commonLayout = {
            paper_bgcolor: 'rgba(0,0,0,0)',
            plot_bgcolor: 'rgba(0,0,0,0)',
            font: { color: '#94a3b8' },
            margin: { t: 40, r: 20, l: 60, b: 50 },
            xaxis: { gridcolor: '#334155', zerolinecolor: '#475569' },
            yaxis: { gridcolor: '#334155', zerolinecolor: '#475569' },
            showlegend: true,
            legend: { x: 1, xanchor: 'right', y: 1 }
        };

        const tracesZ = [];
        const tracesPV = [];

        Object.keys(models).forEach(key => {
            if(!models[key]) return;
            
            tracesZ.push({
                x: data[key].xP,
                y: data[key].yZ,
                type: 'scatter',
                mode: 'lines',
                name: Models[key].label,
                line: { color: Models[key].color, width: 2 }
            });

            tracesPV.push({
                x: data[key].xV,
                y: data[key].yP,
                type: 'scatter',
                mode: 'lines',
                name: Models[key].label,
                line: { color: Models[key].color, width: 2 }
            });
        });

        const layoutZ = {
            ...commonLayout,
            title: 'Compressibility Factor (Z)',
            xaxis: { ...commonLayout.xaxis, title: 'Pressure (bar)' },
            yaxis: { ...commonLayout.yaxis, title: 'Z = PV / RT' }
        };

        const layoutPV = {
            ...commonLayout,
            title: 'P-V Isotherm',
            xaxis: { ...commonLayout.xaxis, title: 'Molar Volume (L/mol)', type: 'log' }, 
            yaxis: { ...commonLayout.yaxis, title: 'Pressure (bar)' }
        };

        window.Plotly.newPlot(zPlotRef.current, tracesZ, layoutZ, {responsive: true, displayModeBar: false});
        window.Plotly.newPlot(pvPlotRef.current, tracesPV, layoutPV, {responsive: true, displayModeBar: false});

    }, [data, models, plotlyLoaded]);

    return (
        <div className="flex-1 p-4 md:p-8 overflow-y-auto grid grid-cols-1 lg:grid-cols-2 gap-6">
            <div className="bg-slate-800/50 rounded-xl border border-slate-700 p-2 shadow-xl backdrop-blur-sm min-h-[400px] flex flex-col">
                <div ref={zPlotRef} className="flex-1 w-full h-full">
                    {!plotlyLoaded && <div className="flex items-center justify-center h-full text-slate-500">Loading Plotly...</div>}
                </div>
                <div className="px-4 pb-2 text-xs text-slate-500 text-center">
                    Shows deviation from Ideal Gas Law (Z=1). Z &lt; 1 implies attractive forces dominate; Z &gt; 1 implies repulsive forces dominate.
                </div>
            </div>
            <div className="bg-slate-800/50 rounded-xl border border-slate-700 p-2 shadow-xl backdrop-blur-sm min-h-[400px] flex flex-col">
                <div ref={pvPlotRef} className="flex-1 w-full h-full">
                    {!plotlyLoaded && <div className="flex items-center justify-center h-full text-slate-500">Loading Plotly...</div>}
                </div>
                <div className="px-4 pb-2 text-xs text-slate-500 text-center">
                    Pressure vs Volume Isotherm (Log scale on V). Visualize the compression path.
                </div>
            </div>

            <div className="lg:col-span-2 bg-slate-900 rounded-lg p-6 border border-slate-800 mb-8">
                <h3 className="text-lg font-semibold text-slate-200 mb-2">Analysis: {gas.name}</h3>
                <p className="text-slate-400 text-sm leading-relaxed">
                    At <strong>{temp} K</strong> and high pressures, real gas behavior deviates significantly from the Ideal Gas Law.
                    { models.vdw && <span> <strong>Van der Waals</strong> accounts for molecular volume (b) and intermolecular forces (a), predicting a Z-factor dip due to attraction. </span> }
                    { models.pr && <span> <strong>Peng-Robinson</strong> improves accuracy near the critical point and for non-spherical molecules (using ω), often showing a sharper deviation in the liquid/dense-gas region. </span> }
                </p>
            </div>
        </div>
    );
};

// --- Main App Component ---

export default function App() {
    const [gas, setGas] = useState(GAS_PRESETS[1]); // Default CO2
    const [temp, setTemp] = useState(310);
    const [pMax, setPMax] = useState(200);
    
    const [models, setModels] = useState({
        ideal: true,
        vdw: true,
        pr: true
    });

    const toggleModel = (key) => {
        setModels(prev => ({...prev, [key]: !prev[key]}));
    };

    return (
        <div className="flex flex-col md:flex-row h-screen w-full bg-slate-950 text-slate-200 font-sans overflow-hidden">
             {/* Slider CSS Injection for Vite compatibility */}
             <style>{`
                .slider-thumb::-webkit-slider-thumb {
                    -webkit-appearance: none; appearance: none;
                    width: 16px; height: 16px;
                    background: #3b82f6; cursor: pointer; border-radius: 50%;
                }
                .slider-thumb::-moz-range-thumb {
                    width: 16px; height: 16px;
                    background: #3b82f6; cursor: pointer; border-radius: 50%;
                }
            `}</style>
            
            <ControlPanel 
                gas={gas} setGas={setGas}
                temp={temp} setTemp={setTemp}
                pMax={pMax} setPMax={setPMax}
                models={models} toggleModel={toggleModel}
            />
            <Charts 
                gas={gas}
                temp={temp}
                pMax={pMax}
                models={models}
            />
        </div>
    );
}