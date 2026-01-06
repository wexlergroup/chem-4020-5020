import React, { useState, useMemo } from 'react';
import {
  LineChart,
  Line,
  XAxis,
  YAxis,
  CartesianGrid,
  Tooltip,
  ResponsiveContainer,
  ScatterChart,
  Scatter,
  ReferenceLine,
  Label
} from 'recharts';
import { Beaker, Calculator, TrendingUp, BookOpen, CheckCircle, AlertCircle, RefreshCw } from 'lucide-react';

const ClausiusClapeyronStudio = () => {
  // Real Vapor Pressure Data for Ethanol (Source: NIST/Wikipedia)
  // T in Celsius, P in Pascals
  const rawData = [
    { t_c: -31.3, p_pa: 133.3 },  // ~1 mmHg
    { t_c: -2.3, p_pa: 1333.2 },  // ~10 mmHg
    { t_c: 19.0, p_pa: 5332.9 },  // ~40 mmHg
    { t_c: 34.9, p_pa: 13332.2 }, // ~100 mmHg
    { t_c: 48.0, p_pa: 26664.0 }, // ~200 mmHg (Interpolated approx)
    { t_c: 63.5, p_pa: 53329.0 }, // ~400 mmHg
    { t_c: 78.4, p_pa: 101325.0 }, // Normal Boiling Point (760 mmHg)
    { t_c: 97.5, p_pa: 202650.0 }, // ~2 atm
  ];

  const R = 8.314; // J/(mol*K)

  // State
  const [step, setStep] = useState(1);
  const [transformedData, setTransformedData] = useState([]);
  const [showBestFit, setShowBestFit] = useState(false);
  const [selectionRange, setSelectionRange] = useState([0, rawData.length - 1]);
  const [studentSlope, setStudentSlope] = useState('');
  const [studentHvap, setStudentHvap] = useState('');
  const [feedback, setFeedback] = useState(null);

  // Constants for Chart
  const processData = () => {
    const processed = rawData.map((point, index) => {
      const t_k = point.t_c + 273.15;
      const inv_t = 1 / t_k;
      const ln_p = Math.log(point.p_pa);
      return {
        id: index,
        ...point,
        t_k: t_k.toFixed(2),
        inv_t: inv_t, // Keep precision for plotting
        inv_t_display: inv_t.toExponential(4),
        ln_p: ln_p, // Keep precision
        ln_p_display: ln_p.toFixed(3)
      };
    });
    setTransformedData(processed);
    setStep(2);
  };

  // Calculate Linear Regression for the selected range
  const regressionStats = useMemo(() => {
    if (transformedData.length === 0) return { slope: 0, intercept: 0, r2: 0 };

    const subset = transformedData.slice(selectionRange[0], selectionRange[1] + 1);
    const n = subset.length;
    if (n < 2) return { slope: 0, intercept: 0, r2: 0 };

    let sumX = 0, sumY = 0, sumXY = 0, sumXX = 0, sumYY = 0;
    subset.forEach(d => {
      sumX += d.inv_t;
      sumY += d.ln_p;
      sumXY += d.inv_t * d.ln_p;
      sumXX += d.inv_t * d.inv_t;
      sumYY += d.ln_p * d.ln_p;
    });

    const slope = (n * sumXY - sumX * sumY) / (n * sumXX - sumX * sumX);
    const intercept = (sumY - slope * sumX) / n;
    
    // Theoretical Hvap from slope
    const calculatedHvap = -slope * R;

    return { slope, intercept, count: n, calculatedHvap };
  }, [transformedData, selectionRange]);

  // Generate line points for the chart based on regression
  const lineData = useMemo(() => {
    if (!transformedData.length) return [];
    const subset = transformedData.slice(selectionRange[0], selectionRange[1] + 1);
    if (subset.length < 2) return [];

    const minX = subset[subset.length - 1].inv_t;
    const maxX = subset[0].inv_t;
    
    return [
      { inv_t: minX, ln_p: regressionStats.slope * minX + regressionStats.intercept },
      { inv_t: maxX, ln_p: regressionStats.slope * maxX + regressionStats.intercept }
    ];
  }, [regressionStats, transformedData, selectionRange]);


  const checkAnswer = () => {
    const userHvap = parseFloat(studentHvap);
    const actualHvap = regressionStats.calculatedHvap / 1000; // Convert to kJ/mol for comparison
    
    if (isNaN(userHvap)) {
      setFeedback({ type: 'error', msg: "Please enter a numeric value." });
      return;
    }

    const error = Math.abs((userHvap - actualHvap) / actualHvap) * 100;

    if (error < 5) {
      setFeedback({ type: 'success', msg: `Excellent! Your value of ${userHvap} kJ/mol is within ${error.toFixed(1)}% of the regression value (${actualHvap.toFixed(2)} kJ/mol).` });
      if (step < 4) setStep(4);
    } else {
      setFeedback({ type: 'warning', msg: `Not quite. Based on the current slope of ${regressionStats.slope.toFixed(0)}, the result should be closer to ${actualHvap.toFixed(2)} kJ/mol. Check your sign or units.` });
    }
  };

  return (
    <div className="max-w-4xl mx-auto p-4 bg-gray-50 min-h-screen font-sans text-gray-800">
      
      {/* Header */}
      <header className="bg-white border-b border-gray-200 p-6 rounded-t-xl shadow-sm mb-6">
        <div className="flex items-center gap-3 mb-2">
          <div className="p-2 bg-blue-100 rounded-lg text-blue-600">
            <Beaker size={24} />
          </div>
          <h1 className="text-2xl font-bold text-gray-900">Clausius-Clapeyron Studio</h1>
        </div>
        <p className="text-gray-600">
          Investigate the relationship between Vapor Pressure and Temperature for <strong>Ethanol</strong> to determine its Enthalpy of Vaporization (<span className="font-serif italic">ΔH<sub>vap</sub></span>).
        </p>
      </header>

      {/* Progress Stepper */}
      <div className="flex gap-2 mb-8 px-2">
        {[1, 2, 3, 4].map((s) => (
          <button
            key={s}
            onClick={() => step >= s ? setStep(s) : null}
            className={`flex-1 py-2 text-sm font-medium rounded-md transition-colors
              ${step === s 
                ? 'bg-blue-600 text-white shadow-md' 
                : step > s 
                  ? 'bg-blue-100 text-blue-800 hover:bg-blue-200' 
                  : 'bg-gray-200 text-gray-400 cursor-not-allowed'}`}
          >
            Phase {s}: {['Data', 'Transform', 'Analysis', 'Interpretation'][s-1]}
          </button>
        ))}
      </div>

      <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
        
        {/* Left Panel: Controls & Context */}
        <div className="lg:col-span-1 space-y-6">
          
          {/* Phase 1: Raw Data View */}
          {step === 1 && (
            <div className="bg-white p-6 rounded-xl shadow-sm border border-gray-200 animate-fade-in">
              <div className="flex items-center gap-2 mb-4 text-blue-700">
                <BookOpen size={20} />
                <h3 className="font-semibold">Raw Experimental Data</h3>
              </div>
              <p className="text-sm text-gray-600 mb-4">
                Here is vapor pressure data for Ethanol collected at various temperatures. To analyze this using linear regression, we need to linearize the Clausius-Clapeyron equation.
              </p>
              <div className="bg-gray-50 p-4 rounded-lg border border-gray-200 mb-6">
                <code className="text-xs text-gray-800 block text-center mb-2">
                  ln(P) = (-ΔH<sub>vap</sub>/R) · (1/T) + C
                </code>
                <p className="text-xs text-center text-gray-500">
                  This resembles y = mx + b
                </p>
              </div>
              <button 
                onClick={processData}
                className="w-full py-3 bg-blue-600 hover:bg-blue-700 text-white rounded-lg font-medium flex items-center justify-center gap-2 transition-all"
              >
                <RefreshCw size={18} />
                Transform Data
              </button>
            </div>
          )}

          {/* Phase 2: Table View */}
          {step >= 2 && (
            <div className="bg-white p-4 rounded-xl shadow-sm border border-gray-200 overflow-hidden">
              <h3 className="font-semibold text-gray-700 mb-3 flex items-center justify-between">
                <span>Data Table</span>
                <span className="text-xs font-normal text-gray-400">Ethanol</span>
              </h3>
              <div className="overflow-x-auto">
                <table className="w-full text-sm text-left">
                  <thead className="bg-gray-50 text-gray-500 font-medium border-b">
                    <tr>
                      <th className="p-2">T (°C)</th>
                      <th className="p-2">1/T (K⁻¹)</th>
                      <th className="p-2">ln P</th>
                    </tr>
                  </thead>
                  <tbody className="divide-y divide-gray-100">
                    {transformedData.map((row, i) => (
                      <tr 
                        key={i} 
                        className={`hover:bg-blue-50 cursor-pointer transition-colors ${
                          (i >= selectionRange[0] && i <= selectionRange[1]) ? 'bg-blue-50/50' : 'opacity-50'
                        }`}
                        onClick={() => {
                          // Simple range selection logic for demo
                          if (i < selectionRange[0]) setSelectionRange([i, selectionRange[1]]);
                          else if (i > selectionRange[1]) setSelectionRange([selectionRange[0], i]);
                          else setSelectionRange([i, i]); // Reset if clicking inside
                        }}
                      >
                        <td className="p-2 text-gray-900 font-medium">{row.t_c}</td>
                        <td className="p-2 font-mono text-xs text-gray-600">{row.inv_t_display}</td>
                        <td className="p-2 font-mono text-xs text-blue-600">{row.ln_p_display}</td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
              <p className="text-xs text-gray-400 mt-2 italic text-center">
                * Rows highlighted are included in regression
              </p>

              {step === 2 && (
                <button 
                  onClick={() => setStep(3)}
                  className="w-full mt-4 py-2 bg-blue-600 hover:bg-blue-700 text-white rounded-lg font-medium transition-colors flex items-center justify-center gap-2"
                >
                  <TrendingUp size={16} />
                  Proceed to Analysis
                </button>
              )}
            </div>
          )}

          {/* Phase 3 & 4: Calculation Tools */}
          {step >= 3 && (
            <div className="bg-white p-6 rounded-xl shadow-sm border border-gray-200 animate-fade-in">
              <div className="flex items-center gap-2 mb-4 text-purple-700">
                <Calculator size={20} />
                <h3 className="font-semibold">Analysis</h3>
              </div>
              
              <div className="space-y-4">
                <div>
                  <label className="block text-xs font-medium text-gray-500 uppercase mb-1">Slope (m)</label>
                  <div className="text-lg font-mono bg-gray-100 p-2 rounded border border-gray-200 text-gray-800">
                    {regressionStats.slope.toFixed(1)} <span className="text-xs text-gray-400">K</span>
                  </div>
                </div>

                <div>
                  <label className="block text-xs font-medium text-gray-500 uppercase mb-1">Calculate ΔH<sub>vap</sub> (kJ/mol)</label>
                  <div className="flex gap-2">
                    <input 
                      type="number" 
                      placeholder="Enter value..."
                      className="w-full p-2 border border-gray-300 rounded focus:ring-2 focus:ring-purple-500 outline-none"
                      value={studentHvap}
                      onChange={(e) => setStudentHvap(e.target.value)}
                    />
                    <button 
                      onClick={checkAnswer}
                      className="px-4 bg-purple-600 text-white rounded hover:bg-purple-700"
                    >
                      Check
                    </button>
                  </div>
                  <p className="text-xs text-gray-400 mt-1">Hint: Slope = -ΔH<sub>vap</sub> / R</p>
                </div>

                {feedback && (
                  <div className={`p-3 rounded-lg text-sm flex items-start gap-2 ${
                    feedback.type === 'success' ? 'bg-green-50 text-green-800 border border-green-200' : 'bg-amber-50 text-amber-800 border border-amber-200'
                  }`}>
                    {feedback.type === 'success' ? <CheckCircle size={16} className="mt-0.5" /> : <AlertCircle size={16} className="mt-0.5" />}
                    {feedback.msg}
                  </div>
                )}
              </div>
            </div>
          )}

        </div>

        {/* Right Panel: Visualization */}
        <div className="lg:col-span-2 space-y-6">
          
          {step >= 2 ? (
            <div className="bg-white p-6 rounded-xl shadow-sm border border-gray-200 h-[500px] flex flex-col">
               <div className="flex justify-between items-center mb-4">
                 <h2 className="text-lg font-semibold flex items-center gap-2">
                   <TrendingUp size={20} className="text-gray-400"/>
                   Clausius-Clapeyron Plot
                 </h2>
                 <div className="flex items-center gap-3">
                   <label className="flex items-center gap-2 text-sm text-gray-600 cursor-pointer select-none">
                     <input 
                       type="checkbox" 
                       checked={showBestFit} 
                       onChange={() => setShowBestFit(!showBestFit)}
                       className="rounded text-blue-600 focus:ring-blue-500"
                     />
                     Show Regression Line
                   </label>
                 </div>
               </div>

               <div className="flex-1 w-full min-h-0">
                <ResponsiveContainer width="100%" height="100%">
                  <ScatterChart margin={{ top: 20, right: 30, bottom: 40, left: 20 }}>
                    <CartesianGrid strokeDasharray="3 3" stroke="#e5e7eb" />
                    <XAxis 
                      type="number" 
                      dataKey="inv_t" 
                      name="1/T" 
                      unit=" K⁻¹" 
                      domain={['auto', 'auto']}
                      tickFormatter={(val) => val.toExponential(2)}
                      reversed={true} // Convention: High Temp (Low 1/T) on right? Actually standard is low 1/T on left (High T). But often chemists plot 1/T increasing to right. Let's stick to standard math axis (increasing right).
                      label={{ value: 'Inverse Temperature (1/T)', position: 'bottom', offset: 25 }}
                    />
                    <YAxis 
                      type="number" 
                      dataKey="ln_p" 
                      name="ln(P)" 
                      domain={['auto', 'auto']}
                      label={{ value: 'ln(Pressure)', angle: -90, position: 'left', offset: 10 }}
                    />
                    <Tooltip 
                      cursor={{ strokeDasharray: '3 3' }}
                      content={({ active, payload }) => {
                        if (active && payload && payload.length) {
                          const data = payload[0].payload;
                          return (
                            <div className="bg-white p-3 border border-gray-200 shadow-lg rounded-lg text-sm">
                              <p className="font-semibold mb-1">T: {data.t_c}°C</p>
                              <p>1/T: {data.inv_t.toExponential(4)} K⁻¹</p>
                              <p>ln(P): {data.ln_p.toFixed(3)}</p>
                              <p className="text-gray-400 text-xs mt-1">P: {(data.p_pa/1000).toFixed(1)} kPa</p>
                            </div>
                          );
                        }
                        return null;
                      }}
                    />
                    
                    {/* The Data Points */}
                    <Scatter 
                      name="Ethanol Data" 
                      data={transformedData} 
                      fill="#3b82f6" 
                      shape="circle"
                    />

                    {/* The Regression Line */}
                    {showBestFit && lineData.length > 0 && (
                      <Scatter 
                        data={lineData} 
                        line={{ stroke: '#9333ea', strokeWidth: 2 }} 
                        shape={() => null} 
                        legendType="none" 
                      />
                    )}
                  </ScatterChart>
                </ResponsiveContainer>
               </div>
               
               {/* Interactive Range Slider (Simplified representation) */}
               {step >= 3 && (
                 <div className="mt-4 px-4 py-3 bg-gray-50 rounded-lg border border-gray-200">
                   <div className="flex justify-between text-xs text-gray-500 mb-2">
                     <span>High Temp (Low 1/T)</span>
                     <span className="font-semibold text-purple-700">Adjust range to see slope changes</span>
                     <span>Low Temp (High 1/T)</span>
                   </div>
                   <input 
                      type="range"
                      min="0"
                      max={transformedData.length - 1}
                      value={selectionRange[1]}
                      onChange={(e) => {
                        const val = parseInt(e.target.value);
                        if(val > selectionRange[0]) setSelectionRange([selectionRange[0], val]);
                      }}
                      className="w-full h-2 bg-gray-200 rounded-lg appearance-none cursor-pointer accent-purple-600"
                   />
                   <div className="flex justify-between text-xs mt-1 font-mono">
                     <span>{transformedData[selectionRange[1]]?.t_c}°C</span>
                     <span>Current Range</span>
                     <span>{transformedData[selectionRange[0]]?.t_c}°C</span>
                   </div>
                 </div>
               )}
            </div>
          ) : (
            <div className="h-[500px] flex items-center justify-center bg-gray-100 rounded-xl border-2 border-dashed border-gray-300 text-gray-400">
              Transform data to view plot
            </div>
          )}

          {/* Phase 4: Interpretation Panel */}
          {step === 4 && (
            <div className="bg-blue-50 p-6 rounded-xl border border-blue-200 animate-fade-in">
              <h3 className="font-bold text-blue-900 mb-3">Interpretation & Linearity Limits</h3>
              <p className="text-sm text-blue-800 mb-4">
                Literature value for Ethanol ΔH<sub>vap</sub> at boiling point is typically quoted around <strong>38.6 kJ/mol</strong>. 
                However, at 25°C, it's closer to <strong>42.3 kJ/mol</strong>.
              </p>
              <div className="space-y-3">
                <div className="flex gap-3 items-start">
                  <div className="bg-white p-2 rounded shadow-sm shrink-0 font-bold text-blue-600">Q1</div>
                  <p className="text-sm text-gray-700">
                    Try adjusting the slider to select only the <strong>low temperature points</strong> (left side of slider, high 1/T). Does the slope become steeper (more negative) or flatter? What does this imply about ΔH<sub>vap</sub> at lower temperatures?
                  </p>
                </div>
                <div className="flex gap-3 items-start">
                  <div className="bg-white p-2 rounded shadow-sm shrink-0 font-bold text-blue-600">Q2</div>
                  <p className="text-sm text-gray-700">
                    Why might the line curve slightly at the very highest temperatures (nearer to critical point)? 
                    <span className="text-gray-500 italic ml-1">(Hint: The assumption that vapor is an Ideal Gas breaks down).</span>
                  </p>
                </div>
              </div>
            </div>
          )}

        </div>
      </div>
    </div>
  );
};

export default ClausiusClapeyronStudio;