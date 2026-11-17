/*
  Water: the three coexistence curves and the density of each phase.

  Coexistence curves (pressures in Pa): the IAPWS equations for the vapor pressure (SR1-86(1992)) and for the melting
  and sublimation pressures of ice Ih (R14-08(2011)). Densities (kg/m3): the liquid and the gas from regions 1, 2,
  and 3 of IAPWS-IF97, the two fluid phases on the liquid-gas curve from the IAPWS equations for the saturated
  densities (SR1-86(1992)), and ice Ih from a fit of the IAPWS equation of state of 2006 (src/if97.js).
*/
import { R1, R2, R3, N3, ICE } from './if97.js';

export const M = 18.015268e-3;                        // kg/mol
export const NA = 6.02214076e23;                      // 1/mol
export const TRIPLE = { T: 273.16, P: 611.657 };      // K, Pa
export const CRITICAL = { T: 647.096, P: 22.064e6, rho: 322 };   // K, Pa, kg/m3
export const ICE_III = 251.165;                       // K, where the melting curve of ice Ih ends (at 2086 bar)
const R = 461.526;                                    // J/(kg K), specific gas constant of water

const series = (a, e, x) => a.reduce((sum, ai, i) => sum + ai * Math.pow(x, e[i]), 0);

// Liquid-gas curve, from the triple point to the critical point
export function saturationPressure(T) {
  const tau = 1 - T / CRITICAL.T;
  const sum = series([-7.85951783, 1.84408259, -11.7866497, 22.6807411, -15.9618719, 1.80122502], [1, 1.5, 3, 3.5, 4, 7.5], tau);
  return CRITICAL.P * Math.exp((CRITICAL.T / T) * sum);
}
// Solid-gas curve, below the triple point
export function sublimationPressure(T) {
  const theta = T / TRIPLE.T;
  const sum = series([-21.2144006, 27.3203819, -6.10598130], [0.00333333333, 1.20666667, 1.70333333], theta);
  return TRIPLE.P * Math.exp(sum / theta);
}
// Solid-liquid curve of ice Ih. The pressure rises as the temperature falls.
export function meltingPressure(T) {
  const theta = T / TRIPLE.T;
  const a = [1195393.37, 80818.3159, 3338.26860], b = [3.0, 25.75, 103.75];
  return TRIPLE.P * (1 + a.reduce((sum, ai, i) => sum + ai * (1 - Math.pow(theta, b[i])), 0));
}

function solve(f, lo, hi) {                           // the root of f between lo and hi, by bisection
  const rising = f(hi) > f(lo);
  for (let i = 0; i < 60; i++) {
    const mid = 0.5 * (lo + hi);
    if ((f(mid) > 0) === rising) hi = mid; else lo = mid;
  }
  return 0.5 * (lo + hi);
}
// The temperature on each curve at a pressure
export const saturationTemperature = (P) => solve((T) => saturationPressure(T) - P, TRIPLE.T, CRITICAL.T);
export const sublimationTemperature = (P) => solve((T) => sublimationPressure(T) - P, 150, TRIPLE.T);
export const meltingTemperature = (P) => solve((T) => meltingPressure(T) - P, ICE_III, TRIPLE.T);

// The stable phase at a temperature and a pressure. Above the critical temperature and the critical pressure the
// fluid is called supercritical; the two limits are a convention and not coexistence curves.
export function phaseAt(T, P) {
  if (T >= CRITICAL.T) return P >= CRITICAL.P ? 'supercritical' : 'gas';
  if (T < TRIPLE.T) {
    if (P < sublimationPressure(T)) return 'gas';
    return T >= ICE_III && P > meltingPressure(T) ? 'liquid' : 'solid';
  }
  return P > saturationPressure(T) ? 'liquid' : 'gas';
}

// IAPWS-IF97. Pressures in MPa inside these functions.
function v1(T, p) {                                   // region 1, specific volume (m3/kg)
  const pi = p / 16.53, tau = 1386 / T;
  let g = 0;
  for (const [I, J, n] of R1) g -= n * I * Math.pow(7.1 - pi, I - 1) * Math.pow(tau - 1.222, J);
  return pi * g * R * T / (p * 1e6);
}
function v2(T, p) {                                   // region 2
  const tau = 540 / T;
  let g = 1 / p;
  for (const [I, J, n] of R2) g += n * I * Math.pow(p, I - 1) * Math.pow(tau - 0.5, J);
  return p * g * R * T / (p * 1e6);
}
function p3(rho, T) {                                 // region 3, pressure (MPa) at a density
  const delta = rho / CRITICAL.rho, tau = CRITICAL.T / T;
  let f = N3 / delta;
  for (const [I, J, n] of R3) f += n * I * Math.pow(delta, I - 1) * Math.pow(tau, J);
  return rho * R * T * delta * f / 1e6;
}
const p23 = (T) => 0.34805185628969e3 - 0.11671859879975e1 * T + 0.10192970039326e-2 * T * T;   // between regions 2 and 3

// Densities of the liquid and of the gas that coexist at a temperature
export function saturatedLiquidDensity(T) {
  const x = Math.pow(1 - T / CRITICAL.T, 1 / 3);
  return CRITICAL.rho * (1 + series([1.99274064, 1.09965342, -0.510839303, -1.75493479, -45.5170352, -6.74694450e5], [1, 2, 5, 16, 43, 110], x));
}
export function saturatedGasDensity(T) {
  const x = Math.pow(1 - T / CRITICAL.T, 1 / 6);
  return CRITICAL.rho * Math.exp(series([-2.03150240, -2.68302940, -5.38626492, -17.2991605, -44.7586581, -63.9201063], [2, 4, 8, 18, 37, 71], x));
}

function fluidDensity(T, P, liquid) {                 // P in Pa; liquid: which side of the liquid-gas curve
  const p = P / 1e6;
  if (T <= 623.15) return 1 / (liquid ? v1(T, p) : v2(T, p));
  if (p < p23(T)) return 1 / v2(T, p);
  const f = (rho) => p3(rho, T) - p;
  if (T >= CRITICAL.T) return solve(f, 20, 850);                  // region 3 holds up to about 850 kg/m3
  const edge = liquid ? saturatedLiquidDensity(T) : saturatedGasDensity(T);
  if (liquid) return f(edge) >= 0 ? edge : solve(f, edge, 850);
  return f(edge) <= 0 ? edge : solve(f, 20, edge);
}
function iceDensity(T, P) {
  const x = T - TRIPLE.T, p = P / 1e6, c = ICE;
  return 1 / (c[0] + c[1] * x + c[2] * x * x + c[3] * p + c[4] * p * x + c[5] * p * p);
}

// Density (kg/m3) of a phase at a temperature and a pressure
export function density(phase, T, P) {
  if (phase === 'solid') return iceDensity(T, P);
  return fluidDensity(T, P, phase === 'liquid');
}

/*
  The phases present at a state. `on` names the place of the state when it lies on a coexistence curve:
  'lg', 'sg', 'sl', 'triple', or 'critical'; otherwise it is null and one phase is present.
  Each phase: { name, rho (kg/m3), v (molar volume, m3/mol), number (molecules per nm3), spacing (nm) }.
*/
export function phases(T, P, on) {
  let list;
  if (on === 'lg') list = [['liquid', saturatedLiquidDensity(T)], ['gas', saturatedGasDensity(T)]];
  else if (on === 'sg') list = [['solid', iceDensity(T, P)], ['gas', fluidDensity(T, P, false)]];
  else if (on === 'sl') list = [['solid', iceDensity(T, P)], ['liquid', fluidDensity(T, P, true)]];
  else if (on === 'triple') list = [['solid', iceDensity(T, P)], ['liquid', saturatedLiquidDensity(T)], ['gas', saturatedGasDensity(T)]];
  else if (on === 'critical') list = [['critical', CRITICAL.rho]];
  else { const name = phaseAt(T, P); list = [[name, density(name, T, P)]]; }
  return list.map(([name, rho]) => {
    const v = M / rho, number = rho / M * NA * 1e-27;
    return { name, rho, v, number, spacing: Math.pow(1 / number, 1 / 3) };
  });
}
