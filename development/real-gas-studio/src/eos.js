/*
  Equations of state of the Real Gas studio (Section 1.4). Units: L, bar, mol, K.

      ideal gas        P = R T / V_m
      van der Waals    P = R T / (V_m - b_m) - a_m / V_m^2
      Peng-Robinson    P = R T / (V_m - b) - a alpha(T) / (V_m^2 + 2 b V_m - b^2)

  The constants of both models follow from the critical temperature and pressure of the gas. Peng-Robinson
  also uses the acentric factor omega of the gas.

  Below the critical temperature, a model has three volumes at one pressure over a range of pressures. The
  state with the lower fugacity is the stable one. Liquid and vapor coexist at the pressure at which the two
  fugacities are equal, which is the tie line of the Maxwell construction.
*/
export const R = 0.08314;                             // L bar / (mol K), as in the notes

export const GASES = [
  { key: 'co2', name: 'Carbon dioxide', Tc: 304.13, Pc: 73.77, omega: 0.224 },
  { key: 'n2', name: 'Nitrogen', Tc: 126.2, Pc: 33.9, omega: 0.037 },
  { key: 'h2o', name: 'Water', Tc: 647.1, Pc: 220.6, omega: 0.344 },
  { key: 'ch4', name: 'Methane', Tc: 190.56, Pc: 45.99, omega: 0.011 },
  { key: 'he', name: 'Helium', Tc: 5.19, Pc: 2.27, omega: -0.385 },
];

const SQRT2 = Math.SQRT2;

// A model at one temperature: its constants a and b, and the pressure at a molar volume
export function vanDerWaals(gas, T) {
  const a = 27 * (R * gas.Tc) ** 2 / (64 * gas.Pc), b = R * gas.Tc / (8 * gas.Pc);
  return {
    a, b, T,
    pressure: (V) => R * T / (V - b) - a / (V * V),
    slope: (V) => -R * T / (V - b) ** 2 + 2 * a / V ** 3,                 // dP/dV_m
    cubic: (P) => {                                                        // Z^3 + c2 Z^2 + c1 Z + c0 = 0
      const A = a * P / (R * T) ** 2, B = b * P / (R * T);
      return [-(1 + B), A, -A * B];
    },
    lnPhi: (Z, P) => {                                                     // logarithm of the fugacity coefficient
      const A = a * P / (R * T) ** 2, B = b * P / (R * T);
      return Z - 1 - Math.log(Z - B) - A / Z;
    },
  };
}

export function pengRobinson(gas, T) {
  const kappa = 0.37464 + 1.54226 * gas.omega - 0.26992 * gas.omega ** 2;
  const alpha = (1 + kappa * (1 - Math.sqrt(T / gas.Tc))) ** 2;
  const a = 0.45724 * (R * gas.Tc) ** 2 / gas.Pc * alpha, b = 0.07780 * R * gas.Tc / gas.Pc;
  return {
    a, b, T,
    pressure: (V) => R * T / (V - b) - a / (V * V + 2 * b * V - b * b),
    slope: (V) => -R * T / (V - b) ** 2 + 2 * a * (V + b) / (V * V + 2 * b * V - b * b) ** 2,
    cubic: (P) => {
      const A = a * P / (R * T) ** 2, B = b * P / (R * T);
      return [-(1 - B), A - 3 * B * B - 2 * B, -(A * B - B * B - B ** 3)];
    },
    lnPhi: (Z, P) => {
      const A = a * P / (R * T) ** 2, B = b * P / (R * T);
      return Z - 1 - Math.log(Z - B) - A / (2 * SQRT2 * B) * Math.log((Z + (1 + SQRT2) * B) / (Z + (1 - SQRT2) * B));
    },
  };
}

// Real roots of Z^3 + c2 Z^2 + c1 Z + c0 = 0, in rising order
export function realRoots([c2, c1, c0]) {
  const q = (3 * c1 - c2 * c2) / 9, r = (9 * c2 * c1 - 27 * c0 - 2 * c2 ** 3) / 54, d = q ** 3 + r * r;
  let roots;
  if (d > 0 || q >= 0) {
    const root = Math.sqrt(Math.max(d, 0));
    roots = [Math.cbrt(r + root) + Math.cbrt(r - root) - c2 / 3];
  } else {
    const theta = Math.acos(Math.max(-1, Math.min(1, r / Math.sqrt(-(q ** 3))))), m = 2 * Math.sqrt(-q);
    roots = [0, 2, 4].map((k) => m * Math.cos((theta + k * Math.PI) / 3) - c2 / 3);
  }
  // A few steps of Newton's method remove the rounding error of the formulas, which matters for small roots
  roots = roots.map((z) => {
    for (let i = 0; i < 4; i++) {
      const f = ((z + c2) * z + c1) * z + c0, slope = (3 * z + 2 * c2) * z + c1;
      if (Math.abs(slope) < 1e-12) break;
      z -= f / slope;
    }
    return z;
  });
  return roots.sort((x, y) => x - y);
}

function bisect(f, lo, hi) {
  let flo = f(lo);
  for (let i = 0; i < 80; i++) {
    const mid = 0.5 * (lo + hi), fmid = f(mid);
    if ((fmid > 0) === (flo > 0)) { lo = mid; flo = fmid; } else hi = mid;
  }
  return 0.5 * (lo + hi);
}

// Coexistence of liquid and vapor in a model: { P, Vl, Vg }, or null at and above the critical temperature
export function coexistence(model) {
  const { b, T } = model;
  // The loop of the isotherm lies between the two volumes at which dP/dV_m = 0
  const turns = [];
  let before = b * 1.0005, slopeBefore = model.slope(before);
  for (let i = 1; i <= 600; i++) {
    const V = b * 1.0005 * 10 ** (3 * i / 600), s = model.slope(V);       // from b to 1000 b
    if ((s > 0) !== (slopeBefore > 0)) turns.push(bisect(model.slope, before, V));
    before = V; slopeBefore = s;
  }
  if (turns.length < 2) return null;
  const low = Math.max(model.pressure(turns[0]), 0), high = model.pressure(turns[1]);
  if (!(high > low)) return null;
  // At a pressure inside the loop: the volume of the liquid (below the first turn) and of the vapor (above the second)
  const volumes = (P) => [
    bisect((V) => model.pressure(V) - P, b * (1 + 1e-12), turns[0]),
    bisect((V) => model.pressure(V) - P, turns[1], R * T / P + b),
  ];
  const difference = (P) => {
    const [Vl, Vg] = volumes(P);
    return model.lnPhi(P * Vl / (R * T), P) - model.lnPhi(P * Vg / (R * T), P);
  };
  const span = high - low, P = bisect(difference, low + 1e-7 * span, high - 1e-9 * span);
  const [Vl, Vg] = volumes(P);
  if (!(Vg - Vl > 1e-9 * Vg)) return null;
  return { P, Vl, Vg };
}

// Compressibility factor of the stable state at the pressure P; sat is the result of coexistence(model)
export function stableZ(model, P, sat) {
  const roots = realRoots(model.cubic(P)).filter((x) => x > model.b * P / (R * model.T));
  if (roots.length === 1 || !sat) return roots[roots.length - 1];
  return P < sat.P ? roots[roots.length - 1] : roots[0];
}
