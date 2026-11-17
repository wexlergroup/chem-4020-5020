/*
  The Clausius–Clapeyron fit of Section 6.2: a straight line through ln(P_sat / 1 bar) against 1/T, by least squares.
  The slope is -Δh_vap / R (the integrated form with a constant enthalpy of vaporization, Eq. 6.2.17 of the notes).
*/
export const R = 8.314462618;                        // J/(mol K)

// points: [t (°C), P (kPa), ...]; from and to: the first and the last point of the fit (indices)
export function fit(points, from, to) {
  const all = points.map(([t, P]) => ({ t, T: t + 273.15, P, x: 1 / (t + 273.15), y: Math.log(P / 100) }));
  const used = all.slice(from, to + 1), n = used.length;
  const mx = used.reduce((s, p) => s + p.x, 0) / n, my = used.reduce((s, p) => s + p.y, 0) / n;
  let sxy = 0, sxx = 0, syy = 0;
  for (const p of used) { sxy += (p.x - mx) * (p.y - my); sxx += (p.x - mx) ** 2; syy += (p.y - my) ** 2; }
  const slope = sxy / sxx, intercept = my - slope * mx;
  for (const p of all) p.residual = p.y - (slope * p.x + intercept);
  const ssr = used.reduce((s, p) => s + p.residual ** 2, 0);
  return { all, used, n, slope, intercept, r2: syy > 0 ? 1 - ssr / syy : 1, dh: -slope * R / 1000, Tmean: 1 / mx };
}

// The reference enthalpy of vaporization at T (K), by linear interpolation in the table of the liquid
export function referenceAt(reference, T) {
  for (let i = 1; i < reference.length; i++) {
    const [T0, h0] = reference[i - 1], [T1, h1] = reference[i];
    if (T <= T1) return h0 + (h1 - h0) * (T - T0) / (T1 - T0);
  }
  return reference[reference.length - 1][1];
}
