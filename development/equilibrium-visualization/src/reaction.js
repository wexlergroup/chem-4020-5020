/*
  2 NO2(g) <=> N2O4(g), starting from 2 mol of NO2, as in the worked example of Section 7.1.
  The standard reaction enthalpy and entropy are the values at 298.15 K from the NIST-JANAF tables, taken as constant
  over the temperature range of the studio, as in the integrated van't Hoff equation of the section.
*/
export const DH = -57.112;                           // kJ/mol
export const DS = -0.17569;                          // kJ/(mol K)
export const R = 0.008314462618;                     // kJ/(mol K)

export function reaction(T, P) {                     // T in K, P in bar (P° = 1 bar)
  const RT = R * T, dG0 = DH - T * DS, K = Math.exp(-dG0 / RT);
  // G(xi) - G(0) in kJ, Eq. (G-xi-NO2) of the notes with the total pressure P
  const g = (xi) => {
    const n1 = 2 - 2 * xi, n2 = xi, n = 2 - xi;
    let mix = 0;
    if (n1 > 0) mix += n1 * Math.log(n1 / n);
    if (n2 > 0) mix += n2 * Math.log(n2 / n);
    return xi * (dG0 - RT * Math.log(P)) + RT * mix;
  };
  const xiEq = 1 - 1 / Math.sqrt(1 + 4 * K * P);     // root of K P = xi (2 - xi) / (2 - 2 xi)^2
  const at = (xi) => {                               // the state at an extent xi
    const y1 = (2 - 2 * xi) / (2 - xi), y2 = xi / (2 - xi);
    const Q = y2 / (y1 * y1 * P);
    return { xi, y1, y2, Q, slope: RT * Math.log(Q / K), G: g(xi), c1: y1 * P / (0.08314462618 * T) };   // c1: mol/L of NO2
  };
  return { RT, dG0, K, g, xiEq, at };
}
