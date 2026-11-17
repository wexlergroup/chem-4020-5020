/*
  Linear rigid rotor (Section 2.7), in SI units.

      mu        = m_A m_B / (m_A + m_B)                 reduced mass
      I         = mu r^2                                 moment of inertia
      Theta_rot = hbar^2 / (2 k_B I)                     rotational temperature
      B         = h / (8 pi^2 c I)                       rotational constant as a wavenumber
      E_J       = k_B Theta_rot J (J + 1),  g_J = 2J + 1
      q_rot     = (1/sigma) sum_J g_J exp(-E_J / k_B T)  which is close to T / (sigma Theta_rot) for T >> Theta_rot
      p_J       = g_J exp(-E_J / k_B T) / sum
      S_rot     = R [ln q_rot + U_rot / (k_B T)]         molar rotational entropy (Section 4.3)

  The sums run until the terms are negligible. For a molecule of two identical nuclei, the nuclear spins change
  the weights of the levels with even and odd J. The symmetry number sigma = 2 accounts for this at high
  temperature only, and the probabilities p_J here do not include it.
*/
export const H = 6.62607015e-34;                      // J s
export const HBAR = H / (2 * Math.PI);
export const KB = 1.380649e-23;                       // J/K
export const C_CM = 2.99792458e10;                    // cm/s
export const U_KG = 1.66053907e-27;                   // kg
export const R_GAS = 8.314462618;                     // J/(mol K)

// Masses of the most abundant isotopes (u)
const M = { H: 1.00782503, C: 12, N: 14.003074, O: 15.99491462, Cl: 34.96885268, Br: 78.9183376, I: 126.904473 };

// The molecules of the table in Section 2.7, with their rotational constants (cm^-1, NIST Chemistry WebBook)
export const MOLECULES = [
  ['H₂', M.H, M.H, 60.853], ['HCl', M.H, M.Cl, 10.593], ['N₂', M.N, M.N, 1.9982], ['CO', M.C, M.O, 1.9313],
  ['O₂', M.O, M.O, 1.4456], ['Cl₂', M.Cl, M.Cl, 0.24399], ['Br₂', M.Br, M.Br, 0.082107], ['I₂', M.I, M.I, 0.037372],
].map(([name, mA, mB, B]) => {
  const mu = mA * mB / (mA + mB) * U_KG;
  const r = Math.sqrt(H / (8 * Math.PI ** 2 * C_CM * mu * B)) * 1e10;   // bond length in angstrom that gives this B
  return { name, mA, mB, r, same: mA === mB };
});

// A rotor of the masses mA and mB (u) and the bond length r (angstrom) at the temperature T (K)
export function rotor(mA, mB, r, T, sigma) {
  const mu = mA * mB / (mA + mB), I = mu * U_KG * (r * 1e-10) ** 2;
  const theta = HBAR ** 2 / (2 * KB * I);
  const levels = [];
  let sum = 0, energySum = 0;
  for (let J = 0; J < 5000; J++) {
    const g = 2 * J + 1, e = J * (J + 1) * theta, w = g * Math.exp(-e / T);    // e = E_J / k_B in K
    if (J > 2 && w < 1e-14 * sum) break;
    levels.push({ J, e, w });
    sum += w; energySum += w * e;
  }
  for (const level of levels) level.p = level.w / sum;
  const q = sum / sigma, u = energySum / sum;          // u = U_rot / k_B in K
  return {
    mu, I, theta, levels, q, u,
    B: H / (8 * Math.PI ** 2 * C_CM * I),              // cm^-1
    high: T / (sigma * theta),                         // high-temperature value of q_rot
    S: R_GAS * (Math.log(q) + u / T),                  // J/(mol K)
    omega: Math.sqrt(2 * KB * T / I),                  // rad/s, classical root mean square angular velocity
  };
}
