/*
  Ideal monatomic gas from its partition function (Section 2.5), in SI units.

      Lambda  = (h^2 / 2 pi m k_B T)^(1/2)      thermal de Broglie wavelength
      q_trans = V / Lambda^3                    partition function of one atom
      Q       = q_trans^N / N!                  partition function of N atoms
      U       = (3/2) N k_B T
      P       = N k_B T / V
      S       = U / T + k_B ln Q                (Section 4.3)

  The number of atoms is small here, so ln N! is the exact sum and not the Stirling approximation.
*/
export const H = 6.62607015e-34;                      // J s
export const KB = 1.380649e-23;                       // J/K
export const U_KG = 1.66053907e-27;                   // kg, atomic mass unit

export const GASES = [
  { key: 'he', name: 'Helium', mass: 4.0026 },
  { key: 'ne', name: 'Neon', mass: 20.180 },
  { key: 'ar', name: 'Argon', mass: 39.948 },
];

export function lnFactorial(n) {
  let sum = 0;
  for (let i = 2; i <= n; i++) sum += Math.log(i);
  return sum;
}

// N atoms of mass m (in u) in a cube with the edge L (in nm) at the temperature T (in K)
export function gas(N, edge, T, mass) {
  const m = mass * U_KG, V = (edge * 1e-9) ** 3;
  const lambda = Math.sqrt(H * H / (2 * Math.PI * m * KB * T));
  const q = V / lambda ** 3;
  const lnQ = N * Math.log(q) - lnFactorial(N);
  const U = 1.5 * N * KB * T;
  return {
    V, lambda, q, lnQ, U,
    spacing: Math.cbrt(V / N),                        // m
    P: N * KB * T / V / 1e5,                          // bar
    S: U / T + KB * lnQ,                              // J/K
    speed: Math.sqrt(KB * T / m),                     // m/s, root mean square of one component of the velocity
  };
}
