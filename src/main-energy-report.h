// =============================================================================
// =============================================================================
// ============ Energy reports and the calculation of free energies ============
// =============================================================================
// =============================================================================
class IET_Report {
  public:
  // fundamental parameters
    double mass, mass_mol, density;
  public:
  // basic energies
    double lj, coulsr, coullr, entropy, N, N0, dN, Ng, dNg;
    double Uef0, Uef1; // electric field energy
  // free energy
    double excess_chem[3];  // [0]: GF, [1]: RISM, [2]: hybrid
  public:
    void operator += (IET_Report & o){
      // basic energies
        lj += o.lj; coulsr += o.coulsr; coullr += o.coullr;
        entropy += o.entropy * o.mass / o.mass_mol;
        N += o.N * o.mass / o.mass_mol; N0 += o.N0 * o.mass / o.mass_mol; dN = N - N0;
        Ng += o.Ng * o.mass / o.mass_mol; dNg = Ng - N0;
        Uef0 += o.Uef0;
        Uef1 += o.Uef1;
      // free energy
        for (int i=0; i<sizeof(excess_chem)/sizeof(excess_chem[0]); i++) excess_chem[i] += o.excess_chem[i];
    }
    void operator *= (double scaling){
        mass *= scaling; mass_mol *= scaling; density *= scaling;
        lj *= scaling; coulsr *= scaling; coullr *= scaling; entropy *= scaling;
        N *= scaling; N0 *= scaling; Ng *= scaling; dNg *= scaling;
        Uef0 *= scaling; Uef1 *= scaling;
        for (int i=0; i<sizeof(excess_chem)/sizeof(excess_chem[0]); i++) excess_chem[i] *= scaling;
    }
};
