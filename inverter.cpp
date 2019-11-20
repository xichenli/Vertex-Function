#include"task.h"

/*! \mainpage Fully Irreducible Vertex Inversion Code
 *
 * \section intro_sec Introduction
 *
 * This code takes the result of a continuous-time auxiliary field simulation file and extracts several diagrammatic quantities from it:
 *
 * - Single-particle quantities:
 *  -# The Green's function \f$G\f$
 *  -# The self-energy \f$\Sigma\f$
 * - Two-particle quantities
 *  -# The two-particle Green's function \f$G_4\f$, a.k.a. the four-point function.
 *  -# The generalized susceptibility \f$\chi\f$
 *  -# The bare susceptibility \f$\chi_0\f$
 *  -# The reducible vertex \f$\Gamma\f$
 *  -# The vertex \f$F\f$. \f$F\f$ contains the connected part of the four-point function and is defined by equation 11,
 \f{align}{ \chi = \chi_0 - \frac{1}{\beta^2} \sum \chi_0 F \chi_0 \f}
 *  -# The fully irreducible vertex \f$\Lambda\f$
 *
 * \section Input Input data
 * Data is read from the Monte Carlo hdf5 file, typically called sim.h5, specified by the parameter sim.
 *
 * \subsection input_static Static quantities
 * The density is stored in the hdf5 file at the path:
 * <CODE>/simulation/results/density_up/mean/value</CODE>
 *
 * \subsection input_single Single particle quantities
 * The two-frequency single-particle Green's function is stored in the hdf5 file at the path:
 *
 * <CODE>"/simulation/results/G_omega_omega_up_re"  <<k<<"_"<<k<<"_times_sign/mean/value </CODE>
 *
 * It is frequency diagonal and guaranteed to have the same error structure as \f$\chi\f$.
 *
 * The more precisely measured single frequency Green's function is stored at:
 *
 * <CODE>"/simulation/results/G_omega_up_re"  <<k<<"/mean/value";</CODE>
 *
 * and read in as well. The two are equivalent.
 *
 * \subsection input_two Two-particle quantities
 * The CT-AUX code measures the two-particle Green's function and stores it into
 *
 * <CODE>/simulation/results/G4_Q_K_Kprime_nu_omega_omegaprime_ph_ud_re_"<<Q<<"_"<<K<<"_"<<Kprime<<"_times_sign/mean/value"</CODE>
 *
 * This object is dependent on three frequencies and three K-vectors.
 * \subsection input_lattice Lattice descriptions
 * Lattice descriptions are taken from the ALPS DMFT framework. It provides momenta and symmetries. Have a look at the impurity solver output for a list of momenta and indices.
 * The vienna convention is different from our convention. We have special functions to map our convention to the vienna convention and write vienna output files.
 *
 * \section Output Output Data
 * The output file prefix is <CODE>vert</CODE>. All output files are generated with this prefix.
 * In the case of the four-site cluster we use the 'Vienna' output format.
 * Output files are:
 * -# <CODE>vert_chi_ph</CODE> Particle hole channel \f$\chi_{ph}\f$
 * -# <CODE>vert_chi_pp</CODE> Particle particle channel \f$\chi_{pp}\f$
 * -# <CODE>vert_chi_dm</CODE> Density and magnetic channel \f$\chi_{d}, \chi_m\f$
 * -# <CODE>vert_chi_st</CODE> Singlet and triplet channel \f$\chi_{s}, \chi_d\f$
 *
 * \section Equations Equations
 * The particle hole and particle particle notation is rather weird and illustrated in the following two plots:
 * \image latex ph_notation.pdf "Vertex in the particle-hole notation" width=10cm
 * \image html ph_notation.jpeg "Vertex in the particle-hole notation" width=10cm
 * \image latex pp_notation.pdf "Vertex in the particle-particle notation" width=10cm
 * \image html pp_notation.jpeg "Vertex in the particle-particle notation" width=10cm
 *
 * The big difference between particle hole and particle particle notation is that the momentum transfer in one notation is Q, and in the other notation is Q-K-Kprime. In other words,
 * \f{align}{ \chi^{ph}(K,K',Q) &= \chi(K,K+Q,K'+Q,K')=\chi^{pp}(K,K',Q+K+K')\\
 *            \chi^{pp}(\tilde K,\tilde K',\tilde Q) &= \chi(K,Q-K',Q-K,K')=\chi^{ph}(\tilde K,\tilde K',\tilde Q-\tilde K - \tilde K') \f}
 * The first particle always arrives with momentum \f$K\f$, and the second one always leaves with momentum \f$K'\f$.
 */


int main(int argc, char**argv) {
  
  //create a task
  task *work_to_do=task::task_factory(argc, argv);
  
  //do the work
  work_to_do->work();

  exit(0);
}
