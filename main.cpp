#ifndef SPACEDIM
#define SPACEDIM 3
#endif

#include <iostream>
#include <math.h>
#include <sstream>
#include "BosonStar.h"
#include "EvolutionVariables.h"
#include "LinearPerturbation.h"
#include "mathutils.h"
#include <iomanip>


using namespace std;

void slice_convergence_test (BSSNSlice& sl, BSSNSlice& sm, BSSNSlice& sh);
void gauss_initialize(BosonStar& boson_star);
void compute_linear_perturbation(BosonStar& boson_star, double A0, double dA, double n_stars);

int main()
{
    //construct initial BS solution
    BosonStar boson_star{};
    boson_star.read_parameters();

    boson_star.solve();
    boson_star.fill_isotropic_arrays();
    boson_star.write_isotropic();
    boson_star.write_field();

    if (boson_star.gaussian_start)
        gauss_initialize(boson_star);

    //main run: either cycles models/ radial oscillation freqs/ dynamical evolution depending on param choices
    if (boson_star.cycle_only)
        boson_star.cycle_models(boson_star.n_stars, boson_star.A0, boson_star.dA);
    else if (boson_star.pert_only)
        compute_linear_perturbation(boson_star, boson_star.A0, boson_star.dA, boson_star.n_stars);
    else
    {
        Spacetime st{};
        st.initialize(boson_star);
        st.slices[0].write_slice();
        st.write_diagnostics();
        st.evolve();
    }

    cout << "Ending..."  << endl;
    return 0;
}

void gauss_initialize(BosonStar& boson_star)
{

    //TODO NEXT: add proper logic for when to add perturbations + maybe start with random mini BS model chi + alpha or similar

    std::ofstream nm_file{"mass_charge.dat"};

    int num_gaussians = 1; //change this to produce a gaussian model cycle to explore 1-parameter families
    double step = 0.0001;

    double& step_var = boson_star.perturb_amp;

    for (int k = 0; k < num_gaussians; k++)
    {
        boson_star.clear_BS();
        boson_star.omega = boson_star.enforced_freq;
        boson_star.add_perturbation(boson_star.perturb_amp, boson_star.perturb_spread, 0.);

        if (!boson_star.fill_given_A(boson_star.omega, 0) )
        {
            //boson_star.solitonic = 0;
            //boson_star.solve();
            //boson_star
            boson_star.default_metric_vars();
        }
        boson_star.fill_isotropic_arrays();
        boson_star.write_isotropic();
        boson_star.write_field();

        Spacetime st_gauss{};
        st_gauss.initialize(boson_star);
        st_gauss.slices[0].write_slice();
        st_gauss.write_diagnostics();

        nm_file << st_gauss.slice_mass(&st_gauss.slices[0]) << "    " << st_gauss.slice_charge(&st_gauss.slices[0]) <<  "    " <<  step_var << endl;

        step_var += step;
        }
}

//computes oscillation frequencies for sequence of models starting at A0 separated by dA
void compute_linear_perturbation(BosonStar& boson_star, double A0, double dA, double n_stars)
{
    LinearPerturbation lp{&boson_star, -2e-5, 0.025, 0.0001, 30.}; //s=0.08, A = 0.06: 0.00003, 0.18, 0.00003

    //lp.rk4_solve(-0.00001, 0.0235);
    //lp.get_best_gamma(-4.5e-5);
    //lp.get_chi_sq();
    //lp.get_chi_sq_newton();

    lp.read_parameters(0);
    lp.pert_cycle(A0, dA, n_stars); //lp.pert_cycle(0.098, 0.0005, 200);
    cout << "Noether charge perturbation is " << lp.get_noether_perturbation() << endl;
    lp.write_pert();
}

/*BSSNSlice sl{}; BSSNSlice sm{}; BSSNSlice sh{};
    sl.read_BS_data(boson_star);

    boson_star.double_resolution();
    boson_star.solve(1);
    boson_star.fill_isotropic_arrays();
    sm.read_BS_data(boson_star);

    boson_star.double_resolution();
    boson_star.solve(1);
    boson_star.fill_isotropic_arrays();
    sh.read_BS_data(boson_star);

    slice_convergence_test(sl, sm, sh);*/

  /*Spacetime stl{}; Spacetime stm{}; Spacetime sth{};
    stl.initialize(boson_star);
    stl.evolve();
    boson_star.double_resolution();
    boson_star.solve(1);
    boson_star.fill_isotropic_arrays();
    stm.initialize(boson_star);
    stm.evolve();

    boson_star.double_resolution();
    boson_star.solve(1);
    boson_star.fill_isotropic_arrays();
    sth.initialize(boson_star);
    sth.evolve();

    slice_convergence_test(stl.slices[stl.slices.size() - 1], stm.slices[stm.slices.size() - 1], sth.slices[sth.slices.size() - 1]);*/

    //LinearPerturbation lp{&boson_star, 0.0, 0.0, 0.0001, 150.}; //s=0.08, A = 0.06: 0.00003, 0.18, 0.00003

    //FieldState f = (FieldState){0.1, 0.2, 0.3, 0.4};
    //PertState p = (PertState){0.1, 0.2, 0.3, 0.4};
    //lp.test_rhs(1.0, f, p, 0.0145, 6.3);

    //lp.rk4_solve(0.0145, 6.3); //0.00003519, 0.3097 // 0.0, 0.30965
    //lp.get_best_gamma(0.000084); //0.000083
    //lp.get_chi_sq();
    //lp.pert_cycle(0.0001, 0.0001, 150); //lp.pert_cycle(0.098, 0.0005, 200);
    //cout << "Noether charge perturbation is " << lp.get_noether_perturbation() << endl;
    //lp.write_pert();
    //lp.write_chi_results();


    //st.slices[st.slices.size() - 1].write_slice();
