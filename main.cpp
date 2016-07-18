#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <alps/gf/tail.hpp>
#include "alps/hdf5.hpp"
#include <tclap/CmdLine.h>
#include<gsl/gsl_sf_psi.h>

typedef alps::gf::omega_k_sigma_gf_with_tail cluster_matsubara_kspace_gf;
typedef alps::gf::two_index_gf<double, alps::gf::momentum_index_mesh, alps::gf::index_mesh> cluster_matsubara_kspace_gf_tail;
typedef alps::gf::matsubara_index freq_index;
typedef alps::gf::momentum_index mom_index;
typedef alps::gf::index spin_index;

void read_iteration(int iteration, const std::string &directory, const std::string &GF_h5_filename, cluster_matsubara_kspace_gf &green,
                    cluster_matsubara_kspace_gf &sigma, std::vector<double> &density);

int main(int argc, char **argv) {
    std::cout << "Calculating Energy from Self-Energy and Matsubara GF" << std::endl;
    std::string directory;
    std::string GF_h5_filename, Sim_h5_filename;
    int niteration, siteration;
    int nsites, nfreq, mu, beta, U;

    try{    //command line parser
        TCLAP::CmdLine cmd("Calculate Energy from Matsubara Self Energies and Greens functions");
        TCLAP::ValueArg<std::string> directory_arg("d", "directory", "directory containing SE, GF (.)", false, ".", "string", cmd);
        TCLAP::ValueArg<std::string> GF_filename_arg("G", "greens_h5", "Greens function h5 base filename (greens)", false, "greens_h5_", "string", cmd);
        TCLAP::ValueArg<std::string> h5_filename_arg("x", "hdf5", "hdf5 file with parameters (sim.h5)", false, "sim.h5", "string", cmd);
        TCLAP::ValueArg<int> niteration_arg("n", "niteration", "number of iterations (0)", false, 0, "int", cmd);
        TCLAP::ValueArg<int> siteration_arg("s", "siteration", "starting iteration for energy calculation (0)", false, 0, "int", cmd);

        cmd.parse( argc, argv );

        directory = directory_arg.getValue();
        GF_h5_filename = GF_filename_arg.getValue();
        Sim_h5_filename = h5_filename_arg.getValue();
        niteration = niteration_arg.getValue();
        siteration = siteration_arg.getValue();

        std::cout<<"Directory = "<<directory<<std::endl;
        std::cout<<"GF_filename = "<<GF_h5_filename<<std::endl;
        std::cout<<"h5_filename = "<<Sim_h5_filename<<std::endl;
        std::cout<<"niteration = "<<niteration<<std::endl;
        std::cout<<"siteration = "<<siteration<<std::endl;
        //help, directory, nfreq, nsite, mu, beta, U, selfenergy, gf, niteration, siteration
        //should get nfreq, nsite, mu, beta, and U from h5 file
        //also, self energy and gf file names optionally, default to G_omega,
        // niteration, siteration
        //directory should be optional, default to current directory
        //for now, assume we are working with a single iteration if  niteration=1

    }
    catch(TCLAP::ArgException &e) {
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; exit(1);
    }

    int num_iterations = niteration - siteration + 1;

    //Get parameters from simulation h5 file
    std::stringstream h5_file_fullpath; h5_file_fullpath << directory<< "/" <<Sim_h5_filename;
    alps::hdf5::archive ar(h5_file_fullpath.str().c_str(), alps::hdf5::archive::READ);
    ar.read("/parameters/dca.SITES", nsites);
    ar.read("/parameters/NMATSUBARA", nfreq);
    ar.read("/parameters/U", U);
    ar.read("/parameters/MU", mu);
    ar.read("/parameters/BETA", beta);

    std::cout<<"nsites = "<<nsites<<std::endl;
    std::cout<<"nfreq = "<<nfreq<<std::endl;
    std::cout<<"U = "<<U<<std::endl;
    std::cout<<"mu = "<<mu<<std::endl;
    std::cout<<"beta = "<<beta<<std::endl;
    ar.close();

    cluster_matsubara_kspace_gf dummy_gf(alps::gf::omega_k_sigma_gf(
            alps::gf::matsubara_positive_mesh(beta, nfreq),
            alps::gf::momentum_index_mesh(9,2),
            alps::gf::index_mesh(2))
    );
    std::vector <cluster_matsubara_kspace_gf> green(num_iterations, dummy_gf);
    std::vector <cluster_matsubara_kspace_gf> sigma(num_iterations, dummy_gf);
    std::vector<std::vector<double> > density(num_iterations, std::vector<double>(nsites*2, 0.5));

    /*
     * Read in the Greens function and self-energies
     */
    for (int i = siteration; i < niteration + 1; ++i) {
        read_iteration(i, directory, GF_h5_filename, green[i-siteration], sigma[i-siteration], density[i-siteration]);
    }

    return 0;
}

void read_iteration(int iteration, const std::string &directory, const std::string &GF_h5_filename, cluster_matsubara_kspace_gf &green,
                    cluster_matsubara_kspace_gf &sigma, std::vector<double> &density){
    std::stringstream g_filename; g_filename << directory<< "/" <<GF_h5_filename<<iteration<<".h5";
    alps::hdf5::archive ar(g_filename.str().c_str(), alps::hdf5::archive::READ);

    green.load(ar, "/G_omega");  
    sigma.load(ar, "/Sigma_omega");
    std::vector <cluster_matsubara_kspace_gf_tail> tail;
    tail=green.tail();
    std::cout<<"Num coeffs = "<<tail.size()<<std::endl;
    std::cout<<"Min Order = "<<green.min_tail_order()<<std::endl;
    std::cout<<"Max Order = "<<green.max_tail_order()<<std::endl;

    for(int i=0; i<tail.size(); i++){
        cluster_matsubara_kspace_gf_tail temp = tail[i];
        std::cout<<"Tail "<< i<<std::endl;
        std::cout<<"Mesh1 " <<temp.mesh1().extent()<<std::endl;
        std::cout<<"Mesh2 " <<temp.mesh2().extent()<<std::endl;
        std::cout<<"Mesh1 kind "<<temp.mesh1().kind()<<std::endl;
        for(alps::gf::momentum_index mom(0); mom<temp.mesh1().extent(); mom++){
            for(int dim = 0; dim<temp.mesh1().dimension(); dim++){
                std::cout<<temp.mesh1().points()[mom()][dim]<<"\t";
            }

            for(alps::gf::index spin(0); spin<temp.mesh2().extent(); spin++){
                std::cout<<temp(mom, spin)<<"\t";
            }
            std::cout<<std::endl;
        }

    }
}

void compute_energy(double &kin, double &pot, int nfreq, double  beta, double U, int nsite,
                    cluster_matsubara_kspace_gf &green, cluster_matsubara_kspace_gf &sigma){
    double energy_sigmag_hifreq_1=0, energy_sigmag_hifreq_2=0;
    double hifreq_integral_square=beta*beta*gsl_sf_psi_n (1, 0.5+nfreq)/(4.*M_PI*M_PI); //this takes care of the 1/omega^2 part
    std::vector<double> ebuf(nfreq,0.);
    std::vector<double> I1(nfreq,0.);
    std::vector<double> I2(nfreq,0.);
    double I1tot=0, I2tot=0;
    for(freq_index w(0);w<nfreq;++w){
        double wn((2.*w()+1)*M_PI/beta);
        for(spin_index f(0);f<2;++f){
            for(mom_index p(0);p<nsite;++p){
                I1[w()]+=2./(beta*nsite)*(-wn*green(w,p,f).imag()-1.);
                I2[w()]+=2./(beta*nsite)*(-sigma(w,p,f).imag()*green(w,p,f).imag() + sigma(w,p,f).real()*green(w,p,f).real());
            }
        }
        ebuf[w()]=I1[w()]+I2[w()];
        I1tot+=I1[w()];
        I2tot+=I2[w()];
    }
    //Get tails
    cluster_matsubara_kspace_gf_tail green_c1 = green.tail(1);
    cluster_matsubara_kspace_gf_tail green_c2 = green.tail(2);
    cluster_matsubara_kspace_gf_tail green_c3 = green.tail(3);

    cluster_matsubara_kspace_gf_tail sigma_c0 = sigma.tail(0);
    cluster_matsubara_kspace_gf_tail sigma_c1 = sigma.tail(1);

    //hifreq of sigmag
    {
        for(mom_index i(0);i<nsite;++i){
            for(spin_index f(0);f<2;++f){
                energy_sigmag_hifreq_1 -= 2./(beta*nsite)*(green_c1(i,f)*sigma_c1(i,f))*hifreq_integral_square;
                energy_sigmag_hifreq_2 -= 2./(beta*nsite)*(green_c2(i,f)*sigma_c0(i,f))*hifreq_integral_square;
            }
        }
    }
    //hifreq of G
    double energy_gomega_hifreq=0.;
    {
        for(mom_index p(0);p<nsite;++p){
            for(spin_index f(0);f<2;++f){
                energy_gomega_hifreq -= 2./(beta*nsite)*(green_c3(p,f))*hifreq_integral_square;
            }
        }
    }

    I2tot += (energy_sigmag_hifreq_1+energy_sigmag_hifreq_2);
    I1tot += energy_gomega_hifreq;

    kin = I1tot - I2tot;
    pot = U/4. + I2tot/2.;
}