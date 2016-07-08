#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <alps/gf/tail.hpp>
#include "alps/hdf5.hpp"
#include <tclap/CmdLine.h>

typedef alps::gf::omega_k_sigma_gf_with_tail cluster_matsubara_kspace_gf;

int main(int argc, char **argv) {
    std::cout << "Calculating Energy from Self-Energy and Matsubara GF" << std::endl;
    std::string directory;
    std::string GF_h5_filename, Sim_h5_filename;
    int niteration, siteration;
    int nsites, nfreq, mu, beta, U;

    try{    //command line parser
        TCLAP::CmdLine cmd("Calculate Energy from Matsubara Self Energies and Greens functions");
        TCLAP::ValueArg<std::string> directory_arg("d", "directory", "directory containing SE, GF (.)", false, ".", "string", cmd);
        TCLAP::ValueArg<std::string> GF_filename_arg("G", "greens", "Greens function h5 base filename (greens)", false, "greens", "string", cmd);
        TCLAP::ValueArg<std::string> h5_filename_arg("x", "hdf5", "hdf5 file with parameters (sim.h5)", false, "sim.h5", "string", cmd);
        TCLAP::ValueArg<int> niteration_arg("n", "niteration", "number of iterations (1)", false, 1, "int", cmd);
        TCLAP::ValueArg<int> siteration_arg("s", "siteration", "starting iteration for energy calculation (1)", false, 1, "int", cmd);

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
    if(num_iterations == 1){

    }
    else {
        for (int i = siteration; i < niteration + 1; ++i) {

        }
    }
    return 0;
}

void read_iteration(int iteration, const std::string &directory, const std::string &GF_h5_filename, cluster_matsubara_kspace_gf &green,
                    cluster_matsubara_kspace_gf &sigma, std::vector<double> &density){
    std::stringstream g_filename; g_filename << directory<< "/" <<GF_h5_filename<<"_"<<iteration;
    alps::hdf5::archive ar(g_filename.str().c_str(), alps::hdf5::archive::READ);

    green.load(ar, "/Greens");  //TODO I think this loads the high frequency tail as well? Check
    sigma.load(ar, "/Sigma");
}