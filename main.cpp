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
void print_tails(cluster_matsubara_kspace_gf green);
void calculate_density(int nfreq, int nsite, double beta, cluster_matsubara_kspace_gf &green,
                       cluster_matsubara_kspace_gf &sigma, std::vector<double> & density);
void calculate_green_sigma_ave(int M, std::vector<cluster_matsubara_kspace_gf> greens, std::vector<cluster_matsubara_kspace_gf> sigmas,
                               cluster_matsubara_kspace_gf &green_ave, cluster_matsubara_kspace_gf &sigma_ave, bool partial, int partial_index, bool info);
void compute_energy(double &kin, double &pot, int nfreq, double  beta, double U, int nsite,
                    cluster_matsubara_kspace_gf &green, cluster_matsubara_kspace_gf &sigma, bool info);
void print_head_gf(cluster_matsubara_kspace_gf green);

int main(int argc, char **argv) {
    std::string directory;
    std::string GF_h5_filename, Sim_h5_filename;
    int niteration, siteration;
    int nsites, nfreq, mu, beta, U;
    bool info;

    try{    //command line parser
        TCLAP::CmdLine cmd("Calculate Energy from Matsubara Self Energies and Greens functions");
        TCLAP::ValueArg<std::string> directory_arg("d", "directory", "directory containing SE, GF (.)", false, ".", "string", cmd);
        TCLAP::ValueArg<std::string> GF_filename_arg("G", "greens_h5", "Greens function h5 base filename (greens)", false, "Greens_", "string", cmd);
        TCLAP::ValueArg<std::string> h5_filename_arg("x", "hdf5", "hdf5 file with parameters (sim.h5)", false, "sim.h5", "string", cmd);
        TCLAP::ValueArg<int> niteration_arg("n", "niteration", "number of iterations (0)", false, 0, "int", cmd);
        TCLAP::ValueArg<int> siteration_arg("s", "siteration", "starting iteration for energy calculation (0)", false, 0, "int", cmd);
        TCLAP::ValueArg<bool> info_arg("i", "info", "output info about the calculation", false, false, "bool", cmd);

        cmd.parse( argc, argv );

        directory = directory_arg.getValue();
        GF_h5_filename = GF_filename_arg.getValue();
        Sim_h5_filename = h5_filename_arg.getValue();
        niteration = niteration_arg.getValue();
        siteration = siteration_arg.getValue();
        info = info_arg.getValue();

        if(info) {
            std::cout << "Calculating Energy from Self-Energy and Matsubara GF" << std::endl;
            std::cout << "Directory = " << directory << std::endl;
            std::cout << "GF_filename = " << GF_h5_filename << std::endl;
            std::cout << "h5_filename = " << Sim_h5_filename << std::endl;
            std::cout << "niteration = " << niteration << std::endl;
            std::cout << "siteration = " << siteration << std::endl;
        }
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
    if(info) {
        std::cout << "Num Iterations = " << num_iterations << std::endl;
    }
    //Get parameters from simulation h5 file
    std::stringstream h5_file_fullpath; h5_file_fullpath << directory<< "/" <<Sim_h5_filename;
    alps::hdf5::archive ar(h5_file_fullpath.str().c_str(), alps::hdf5::archive::READ);
    ar.read("/parameters/dca.SITES", nsites);
    ar.read("/parameters/NMATSUBARA", nfreq);
    ar.read("/parameters/U", U);
    ar.read("/parameters/MU", mu);
    ar.read("/parameters/BETA", beta);

    if(info) {
        std::cout << "nsites = " << nsites << std::endl;
        std::cout << "nfreq = " << nfreq << std::endl;
        std::cout << "U = " << U << std::endl;
        std::cout << "mu = " << mu << std::endl;
        std::cout << "beta = " << beta << std::endl;
    }
    ar.close();

    //Dummy gf just so we can initialize vector
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

    /*
     * Prepare the jackknife quantities
     */
    cluster_matsubara_kspace_gf template_gf(green[0]);    //Copy one of the read GF, get correct meshes
    template_gf.initialize();  //Set entries to zero
    std::vector<cluster_matsubara_kspace_gf> sigma_i(num_iterations, template_gf);
    cluster_matsubara_kspace_gf sigma_zero(template_gf);
    cluster_matsubara_kspace_gf sigma_bar(template_gf);

    std::vector<cluster_matsubara_kspace_gf> green_i(num_iterations, template_gf);
    cluster_matsubara_kspace_gf green_zero(template_gf);
    cluster_matsubara_kspace_gf green_bar(template_gf);

    /*
     * Calculate each of the needed jackknife quantities
     */
    calculate_green_sigma_ave(num_iterations, green, sigma, green_zero, sigma_zero, false, 0, info);
    for(int z=0; z<num_iterations; ++z){
        calculate_green_sigma_ave(num_iterations, green, sigma, green_i[z], sigma_i[z], true, z, info);
    }
    calculate_green_sigma_ave(num_iterations, green_i, sigma_i, green_bar, sigma_bar, false, 0, info);

    if(info) {
        std::cout<<"G_0"<<std::endl;
        print_head_gf(green[0]);
        std::cout<<"Sigma_0"<<std::endl;
        print_head_gf(sigma[0]);
        std::cout<<"Sigma_zero"<<std::endl;
        print_head_gf(sigma_zero);
        std::cout<<"Sigma_i0"<<std::endl;
        print_head_gf(sigma_i[0]);
        std::cout<<"Sigma_bar"<<std::endl;
        print_head_gf(sigma_bar);
    }

    /*
     * Compute the potential and kinetic energies
     */
    std::vector<double> ekin_i(num_iterations);
    std::vector<double> epot_i(num_iterations);
    double ekin_zero, ekin_bar;
    double epot_zero, epot_bar;
    for(int z=0;z<num_iterations;++z){
        compute_energy(ekin_i[z], epot_i[z], nfreq,  beta, U, nsites, green_i[z], sigma_i[z], info);
    }
    compute_energy(ekin_zero, epot_zero, nfreq,  beta, U, nsites, green_zero, sigma_zero, info);
    compute_energy(ekin_bar, epot_bar, nfreq,  beta, U, nsites, green_bar, sigma_bar, info);

    /***********************************/
    //Jackknife mean and error
    /***********************************/
    double ekin_mean=ekin_zero-(num_iterations-1.)*(ekin_bar-ekin_zero);
    double epot_mean=epot_zero-(num_iterations-1.)*(epot_bar-epot_zero);
    double ekin_isquare=0.;
    double epot_isquare=0.;
    for(int z=0;z<num_iterations;++z){
        ekin_isquare+=ekin_i[z]*ekin_i[z];
        epot_isquare+=epot_i[z]*epot_i[z];
    }
    //std::cout<< ekin_isquare<< " "<< M<< " " << ekin_bar<< " "<<ekin_isquare/M-ekin_bar*ekin_bar<< std::endl;
    //std::cout<< epot_isquare<< " "<< M<< " " << epot_bar<< " "<<epot_isquare/M-epot_bar*epot_bar<< std::endl;
    double ekin_error=std::sqrt(num_iterations-1.)*std::sqrt(fabs(ekin_isquare/num_iterations-ekin_bar*ekin_bar));
    double epot_error=std::sqrt(num_iterations-1.)*std::sqrt(fabs(epot_isquare/num_iterations-epot_bar*epot_bar));

    std::cout<<1.0/nfreq<<" "<<ekin_mean<<" "<<ekin_error<<" "<<epot_mean<<" "<<epot_error<<" "<<ekin_mean+epot_mean<<" "<<nsites<<std::endl;


    return 0;
}

void calculate_green_sigma_ave(int M, std::vector<cluster_matsubara_kspace_gf> greens, std::vector<cluster_matsubara_kspace_gf> sigmas,
                    cluster_matsubara_kspace_gf &green_ave, cluster_matsubara_kspace_gf &sigma_ave, bool partial, int partial_index, bool info){
    //std::cout<<"Gathering tails"<<std::endl;
    //Get tails
    cluster_matsubara_kspace_gf_tail green_c1 = greens[0].tail(1);
    cluster_matsubara_kspace_gf_tail green_c2 = greens[0].tail(2);
    cluster_matsubara_kspace_gf_tail green_c3 = greens[0].tail(3);
    green_c1.initialize();
    green_c2.initialize();
    green_c3.initialize();

    cluster_matsubara_kspace_gf_tail sigma_c0 = sigmas[0].tail(0);
    cluster_matsubara_kspace_gf_tail sigma_c1 = sigmas[0].tail(1);
    sigma_c0.initialize();
    sigma_c1.initialize();

    green_ave.initialize();
    sigma_ave.initialize();

    double N = (partial ? M-1 : M);
    if(int(N) == 0){
        if(info) {
            std::cout << "Cannot calculate partial for only one iteration, returning zero" << std::endl;
            std::cout << "Setting blank tails" << std::endl;
        }
        green_ave.set_tail(1, green_c1);    sigma_ave.set_tail(0, sigma_c0);
        green_ave.set_tail(2, green_c2);    sigma_ave.set_tail(1, sigma_c1);
        green_ave.set_tail(3, green_c3);
        return;
    }
    //std::cout<<"Calculating aves"<<std::endl;
    for(int z=0; z<M; ++z){
        if(partial && partial_index == z) {
            continue;
        }
        green_ave += greens[z];
        sigma_ave += sigmas[z];

        green_c1 += greens[z].tail(1);
        green_c2 += greens[z].tail(2);
        green_c3 += greens[z].tail(3);

        sigma_c0 += sigmas[z].tail(0);
        sigma_c1 += sigmas[z].tail(1);
    }

    green_ave /= N;
    sigma_ave /= N;

    green_c1 /= N;  sigma_c0 /= N;
    green_c2 /= N;  sigma_c1 /= N;
    green_c3 /= N;

    //std::cout<<"Setting tails"<<std::endl;
    green_ave.set_tail(1, green_c1);    sigma_ave.set_tail(0, sigma_c0);
    green_ave.set_tail(2, green_c2);    sigma_ave.set_tail(1, sigma_c1);
    green_ave.set_tail(3, green_c3);

}

void print_head_gf(cluster_matsubara_kspace_gf green){
    freq_index freq(0);
    for(mom_index k(0); k<green.mesh2().extent(); ++k){
        for(spin_index s(0); s<green.mesh3().extent(); ++s){
            std::cout << green(freq, k, s)<<"  ";
        }
    }
    std::cout<<std::endl;
}

void print_tails(cluster_matsubara_kspace_gf green){
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

void read_iteration(int iteration, const std::string &directory, const std::string &GF_h5_filename, cluster_matsubara_kspace_gf &green,
                    cluster_matsubara_kspace_gf &sigma, std::vector<double> &density){
    std::stringstream g_filename; g_filename << directory<< "/" <<GF_h5_filename<<iteration<<".h5";
    alps::hdf5::archive ar(g_filename.str().c_str(), alps::hdf5::archive::READ);

    green.load(ar, "/G_komega");
    sigma.load(ar, "/Sigma_komega");
}

void calculate_density(int nfreq, int nsite, double beta, cluster_matsubara_kspace_gf &green,
                       cluster_matsubara_kspace_gf &sigma, std::vector<double> & density){
    double hifreq_integral_square=beta*beta*gsl_sf_psi_n (1, 0.5+nfreq)/(4.*M_PI*M_PI); //this takes care of the 1/omega^2 part
    cluster_matsubara_kspace_gf_tail green_c2 = green.tail(2);
    for(mom_index i(0);i<nsite;++i){
        for(spin_index f(0);f<2;++f){
            for(freq_index w(0);w<nfreq;++w){
                density[f()*nsite+i()]+=2.*green(w,i,f).real()/beta;
            }
            density[f()*nsite+i()]-=2.*hifreq_integral_square*green_c2(i,f)/beta;
        }
    }
    //real space density for Hartree term
    double density_up=0., density_dn=0.;
    for(int i=0;i<nsite;++i){
        density_up+=density[      i]/(double)nsite;
        density_dn+=density[nsite+i]/(double)nsite;
    }

    std::cout<<"Density Up = "<<density_up<<std::endl;
    std::cout<<"Density Down = "<<density_dn<<std::endl;
}

void compute_energy(double &kin, double &pot, int nfreq, double  beta, double U, int nsite,
                    cluster_matsubara_kspace_gf &green, cluster_matsubara_kspace_gf &sigma, bool info){
    double energy_sigmag_hifreq_1=0, energy_sigmag_hifreq_2=0;
    double hifreq_integral_square=beta*beta*gsl_sf_psi_n (1, 0.5+nfreq)/(4.*M_PI*M_PI); //this takes care of the 1/omega^2 part
    std::vector<double> ebuf(nfreq,0.);
    std::vector<double> I1(nfreq,0.);
    std::vector<double> I2(nfreq,0.);
    std::vector<double> I3(nfreq,0.);
    double I1tot=0, I2tot=0, I3tot=0;
    for(freq_index w(0);w<nfreq;++w){
        double wn((2.*w()+1)*M_PI/beta);
        for(spin_index f(0);f<2;++f){
            for(mom_index p(0);p<nsite;++p){
                I1[w()] += 2./(beta*nsite)*(-wn*green(w,p,f).imag()-1.);
                I2[w()] += 2./(beta*nsite)*(-sigma(w,p,f).imag()*green(w,p,f).imag() + sigma(w,p,f).real()*green(w,p,f).real());
                I3[w()] += 2./(beta*nsite)*green(w,p,f).real();
            }
        }
        ebuf[w()]=I1[w()]+I2[w()];
        I1tot+=I1[w()];
        I2tot+=I2[w()];
        I3tot+=I3[w()];
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
    double energy_gomega_hifreq_3=0.;
    if(info) {
        std::cout << "G.c3, highfreq_int=" << hifreq_integral_square << std::endl;
    }
    {
        for(mom_index p(0);p<nsite;++p){
            for(spin_index f(0);f<2;++f){
                energy_gomega_hifreq -= 2./(beta*nsite)*(green_c3(p,f))*hifreq_integral_square;
                energy_gomega_hifreq_3 -= 2./(beta*nsite)*green_c2(p,f)*hifreq_integral_square;
                //std::cout<<green_c3(p,f)<<std::endl;
            }
        }
    }
    if(info) {
        std::cout << "energy_gomega_hifreq=" << energy_gomega_hifreq << std::endl;
        std::cout << "I1tot=" << I1tot << ",\tI2tot=" << I2tot << std::endl;
    }
    I3tot += energy_gomega_hifreq_3;
    I2tot += (energy_sigmag_hifreq_1+energy_sigmag_hifreq_2);
    I1tot += energy_gomega_hifreq;

    kin = I1tot - I2tot;
    pot = U/4. + I2tot/2. - U*I3tot/4.; //Check the coeff on I3tot term

    if(info) {
        std::cout << "I1tot=" << I1tot << ",\tI2tot=" << I2tot << std::endl;
        std::cout << I3tot << std::endl;
        std::cout << U << std::endl;
        std::cout << "kin=" << kin << ",\tpot=" << pot << std::endl;
    }
}