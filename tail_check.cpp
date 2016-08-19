//
// Created by oryx on 8/17/16.
//
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <alps/gf/tail.hpp>
#include "alps/hdf5.hpp"
#include <tclap/CmdLine.h>
#include<gsl/gsl_sf_psi.h>
#include <math.h>

typedef alps::gf::omega_k_sigma_gf_with_tail cluster_matsubara_kspace_gf;
typedef alps::gf::two_index_gf<double, alps::gf::momentum_index_mesh, alps::gf::index_mesh> cluster_matsubara_kspace_gf_tail;
typedef alps::gf::matsubara_index freq_index;
typedef alps::gf::momentum_index mom_index;
typedef alps::gf::index spin_index;

int main(int argc, char **argv){
    std::cout<<"Hi! I check tails!"<<std::endl;

    std::string directory;
    std::string GF_h5_filename, Sim_h5_filename;
    std::string GFmTail_output, SEmTail_output;
    std::string GF_output, SE_output;
    int nsites, nfreq, mu, beta, U;
    bool info;

    try{    //command line parser
        TCLAP::CmdLine cmd("Subtract tail from Greens function, Self-Energy to check validity");
        TCLAP::ValueArg<std::string> directory_arg("d", "directory", "directory containing SE, GF (.)", false, ".", "string", cmd);
        TCLAP::ValueArg<std::string> GF_filename_arg("G", "greens_h5", "Greens function h5 filename (greens)", false, "Greens", "string", cmd);
        TCLAP::ValueArg<std::string> h5_filename_arg("x", "hdf5", "hdf5 file with parameters (sim.h5)", false, "sim.h5", "string", cmd);
        TCLAP::ValueArg<std::string> GFmTail_output_arg("g", "greensmtail_out", "output file with G-Gtail", false, "GmTail.dat", "string", cmd);
        TCLAP::ValueArg<std::string> SEmTail_output_arg("s", "selfenergymtail_out", "output file with SE-SEtail", false, "SEmTail.dat", "string", cmd);
        TCLAP::ValueArg<std::string> GF_output_arg("j", "greens_out", "output file with G", false, "G.dat", "string", cmd);
        TCLAP::ValueArg<std::string> SE_output_arg("t", "selfenergy_out", "output file with SE", false, "SE.dat", "string", cmd);
        TCLAP::ValueArg<bool> info_arg("i", "info", "output info about the calculation", false, false, "bool", cmd);

        cmd.parse( argc, argv );

        directory = directory_arg.getValue();
        GF_h5_filename = GF_filename_arg.getValue();
        Sim_h5_filename = h5_filename_arg.getValue();
        GFmTail_output = GFmTail_output_arg.getValue();
        SEmTail_output = SEmTail_output_arg.getValue();
        GF_output = GF_output_arg.getValue();
        SE_output = SE_output_arg.getValue();
        info = info_arg.getValue();

        if(info) {
            std::cout << "Subtracting tails from GF and SE" << std::endl;
            std::cout << "Directory = " << directory << std::endl;
            std::cout << "GF_filename = " << GF_h5_filename << std::endl;
            std::cout << "h5_filename = " << Sim_h5_filename << std::endl;
            std::cout << "GFmTail output = " << GFmTail_output << std::endl;
            std::cout << "SEmTail output = " << SEmTail_output << std::endl;
            std::cout << "GF output = " << GF_output << std::endl;
            std::cout << "SE output = " << SE_output << std::endl;

        }

    }
    catch(TCLAP::ArgException &e) {
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; exit(1);
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

    cluster_matsubara_kspace_gf green(alps::gf::omega_k_sigma_gf(
            alps::gf::matsubara_positive_mesh(beta, nfreq),
            alps::gf::momentum_index_mesh(9,2),
            alps::gf::index_mesh(2))
    );

    cluster_matsubara_kspace_gf sigma(alps::gf::omega_k_sigma_gf(
            alps::gf::matsubara_positive_mesh(beta, nfreq),
            alps::gf::momentum_index_mesh(9,2),
            alps::gf::index_mesh(2))
    );

    alps::hdf5::archive G_ar(GF_h5_filename.c_str(), alps::hdf5::archive::READ);
    green.load(G_ar, "/G_komega");
    sigma.load(G_ar, "/Sigma_komega");
    G_ar.close();

    //Get tails
    cluster_matsubara_kspace_gf_tail green_c1 = green.tail(1);
    cluster_matsubara_kspace_gf_tail green_c2 = green.tail(2);
    cluster_matsubara_kspace_gf_tail green_c3 = green.tail(3);

    cluster_matsubara_kspace_gf_tail sigma_c0 = sigma.tail(0);
    cluster_matsubara_kspace_gf_tail sigma_c1 = sigma.tail(1);

    std::ofstream G_out;
    std::ofstream SE_out;
    std::ofstream GmTail_out;
    std::ofstream SEmTail_out;
    G_out.open(GF_output);
    SE_out.open(SE_output);
    GmTail_out.open(GFmTail_output);
    SEmTail_out.open(SEmTail_output);

    G_out<<"#Freq\t";
    SE_out<<"#Freq\t";
    GmTail_out<<"#Freq\t";
    SEmTail_out<<"#Freq\t";
    for(mom_index k(0); k<nsites; k++){
        double k0 = green.mesh2().points()[k()][0];
        double k1 = green.mesh2().points()[k()][1];
        G_out<<"Re(k="<<k0<<","<<k1<<",up)\tIm(k="<<k0<<","<<k1<<",up)\tRe(k="<<k0<<","<<k1<<",dn)\tIm(k="<<k0<<","<<k1<<",dn)\t";
        SE_out<<"Re(k="<<k0<<","<<k1<<",up)\tIm(k="<<k0<<","<<k1<<",up)\tRe(k="<<k0<<","<<k1<<",dn)\tIm(k="<<k0<<","<<k1<<",dn)\t";
        GmTail_out<<"Re(k="<<k0<<","<<k1<<",up)\tIm(k="<<k0<<","<<k1<<",up)\tRe(k="<<k0<<","<<k1<<",dn)\tIm(k="<<k0<<","<<k1<<",dn)\t";
        SEmTail_out<<"Re(k="<<k0<<","<<k1<<",up)\tIm(k="<<k0<<","<<k1<<",up)\tRe(k="<<k0<<","<<k1<<",dn)\tIm(k="<<k0<<","<<k1<<",dn)\t";
    }
    G_out<<std::endl;
    SE_out<<std::endl;
    GmTail_out<<std::endl;
    SEmTail_out<<std::endl;
    for(freq_index f(0); f<nfreq; f++){
        double freq = green.mesh1().points()[f()];
        G_out<<freq<<"\t";
        SE_out<<freq<<"\t";
        GmTail_out<<std::log(freq)<<"\t";
        SEmTail_out<<std::log(freq)<<"\t";
        std::complex<double> ifreq(0.0, freq);
        //std::cout<<ifreq<<std::endl;

        for(mom_index k(0); k<nsites; k++){
           for(spin_index s(0); s<2; s++){
               if(f()==0){
                   std::cout<<"k = ("<<green.mesh2().points()[k()][0]<<", "<<green.mesh2().points()[k()][1]<<"), spin = "<<s()<<std::endl;
                   std::cout<<green_c1(k,s) <<"\t"<< green_c2(k,s) <<"\t"<< green_c3(k,s)<<"\t"<<(U*0.5 - 2*(std::cos(green.mesh2().points()[k()][0])+std::cos(green.mesh2().points()[k()][1])) - mu)<<std::endl;
               }
               if(f()==nfreq-1){
                   std::cout<<"k = ("<<green.mesh2().points()[k()][0]<<", "<<green.mesh2().points()[k()][1]<<"), spin = "<<s()<<std::endl;
                   std::cout<<"Fit value of c2(k) = "<<green(f,k,s).real()*ifreq*ifreq<<std::endl;
                   std::cout<<"Attached value = "<<green_c2(k,s)<<std::endl;
                   std::cout<<"Theory value = "<<(U*0.5 - 2*(std::cos(green.mesh2().points()[k()][0])+std::cos(green.mesh2().points()[k()][1])) - (mu+U/2.))<<std::endl;
               }
               std::complex<double> g_val = green(f,k,s);
               std::complex<double> gtail_val = green_c1(k,s)/ifreq + green_c2(k,s)/(ifreq*ifreq) + green_c3(k,s)/(ifreq*ifreq*ifreq);
               //std::complex<double> gtail_val = green_c1(k,s)/ifreq +
               //        (U*0.5 - 2*(std::cos(green.mesh2().points()[k()][0])+std::cos(green.mesh2().points()[k()][1])) - (mu+U/2))/(ifreq*ifreq) +
               //        green_c3(k,s)/(ifreq*ifreq*ifreq);
               std::complex<double> se_val = sigma(f,k,s);
               std::complex<double> setail_val = sigma_c0(k,s) + sigma_c1(k,s)/ifreq;
               G_out<<g_val.real()<<"\t"<<g_val.imag()<<"\t";
               SE_out<<se_val.real()<<"\t"<<se_val.imag()<<"\t";
               GmTail_out<<std::log(std::abs((g_val-gtail_val).real()))<<"\t"<<std::log(std::abs((g_val-gtail_val).imag()))<<"\t";
               SEmTail_out<<std::log(std::abs((se_val-setail_val).real()))<<"\t"<<std::log(std::abs((se_val-setail_val).imag()))<<"\t";
           }
        }

        G_out<<std::endl;
        SE_out<<std::endl;
        GmTail_out<<std::endl;
        SEmTail_out<<std::endl;
    }
    G_out.close();
    SE_out.close();
    GmTail_out.close();
    SEmTail_out.close();

    return 0;
}
