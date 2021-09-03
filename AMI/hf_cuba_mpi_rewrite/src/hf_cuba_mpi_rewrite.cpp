//=======================================================================
// Copyright 2018 JPF LeBlanc
//=======================================================================
#include "amigraph.hpp"
#include <ctime>
#include <Eigen/Dense>
#include "cuba.h"
#include "mpi.h"
#include <string> 
#include <sstream> 
#include <fstream>

///// Parameters


double global_lambda;
double global_rs;
double global_esquare;
double kc=5.0; // cutoff in integration 

int global_integral=1;// use VEGAS by default 
int sigma_max=1;
int bubble_max=2;

double ereg=1e-8;

int AMI_REDUCE_TRIES=40;

bool global_hf=true;// don't need to change this 

int global_seed=0;
double global_epsrel=1e-5;
double global_epsabs=1e-10;
int global_max, global_min, global_maxeval;
double global_nstart=0.01;
int global_use_bare=0;
int global_use_groups=0;
int ZERO_EXTERNAL_Q=0;



////


static int ami_integrand(const int *ndim, const double xx[],
  const int *ncomp, double ff[], void *userdata);
	
	static int group_ami_integrand(const int *ndim, const double xx[],
  const int *ncomp, double ff[], void *userdata);
 
// Assigns V product - pilfered from libamigraph
double get_V(NewAmiCalc::internal_state &state, NewAmiCalc::solution_set &AMI,NewAmiCalc::ext_vars &external, double lambda);
  void get_q_list(std::vector<double> &qs, NewAmiCalc::solution_set &AMI, NewAmiCalc::internal_state &state, NewAmiCalc::ext_vars &external);

//int seed=0; // but this can be set to anything
auto now =std::chrono::high_resolution_clock::now();
auto seed = std::chrono::duration_cast<std::chrono::microseconds>(now.time_since_epoch()).count();

AmiGraph graph(AmiBase::Pi_phud, seed);


// std::vector< double > integrator(NewAmiCalc::solution_set &sol, NewAmiCalc::ext_vars &ext, int nsamples, int seed,std::string SFILE);

void integrator(NewAmiCalc::solution_set &sol, NewAmiCalc::ext_vars &ext, int nsamples, int seed, int eff_ord, int ct, int sigct,int num, int extnum, std::vector<double> &out);

void group_integrator(NewAmiCalc::solution_set_vec_t &ssamiv, NewAmiCalc::ext_vars &ext, int nsamples, int seed, int eff_ord, int ct, int sigct,int num, int extnum, std::vector<double> &out);

// NewAmiCalc::ext_vars &ext, int nsamples, int seed, const char* STATEFILE);

// group_integrator(GG_AMI_MATRIX[ord][num], extern_list[extvar], nsamples, intseed, ord+1,GG_AMI_MATRIX[ord][num].ct_count_ ,AMI_MATRIX[ord][num].sigma_ct_count_, num, extvar, this_result);


void load_settings();
void load_hf();


struct evaluation_vector_set{
evaluation_vector_set(){}

evaluation_vector_set(std::vector<NewAmiCalc::solution_set> s,	NewAmiCalc::ext_vars ext){
ext_vars_=ext;
sol_=s;	
}
	
	
NewAmiCalc::ext_vars ext_vars_;
std::vector<NewAmiCalc::solution_set> sol_;
	
};

void zero_external(AmiGraph::graph_t &g, int index);



// NOte, 3rd order translation is ( 4,2,1,3,5,6)

int main(int argc, char *argv[])
{
  
int mpi_rank;

int comm_size;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);	
		
		
graph.mpi_rank=mpi_rank;
graph.mpi_size=comm_size;
  
	
  
  int max=2;
  int min=0;
  int maxeval=10000;
  int intseed=0;
  
  bool use_groups=false;//true;
	bool use_bare=false;
  
  
  if(argc==1) {
        printf("\nNo Extra Command Line Argument Passed Other Than Program Name");
				printf("Will try loading parameters from paramters.param \nNo Extra Command Line Argument Passed Other Than Program Name");		
load_settings();	

max=global_max;
min=global_min;
maxeval=global_maxeval;
intseed=global_seed;
use_bare=bool(global_use_bare);

				// exit(0);
	}
    if(argc==2) 
    { 
        max=atoi(argv[1]);
    }  
	if(argc==3) 
    { 
        min=atoi(argv[1]);
		max=atoi(argv[2]);
    } 
	if(argc==4) 
    { 
        min=atoi(argv[1]);
		max=atoi(argv[2]);
		maxeval=atoi(argv[3]);
				
    } 
	if(argc==5) 
    { 
        min=atoi(argv[1]);
		max=atoi(argv[2]);
		maxeval=atoi(argv[3]);
		intseed=atoi(argv[4]);
		
    } 
	if(argc==6) 
    { 
        min=atoi(argv[1]);
		max=atoi(argv[2]);
		maxeval=atoi(argv[3]);
		intseed=atoi(argv[4]);
		global_integral=atoi(argv[5]);
    } 
		
		if(argc==7) 
    { 
        min=atoi(argv[1]);
		max=atoi(argv[2]);
		maxeval=atoi(argv[3]);
		intseed=atoi(argv[4]);
		global_integral=atoi(argv[5]);
		use_bare=atoi(argv[6]);
    } 
		
		
	// if(argc==7) 
    // { 
        // min=atoi(argv[1]);
		// max=atoi(argv[2]);
		// maxeval=atoi(argv[3]);
		// intseed=atoi(argv[4]);
		// global_integral=atoi(argv[5]);
		// use_groups=(bool)atoi(argv[6]);
    // } // TODO the grouped case hangs the MPI process if the number of groups to be evaluated is too small 

graph.ami.amibase.drop_bosonic_diverge=false;
graph.ami.amibase.drop_matsubara_poles=true;
// graph.ami.amibase.is_real_external=false;//true;

/* std::cout<<"M=0 "<< graph.ami.amibase.fermi_bose(0,1.0,20.0,0.3)<<std::endl;
std::cout<<"M=1 "<< graph.ami.amibase.fermi_bose(1,1.0,20.0,0.3)<<std::endl;
std::cout<<"M=2 at positive E"<< graph.ami.amibase.fermi_bose(2,1.0,20.0,0.3)<<std::endl;
std::cout<<"M=2 at negative E"<< graph.ami.amibase.fermi_bose(2,1.0,20.0,-0.3)<<std::endl;
std::cout<<"M=3 "<< graph.ami.amibase.fermi_bose(3,1.0,20.0,0.3)<<std::endl;
std::cout<<"M=4 "<< graph.ami.amibase.fermi_bose(4,1.0,20.0,0.3)<<std::endl;

std::cout<<"M=0 "<< graph.ami.amibase.fermi_bose(0,-1.0,20.0,0.3)<<std::endl;
std::cout<<"M=1 "<< graph.ami.amibase.fermi_bose(1,-1.0,20.0,0.3)<<std::endl;
std::cout<<"M=2 "<< graph.ami.amibase.fermi_bose(2,-1.0,20.0,0.3)<<std::endl;
std::cout<<"M=3 "<< graph.ami.amibase.fermi_bose(3,-1.0,20.0,0.3)<<std::endl;
std::cout<<"M=4 "<< graph.ami.amibase.fermi_bose(4,-1.0,20.0,0.3)<<std::endl;
 */

load_hf();
/* 
/////////////
					// read lambda rs etc
						
					std::ifstream infile_stream;
					infile_stream.open("input_hf.dat");

					if(infile_stream.fail()) // checks to see if file opended 
							{ 
							throw std::runtime_error("Could not open input_hf.dat file");
							} 	
						
					std::string line;
					std::getline(infile_stream,line);	
					
					std::stringstream ss;
					std::getline(infile_stream,line);
					
					ss<< line;
					ss>>global_rs;
					std::getline(infile_stream,line);
					
					std::stringstream ss2;
					ss2<<line;
					double kappa;
					
					ss2>>kappa;
					global_lambda=kappa*kappa;
					global_esquare=global_rs*pow(32.0/(9*M_PI), 1.0/3.0);

// std::cout<<kappa<<" "<<test<<" "<<test2<<std::endl; */
std::cout<<"Running for lambda="<<global_lambda<<" rs="<<global_rs<<" esquare="<<global_esquare<<std::endl;
 

///////////// 



std::ofstream checkfile;
checkfile.open("fermi_der.dat");

std::ofstream partfile;
partfile.open("partial.dat",std::ofstream::out | std::ofstream::app);
// file<<"# order num ext_index eff_ord ct_count sigma_ct_count ext_freq Im_ext_freq Beta reMu imMu kx ky kz Re err_Re Im err_Im"<<std::endl;
for(double x=-2; x<2; x+=0.005){
	
checkfile<<x<<" "<<	std::real(graph.ami.amibase.fermi_bose(2,1.0,20.0,x))<<std::endl;
	
	
}
checkfile.close();


//true;//false;//true;//false;//true;
std::cout<<min<<" "<<max<<" "<<maxeval<<std::endl;

  

//int seed=0;

std::cout<< seed <<" " <<std::endl;
// graph=new AmiGraph(AmiCalc::Pi_phud, seed);//g(AmiCalc::Pi,seed);//g(AmiCalc::Hartree,seed);



AmiGraph::vertex_vector_list_t bubble_vertex_list;
AmiGraph::edge_vector_list_t bubble_edge_list;
AmiGraph::edge_vector_list_t legs_edge_list;
graph.bubble_finder(graph.current_graph, bubble_vertex_list,  bubble_edge_list, legs_edge_list);

AmiGraph::edge_t one, two;




	
		std::cout<<"Reading external parameters from ext_vars.dat"<<std::endl;
		NewAmiCalc::external_variable_list extern_list;

		std::string infile("ext_vars.dat");
		graph.ami.read_external(infile, extern_list);	
		std::cout<<"external parameters read on rank "<<mpi_rank<<" are"<<std::endl;
for(int i=0; i<extern_list.size();i++){
		std::cout<<extern_list[i].BETA_<<" "<<extern_list[i].MU_<<" "<< extern_list[i].H_<<" "<<extern_list[i].KDIM_<<" "<<extern_list[i].external_k_list_[0][0]<<" "<<extern_list[i].external_k_list_[0][1] <<" "<<extern_list[i].external_k_list_[0][2]<<" "<<extern_list[i].external_freq_[0]<<std::endl;
		}	
	
std::vector< std::vector< AmiGraph::graph_t>> graph_matrix;
// graph.generate_bubble_graphs(graph_matrix,max,mpi_rank);

// bool success;


// graph.label_graphs(graph_matrix,1,max);
// // graph.min_ext_label_counts(graph_matrix);
// graph.reduce_gm_rf(graph_matrix,mpi_rank);
// graph.reduce_gm_tp(graph_matrix, mpi_rank, 1);
// //graph.reduce_gm_oneleg(graph_matrix, mpi_rank);
// graph.reduce_gm_1PBose(graph_matrix, mpi_rank,1);


// graph.reduce_gm_ladder(graph_matrix, mpi_rank,1);

AmiGraph::gg_matrix_t ggm;

			// ggm.resize(graph_matrix.size());
			// graph.read_ggmp("../../coulomb_graphs/ggm_coulomb_irr_nofock_RC2/",ggm, max);
			graph.read_ggmp("../../coulomb_graphs/ggm_coulomb_irr/",ggm, max);
			// graph.read_ggmp("../../coulomb_graphs/ggm_coulomb_irr_grouped/",ggm, max);
			// graph.read_ggmp("../../coulomb_graphs/ggm_coulomb_o2_ladder/",ggm, max);
			std::cout<<"Completed read"<<std::endl;


			// graph.label_graphs_sys(graph_matrix,0,max);
			// graph.gm_to_ggm(graph_matrix, ggm);
			graph.mpi_print_ggm(ggm,mpi_rank);
			
			// for(int ord=0; ord<ggm.size(); ord++){
				// for(int group=0; group<ggm[ord].size(); group++){
					// for(int num=0; num<ggm[ord][group].graph_vec.size(); num++){
						
						// std::cout<<"CT count is "<<ggm[ord][group].graph_vec[0][boost::graph_bundle].ct_count<<std::endl;
						
					// }
				// }
				
			// }

			// std::cout<<"Filtering out fock diagrams"<<std::endl;
			// graph.ggm_filter_fock(ggm);
			// graph.mpi_print_ggm(ggm,mpi_rank);

			// return 0;


			// purge anything lower than min order from ggm 
			for(int i=0; i<min; i++){
				ggm[i].clear();
			}



			graph.mpi_print_ggm(ggm, mpi_rank);

			graph.ggm_label(ggm, 0);

			std::vector< AmiGraph::graph_t > ct_vec;
			AmiGraph::gg_matrix_t ct_ggm;
			ct_ggm.resize(ggm.size());
			for(int i=0;i<ct_ggm.size(); i++){

			ct_ggm[i].resize(1);	
				
			}
			// ct_groups.resize(max);

			// graph.print_all_edge_info(ggm[2][0].graph_vec[0]);
			for(int m=0; m< max; m++){
				
				// AmiGraph::graph_group ct_group; 
				for(int g=0; g<ggm[m].size(); g++){
					for(int n=0; n< ggm[m][g].graph_vec.size(); n++){
						std::vector< AmiGraph::graph_t > ct_temp;
						
						// std::cout<<"On "<<m<<" "<< g<<" "<<n<<std::endl;
			graph.generate_bubble_ct( ggm[m][g].graph_vec[n] , ct_temp);
			// std::cout<<"done"<<std::endl;
			///////*******Sorting of bubble ct is suspect *********//////////
			
			// std::cout<<"CT temp size "<<ct_temp.size()<<std::endl;
			for(int i=0; i<ct_temp.size(); i++){
				// std::cout<<"Alpha size and ct_count are "<<graph.alpha_size(ct_temp[i])<<" "<<ct_temp[i][boost::graph_bundle].ct_count<<std::endl;
			int rel_ord=graph.alpha_size(ct_temp[i])-2;  //         -2;
			//lets insert the ct into the ct_ggm in the appropriate spot 
			ct_ggm[rel_ord][0].graph_vec.push_back(ct_temp[i]);
			}
			///////****************//////////

			ct_vec.insert(ct_vec.end(), ct_temp.begin(),ct_temp.end());
			// std::cout<<"Total size is "<< ct_vec.size()<<std::endl;
			// ct_temp.clear();
			}}


			}


			// for(int m=0; m< ct_ggm.size();m++){
			// for(int i=0; i< ct_ggm[m].size();i++){
			// for(int j=0; j< ct_ggm[m][i].graph_vec.size(); j++){
				// std::cout<<"Printing graph "<<m<<" "<<i<<" "<<j<<" with ct count "<<ct_ggm[m][i].graph_vec[j][boost::graph_bundle].ct_count<< std::endl;
			// graph.print_all_edge_info(ct_ggm[m][i].graph_vec[j]);
			 // graph.ami.print_final(m+1, ggm[m][i].ss_vec[j].R_, ggm[m][i].ss_vec[j].P_, ggm[m][i].ss_vec[j].S_);
			// }
			// }} 

			// return 0;



			std::cout<<"Returned counter terms "<< ct_vec.size()<<std::endl;
			/* 
			for(int i=0; i< ct_vec.size(); i++){
			std::cout<<"Counter term "<< i<<std::endl;
			std::cout<<"Has "<< graph.count_fermi_loops(ct_vec[i])<<" loops"<<std::endl;;
			// graph.print_all_edge_info(ct_vec[i]);
			} */

			std::cout<<"ct_ggm contains "<<std::endl;
			graph.mpi_print_ggm(ct_ggm, mpi_rank);

			// std::cout<<"ct_ggm without fock contains "<<std::endl;
			// graph.ggm_filter_fock(ct_ggm);
			// graph.mpi_print_ggm(ct_ggm,mpi_rank);

			// return 0;

			std::cout<<"ggm contains "<<std::endl;
			graph.mpi_print_ggm(ggm, mpi_rank);


			std::vector< std::vector< AmiGraph::graph_t>> temp_gm;
			graph.ggm_to_gm(ct_ggm,temp_gm,0);
			graph.gm_to_ggm(temp_gm,ct_ggm);

			graph.merge_ggm(ggm,ct_ggm);
			std::cout<<"AFTER merge ggm contains "<<std::endl;
			graph.mpi_print_ggm(ggm, mpi_rank);


			// working on sigma CT 

			std::vector< AmiGraph::graph_t > sig_ct_vec;
			AmiGraph::gg_matrix_t sig_ct_ggm;
			sig_ct_ggm.resize(ggm.size());
			for(int i=0;i<ct_ggm.size(); i++){

			sig_ct_ggm[i].resize(1);	
				
			}
			for(int m=0; m< max; m++){
				
				// AmiGraph::graph_group ct_group; 
				for(int g=0; g<ggm[m].size(); g++){
					for(int n=0; n< ggm[m][g].graph_vec.size(); n++){
						std::vector< AmiGraph::graph_t > ct_temp;
						
						// std::cout<<"On "<<m<<" "<< g<<" "<<n<<std::endl;
						// graph.print_all_edge_info(ggm[m][g].graph_vec[n]);


			graph.generate_sigma_ct( ggm[m][g].graph_vec[n] , ct_temp, sigma_max);


			// std::cout<<"Generated sigma counter terms "<<std::endl;
			// std::cout<<"sigma ct generated with size "<<ct_temp.size()<<std::endl;
			// for(int temp=0; temp<ct_temp.size(); temp++){
			// std::cout<<"CT "<<temp<<std::endl;
			// graph.print_all_edge_info(ct_temp[temp]);
				
			// }



			// std::cout<<"sigma ct generated with size "<<ct_temp.size()<<std::endl;

			for(int i=0; i<ct_temp.size(); i++){
			int rel_ord=graph.alpha_size(ct_temp[i])-2;

			// std::cout<<i<<" "<<rel_ord<<std::endl;
			//lets insert the ct into the ct_ggm in the appropriate spot 
			sig_ct_ggm[rel_ord][0].graph_vec.push_back(ct_temp[i]);
			}

			sig_ct_vec.insert(sig_ct_vec.end(), ct_temp.begin(),ct_temp.end());
			// std::cout<<"Total size is "<< sig_ct_vec.size()<<std::endl;
			// ct_temp.clear();
			}}


			}

			std::cout<<"sigma_ct_ggm contains "<<std::endl;
			graph.mpi_print_ggm(sig_ct_ggm, mpi_rank);
			// for(int ord=0; ord< sig_ct_ggm.size(); ord++){
			// for(int i=0; i< sig_ct_ggm[ord][0].graph_vec.size(); i++){
			// std::cout<<"On sig ct "<<ord<<" "<<i<<std::endl;
			// graph.print_all_edge_info(sig_ct_ggm[ord][0].graph_vec[i]);
			// std::cout<<std::endl<<std::endl;

			// }
			// }
			// graph.save_ct_eff_orders(sig_ct_ggm);

			// return 0;

			graph.ggm_to_gm(sig_ct_ggm,temp_gm,0);
			graph.gm_to_ggm(temp_gm,sig_ct_ggm);
			// std::cout<<"Testing line 1"<<std::endl;
			// graph.print_all_edge_info(sig_ct_ggm[1][0].graph_vec[0]);

			graph.merge_ggm(ggm,sig_ct_ggm);
			// ggm=sig_ct_ggm;
			// std::cout<<"Testing line 3"<<std::endl;
			// graph.print_all_edge_info(ggm[1][0].graph_vec[0]);
			// ggm[1].clear();
			// ggm[1].push_back(sig_ct_ggm[1][0]);
			// std::cout<<"Testing line 4"<<std::endl;
			// graph.print_all_edge_info(ggm[1][0].graph_vec[0]);

			std::cout<<"AFTER sigma ct merge ggm contains "<<std::endl;
			graph.mpi_print_ggm(ggm, mpi_rank);


			std::cout<<"Filtering out fock diagrams"<<std::endl;
			graph.ggm_filter_fock(ggm);
			graph.mpi_print_ggm(ggm,mpi_rank);
			
			// graph.ggm_filter_hfspinsusc(ggm);
			// graph.mpi_print_ggm(ggm,mpi_rank);
			
			
			std::cout<<"Filtering out diagrams based on max effective order"<<std::endl;
			graph.ggm_filter_max_eff_order(ggm,max,bubble_max,sigma_max);
			graph.mpi_print_ggm(ggm,mpi_rank);


			// return 0;



			// graph.merge_ggm(ggm,ct_ggm);
			std::cout<<"AFTER merge ggm contains "<<std::endl;
			graph.mpi_print_ggm(ggm, mpi_rank);

			// std::cout<<"Testing line 2"<<std::endl;
			// graph.print_all_edge_info(ggm[1][0].graph_vec[0]);

			std::cout<<"Grouping by effective orders "<<std::endl;
			graph.ggm_group_eff_order_CT(ggm, bubble_max, sigma_max);
			graph.mpi_print_ggm(ggm, mpi_rank);

			graph.save_ct_eff_orders(ggm);

/* 
			for(int ord=0; ord<ggm.size(); ord++){
				for(int group=0; group<ggm[ord].size(); group++){
					for(int num=0; num<ggm[ord][group].graph_vec.size(); num++){
						
						std::cout<<"AFTER merge CT count is "<<ggm[ord][group].graph_vec[0][boost::graph_bundle].ct_count<<std::endl;
						
					}
				}
				
			}
 */

			// return 0;

			/// ZERO EXTERNAL PARAMETERS
			// for Q=0 and Omega=0 set externals to zero 
			/* 
			for(int m=0; m< ggm.size();m++){
			for(int i=0; i< ggm[m].size();i++){
			for(int j=0; j< ggm[m][i].graph_vec.size(); j++){
				// std::cout<<"-----"<<std::endl;
				// std::cout<<"Zeroing external index="<<m+2<<std::endl;
				// graph.print_all_edge_info(ggm[m][i].graph_vec[j]);
				zero_external(ggm[m][i].graph_vec[j], m+1);
				// std::cout<<std::endl;
				// graph.print_all_edge_info(ggm[m][i].graph_vec[j]);
				// std::cout<<"-----"<<std::endl;
			}
			}}
				*/

			//////

			// return 0;


			// graph.print_all_edge_info(ggm[2][0].graph_vec[0]);

			// return 0;


			// if(all){

			// std::vector< std::vector< AmiGraph::graph_t>> temp_graph_matrix;
			// graph.ggm_to_gm(ggm,temp_graph_matrix, 0);
			// graph.gm_to_ggm_all(temp_graph_matrix, ggm,0);	
				
			// }


			// graph.mpi_print_ggm(ggm, mpi_rank);

if(ZERO_EXTERNAL_Q){
std::cout<<"Zeroing external Q"<<std::endl;

for(int m=0; m< ggm.size();m++){
for(int i=0; i< ggm[m].size();i++){
for(int j=0; j< ggm[m][i].graph_vec.size(); j++){
	// std::cout<<"-----"<<std::endl;
	// std::cout<<"Zeroing external Q"<<std::endl;
	// graph.print_all_edge_info(ggm[m][i].graph_vec[j]);
	graph.zero_external_Q(ggm[m][i].graph_vec[j]);
	// std::cout<<std::endl;
	// graph.print_all_edge_info(ggm[m][i].graph_vec[j]);
	// std::cout<<"-----"<<std::endl;
}
}}	
	
}

			graph.ggm_construct_ami_sol(ggm, ereg, mpi_rank);
			
			

			// graph.print_all_edge_info(ggm[2][2].graph_vec[0]);
			// graph.ami.print_final(3, ggm[2][2].ss_vec[0].R_, ggm[2][2].ss_vec[0].P_, ggm[2][2].ss_vec[0].S_);
			// std::cout<<"Sol size is "<<ggm[2][2].ss_vec[0].R_[3].size()<<std::endl;

			bool success=false;


			// for(int n=0; n< 100; n++){

			// int first=graph.random_int(0,2);
			// int second=graph.random_int(0,2);
			// graph.swap_alphas(ggm[2][2].graph_vec[0], 0,2);


			// std::cout<<"After Swap"<<std::endl;
			// graph.print_all_edge_info(ggm[2][2].graph_vec[0]);

			// graph.gg_construct_ami_sol(ggm[2][2], 1e-8);
			// std::cout<<"Reconstructed Sol size is "<<ggm[2][2].ss_vec[0].R_[3].size()<<std::endl;
			// }
			 // graph.ami.print_final(3, ggm[2][2].ss_vec[0].R_, ggm[2][2].ss_vec[0].P_, ggm[2][2].ss_vec[0].S_);

			// return 0;



			std::ofstream sizefile;
			sizefile.open("initial_sol_sizes.dat");
			int total=0;
			for(int m=0; m< ggm.size();m++){
			for(int i=0; i< ggm[m].size();i++){
			for(int j=0; j< ggm[m][i].graph_vec.size(); j++){
				// std::cout<<"Printing graph "<<m<<" "<<i<<" "<<j<<std::endl;
			// graph.print_all_edge_info(ggm[m][i].graph_vec[j]);
			// graph.ami.print_final(m+1, ggm[m][i].ss_vec[j].R_, ggm[m][i].ss_vec[j].P_, ggm[m][i].ss_vec[j].S_);
			sizefile<<m<<" "<<i<<" "<<j<<" "<<ggm[m][i].ss_vec[j].R_[m+1].size()<<std::endl;

			total+=ggm[m][i].ss_vec[j].R_[m+1].size();

			}
			}}

			sizefile.close();

int after_total=0;
if(mpi_rank==0){
	
	// std::cout<<"Before"<<std::endl;
	// for(int m=0; m< ggm.size();m++){
			// for(int i=0; i< ggm[m].size();i++){
			// for(int j=0; j< ggm[m][i].ss_vec.size(); j++){
				// std::cout<<ggm[m][i].ss_vec[j].ct_count_<<" "<<ggm[m][i].ss_vec[j].sigma_ct_count_<<std::endl;
				
	// }}}
	
	
	
	
			graph.ggm_reduce_ami_terms(ggm, ereg, mpi_rank, AMI_REDUCE_TRIES);

// std::cout<<"After"<<std::endl;
	// for(int m=0; m< ggm.size();m++){
			// for(int i=0; i< ggm[m].size();i++){
			// for(int j=0; j< ggm[m][i].ss_vec.size(); j++){
				// std::cout<<ggm[m][i].ss_vec[j].ct_count_<<" "<<ggm[m][i].ss_vec[j].sigma_ct_count_<<std::endl;
				
	// }}}
	// exit(0);
	


			std::ofstream aftersizefile;
			aftersizefile.open("after_sol_sizes.dat");
			
			for(int m=0; m< ggm.size();m++){
			for(int i=0; i< ggm[m].size();i++){
			for(int j=0; j< ggm[m][i].graph_vec.size(); j++){
				// std::cout<<"Printing graph "<<m<<" "<<i<<" "<<j<<std::endl;
			// graph.print_all_edge_info(ggm[m][i].graph_vec[j]);
			// graph.ami.print_final(m+1, ggm[m][i].ss_vec[j].R_, ggm[m][i].ss_vec[j].P_, ggm[m][i].ss_vec[j].S_);
			aftersizefile<<m<<" "<<i<<" "<<j<<" "<<ggm[m][i].ss_vec[j].R_[m+1].size()<<std::endl;

			after_total+=ggm[m][i].ss_vec[j].R_[m+1].size();

			}
			}}

			std::cout<<"Solution sizes went from "<<total<<" to "<<after_total<<std::endl;
			aftersizefile.close();
			
			graph.ggm_optimize_ami(ggm,mpi_rank);
}




graph.broadcast_ggm(ggm,0); // note, this does not communicate the graphs, only the solutions



			// for(int ord=0; ord<ggm.size(); ord++){
				// for(int group=0; group<ggm[ord].size(); group++){
					// for(int num=0; num<ggm[ord][group].graph_vec.size(); num++){
						
						// std::cout<<"AFTER broadcast CT count is "<<ggm[ord][group].graph_vec[0][boost::graph_bundle].ct_count<<std::endl;
						
					// }
				// }
				
			// }
/* 
for(int m=0; m< ggm.size();m++){
			for(int i=0; i< ggm[m].size();i++){
			for(int j=0; j< ggm[m][i].graph_vec.size(); j++){
std::cout<<"This is rank "<<mpi_rank <<" with unique "<<ggm[m][i].ss_vec[j].Unique.size()<<std::endl;

			}
			}
} */

// hf stuff

graph.current_state.disp_=AmiBase::hf;
if(use_bare){
	graph.current_state.disp_=AmiBase::fp;
  global_hf=false;	
}




std::cout<<"Reading ptoi.dat, pgrid.dat and Sigma_HF.dat"<<std::endl;
graph.ami.read_hf("pgrid.dat", "ptoi.dat", "Sigma_HF.dat");
std::cout<<"Completed on rank "<< mpi_rank<<std::endl;

std::cout<<"Mu is "<< graph.ami.hf_mu<<std::endl;

if(mpi_rank==0){
	std::cout<<"Dumping epsilon to file"<<std::endl;
std::ofstream efile;
efile.open("eps_dump.dat");


for(double q=0; q< 16; q=q+.05){

efile <<q<<" "<< graph.ami.hf_energy(q)<<std::endl;	
	
	
}

efile.close();
}

//



NewAmiCalc::solution_set_matrix_t AMI_MATRIX;
NewAmiCalc::gg_solution_set_matrix_t GG_AMI_MATRIX;

if (! bool(global_use_groups)){
std::cout<<"ggm to amim"<<std::endl;
graph.ggm_to_amim(ggm, AMI_MATRIX);

}else{
	
	// graph.mpi_print_ggm(ggm, mpi_rank);
	
// graph.mpi_print_ggm(ggm, mpi_rank);
std::cout<<"ggm to gg_amim"<<std::endl;
graph.ggm_to_ggamim(ggm,GG_AMI_MATRIX,0);

} 


/* 
std::cout<<"Testing amim"<<std::endl;
for(int i=0; i<AMI_MATRIX.size(); i++){
	for(int j=0; j< AMI_MATRIX[i].size(); j++){
		// for(int k=0; k< AMI_MATRIX[i][j].size(); k++){
		std::cout<<i<<" "<<j<<" "<<" "<<AMI_MATRIX[i][j].ct_count_	<<" "<<AMI_MATRIX[i][j].sigma_ct_count_	<<std::endl;
		// }
	}	
} */

// exit(0);

	

/* std::cout<<"Dimensions of GG are "<<std::endl;
for(int i=0; i<GG_AMI_MATRIX.size(); i++){
std::cout<<i<<" "<<GG_AMI_MATRIX[i].size()<<" on rank "<< mpi_rank<<std::endl;	

for(int j=0; j< GG_AMI_MATRIX[i].size(); j++){
	
	// std::cout<<"on "<<j<< " with size "<< GG_AMI_MATRIX[i][j].size()<<std::endl;
	
	// for(int k=0; k< GG_AMI_MATRIX[i][j].size();k++){
		
		// std::cout<<"last "<<k<<std::endl;
		
	// }
	
}
	
} */

std::cout<<"Completed graph content on rank "<< mpi_rank<<" of "<< comm_size<<std::endl;

std::vector< std::vector< std::vector< std::vector<double> > > > result_matrix, g_result_matrix;
// resize vectors
double last, g_last;

result_matrix.resize(max+1);
g_result_matrix.resize(max+1);


int ord_size;
if(!bool(global_use_groups)){
ord_size=AMI_MATRIX.size();
}else{
ord_size=GG_AMI_MATRIX.size();	
}

for(int ord=0; ord<ord_size; ord++){
	
	std::cout<<ord<<std::endl;

int group_size;
if(!bool(global_use_groups)){
group_size=AMI_MATRIX[ord].size();
}else{
group_size=GG_AMI_MATRIX[ord].size();	
}
	
	// std::cout<<"AMI_MATRIX of ord and size "<<ord<<" "<<AMI_MATRIX[ord].size() <<std::endl;
result_matrix[ord].resize(group_size);
g_result_matrix[ord].resize(group_size);
	for(int i=0; i< result_matrix[ord].size(); i++){
	result_matrix[ord][i].resize(extern_list.size());
	g_result_matrix[ord][i].resize(extern_list.size());
	for (int j=0; j< result_matrix[ord][i].size(); j++){
	result_matrix[ord][i][j].resize(4,0);	
	g_result_matrix[ord][i][j].resize(4,0);	
	}
		
	}
	
}
// std::cout<<"result_matrix[0].size is "<<result_matrix[0].size()<<std::endl;

std::cout<<"Finished setting up result matrix on rank "<< mpi_rank<<" of "<< comm_size<<std::endl;




// create a linearized index vector
std::vector<int> index(4);
std::vector< std::vector<int>> index_vec;
int count=0;

// for (int ord=0; ord< max+1; ord++){
for (int ord=0; ord< max; ord++){
int ami_size;

if(!bool(global_use_groups)){
ami_size=AMI_MATRIX[ord].size();
 }else{
ami_size=GG_AMI_MATRIX[ord].size();	
// std::cout<<"GG_AMI_MATRIX is size "<<ami_size<<" on rank "<<mpi_rank <<" at order "<<ord<<std::endl;
} 	

int extern_size=extern_list.size();
	for (int i=0; i< ami_size; i++){
	// #pragma omp parallel for 
		for (int j=0; j<extern_size; j++){
			// std::cout<<" ord is "<<ord<<std::endl;
		index[0]=ord;
		index[1]=i;
		index[2]=j;
		index[3]=count;
		count++;
		

int eff_ord;
if(!bool(global_use_groups)){
	eff_ord=ord+1+AMI_MATRIX[ord][i].ct_count_+2*AMI_MATRIX[ord][i].sigma_ct_count_;
}else{
	eff_ord=ord+1+GG_AMI_MATRIX[ord][i][0].ct_count_+2*GG_AMI_MATRIX[ord][i][0].sigma_ct_count_;// assumes groups contain same values 
}

if(eff_ord<= max){ index_vec.push_back(index);}		
		
		// index_vec.push_back(index);
			
			
		}
	}
}


// std::cout<<"Printing linearized index "<<std::endl;
// for(int i=0; i<index_vec.size(); i++){
	
// std::cout<<"On rank "<<mpi_rank<<" with i="<<i<<" and entries "<< index_vec[i][0]<<" "<<index_vec[i][1]<<" "<<index_vec[i][2]<<" "<<index_vec[i][3]<<" "<<std::endl;	
	
	
// }

// MPI_Finalize();
// return 0;

std::cout<<"Starting loop on rank "<< mpi_rank<<std::endl;


for(int i=0; i<index_vec.size(); i++){
// std::cout<<"Now on i"<<i<<std::endl;
// std::cout<<mpi_rank<<" "<<comm_size<<" "<<i%comm_size<<std::endl;

			
if (i%comm_size != mpi_rank) continue;			
std::cout<<"This is rank "<<mpi_rank <<" on i="<<i<<std::endl;
		//printf("ord = %d, i = %d, j= %d, threadId = %d \n",ord, i, j, omp_get_thread_num());
int ord=index_vec[i][0];
int num=index_vec[i][1];
int extvar=index_vec[i][2];

	std::cout<<"This is rank "<< mpi_rank <<" working on ord="<<ord<<", num="<< num <<", and extvar="<<extvar<<" using seed "<< intseed<<std::endl;
				
		// std::stringstream state_filename;
		// state_filename << ord+1;
		// state_filename<<AMI_MATRIX[ord][num].ct_count_<<AMI_MATRIX[ord][num].sigma_ct_count_;
		// state_filename << "_"<<extvar<<".state";
	
// if(ord!=2){continue;}
// if(num!=6){continue;}
	
		int nsamples=maxeval;
		
			if(!bool(global_use_groups)){
			std::vector<double> this_result;
			integrator(AMI_MATRIX[ord][num], extern_list[extvar], nsamples, intseed, ord+1,AMI_MATRIX[ord][num].ct_count_ ,AMI_MATRIX[ord][num].sigma_ct_count_, num, extvar, this_result);
		  result_matrix[ord][num][extvar]=this_result;
			}else{
				
				std::vector<double> this_result;
			group_integrator(GG_AMI_MATRIX[ord][num], extern_list[extvar], nsamples, intseed, ord+1,GG_AMI_MATRIX[ord][num][0].ct_count_ ,GG_AMI_MATRIX[ord][num][0].sigma_ct_count_, num, extvar, this_result);
		  result_matrix[ord][num][extvar]=this_result;
				
				
			}
		
		
} 

std::cout<<"This is rank "<< mpi_rank<<" waiting to reduce "<<std::endl;

MPI_Barrier(MPI_COMM_WORLD);
 

for (int ord=0; ord< max+1; ord++){
	for(int i=0; i< g_result_matrix[ord].size(); i++){
		for(int j=0; j< g_result_matrix[ord][i].size(); j++){
			// MPI_Reduce(&last, &g_last, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			for(int m=0; m<4; m++){
				// std::cout<<"Reducing "<<ord<< " "<<i<<" "<< j<<" "<<m<<std::endl;
		MPI_Reduce(&result_matrix[ord][i][j][m], &g_result_matrix[ord][i][j][m], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			}
		}
	}
}
 



if(mpi_rank==0){
	std::cout<<"This is master rank "<<mpi_rank<<" finalizing result."<<std::endl;
std::ofstream file;
file.open("outfile.dat");
file<<"# order num ext_index eff_ord ct_count sigma_ct_count ext_freq Im_ext_freq Beta reMu imMu kx ky kz Re err_Re Im err_Im"<<std::endl;

for (int ord=0; ord< max+1; ord++){
	std::cout<<"Writing order "<<ord<<std::endl;
	for(int i=0; i< g_result_matrix[ord].size(); i++){
		// std::cout<<"Writing graph "<<i<<std::endl;
		
		for(int j=0; j< g_result_matrix[ord][i].size(); j++){
			// std::cout<<"Writing extern "<<j<<std::endl;

file<< ord<<" "<<i<<" "<<j<<" ";

if(!global_use_groups){
file<<ord+1<<" "<<AMI_MATRIX[ord][i].ct_count_<<" "<<AMI_MATRIX[ord][i].sigma_ct_count_<<" ";
}else{
file<< ord+1<<" "<<GG_AMI_MATRIX[ord][i][0].ct_count_	<<" "<<GG_AMI_MATRIX[ord][i][0].sigma_ct_count_	<<" ";
	
}

file<< extern_list[j].external_freq_[0].real()<<" "<<extern_list[j].external_freq_[0].imag()<<" ";
file<<extern_list[j].BETA_<<" "<< extern_list[j].MU_.real()<<" "<<extern_list[j].MU_.imag()<<" ";
file<<extern_list[j].external_k_list_[0][0]<<" "<<extern_list[j].external_k_list_[0][1]<<" "<<extern_list[j].external_k_list_[0][2]<<" ";
file<< g_result_matrix[ord][i][j][0]<<" "<<g_result_matrix[ord][i][j][1] <<" "<<g_result_matrix[ord][i][j][2]<<" "<<g_result_matrix[ord][i][j][3]<<std::endl;

		}

	}
}
std::cout<<"Closing file "<<std::endl;

file.close();
std::cout<<"File Closed on Master "<<std::endl;
} 
 // exit(0);



MPI_Finalize();

return 0;
}


#define NCOMP 2
#define NVEC 1
// #define EPSREL 1e-5
// #define EPSABS 1e-10
// #define VERBOSE 0
#define FLAG 16
#define LAST 4
#define SEED 0
// #define MINEVAL 2000
//#define NSTART 10000
//#define NINCREASE 500
//#define NBATCH 5000
#define GRIDNO 0

// #define STATEFILE NULL
#define SPIN NULL

#define NNEW 500
#define NMIN 2
#define FLATNESS 1

#define KEY1 -2470 //47//47
#define KEY2 -2100 // 1
#define KEY3 1
#define MAXPASS 5
#define BORDER 0.001
#define MAXCHISQ 10// 10
#define MINDEVIATION .0125
#define NGIVEN 0

#define NEXTRA 0

#define KEY 0

// #define PseudoRandom True



void group_integrator(NewAmiCalc::solution_set_vec_t &ssamiv, NewAmiCalc::ext_vars &ext, int nsamples, int seed, int eff_ord, int ct, int sigct,int num, int extnum, std::vector<double> &out){
	
	
	out.clear();

	
	std::stringstream state_filename;
		state_filename << "statedir/"<<eff_ord;
		state_filename<<ct<<sigct<<"_n_"<<num<<"_"<<extnum<<".state";
	
	const std::string tmp=state_filename.str();
	const char* vegas_state=tmp.c_str();	
	

 std::cout<<"Statefile in integrator  is "<<vegas_state<<std::endl;
	
	
	std::vector<double> output;
// std::complex<double> output;
std::cout<<ssamiv.size()<<std::endl;


	
double EPSREL=global_epsrel;
double EPSABS=global_epsabs;

int NDIM=ext.KDIM_*ssamiv[0].ami_parms_.N_INT_;
std::cout<<"NDIM is "<< NDIM<<std::endl;
int MINEVAL=50000*NDIM;
int NBATCH=500*NDIM;
int NINCREASE=0;//250*NDIM;

int MAXEVAL=nsamples*NDIM;
int NSTART=global_nstart*MAXEVAL;
// int NSTART=MAXEVAL;
int LDXGIVEN=NDIM;


int comp, nregions, neval, fail;
double integral[NCOMP], error[NCOMP], prob[NCOMP];


evaluation_vector_set eval_v(ssamiv,ext);




//////////////////////////
// double complex saved;

if(global_integral==0){
  printf("\n------------------- Divonne test -------------------\n");

  Divonne(NDIM, NCOMP, group_ami_integrand, (void*)&eval_v, NVEC,
    EPSREL, EPSABS, FLAG, seed,
    MINEVAL, MAXEVAL, KEY1, KEY2, KEY3, MAXPASS,
    BORDER, MAXCHISQ, MINDEVIATION,
    NGIVEN, LDXGIVEN, NULL, NEXTRA, NULL,
    vegas_state, SPIN,
    &nregions, &neval, &fail, integral, error, prob);
	
	

  printf("DIVONNE RESULT:\tnregions %d\tneval %d\tfail %d\n",  nregions, neval, fail);
  for( comp = 0; comp < NCOMP; ++comp ){
    printf("DIVONNE RESULT:\t%.6e +- %.6e\tp = %.3f\n",      (double)integral[comp], (double)error[comp], (double)prob[comp]);
	  
	  out.push_back((double)integral[comp]);
	  out.push_back((double)error[comp]);
  }
}
  
if(global_integral==1){  
   printf("-------------------- Vegas test --------------------\n");

  Vegas(NDIM, NCOMP, group_ami_integrand, (void*)&eval_v, NVEC,
    EPSREL, EPSABS, FLAG, seed,
    MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
    GRIDNO, vegas_state, SPIN,
    &neval, &fail, integral, error, prob);

  printf("VEGAS RESULT:\tneval %d\tfail %d\n",
    neval, fail);
  for( comp = 0; comp < NCOMP; ++comp ){
    printf("VEGAS RESULT:\t%.6e +- %.6e\tp = %.3f\n",
      (double)integral[comp], (double)error[comp], (double)prob[comp]);
	out.push_back((double)integral[comp]);
	  out.push_back((double)error[comp]);
  }
}

 
 if(global_integral==2){
 printf("\n-------------------- Cuhre test --------------------\n");

  Cuhre(NDIM, NCOMP,  group_ami_integrand, (void*)&eval_v, NVEC,
    EPSREL, EPSABS, FLAG,
    MINEVAL, MAXEVAL, KEY,
    vegas_state, SPIN,
    &nregions, &neval, &fail, integral, error, prob);

  printf("CUHRE RESULT:\tnregions %d\tneval %d\tfail %d\n",
    nregions, neval, fail);
  for( comp = 0; comp < NCOMP; ++comp ){
    printf("CUHRE RESULT:\t%.6e +- %.6e\tp = %.3f\n",
      (double)integral[comp], (double)error[comp], (double)prob[comp]);


	out.push_back((double)integral[comp]);
	  out.push_back((double)error[comp]);

  }
}




// return output ;///M_PI;
	
	
}

void integrator(NewAmiCalc::solution_set &sol, NewAmiCalc::ext_vars &ext, int nsamples, int seed, int eff_ord, int ct, int sigct, int num, int extnum, std::vector<double> &out){
	
	out.clear();

	
	std::stringstream state_filename;
		state_filename << "statedir/"<<eff_ord;
		state_filename<<ct<<sigct<<"_n_"<<num<<"_"<<extnum<<".state";
	
	const std::string tmp=state_filename.str();
	const char* vegas_state=tmp.c_str();	
	

 std::cout<<"Statefile in integrator  is "<<vegas_state<<std::endl;
 // exit(0);

// seed=0;


	
double EPSREL=global_epsrel;
double EPSABS=global_epsabs;
	
	std::vector<double> output;
// std::complex<double> output;	
int NDIM=ext.KDIM_*sol.ami_parms_.N_INT_;
std::cout<<"NDIM is "<< NDIM<<std::endl;
int MINEVAL=50000*NDIM;

//2000*NDIM;
int NBATCH=500*NDIM;
// int NINCREASE=250*NDIM;
int NINCREASE=0;//250*NDIM;

// int KEY1=1000*NDIM;
// int KEY2=200*NDIM;
// int NCOMP=2;
//#define USERDATA NULL
// int NVEC=1;
// double EPSREL=1e-4;//5e-5;
// double EPSABS=1e-7;//5e-5;
// int VERBOSE=0;
// int LAST=0;//4
// int SEED=0;//(unsigned int)seedvar;//0
// int MINEVAL=000;
int MAXEVAL=nsamples*NDIM;//200000;
int NSTART=global_nstart*MAXEVAL;
// int NSTART=MAXEVAL;
// int NSTART=1000;
// int NINCREASE=500;
// int NBATCH=1000;
// int GRIDNO=-1;
// char* STATEFILE=NULL;

// int NNEW=1000;
// double FLATNESS=1.0/2.0;//25.

// int KEY1=3000;
// int KEY2=1800;
// int KEY3=1;
// int MAXPASS=10;//5
// double BORDER=1e-6;
// double MAXCHISQ=4;
// double MINDEVIATION=.52500;
// int NGIVEN=0;
int LDXGIVEN=NDIM;
// int NEXTRA= 0;

// int KEY=0;  

// #define STATEFILE NULL
// #define SPIN NULL

int comp, nregions, neval, fail;
double integral[NCOMP], error[NCOMP], prob[NCOMP];

NewAmiCalc::evaluation_set eval(sol,ext);




//////////////////////////
// double complex saved;

if(global_integral==0){
  printf("\n------------------- Divonne test -------------------\n");

  Divonne(NDIM, NCOMP, ami_integrand, (void*)&eval, NVEC,
    EPSREL, EPSABS, FLAG, seed,
    MINEVAL, MAXEVAL, KEY1, KEY2, KEY3, MAXPASS,
    BORDER, MAXCHISQ, MINDEVIATION,
    NGIVEN, LDXGIVEN, NULL, NEXTRA, NULL,
    vegas_state, SPIN,
    &nregions, &neval, &fail, integral, error, prob);
	
	

  printf("DIVONNE RESULT:\tnregions %d\tneval %d\tfail %d\n",  nregions, neval, fail);
  for( comp = 0; comp < NCOMP; ++comp ){
    printf("DIVONNE RESULT:\t%.6e +- %.6e\tp = %.3f\n",      (double)integral[comp], (double)error[comp], (double)prob[comp]);
	  
	  out.push_back((double)integral[comp]);
	  out.push_back((double)error[comp]);
  }
}
  
if(global_integral==1){  
   printf("-------------------- Vegas test --------------------\n");

  Vegas(NDIM, NCOMP, ami_integrand, (void*)&eval, NVEC,
    EPSREL, EPSABS, FLAG, seed,
    MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
    GRIDNO, vegas_state, SPIN,
    &neval, &fail, integral, error, prob);

  printf("VEGAS RESULT:\tneval %d\tfail %d\n",
    neval, fail);
  for( comp = 0; comp < NCOMP; ++comp ){
    printf("VEGAS RESULT:\t%.6e +- %.6e\tp = %.3f\n",
      (double)integral[comp], (double)error[comp], (double)prob[comp]);
	out.push_back((double)integral[comp]);
	  out.push_back((double)error[comp]);
  }
}

 
 if(global_integral==2){
 printf("\n-------------------- Cuhre test --------------------\n");

  Cuhre(NDIM, NCOMP,  ami_integrand, (void*)&eval, NVEC,
    EPSREL, EPSABS, FLAG,
    MINEVAL, MAXEVAL, KEY,
    vegas_state, SPIN,
    &nregions, &neval, &fail, integral, error, prob);

  printf("CUHRE RESULT:\tnregions %d\tneval %d\tfail %d\n",
    nregions, neval, fail);
  for( comp = 0; comp < NCOMP; ++comp ){
    printf("CUHRE RESULT:\t%.6e +- %.6e\tp = %.3f\n",
      (double)integral[comp], (double)error[comp], (double)prob[comp]);


	out.push_back((double)integral[comp]);
	  out.push_back((double)error[comp]);

  }
}


// std::cout<<"Exiting integrator with output" <<std::endl;

// return output ;///M_PI;
	
	
}


static int group_ami_integrand(const int *ndim, const double xx[],
  const int *ncomp, double ff[], void *userdata){
		
double kmin,kmax;
double ktemp;
		
evaluation_vector_set *eval;	 
std::vector<NewAmiCalc::solution_set> sol_v;
NewAmiCalc::ext_vars ext ;

eval=(evaluation_vector_set *) userdata;

sol_v=eval->sol_;
ext=eval->ext_vars_;

NewAmiCalc::k_vector_t k;
NewAmiCalc::k_vect_list_t k_list;
		
		
	
int kdim=ext.KDIM_;
if (kdim !=3){
	throw std::runtime_error("Only works for kdim=3");
}

int gsize=sol_v[0].R0_.size();

double magk=std::sqrt(std::pow(ext.external_k_list_[0][0],2)+std::pow(ext.external_k_list_[0][1],2)+std::pow(ext.external_k_list_[0][2],2));
double maxf=std::abs(ext.external_freq_[0]);

kmin=0; 
kmax=kc;	
		
		

double kx,ky,kz;
double mag,theta,phi;

double Jk=1;
int j=0;
do{


mag=xx[j]*kmax;
theta=xx[j+1]*M_PI;
phi=xx[j+2]*2.0*M_PI;
Jk=Jk*std::pow(mag,2)*std::sin(theta)*kmax*2.0*std::pow(M_PI,2);

kx=mag*std::sin(theta)*std::cos(phi);
ky=mag*std::sin(theta)*std::sin(phi);
kz=mag*std::cos(theta);

k.push_back(kx);
k.push_back(ky);
k.push_back(kz);

k_list.push_back(k);
k.clear();	
	
	
j=j+kdim;
}while(j< *ndim);

	
NewAmiCalc::internal_state state(k_list.size(), k_list[0].size());


state.t_list_.clear();
state.t_list_.resize(gsize, 1);	
state.mink_=kmin;
state.maxk_=kmax;

if(global_hf){
state.disp_=AmiBase::hf;
}else{
	state.disp_=AmiBase::fp;
}

state.internal_k_list_=k_list;


std::complex<double> result_sum(0,0);

for(int s=0; s< sol_v.size(); s++){
	
NewAmiCalc::solution_set sol;


sol=sol_v[s];

AmiBase::ami_vars vars=graph.ami.construct_ami_vars(sol.R0_,sol.prefactor_, state, ext);


double norm= 1.0/std::pow(2*M_PI,state.dim_*state.internal_k_list_.size());
double V=get_V(state, sol, ext,  global_lambda); //1.44); // lambda should not be hard-coded
double sq=std::pow(2, sol.loops_);
double ct=std::pow(-global_lambda/4.0/M_PI/global_esquare, sol.ct_count_); 

std::complex<double> ami_result=graph.ami.amibase.evaluate(sol.ami_parms_, sol.R_, sol.P_, sol.S_,  vars, sol.Unique, sol.Rref, sol.Eval_list );

std::complex<double> calc_result=Jk*norm*V*sq*ct*ami_result;


if(abs(calc_result.real())>1e8 || abs(calc_result.imag())>1e8){ 

calc_result=(0,0);

}


if(std::isnan(calc_result.real()) || std::isnan(calc_result.imag())){
		std::cout<<"Nan ignored"<<std::endl;
}

result_sum+=calc_result;


}

ff[0]=result_sum.real();//*factor;
ff[1]=result_sum.imag();//*factor;
	  
	  
	return 0;  		
		
	}


static int ami_integrand(const int *ndim, const double xx[],
  const int *ncomp, double ff[], void *userdata){
	  
// std::cout<<"In integrand function"<<std::endl;		
		
// AmiCalc ami;	  

double kmin,kmax;
double ktemp;

// hard-coded integration range
// kmin=-1.00+.001;
// kmax=1.0+.001;

NewAmiCalc::evaluation_set *eval;	 
NewAmiCalc::solution_set sol;
NewAmiCalc::ext_vars ext ;

eval=(NewAmiCalc::evaluation_set *) userdata;

sol=eval->sol_;
ext=eval->ext_vars_;

NewAmiCalc::k_vector_t k;
NewAmiCalc::k_vect_list_t k_list;

int kdim=ext.KDIM_;
if (kdim !=3){
	throw std::runtime_error("Only works for kdim=3");
}

int gsize=sol.R0_.size();

double magk=std::sqrt(std::pow(ext.external_k_list_[0][0],2)+std::pow(ext.external_k_list_[0][1],2)+std::pow(ext.external_k_list_[0][2],2));
double maxf=std::abs(ext.external_freq_[0]);

kmin=0; //-2.500-magk*5.0 - maxf+.001;//-6.0+.001;//-1.00-magk*5.0 - maxf+.001;
kmax=kc;//3.001; //2.50+magk*5.0+maxf+.001;//6.0+.001;//1.0+magk*5.0+maxf+.001;

// ext.external_k_list_[0][0]	<<" "<< ext.external_k_list_[0][1]<<" "<< ext.external_k_list_[0][2]<<std::endl


// kmax=ext.maxk_;
// kmin=ext.mink_;
// std::cout<<"Kdim is "<<kdim<<std::endl;


// first, translate the array of xx[i] into kx,ky sets 
// tehn, construct a state 

// new function - takes the array of xx[] and returns the state

double kx,ky,kz;
double mag,theta,phi;

double Jk=1;
int j=0;
do{

mag=xx[j]*kmax;
theta=xx[j+1]*M_PI;
phi=xx[j+2]*2.0*M_PI;
Jk=Jk*std::pow(mag,2)*std::sin(theta)*kmax*2.0*std::pow(M_PI,2);

kx=mag*std::sin(theta)*std::cos(phi);
ky=mag*std::sin(theta)*std::sin(phi);
kz=mag*std::cos(theta);

k.push_back(kx);
k.push_back(ky);
k.push_back(kz);

// if(mag>2){
// std::cout<<xx[j]<<" "<<xx[j+1]<<" "<<xx[j+2]<<std::endl;

// std::cout<<mag<<" "<<theta<<" "<<phi<<std::endl;
// std::cout<<"kxkykz on j with jk "<< kx<<" "<< ky<<" "<<kz<<" "<<j<<" "<<Jk<<std::endl;	
// }


k_list.push_back(k);
k.clear();	
	
	
	
// for(int i=j; i<kdim+j; i++){
	// ktemp=kmin+(kmax-kmin)*xx[i];
	// std::cout<<"ktemp and xx are "<<ktemp<<" "<<xx[i]<<" "<<i<<std::endl;
	// k.push_back(ktemp);
	
// }
// k_list.push_back(k);
// k.clear();
j=j+kdim;
}while(j< *ndim);

// double Jk=kmax-kmin;

// std::cout<<k_list.size()<<std::endl;
// std::cout<<k_list[0].size()<<std::endl;

// AmiCalc::ami_vars_list vars_list;
// graph.ami.construct_ami_vars_list(AMI_MATRIX[ord][num].R0_, graph.current_state, extern_list, vars_list);

NewAmiCalc::internal_state state(k_list.size(), k_list[0].size());
// reset all hoppings to 1
state.t_list_.clear();
state.t_list_.resize(gsize, 1);	
state.mink_=kmin;
state.maxk_=kmax;

if(global_hf){
state.disp_=AmiBase::hf;
}else{
	state.disp_=AmiBase::fp;
}
// std::cout<< "hopping is size "<<state.t_list_.size()<<std::endl;


state.internal_k_list_=k_list;
// state.prefactor_=sol.prefactor_; // TODO: depricated

// graph.print_fullstate(state);
// std::cout<<"Prefactor is "<<sol.prefactor_<<std::endl;
AmiBase::ami_vars vars=graph.ami.construct_ami_vars(sol.R0_,sol.prefactor_, state, ext);

// double norm=std::pow(state.maxk_-state.mink_, state.dim_*state.internal_k_list_.size())/std::pow(2*M_PI,state.dim_*state.internal_k_list_.size());
double norm= 1.0/std::pow(2*M_PI,state.dim_*state.internal_k_list_.size());
double V=get_V(state, sol, ext, global_lambda);// 1.44); // lambda should not be hard-coded
double sq=std::pow(2, sol.loops_);
double ct=std::pow(-global_lambda/4.0/M_PI/global_esquare, sol.ct_count_); 

// std::cout<<"norm V and sq "<< norm<<" "<< V<<" "<<sq<<std::endl;
// std::cout<<"norm V sq ct "<< norm<<" "<< V<<" "<<sq<<" "<<ct<<std::endl;

// std::cout<<"EVAL-------"<<std::endl;
std::complex<double> ami_result=graph.ami.amibase.evaluate(sol.ami_parms_, sol.R_, sol.P_, sol.S_,  vars, sol.Unique, sol.Rref, sol.Eval_list );

std::complex<double> calc_result=Jk*norm*V*sq*ct*ami_result;
//ami_result*mag*mag*sin(theta)*2.0/std::pow(2.0*M_PI,3)*2.0*M_PI*M_PI*kc*2.0*M_PI*M_PI;;//Jk*norm*V*sq*ct*ami_result;
// mag*mag*sin(theta)*4.0/std::pow(2.0*M_PI,3)*2.0*M_PI*M_PI*kc*2.0*M_PI*M_PI;

// std::cout<<"TESTING-------"<<std::endl;
// double kf=1.919158;
// std::complex<double> this_E= mag*mag-kf*kf;  //graph.ami.hf_energy(mag);
// double beta=25.0/kf/kf;
// std::complex<double> test_result=beta*beta*(-std::exp(-beta*this_E)+1.0)/( 1.0+std::exp(beta*this_E) )/( 1.0+std::exp(-beta*this_E) )/( 1.0+std::exp(-beta*this_E) )/2.0;
// std::complex<double> test_result=graph.ami.amibase.fermi_bose(2,1.0,25.0,this_E);
// std::complex<double> test_100_result=graph.ami.amibase.fermi_bose(1,1.0,25.0,this_E);

// std::cout<<"Integrand gave "<<ami_result<<std::endl;
// std::cout<<"Test function gave "<<test_result<<std::endl;
// std::cout<<"jk norm V sq ct "<<Jk<<" "<< norm<<" "<< V<<" "<<sq<<" "<<ct<<std::endl;
// std::cout<<"Actual output should be "<<calc_result<<std::endl;
// std::cout<<"Actual Test output should be "<<-0.5*Jk*norm*V*sq*ct*test_result<<std::endl;

if(abs(calc_result.real())>1e8 || abs(calc_result.imag())>1e8){ 
// std::cout<<"overflow at "<<std::endl;
// graph.print_fullstate(state);
// throw std::runtime_error("Overflow detected, exiting");
// graph.print_fullstate(state);
// calc_result=(0,0);
// std::cout<<"Cutoff triggered"<<std::endl;
}

// double factor=std::pow(Jk, *ndim);
// double factor=10000;
if(std::isnan(calc_result.real()) || std::isnan(calc_result.imag())){
	
	// std::cout<<"nan detected"<<std::endl;
	// graph.print_fullstate(state);
	// throw std::runtime_error("Nan detected");
	calc_result=(0,0);
	// std::cout<<"Nan ignored"<<std::endl;
}

ff[0]=calc_result.real();//(-Jk*norm*V*sq*ct*test_100_result).real();//calc_result.real();//*factor;
ff[1]=calc_result.imag();
// (4.0)*(-1.0)*(mag*mag*std::sin(theta))*(2.0*M_PI*M_PI*kc)/(kf/(2.0*M_PI*M_PI))/std::pow(2.0*M_PI,3)*test_result.real();
//(mag*mag*sin(theta)*test_result).real()*4.0/std::pow(2.0*M_PI,3)*2.0*M_PI*M_PI*kc*2.0*M_PI*M_PI/1.919158/1.919158;//calc_result.imag();//(-0.5*Jk*norm*V*sq*ct*test_result).real();//calc_result.imag();//*factor;

// std::cout<<"Exiting integrand "<<ff[0]<< " "<<ff[1]<<std::endl;	  
	  
	return 0;  
  }



double get_V(NewAmiCalc::internal_state &state, NewAmiCalc::solution_set &AMI, NewAmiCalc::ext_vars &external, double lambda){
	
double result=1;	

std::vector<double> q_list;

// std::cout<<"Getting q's on rank "<<mpi_rank<<std::endl;

get_q_list(q_list,AMI, state, external);	
	
for(int m=0; m< q_list.size();m++){

// std::cout<<"q_"<<m<<" = "<<q_list[m]<<std::endl;
result=result/(q_list[m] +lambda); 
}	
// result=result*pow(8.0*M_PI, q_list.size());
// 4 pi or 8pi depending on dispersion 
// double rs= 2.0;// hard-coded rs   parameters["RS"];
double esquare=global_rs*pow(32.0/(9*M_PI), 1.0/3.0);
result=result*pow(4.0*M_PI*esquare, q_list.size());
	

return result;	
}


void get_q_list(std::vector<double> &qs, NewAmiCalc::solution_set &AMI, NewAmiCalc::internal_state &state, NewAmiCalc::ext_vars &external){
	
	qs.clear();
	
	
NewAmiCalc::k_vect_list_t k_list;



k_list=state.internal_k_list_;
k_list.push_back(external.external_k_list_[0]);
// std::cout<<"Printing k list"<<std::endl;
// for(int i=0; i< k_list.size(); i++){
// for(int m=0; m< k_list[i].size(); m++){	
	// std::cout<<k_list[i][m]<<std::endl;
// }
// }
// std::cout<<"Printing k list"<<std::endl;


// std::cout<<"Bose alphas size is "<< AMI.bose_alphas_.size()<<std::endl;

for(int i=0; i< AMI.bose_alphas_.size(); i++){

std::vector<double> q_vec(state.dim_,0.0);	

double sum=0;
if(AMI.bose_alphas_[i].size()!=k_list.size()){
	throw std::runtime_error("alpha size and momentum list sizes don't match exiting");
}
for(int j=0; j< AMI.bose_alphas_[i].size(); j++){

for(int m=0; m<state.dim_; m++){
	
	// std::cout<<"k and alpha "<< k_list[j][m]<<" "<<AMI.bose_alphas_[i][j]<<std::endl;

q_vec[m]+=k_list[j][m]*AMI.bose_alphas_[i][j];

}
}



double temp=0;
for(int m=0; m<state.dim_; m++){


// std::cout<<q_vec[m]<<" ";

temp+=q_vec[m]*q_vec[m];
}
// std::cout<<std::endl;

qs.push_back(temp);

}

if(qs.size()!= AMI.bose_alphas_.size()){
	throw std::runtime_error("Didn't get all the q values");
}
	
}


void zero_external(AmiGraph::graph_t &g, int index){
	
boost::graph_traits<AmiGraph::graph_t>::edge_iterator ei, ei_end;

for (boost::tie(ei,ei_end)=edges(g); ei!=ei_end; ++ei){
	

g[*ei].g_struct_.alpha_[index]=0;	
	
	
}
	
graph.fix_epsilons(g);	
	
}




void load_hf(){
	
	
/////////////
					// read lambda rs etc
						
					std::ifstream infile_stream;
					infile_stream.open("input_hf.dat");

					if(infile_stream.fail()) // checks to see if file opended 
							{ 
							throw std::runtime_error("Could not open input_hf.dat file");
							} 	
						
					std::string line;
					std::getline(infile_stream,line);	
					
					std::stringstream ss;
					std::getline(infile_stream,line);
					
					ss<< line;
					ss>>global_rs;
					std::getline(infile_stream,line);
					
					std::stringstream ss2;
					ss2<<line;
					double kappa;
					
					ss2>>kappa;
					global_lambda=kappa*kappa;
					global_esquare=global_rs*pow(32.0/(9*M_PI), 1.0/3.0);

// std::cout<<kappa<<" "<<test<<" "<<test2<<std::endl;
std::cout<<"Running for lambda="<<global_lambda<<" rs="<<global_rs<<" esquare="<<global_esquare<<std::endl;
 

///////////// 
	
	
	
}


void load_settings(){
	
	
/////////////
					// read lambda rs etc
						
					std::ifstream infile_stream;
					infile_stream.open("parameters.param");

					if(infile_stream.fail()) // checks to see if file opended 
							{ 
							throw std::runtime_error("Could not open parameters.param file");
							} 	
						
					std::string line;
					std::stringstream ss;
					std::getline(infile_stream,line);
					
					std::string junk;
					
					ss<< line;
					ss>>junk>>global_min;
					std::stringstream().swap(ss);
					
					std::getline(infile_stream,line);
					
					ss<< line;
					ss>>junk>>global_max;
					std::stringstream().swap(ss);
					
					std::getline(infile_stream,line);
					
					ss<< line;
					ss>>junk>>global_maxeval;
					std::stringstream().swap(ss);
					
					std::getline(infile_stream,line);
					
					ss<< line;
					ss>>junk>>global_seed;
					std::stringstream().swap(ss);
					
					std::getline(infile_stream,line);
					
					ss<< line;
					ss>>junk>>global_epsabs;
					std::stringstream().swap(ss);
					
					std::getline(infile_stream,line);
					
					ss<< line;
					ss>>junk>>global_epsrel;
					std::stringstream().swap(ss);
					
					std::getline(infile_stream,line);
					
					ss<< line;
					ss>>junk>>global_nstart;
					std::stringstream().swap(ss);
					
					std::getline(infile_stream,line);
					
					ss<< line;
					ss>>junk>>global_use_bare;
					std::stringstream().swap(ss);
					
					std::getline(infile_stream,line);
					
					ss<< line;
					ss>>junk>>kc;
					std::stringstream().swap(ss);
					
					std::getline(infile_stream,line);
					
					ss<< line;
					ss>>junk>>global_use_groups;
					std::stringstream().swap(ss);
					
					std::getline(infile_stream,line);
					
					ss<< line;
					ss>>junk>>ZERO_EXTERNAL_Q;
					std::stringstream().swap(ss);
					
					


std::cout<<"Running for parameters:"<<std::endl;
std::cout<<"min: "<<global_min<<std::endl;
std::cout<<"max: "<<global_max<<std::endl;
std::cout<<"maxeval: "<<global_maxeval<<std::endl;
std::cout<<"seed: "<<global_seed<<std::endl;
std::cout<<"epsabs: "<<global_epsabs<<std::endl;
std::cout<<"epsrel: "<<global_epsrel<<std::endl;
std::cout<<"nstart: "<<global_nstart<<std::endl;
std::cout<<"use_bare: "<<global_use_bare<<std::endl;
std::cout<<"kc: "<<kc<<std::endl;
std::cout<<"use_groups: "<<global_use_groups<<std::endl;
std::cout<<"ZERO_EXTERNAL_Q: "<< ZERO_EXTERNAL_Q <<std::endl;
std::cout<<std::endl;
 

///////////// 
	
	
	
}

