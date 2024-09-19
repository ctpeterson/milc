/*

Rewrite of update_onemass.c that is slightly more clear and 
up-to-date.

(Rewrite) author: Curtis Taylor Peterson

*/

// definitions files and prototypes
#define CONTROL
#include "ks_imp_includes.h"
#include "lattice_qdp.h"

#ifdef HAVE_QUDA
#include <quda_milc_interface.h>
#endif

#ifdef HAVE_QPHIX
#include "../include/generic_qphix.h"
#endif

#ifdef HAVE_QIO
#include <qio.h>
#include "../include/io_scidac.h"
#endif

#include "../generic/generic_includes.h"
#include "../include/openmp_defs.h"

#define special_alloc malloc
#define special_free free

// Input gauge field header
EXTERN gauge_header start_lat_hdr;

int main(int argc, char**argv){
    int meascount,traj_done;
    int prompt,s_iters,avs_iters,avbcorr_iters;
    double dtime, dclock();
  
    // Function for measuring plaquette
    void plaq(){
        double ss_plaquette,st_plaquette;
        su3_matrix *tempmat1,*links;

        tempmat1 = (su3_matrix *)special_alloc(sites_on_node*sizeof(su3_matrix));
        if(tempmat1 == NULL){
            printf("plaq: Can't malloc temporary\n");
            terminate(1);
        }
        links = create_G_from_site();
        d_plaquette(&ss_plaquette, &st_plaquette);
        #if (MILC_PRECISION == 1)
            node0_printf("PLAQ:\t%f\t%f\n",ss_plaquette,st_plaquette);
        #else
            node0_printf("PLAQ:\t%.16f\t%.16f\n",ss_plaquette,st_plaquette);
        #endif
        destroy_G(links);
        special_free(tempmat1);
        if (this_node == 0){fflush(stdout);}
    }

    // Function for measuring Polyakov loop
    void poly(){
        complex plp;
        plp = ploop();
        node0_printf("P_LOOP:\t%e\t%e\n",plp.real,plp.imag);
    }

    // Initialize simulation
    initialize_machine(&argc,&argv);
    if(remap_stdio_from_args(argc, argv) == 1)terminate(1);
    g_sync();
    prompt = setup();

    // Run HMC, measure, IO
    while(readin(prompt) == 0){
        // Warmup loop
        dtime = -dclock();
        for(traj_done=0; traj_done < warms; traj_done++){update();}
        node0_printf("WARMUPS COMPLETED\n"); fflush(stdout);

        // Production loop
        meascount = 0;
        avs_iters = avbcorr_iters = 0;
        for(traj_done=0; traj_done < trajecs; traj_done++){
            // HMC update
            s_iters = update();

            // Measurements
            if((traj_done % propinterval) == (propinterval - 1)){
                // Gauge measurements
                plaq(); poly();

                // Fermion measurements
                restore_fermion_links_from_site(fn_links, prec_pbp);
                #ifdef ONEMASS
                    f_meas_imp( 
                        npbp_reps_in, prec_pbp, F_OFFSET(phi),
                        F_OFFSET(xxx), mass, 0, fn_links
                    );
                #else
                    f_meas_imp( 
                        npbp_reps_in, prec_pbp, F_OFFSET(phi1),
                        F_OFFSET(xxx1), mass1, 0, fn_links
                    );
                    f_meas_imp( 
                        npbp_reps_in, prec_pbp, F_OFFSET(phi2),
                        F_OFFSET(xxx2), mass2, 0, fn_links
                    );
                #endif
                avs_iters += s_iters;
	            ++meascount; fflush(stdout);
            }
            // Cleanup & printout
            node0_printf("RUNNING COMPLETED\n"); fflush(stdout);
            if(meascount > 0){
                node0_printf("average cg iters for step= %e\n",
		        (double)avs_iters/meascount);
            }
            node0_printf("Time = %e seconds\n",dtime);
            node0_printf("total_iters = %d\n",total_iters);
        }
        // Additional cleanup & IO
        fflush(stdout);
        if(saveflag != FORGET){
            rephase(OFF);
            save_lattice( saveflag, savefile, stringLFN );
            rephase(ON);
        }
    }
}