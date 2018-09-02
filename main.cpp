/* 
 * File:   main.cpp
 * Author: Kapil Thakkar
 *
 */

#include <cstdlib>
#include <cstdio>
#include <ctime>
// #include <chrono>
#include <time.h>
#include "SessionOrganizer.h"

// #include<ctime>
using namespace std;

/*
 * 
 */
int main ( int argc, char** argv )
{
    // Parse the input.
    // auto wcts = std::chrono::system_clock::now();

	time_t mainstart;
	time(& mainstart);
	
    if ( argc < 3 )
    {   
        cout<<"Missing arguments\n";
        cout<<"Correct format : \n";
        cout << "./main <input_filename> <output_filename>";
        exit ( 0 );
    }
    string inputfilename ( argv[1] );
    
    // Initialize the conference organizer.
    SessionOrganizer *organizer  = new SessionOrganizer( inputfilename );

    // Organize the papers into tracks based on similarity.
    // organizer->organizePapers ( );


    // organizer->printSessionOrganiser ( argv[2]);

    // // Score the organization against the gold standard.
    // double score = organizer->scoreOrganization ( );
    // cout<< "score_initial:"<<score<<endl;


    // organizer->organizePapers2 ();


    // organizer->printSessionOrganiser ( argv[2]);

    // // Score the organization against the gold standard.
    // double score2 = organizer->scoreOrganization ( );
    // cout<< "score_final_without_random_restart :"<<score2<<endl;


    // organizer->randomState();
    // organizer->randomState();
    // organizer->printSessionOrganiser ( argv[2]);

    // random restarts r times.
    
    double score2 = 0;
    int r=1;
    // if(organizer->conference->getParallelTracks*organizer->conference->getSessionsInTrack*organizer->conference->getPapersInSession > 100) r=1;
    if(organizer->parallelTracks*organizer->sessionsInTrack*organizer->papersInSession > 24) r=1;
    else r = 2;

    for(int i =0;i<1;i++){
    	srand(time(0)+8*(i+9));
    	int seed = rand();
    	int f = rand();
    	seed =seed%f;
    	organizer->randomState(seed);
    	double score = organizer->scoreOrganization ();
	    // cout<< "score_initial:"<<score<<endl;
    	// std::chrono::duration<double> wctduration = (std::chrono::system_clock::now() - wcts);
    	// double t1 =wctduration.count();
    	time_t mainmid;
		time(& mainmid);
		double t1 = mainmid - mainstart;
		if(r==2)organizer->organizePapers2 (t1);
		organizer->organizePapers3 (t1);
	    double score_r = organizer->scoreOrganization ( );
	    cout<<"random score: "<<score_r<<endl;
	    if(score_r>score2) {
			organizer->printSessionOrganiser ( argv[2]);
			score2 = score_r;	    	
	    }
	    cout<< "score_final_with_random_restart :"<<score2<<endl;


    }

 //    double final_check = organizer->scoreOrganization ( );
	// cout<< "final check : "<< final_check<<endl;
    time_t endtime;
        time(& endtime);
        double tf = endtime - mainstart;
        cout<<"total time elspsed : "<< tf<<"\n";
    // std::chrono::duration<double> wctduration = (std::chrono::system_clock::now() - wcts);
	// std::cout << "Finished in " << wctduration.count() << " seconds [Wall Clock]" << std::endl;
    return 0;
}

