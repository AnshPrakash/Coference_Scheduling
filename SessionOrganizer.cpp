
/* 
 * File:   SessionOrganizer.cpp
 * Author: Kapil Thakkar
 * 
 */

#include "SessionOrganizer.h"
#include "Util.h"
#include <ctime>
// #include <chrono>

#include <cstdlib>

#include <random>




SessionOrganizer::SessionOrganizer ( )
{
    parallelTracks = 0;
    papersInSession = 0;
    sessionsInTrack = 0;
    processingTimeInMinutes = 0;
    tradeoffCoefficient = 1.0;
}

SessionOrganizer::SessionOrganizer ( string filename )
{
    readInInputFile ( filename );
    conference = new Conference ( parallelTracks, sessionsInTrack, papersInSession );
}

void SessionOrganizer::organizePapers ( )
{
    int paperCounter = 0;
    for ( int i = 0; i < conference->getSessionsInTrack ( ); i++ )
    {
        for ( int j = 0; j < conference->getParallelTracks ( ); j++ )
        {
            for ( int k = 0; k < conference->getPapersInSession ( ); k++ )
            {
                conference->setPaper ( j, i, k, paperCounter );
                paperCounter++;
            }
        }
    }
}


void SessionOrganizer::organizePapers2 (double t1)
{
    // int paperCounter = 0;
    int flag =0;int max_paper;int min_paper;
    time_t funcstart;
    time(& funcstart);
    do
    {   flag=0;
        int maxi=0,maxj=0,maxk=0;
        int mini=0,minj=0,mink=0;
        double score_initial = scoreOrganization();
        double score_max = score_initial; 
        cout<<"initial score in this state is this "<<score_max<<endl;

        time_t funcmid;
        time(& funcmid);
        if(funcmid-funcstart > processingTimeInMinutes*60.0 - processingTimeInMinutes*0.6 - t1 ){cout<<"time : "<<funcmid-funcstart+t1;   return;}

        for ( int i = 0; i < conference->getSessionsInTrack ( ); i++ )
        {
            for ( int j = 0; j < conference->getParallelTracks ( ); j++ )
            {
                for ( int k = 0; k < conference->getPapersInSession ( ); k++ )
                {   
                    for ( int pi = i; pi < conference->getSessionsInTrack ( ); pi++ )
                    {                
                        for ( int pj = 0; pj < conference->getParallelTracks ( ); pj++ )
                        {
                            for ( int pk = 0; pk < conference->getPapersInSession ( ); pk++ )
                            {        

                                if(pj==j&&pi==i)continue;
                                int paper_initial = conference->getTrack(j).getSession ( i ).getPaper ( k );
                                int paper_final = conference->getTrack(pj).getSession ( pi ).getPaper ( pk );
                                // conference->setPaper ( j, i, k, paper_final );
                                // conference->setPaper ( pj, pi, pk, paper_initial);

                                double sf = fast_score(score_initial,i,j,k,pi,pj,pk,tradeoffCoefficient);
                                // double sf=fast_score(conference);
                                // cout<<"score of this neighbour "<<sf<<endl;
                                if(score_max < sf) {
                                    // cout<<"max till now "<<sf<<endl;
                                    score_max = sf; maxj = pj;maxi = pi;maxk = pk;flag =1;minj =j;mini=i;mink =k;
                                    min_paper = paper_initial;max_paper = paper_final;
                                }  
                                // conference->setPaper ( j, i, k, paper_initial);
                                // conference->setPaper ( pj, pi, pk, paper_final);

                            }
                        }
                    }
                }
            }
        }
        int paper_max = conference->getTrack(maxj).getSession ( maxi ).getPaper ( maxk );
        if(flag == 1) {
            // cout<<"recursing on max neighbour "<<endl;
            conference->setPaper ( maxj, maxi, maxk, min_paper);
            conference->setPaper ( minj, mini, mink, max_paper);
            // return;
            // organizePapers2 ();
        }
        else{
            break;
        }
    } while (flag==1);    
}









void SessionOrganizer::organizePapers3 (double t1)
{
    // int paperCounter = 0;
    int flag =0;int max_paper;int min_paper;
    // std::chrono::duration<double> wctduration;
    // auto t2 = std::chrono::system_clock::now();
    time_t funcstart;
    time(& funcstart);

    do
    {   flag=0;
        int maxi=0,maxj=0,maxk=0;
        int mini=0,minj=0,mink=0;
        double score_initial = scoreOrganization();
        double score_max = score_initial; 
        // cout<<"initial score in this state is this "<<score_max<<endl;

        
        // wctduration = (std::chrono::system_clock::now() - t2 );
        // if(wctduration.count()>=60)return;
        time_t funcmid;
        time(& funcmid);
        // if(processingTimeInMinutes*60< 5)
        if(funcmid-funcstart > processingTimeInMinutes*60.0 - 5 - t1 ){cout<<"time : "<<funcmid-funcstart+t1;   return;}
    


        std::mt19937 rng;
        rng.seed(std::random_device()());
        std::uniform_int_distribution<std::mt19937::result_type> dist6(0,100000);
    




        int t_size=conference->getSessionsInTrack();
        int p_size=conference->getParallelTracks();
        int k_size=conference->getPapersInSession();
        int seed_c=0;
        for ( int i = 0; i < conference->getSessionsInTrack ( ); i++ )
        {
            for ( int j = 0; j < conference->getParallelTracks ( ); j++ )
            {
                for ( int k = 0; k < conference->getPapersInSession ( ); k++ )
                {   
                    for ( int pi = i; pi < conference->getSessionsInTrack ( ); pi++ )
                    {                
                        for ( int pj = 0; pj < conference->getParallelTracks ( ); pj++ )
                        {
                            for ( int pk = 0; pk < conference->getPapersInSession ( ); pk++ )
                            {        
                                srand(seed_c+pk);
                                int t1=dist6(rng)%t_size;
                                int t2=dist6(rng)%t_size;
                                int p1=dist6(rng)%p_size;
                                int p2=dist6(rng)%p_size;
                                int k1=dist6(rng)%k_size;
                                int k2=dist6(rng)%k_size;

                                // int t1=rand()%t_size;
                                // int t2=rand()%t_size;
                                // int p1=rand()%p_size;
                                // int p2=rand()%p_size;
                                // int k1=rand()%k_size;
                                // int k2=rand()%k_size;

                                

                                if(p2==p1&&t1==t2)continue;
                                int paper_initial = conference->getTrack(p1).getSession ( t1 ).getPaper ( k1 );
                                int paper_final = conference->getTrack(p2).getSession ( t2 ).getPaper ( k2 );

                                // int paper_initial = conference->getTrack(j).getSession ( i ).getPaper ( k );
                                // int paper_final = conference->getTrack(pj).getSession ( pi ).getPaper ( pk );
                                // conference->setPaper ( j, i, k, paper_final );
                                // conference->setPaper ( pj, pi, pk, paper_initial);
                                double sf = fast_score(score_initial,t1,p1,k1,t2,p2,k2,tradeoffCoefficient);

                                // double sf = fast_score(score_initial,i,j,k,pi,pj,pk,tradeoffCoefficient);
                                // double sf=fast_score(conference);
                                // cout<<"score of this neighbour "<<sf<<endl;
                                if(score_max < sf) {
                                    // cout<<"max till now "<<sf<<endl;
                                    // score_max = sf; maxj = pj;maxi = pi;maxk = pk;flag =1;minj =j;mini=i;mink =k;
                                    score_max = sf; maxj = p2;maxi = t2;maxk = k2;flag =1;minj =p1;mini=t1;mink =k1;

                                    min_paper = paper_initial;max_paper = paper_final;
                                    break;
                                }  
                                // conference->setPaper ( j, i, k, paper_initial);
                                // conference->setPaper ( pj, pi, pk, paper_final);


                            }
                            if(flag==1)break;
                        }
                        if(flag==1)break;
                    }
                    if(flag==1)break;
                }
                if(flag==1)break;
            }
            if(flag==1)break;
        }
        int paper_max = conference->getTrack(maxj).getSession ( maxi ).getPaper ( maxk );
        if(flag == 1) {
            // cout<<"recursing on max neighbour "<<endl;
            conference->setPaper ( maxj, maxi, maxk, min_paper);
            conference->setPaper ( minj, mini, mink, max_paper);
            // return;
            // organizePapers2 ();
        }
        else{
            break;
        }
    } while (flag==1);    
}














double SessionOrganizer::fast_score(double curr_state_goodness, int t1,int p1,int k1, int t2,int p2, int k2,int c_trade_off){
        double goodness=0;
        // int[][] temp=curr_state;
        int t_size=conference->getSessionsInTrack();
        int p_size=conference->getParallelTracks();
        int k_size=conference->getPapersInSession();
        double sub=0;
        sub+=sim_in_a_session(t1,p1,k1,-1,-1,-1,true);
        sub+=sim_in_a_session(t2,p2,k2,-1,-1,-1,true);
        if(t1==t2){
            sub+=c_trade_off*sim_in_a_session(t1,p1,-1,t2,p2,k2,false);
            sub+=c_trade_off*sim_in_a_session(t2,p2,-1,t1,p1,k1,false);
        }
        else{
            for (int i = 0; i < p_size; ++i){
                if(i!=p1){
                    sub+=c_trade_off*sim_in_a_session(t1,i,-1,t1,p1,k1,false);
                }
            }
            for (int i = 0; i < p_size; ++i){
                if(i!=p2)
                    sub+=c_trade_off*sim_in_a_session(t2,i,-1,t2,p2,k2,false);
            }
               
        }
        int paper_initial = conference->getTrack(p1).getSession (t1).getPaper ( k1 );
        int paper_final = conference->getTrack(p2).getSession ( t2 ).getPaper ( k2 );
        conference->setPaper(p1,t1,k1,paper_final);
        conference->setPaper(p2,t2,k2,paper_initial);
        double add=0;
        add+=sim_in_a_session(t1,p1,k1,-1,-1,-1,true);
        add+=sim_in_a_session(t2,p2,k2,-1,-1,-1,true);
        if(t1==t2){
            add+=c_trade_off*sim_in_a_session(t1,p1,-1,t2,p2,k2,false);
            add+=c_trade_off*sim_in_a_session(t2,p2,-1,t1,p1,k1,false);
        }
        else{
            for (int i = 0; i < p_size; ++i)
            {
                if(i!=p1){
                    add+=c_trade_off*sim_in_a_session(t1,i,-1,t1,p1,k1,false);
                }
            }
            for (int i = 0; i < p_size; ++i){
                if(i!=p2)
                    add+=c_trade_off*sim_in_a_session(t2,i,-1,t2,p2,k2,false);
            }
        }

        goodness=curr_state_goodness+add-sub;
        conference->setPaper(p1,t1,k1,paper_initial);
        conference->setPaper(p2,t2,k2,paper_final);


        // int tmp=temp[p1][k1];
        // temp[p1][k1]=temp[p2][k2];
        // temp[p2][k2]=tmp;
        // goodness=curr_state_goodness-sim_in_a_session(curr_state,p1,k1,-1,-1,true)-sim_in_a_session(curr_state,p2,k2,-1,-1,true)-c_trade_off*sim_in_a_session(curr_state,p2,-1,p1,k1,false)-c_trade_off*sim_in_a_session(curr_state,p1,-1,p2,k2,false)+sim_in_a_session(temp,p1,k1,-1,-1,true)+sim_in_a_session(temp,p2,k2,-1,-1,true)+c_trade_off*sim_in_a_session(temp,p2,-1,p1,k1,false)+c_trade_off*sim_in_a_session(temp,p1,-1,p2,k2,false);
        return(goodness);
}

double SessionOrganizer::sim_in_a_session(int t,int p,int k,int tdiff,int pdiff, int kdiff, bool sim){
        int t_size=conference->getSessionsInTrack();
        int p_size=conference->getParallelTracks();
        int k_size=conference->getPapersInSession();
        double similarity=0;double diff=0;
        // int[][] temp=state;
        for(int i =0;i<k_size;i++){
            if(sim){
                if(i!=k)
                    similarity+=1 - distanceMatrix[conference->getTrack(p).getSession (t).getPaper(i)][conference->getTrack(p).getSession (t).getPaper(k)];
                // similarity += 1-dis_bw_papers[temp[p][i]][temp[p][k]];  
            }
            else{
                diff += distanceMatrix[conference->getTrack(p).getSession (t).getPaper(i)][conference->getTrack(pdiff).getSession (tdiff).getPaper(kdiff)];  
            }
        }   
        if(sim) return similarity;
        else return diff;
}

void SessionOrganizer::randomState (int seed)
{
        // int a[t][p][k];
        int t = conference->getSessionsInTrack ( );
        int p = conference->getParallelTracks ( );
        int k = conference->getPapersInSession ( );
        // Random rand=new Random();
        int r_int=0;
        int n = t*p*k;
        int random_arr [n];
        for (int i=0;i<n;i++){
            random_arr[i]=(i);
        }
        for(int i=0;i<n;i++){
            srand(seed+i);
            // int f1=rand();
            // int f = f1 % n;
            int f = i;
            int s1=rand();
            int s = s1 % n;
            int temp=random_arr[f];
            random_arr[f]=random_arr[s];
            random_arr[s]=temp;
        }
        int count=0;
        for(int j=0;j<p;j++){
            for (int i=0;i<t;i++) {
                for (int z=0;z<k;z++) {
                     conference->setPaper(j,i,z,random_arr[count]);
                    count++;
                }
                
            }
        }
}



void SessionOrganizer::readInInputFile ( string filename )
{
    vector<string> lines;
    string line;
    ifstream myfile ( filename.c_str () );
    if ( myfile.is_open ( ) )
    {
        while ( getline ( myfile, line ) )
        {
            //cout<<"Line read:"<<line<<endl;
            lines.push_back ( line );
        }
        myfile.close ( );
    }
    else
    {
        cout << "Unable to open input file";
        exit ( 0 );
    }

    if ( 6 > lines.size ( ) )
    {
        cout << "Not enough information given, check format of input file";
        exit ( 0 );
    }

    processingTimeInMinutes = atof ( lines[0].c_str () );
    papersInSession = atoi ( lines[1].c_str () );
    parallelTracks = atoi ( lines[2].c_str () );
    sessionsInTrack = atoi ( lines[3].c_str () );
    tradeoffCoefficient = atof ( lines[4].c_str () );

    int n = lines.size ( ) - 5;
    double ** tempDistanceMatrix = new double*[n];
    for ( int i = 0; i < n; ++i )
    {
        tempDistanceMatrix[i] = new double[n];
    }


    for ( int i = 0; i < n; i++ )
    {
        string tempLine = lines[ i + 5 ];
        string elements[n];
        splitString ( tempLine, " ", elements, n );

        for ( int j = 0; j < n; j++ )
        {
            tempDistanceMatrix[i][j] = atof ( elements[j].c_str () );
        }
    }
    distanceMatrix = tempDistanceMatrix;

    int numberOfPapers = n;
    int slots = parallelTracks * papersInSession*sessionsInTrack;
    if ( slots != numberOfPapers )
    {
        cout << "More papers than slots available! slots:" << slots << " num papers:" << numberOfPapers << endl;
        exit ( 0 );
    }
}

double** SessionOrganizer::getDistanceMatrix ( )
{
    return distanceMatrix;
}

void SessionOrganizer::printSessionOrganiser ( char * filename)
{
    conference->printConference ( filename);
}

double SessionOrganizer::scoreOrganization ( )
{
    // Sum of pairwise similarities per session.
    double score1 = 0.0;
    for ( int i = 0; i < conference->getParallelTracks ( ); i++ )
    {
        Track tmpTrack = conference->getTrack ( i );
        for ( int j = 0; j < tmpTrack.getNumberOfSessions ( ); j++ )
        {
            Session tmpSession = tmpTrack.getSession ( j );
            for ( int k = 0; k < tmpSession.getNumberOfPapers ( ); k++ )
            {
                int index1 = tmpSession.getPaper ( k );
                for ( int l = k + 1; l < tmpSession.getNumberOfPapers ( ); l++ )
                {
                    int index2 = tmpSession.getPaper ( l );
                    score1 += 1 - distanceMatrix[index1][index2];
                }
            }
        }
    }

    // Sum of distances for competing papers.
    double score2 = 0.0;
    for ( int i = 0; i < conference->getParallelTracks ( ); i++ )
    {
        Track tmpTrack1 = conference->getTrack ( i );
        for ( int j = 0; j < tmpTrack1.getNumberOfSessions ( ); j++ )
        {
            Session tmpSession1 = tmpTrack1.getSession ( j );
            for ( int k = 0; k < tmpSession1.getNumberOfPapers ( ); k++ )
            {
                int index1 = tmpSession1.getPaper ( k );

                // Get competing papers.
                for ( int l = i + 1; l < conference->getParallelTracks ( ); l++ )
                {
                    Track tmpTrack2 = conference->getTrack ( l );
                    Session tmpSession2 = tmpTrack2.getSession ( j );
                    for ( int m = 0; m < tmpSession2.getNumberOfPapers ( ); m++ )
                    {
                        int index2 = tmpSession2.getPaper ( m );
                        score2 += distanceMatrix[index1][index2];
                    }
                }
            }
        }
    }
    double score = score1 + tradeoffCoefficient*score2;
    return score;
}