/*
Simulation of a series of stacked contact processes with n levels.
It also uses a kinetic update scheme (gillespie-style choice of which reactions occur and when).
This stores absolute next reaction times for each lattice site, not waiting times.
Hard-coded Von Neumann neighbourhood (4 neighbours at NSEW)

Terminal command: java StackedCP gridSize levels birthRates deathRate
                                 int      int    double[]   double[]
*/

import java.io.*;
import java.util.*;

public class StackedCP{


    // Modulus for periodic boundary conditions
    static int Mod(int a, int b) { // "a mod b"
        int solution = (a % b + b) % b;
        return solution;
    }//______________________________________________________________________



    static double[] getDeathProps( int levels, int siteLabel, double[] dRates ){
    double[] dProps = new double[ levels ];
    for(int i=0; i<levels; i++){

        if( i < siteLabel ){
            dProps[i] = dRates[i];
        }
        else{ dProps[i] = 0;  }
    }
    return dProps;

    }//______________________________________________________________________



    static double[] getBirthProps( int levels, int siteLabel, double[] bRates, int[] neighbourList ){

    double[] bProps = new double[ levels + 1 ];
    for( int i=0; i<(levels+1); i++ ){

        if( i <= siteLabel ){ bProps[i] = 0; }
        else{
            if( neighbourList[i] == 0 ){ bProps[i] = 0; }
            else{
                bProps[i] = ( bRates[ siteLabel ] / 4.0 ) * (double)neighbourList[i] ;
            }
        }
    }
    return bProps;

    }//_______________________________________________________________________




    static double getWaitingTime(double[] dProps, double[] bProps){

        double propensity = 0;
        for( int i=0; i<dProps.length; i++){
            propensity = propensity + dProps[i];
        }
        for( int j=0; j< bProps.length; j++){
            propensity = propensity + bProps[j];
        }

        double rand = Math.random();
        double tau=0;;
        if(propensity!=0){    
            tau = -(1.0/propensity)* ( Math.log(1-rand) ) ;  //1000*R since t->infinty very fast
        }
        if(propensity==0){                                     //Empty site surrounded by empty sites
            tau = -(1.0E3);                                                                  
        }

    return tau;
    }//______________________________________________________________________



    static int[] getNeighbourList(int levels, Cell[][] lattice, int xSite, int ySite, int GRID_SIZE  ){

        int[] neighbourList = new int[ levels + 1 ];
    
        int b = lattice[Mod(xSite - 1, GRID_SIZE)][ySite].getLabel();
        int c = lattice[Mod(xSite + 1, GRID_SIZE)][ySite].getLabel();
        int d = lattice[xSite][Mod(ySite - 1, GRID_SIZE)].getLabel();
        int e = lattice[xSite][Mod(ySite + 1, GRID_SIZE)].getLabel();

        neighbourList[b]++;
        neighbourList[c]++;
        neighbourList[d]++;
        neighbourList[e]++;

        return neighbourList;

    }//______________________________________________________________________


    static void perform_random_event( Cell[][] lattice, int ii, int jj ){
    /** Choose a random event for cell[ii][jj] and update label **/
        
        // Sum propensities of all possible cell transitions
        double tot_prop = 0;
        for( int i=0; i < lattice[ii][jj].dProps.length; i++){
            tot_prop = tot_prop + lattice[ii][jj].dProps[i];
        }
        for( int j=0; j< lattice[ii][jj].bProps.length; j++){
            tot_prop = tot_prop + lattice[ii][jj].bProps[j];
        }

        double indicator = Math.random() * tot_prop;
        
        if(tot_prop!=0){
            double counter = 0;
            boolean found = false;
            changeLoop:
            for(int z=0; z < lattice[ii][jj].dProps.length; z++){
                counter = counter + lattice[ii][jj].dProps[z];

                if(indicator < counter){
                lattice[ii][jj].setLabel( z );
                    found = true;
                break changeLoop;
                }
            }
            if( !found ){
                changeLoop2:
                for(int zz=0; zz < lattice[ii][jj].bProps.length; zz++){
                counter = counter + lattice[ii][jj].bProps[zz];

                if(indicator < counter){
                    lattice[ii][jj].setLabel(zz);
                    break changeLoop2;
                }
                }
            } 
        }    

    }//______________________________________________________________________




    static void update_times_propensities( Cell[][] lattice, int GRID_SIZE, double curTime, int levels, int ii, int jj, double[] birth, double[] death){
    /** Update propensities and event time of affected cell + neighbours **/

        lattice[ii][jj].neighbourList =  getNeighbourList( levels, lattice, ii, jj, GRID_SIZE );
        lattice[ii][jj].bProps = getBirthProps( levels, lattice[ii][jj].getLabel(), birth, lattice[ii][jj].neighbourList );
        lattice[ii][jj].dProps = getDeathProps( levels, lattice[ii][jj].getLabel(), death );
        lattice[ii][jj].time = curTime + getWaitingTime( lattice[ii][jj].dProps, lattice[ii][jj].bProps );

        lattice[Mod(ii - 1, GRID_SIZE)][jj].neighbourList =  getNeighbourList( levels, lattice, Mod(ii - 1, GRID_SIZE), jj, GRID_SIZE );
        lattice[Mod(ii - 1, GRID_SIZE)][jj].bProps = getBirthProps( levels, lattice[Mod(ii - 1, GRID_SIZE)][jj].getLabel(), birth, lattice[Mod(ii - 1, GRID_SIZE)][jj].neighbourList );
        lattice[Mod(ii - 1, GRID_SIZE)][jj].dProps = getDeathProps( levels, lattice[Mod(ii - 1, GRID_SIZE)][jj].getLabel(), death );
        lattice[Mod(ii - 1, GRID_SIZE)][jj].time = curTime + getWaitingTime( lattice[Mod(ii - 1, GRID_SIZE)][jj].dProps, lattice[Mod(ii - 1, GRID_SIZE)][jj].bProps );

        lattice[Mod(ii + 1, GRID_SIZE)][jj].neighbourList =  getNeighbourList( levels, lattice,Mod(ii + 1, GRID_SIZE) , jj, GRID_SIZE );
        lattice[Mod(ii + 1, GRID_SIZE)][jj].bProps = getBirthProps( levels, lattice[Mod(ii + 1, GRID_SIZE)][jj].getLabel(), birth, lattice[Mod(ii + 1, GRID_SIZE)][jj].neighbourList );
        lattice[Mod(ii + 1, GRID_SIZE)][jj].dProps = getDeathProps( levels, lattice[Mod(ii + 1, GRID_SIZE)][jj].getLabel(), death );
        lattice[Mod(ii + 1, GRID_SIZE)][jj].time = curTime + getWaitingTime( lattice[Mod(ii + 1, GRID_SIZE)][jj].dProps, lattice[Mod(ii + 1, GRID_SIZE)][jj].bProps );

        lattice[ii][Mod(jj - 1, GRID_SIZE)].neighbourList =  getNeighbourList( levels, lattice, ii, Mod(jj - 1, GRID_SIZE), GRID_SIZE );
        lattice[ii][Mod(jj - 1, GRID_SIZE)].bProps = getBirthProps( levels, lattice[ii][Mod(jj - 1, GRID_SIZE)].getLabel(), birth, lattice[ii][Mod(jj - 1, GRID_SIZE)].neighbourList );
        lattice[ii][Mod(jj - 1, GRID_SIZE)].dProps = getDeathProps( levels, lattice[ii][Mod(jj - 1, GRID_SIZE)].getLabel(), death );
        lattice[ii][Mod(jj - 1, GRID_SIZE)].time = curTime + getWaitingTime( lattice[ii][Mod(jj - 1, GRID_SIZE)].dProps, lattice[ii][Mod(jj - 1, GRID_SIZE)].bProps );

        lattice[ii][Mod(jj + 1, GRID_SIZE)].neighbourList =  getNeighbourList( levels, lattice, ii, Mod(jj + 1, GRID_SIZE), GRID_SIZE );
        lattice[ii][Mod(jj + 1, GRID_SIZE)].bProps = getBirthProps( levels, lattice[ii][Mod(jj + 1, GRID_SIZE)].getLabel(), birth, lattice[ii][Mod(jj + 1, GRID_SIZE)].neighbourList );
        lattice[ii][Mod(jj + 1, GRID_SIZE)].dProps = getDeathProps( levels, lattice[ii][Mod(jj + 1, GRID_SIZE)].getLabel(), death );
        lattice[ii][Mod(jj + 1, GRID_SIZE)].time = curTime + getWaitingTime( lattice[ii][Mod(jj + 1, GRID_SIZE)].dProps, lattice[ii][Mod(jj + 1, GRID_SIZE)].bProps );
   
    }//___________________________________________________________________________




    // Method to find density
    static double[] find_densities_uninfected(int levels, Cell[][] lattice, int GRID_SIZE){

    double[] density = new double[ levels+1 ];
    for(int i=0; i<GRID_SIZE; i++){
        for(int j=0; j<GRID_SIZE; j++){

        int type = lattice[i][j].getLabel();
        density[type] = density[type] + 1.0/( (double)(GRID_SIZE*GRID_SIZE) );
        }
    }
    return density;
    }//_______________________________________________________________________










    public static void main(String[] args) throws IOException{

        PrintWriter out = new PrintWriter( new FileWriter( "density_tot.dat") );
        PrintWriter out2 = new PrintWriter( new FileWriter( "density_uninfected.dat" ) );

        int GRID_SIZE = Integer.parseInt( args[0] );
        int levels = Integer.parseInt( args[1] );
        double equalBIRTHS = Double.parseDouble(args[2]);  
        double equalDEATHS = Double.parseDouble(args[3]);

        // Specify this fraction if we want to vary b or d rate by set 
        // proportion per level.
        //double fraction = Double.parseDouble(args[5]);

        Cell[][] lattice = new Cell[GRID_SIZE][GRID_SIZE];
        double[] density_uninfected = new double[ levels + 1];
        double[] forAverageDens = new double[levels+1];

        double curTime=0;


        //-------   Initialise birth and death rates           
        double[] birth = new double[ levels + 1 ];
        for(int i=0; i<birth.length; i++){
            birth[i] = equalBIRTHS;
        }
        birth[ levels ] = 0;
            
        double[] death = new double[ levels ];
        for(int i=0; i<death.length; i++){
            death[i] = equalDEATHS;
        }    
        /*
        double[] death = new double[levels];
        death[0] = equalDEATHS;
        for(int i=1; i<death.length; i++){
            death[i] = death[i-1]*fraction;
        }
        */
        /*
        double[] birth = new double[ levels + 1 ];
        birth[0] = equalBIRTHS;
        for(int i=1; i<birth.length; i++){
            birth[i] = birth[i-1]*(1.0/fraction);
        }
        birth[ levels ] = 0;
        */


        //Initialise lattice to all having maximum level occupancy
        for(int i=0; i<GRID_SIZE; i++){
            for(int j=0; j<GRID_SIZE; j++){
               lattice[i][j] = new Cell( levels );
            }
        }


        //Initialise propensities and waiting times
        for(int i=0; i<GRID_SIZE; i++){
            for(int j=0; j<GRID_SIZE; j++){
                int siteLabel = lattice[i][j].getLabel();

                lattice[i][j].neighbourList =  getNeighbourList( levels, lattice, i, j, GRID_SIZE  );
                lattice[i][j].bProps = getBirthProps( levels, siteLabel, birth, lattice[i][j].neighbourList );
                lattice[i][j].dProps = getDeathProps( levels, siteLabel, death );
                lattice[i][j].time = curTime + getWaitingTime( lattice[i][j].dProps, lattice[i][j].bProps );
            }
        }







        //------ Main code updating grid state ----------
        long steps = 2000000L;
        //long steps = 20000000000L;
        double maxRunTime = 1.0E20;
        long avgCounter=0;
        mainloop:
        for(long k=0; k<steps; k++){


            //find smallest waiting time
            double ooTime = maxRunTime;  int ii=0;   int jj=0;
            for(int i=0; i<GRID_SIZE; i++){
                for(int j=0; j<GRID_SIZE; j++){
                    if(lattice[i][j].time > curTime  &&  lattice[i][j].time < ooTime){
                        ii = i;    jj = j;
                        ooTime = lattice[i][j].time;
                    }
                }
            }
            if(ooTime < curTime || ooTime == maxRunTime ){ // 
                System.out.println("Everything is dead");
                break mainloop;
            }
            
                     
            //Choose and perform event on site [ii][jj]. 
            perform_random_event( lattice, ii, jj );

            // Update time
            curTime =  ooTime;

            //Update state of cell[ii][jj]'s propensities and event time + those of its nearest neighbours
            update_times_propensities( lattice, GRID_SIZE, curTime, levels, ii, jj, birth, death ) ;


            
            // Print data at desired frequency once in 'steady state'
            if( k%7000==0 && k>10000){

                /*
                // Cluster sizes with Hoshen-Kopelman
                ArrayList<Integer> cluster_sizes = new ArrayList<Integer>();
                cluster_sizes = HoshenKopelman.getClusterSizes( lattice, GRID_SIZE, GRID_SIZE, 4 );
                */

                avgCounter ++;
                
                // Get raw state densities (i.e. lattice labels)
                density_uninfected = find_densities_uninfected(levels, lattice, GRID_SIZE);

                // Get actual densities (sum of each level since a '3' implies a '2' is present also)
                double[] density_actual = new double[ density_uninfected.length ];
                for(int v=1;  v<density_uninfected.length;  v++){    
                    for(int w=v; w<density_uninfected.length; w++){
                        density_actual[v] = density_actual[v] + density_uninfected[w]; 
                    }
                }

                // Print time and densities             
                out.print( curTime + " " );            
                for(int a=1; a<density_actual.length; a++){
                    out.print(density_actual[a] + " ");
                } 
                out.println("");

                // Print uninfected densities of each level
                out2.print( curTime + " " );            
                for(int v=1; v<density_uninfected.length; v++){
                    out2.print(density_uninfected[v] + " ");
                } 
                out2.println("");
               

                //cumulative densities
                //for(int i=0; i<forAverageDens.length; i++){
                //    forAverageDens[i] = forAverageDens[i] + actualDensities[i];
                //}
               

            }

        }// k steps
        


        out.close();
        out2.close();




    }//main



















}//class
