import java.util.*;


// Class to hold HoshenKopelman method which can tell us if grid is connected / how many
// disconnected clusters are there / how big are they(?)... etc...
// Use this to find size of all parasite clusters in a lattice snapshot.
// In StackedCP label "0" is empty site, "1" is a host, "2" is first parasite", etc. 


public class HoshenKopelman {


    static int Mod(int a, int b) { // "a mod b"
    /** Modulus for periodic boundary conditions  **/
	    int solution = (a % b + b) % b;
	    return solution;
    }

    
    public static void replaceAll( int[][] HK_grid, int a, int b, int xLength, int yLength ){
    /**  Method to replace all 'a's with 'b's in an array, int[][].  **/	
	    for( int i=0; i<yLength; i++){
	        for( int j=0; j<xLength; j++){
		        if( HK_grid[j][i] == a ){
		            HK_grid[j][i] = b;
		        }
	        }
	    }	
    }
       

    public static ArrayList<Integer> getClusterSizes( Cell[][] cellGrid, int xLength, int yLength, int labelToCount ){
	/** Hoshen-Kopelman algorithm to return a list of cluster sizes for sites equal or higher than labelToCount" **/

        // Make int[][] grid of Cell labels
        int[][] grid = new int[xLength][yLength];
	    for( int i=0; i<yLength; i++){   
	        for( int j=0; j<xLength; j++){
                grid[j][i] = cellGrid[j][i].getLabel() ;
            }
        }

	    // Array to hold temporary labels. All entries initially 0
	    int[][] HK_grid = new int[xLength][yLength];

	    int HK_label = 1;
	    for( int i=0; i<yLength; i++){   
	        for( int j=0; j<xLength; j++){

	        /* %%%%%%% PRINT OUT HK_grid %%%%%%%%%%
	        // System.out.println();
	        // for( int aaa=0; aaa < HK_grid.length; aaa++ ){
	        //     System.out.println( Arrays.toString( HK_grid[aaa] )  );
	        // }
	        */

		    if(  grid[j][i] >= labelToCount  ){    
		        //Count clusters + sizes of all levels equal and higher than labelToCount.
		        // (i.e. for host we want all sites 1,2,3...,m (susceptible + all infected states)

		        if( HK_grid[j][Mod(i-1,yLength)] != 0){
			        HK_grid[j][i] = HK_grid[j][Mod(i-1,yLength)];
			        //Now if left square is different, do union
			    if( HK_grid[Mod(j-1,xLength)][i]!=0  &&  HK_grid[Mod(j-1,xLength)][i] != HK_grid[j][i] ){
			        //Then go through replacing all occurrences of one with the other
			        replaceAll(HK_grid, HK_grid[Mod(j-1,xLength)][i], HK_grid[j][i], xLength, yLength );
			    }
		
		        }
		        else if( HK_grid[Mod(j-1,xLength)][i] != 0 ){
			        HK_grid[j][i] = HK_grid[Mod(j-1,xLength)][i];
			        //Since we are here, the up square must be a zero so no union needed.						
		        }
		        else{
			        HK_grid[j][i] = HK_label;
			        HK_label++;
		        }


		        //Now check right neighbour and do union if necessary
		        if( HK_grid[j][Mod(i+1,yLength)] !=0  && HK_grid[j][Mod(i+1,yLength)] != HK_grid[j][i]  ){
			    replaceAll(HK_grid, HK_grid[j][i], HK_grid[j][Mod(i+1,yLength)], xLength, yLength );
		        }
		        //Now check down neighbour and do union if necessary
		        if( HK_grid[Mod(j+1,xLength)][i] !=0  && HK_grid[Mod(j+1,xLength)][i] != HK_grid[j][i]  ){
		        	replaceAll(HK_grid,  HK_grid[j][i], HK_grid[Mod(j+1,xLength)][i], xLength, yLength );
		        }

		    }
	        }
	    }


	    //%%%%%%% PRINT OUT HK_grid %%%%%%%%%%
	    // System.out.println();
	    // for( int i=0; i < HK_grid.length; i++ ){
	    //     System.out.println( Arrays.toString( HK_grid[i] )  );
	    // }
	    //%%%%%%%%%%%%%%%


	    //Count the number of disconnected domains (i.e. number of distinct labels in HK_grid)
	    ArrayList<Integer> labels = new ArrayList<Integer>();
	    for( int i=0; i<yLength; i++){                       
	        for( int j=0; j<xLength; j++){		
		        if( HK_grid[j][i] != 0){
		            // Check if this label is already in the labels list, add if not
		            if( !labels.contains( HK_grid[j][i] ) ){
			            labels.add( HK_grid[j][i] );
		            }			
		        }		
	        }
	    }


	    //Now for each cluster label, get cluster size:
	    ArrayList<Integer> clusterSizes = new ArrayList<Integer>();
	    for( Integer lab: labels ){
	        int clusterSize = 0;
	        for( int i=0; i<yLength; i++){
		    for( int j=0; j<xLength; j++){
		        if( HK_grid[j][i] != 0  &&  HK_grid[j][i] == lab ){
			    clusterSize++;
		        }
		    }
	        }
	        clusterSizes.add( clusterSize );
	    }

	    return clusterSizes;
    }













}
