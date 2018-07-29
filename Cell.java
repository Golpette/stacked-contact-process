
public class Cell{

    int label;
    double[] bProps; 
    double[] dProps; 
    double time; 
    int[] neighbourList;

    public Cell(int label){
        this.label = label;
    }

    public int getLabel(){
        return this.label;
    }

    public void setLabel( int a ){
        this.label = a;
    }

}
