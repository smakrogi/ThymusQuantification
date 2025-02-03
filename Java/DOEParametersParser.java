/** DOE command line generation for
 *  thymus analysis.
 *  @author: Sokratis Makrogiannis, 3T @ NIA
 **/

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Iterator;
import java.io.IOException;
import java.lang.Math;


// This program reads in a text file with a description of the
// parameters that we want to use in a DOE and then generates several command
// lines as a helper in a perl script.
// the format of the text file is
// $parametername,$parametertype,$minvalue,$maxvalue,$numberofsamples.
public class DOEParametersParser {

    // Class attributes.
    int nParams;
    ArrayList<AlgorithmParameter> parameterInputs;

    DOEParametersParser() {
	this.nParams = 0;
	this.parameterInputs = new ArrayList<AlgorithmParameter>();
    }

    // Read the parameter names and
    // values from comma-separated-values
    // text file.
    public static void main(String[] args) throws IOException, 
						  java.lang.IllegalArgumentException{

        if(args.length < 2) {
	    System.out.println("syntax: DOEParametersParser <in file> <out file>");
	    System.exit(1);
        }
        
        DOEParametersParser parameterParser = new DOEParametersParser();
        BufferedReader reader = null;
        try{
            reader = new BufferedReader(new FileReader(args[0]));
        }
        catch(IOException ex) {
            System.out.println("Input file not found");
	    System.exit(1);
        }

	String paramString;
	String delims = ",";

	// Read each line and compute the parameter sampling.
	while( (paramString= reader.readLine()) != null ) {

	    // Compute parameter ranges and
	    // create the command line.
	    String[] tokens = paramString.split(delims);

	    // Parse strings.
	    AlgorithmParameter tempParameter = new AlgorithmParameter();
	    tempParameter.setName(tokens[0]);
	    tempParameter.setType(tokens[1]);
	    parameterParser.parameterInputs.add( tempParameter );
	    Double minValue = Double.valueOf(tokens[2]);
	    Double maxValue = Double.valueOf(tokens[3]);
	    Double n = Double.valueOf(tokens[4]);
	    Double step = null;

	    // Compute the step for parameter sampling.
	    if(n==0) {
		throw new java.lang.IllegalArgumentException("Number of samples in parameter space must be greater than zero");
	    }
 	    if (n==1)
 		step = (maxValue - minValue);
 	    else if (n>1)
		step = ( (maxValue - minValue) / (n-1) );
	

	    for (int i=0;i<n.intValue();i++) {
		Double tempValue = minValue + i * step;
		parameterParser.parameterInputs.get(parameterParser.nParams).addValue( tempValue );
	    }
 	    // System.out.println("Parameter values");
 	    System.out.println(parameterParser.parameterInputs.get(parameterParser.nParams).toString());
	    parameterParser.nParams++;
	}
	reader.close();

	// Generate the command lines.
	System.out.println("Command line generation");
	ArrayList<String> commandLines = parameterParser.createCommandLines();

	// Write array of strings to text file.
        BufferedWriter writer = null;
        try {
            writer = new BufferedWriter(new FileWriter(args[1]));
        }
        catch(IOException ex) {
            System.out.println("Cannot write output.");
	    System.exit(1);
        }
	Iterator<String> itr = commandLines.iterator();
	while(itr.hasNext()){
		writer.write(itr.next()+"\n");
	}
	writer.close();
	System.out.println("All done");
    }


    // Create the command lines.
    ArrayList<String> createCommandLines() {

	// Use "Cartesian" combination of the parameter strings.
	int count = 0;
	AlgorithmParameter firstParam, secondParam;
        // firstParam = new AlgorithmParameter();
	firstParam = this.parameterInputs.get(0);
	ArrayList<String> firstCommandStrings = null, secondCommandStrings = null;
	firstCommandStrings = createStringExpressions(this.parameterInputs.get(count));
	// System.out.println(this.nParams);
	while( count < (this.nParams-1) ) {
	    secondCommandStrings = createStringExpressions(this.parameterInputs.get(count+1));
	    firstCommandStrings = combineTwoParameters(firstCommandStrings, secondCommandStrings);
	    count++;
	}

// 	System.out.println("Generated parameter combinations");
// 	System.out.println(firstCommandStrings.toString());
	return firstCommandStrings;
    }


    // Generate command line substrings for one parameter.
    ArrayList<String> createStringExpressions(AlgorithmParameter parameterClass){
	ArrayList<String> outputStrings = new ArrayList<String>();
	Iterator<Double> it = parameterClass.values.iterator();
	String tempString;
	while( it.hasNext() ) {
	    if(parameterClass.getType().equalsIgnoreCase("int")) {
		long temp = Math.round(it.next());
		tempString = parameterClass.name + " " + temp; 
	    }
	    else
		tempString = parameterClass.name + " " + String.format("%.4f", it.next());
	    outputStrings.add(tempString);
	}
	return outputStrings;
    }


    // Cartesian product of two vectors.
    ArrayList<String> combineTwoParameters( ArrayList<String> firstParameter,
					    ArrayList<String> secondParameter) {

	Iterator<String> itr1, itr2;
	ArrayList<String> combinedParameters = new ArrayList<String>();
        
        itr1 = firstParameter.iterator();
	while(itr1.hasNext()) {
	    itr2 = secondParameter.iterator();
	    String tempString1 = itr1.next();
	    while(itr2.hasNext()){
		String tempString2 = tempString1 + " " + itr2.next();
 		combinedParameters.add( tempString2 );
		// System.out.println( tempString2 );
	    }
	}
	return combinedParameters;
    }

}



// The class that handles the doe parameters.
class AlgorithmParameter
{
    String name, type;
    ArrayList<Double> values;

    // Default constructor.
    AlgorithmParameter() {
	this.name = null;
	this.type = null; // "float" or "int"
	this.values = new ArrayList<Double>();
    }

    // Initialization using name and value.
    AlgorithmParameter(String paramName, Double paramValue) {
	this.name = paramName;
	this.values = new ArrayList<Double>();
	this.values.add( paramValue );
    }

    // Name attribute.
    public void setName(String paramName) {
	this.name = paramName;
    }

    public String getName() {
	return this.name;
    }

    // Type attribute.
    public void setType(String paramType) {
	this.type = paramType;
    }

    public String getType() {
	return this.type;
    }

    // Add a new value.
    public void addValue( Double paramValue) {
	this.values.add( paramValue );
    }

    // Return all values.
    public ArrayList<Double> getValues(){
	return this.values;
    }

    // Display parameter.
    public String toString(){
	String tempString =  this.getName() + "=" + this.getValues().toString();
	return tempString;
    }

}