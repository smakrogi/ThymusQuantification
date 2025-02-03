/** DOE command line generation for
 *  thymus analysis.
 *  @author: Sokratis Makrogiannis, 3T @ NIA
 **/

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Iterator;
import java.io.BufferedWriter;
import java.io.FileWriter;

public class ParseDOEParameters
{

ParseDOEParameters() {
	this.nParams = 0;
	this.parameterInputs = new ArrayList<AlgorithmParameter>();
}

// Read the parameter names and
// values from comma-separated-values
// text file.
public static void main(String[] args) {

	BufferedReader reader = new BufferedReader(new FileReader(argv[0]));

	String paramString;
	String delims = ",";
	nParams = 0;

	while( (paramString=reader.readLine()) != null ) {

	// Compute parameter ranges and
	// create the command line.
	String[] tokens = paramString.split(delims);

	// Parse string.
	AlgorithmParameter tempParameter = new AlgorithmParameter();
	tempParameter.setName = tokens[0];
	this.parameterInputs.add( tempParameter );
	Double minValue = Double.valueOf(tokens[1].trim()).DoubleValue();
	Double maxValue = Double.valueOf(tokens[2].trim()).DoubleValue();
	Double step = Double.valueOf(tokens[3].trim()).DoubleValue();
	Double n = (maxValue - minValue) / step;

	for (Int i=1;i<Int(n);i++) {
		Double tempValue = minValue + (n-1) * step;
		this.parameterInputs[count].addValue( tempValue );
	}

	nParams++;
	}
reader.close();

// Generate the command lines.
this.createCommandLines();

// Write array of strings to text file.
ArrayList<String> commandLines = createCommandLines();
File outFile = new File(argv[1]);
Writer writer = null;
writer = new BufferedWriter(new FileWriter(outFile));
Iterator<String> itr;
while(itr.hasNext()){
	writer.write(itr.next());
}
writer.close();
}


// Create the command lines.
void ArrayList<String> createCommandLines() {

	// Use "Cartesian" combination of the parameter strings.
	int count = 0;
	AlgorithmParameter firstParam, secondParam;
	firstParam = this.parameterInputs[0];
	ArrayList<String> firstCommandStrings, secondCommandStrings;
	firstCommandStrings = createStringExpressions(this.parameterInputs[count]);

	while( count < this.nParams ) {
		secondCommandStrings = createStringExpressions(this.parameterInputs[count+1]);
		firstCommandStrings = combineTwoParameters(firstCommandStrings, secondCommandStrings);
	}

return firstCommandStrings;
}


// Generate command line substrings for one parameter.
ArrayList<String> createStringExpressions(AlgorithmParameter parameterClass){
	ArrayList<String> outputStrings;
	Iterator<Double> it = parameterClass.values.iterator();

	while( it.hasNext() ) {
		String tempString = parameterClass.name + " " + it.next().toString;
		outputStrings.add(tempString);
	}
	return outputStrings;
}


// Cartesian product of two vectors.
ArrayList<String> combineTwoParameters( ArrayList<String> firstParameter,
						 ArrayList<String> secondParameter) {

Iterator<String> itr1, itr2;

ArrayList<String> combinedParameters = new ArrayList<String>();
while(itr1.hasNext()) {
	while(itr2.hasNext()){
		String tempString = itr1.next() + " " + itr2.next();
 		combinedParameters.add( tempString );
	}
}
return combinedParameters;
}

// Class attributes.
Int nParams;
ArrayList<AlgorithmParameter> parameterInputs;
}


class AlgorithmParameter
{

	AlgorithmParameter() {
		this.name = null;
		this.values = new ArrayList<Double>();
	}


	AlgorithmParameter(String paramName, Double paramValue) {
		this.name = paramName;
		this.values = new ArrayList<Double>();
		this.values.add( paramValue );
	}

	void setName(String paramName) {
		this.name = paramName;
	}

	String getName() {
		return this.name;
	}

	void addValue( Double paramValue) {
		this.values.add( paramValue );
	}

	ArrayList<Double> getValues(){
		return this.values;
	}

	String name;
	ArrayList<Double> values;
}