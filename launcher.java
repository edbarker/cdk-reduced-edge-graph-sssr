package barkerSSSRFinder;

import org.openscience.cdk.*;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IMoleculeSet;

public class launcher {
	
	public static void main(String[] args) throws Exception 
    {
		boolean test;
	
		try 
		{
			SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
			IMolecule mol = sp.parseSmiles("NC1=CC2=C(C=C1)C(=O)C3=C(C=CC=C3)C2=O");
			IMoleculeSet cycles = DefaultChemObjectBuilder.getInstance().newMoleculeSet();
		   
			test=barkerSSSRFinder.findRings(mol,cycles);	   
			
			if(test)
			{
				System.err.println("\nCycles found");
			}
			else
			{
				System.err.println("\nError - molecule failed in barkerSSSRFinder");				
			}
	 	} 
		catch (InvalidSmilesException ise) 
		{
		
		}
		 
    }
}
