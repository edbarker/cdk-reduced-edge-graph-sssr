package BarkerSSSR;

import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.templates.MoleculeFactory;

public class Test {
    
    public static void main(String[] args) {
        IMolecule mol = MoleculeFactory.makeAlphaPinene();
       boolean found = barkerSSSR.findRings(mol);
    }

}
