package barkerSSSRFinder;

import java.util.Iterator;
import java.util.List;

import org.openscience.cdk.*;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.ringsearch.SSSRFinder;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.tools.manipulator.RingSetManipulator;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IMoleculeSet;
import org.openscience.cdk.graph.PathTools;

/*
 * Algorithm to find the smallest set of smallest rings
 * 
 * Fragments a molecule in a particular way that allows cycles to be easily perceived
 * for most molecules. More complex ring systems can also be solved using more advanced heuristics.
 * The fragmentation gives rise to what I call a reduced edge graph.
 * 
 * This is faster than the original SSSR method because it processes simple ring systems far more
 * efficiently. 
 * 
 * written by Ed Barker June 2011
 * 
 */

public class barkerSSSRFinder {

	public static boolean findRings(IAtomContainer originalmol, IMoleculeSet cycles) throws Exception 
	{	
		IAtomContainer copymol = DefaultChemObjectBuilder.getInstance().newMolecule();

		// set up the ids for atoms and bonds - so we can find them later
		int id = 0;              

		for(IAtom atom : originalmol.atoms())
		{ 
			atom.setID(Integer.toString(id));
			id++;
		}
		
		// needed for the checksplitting optimisation
		id=0;
		for(IBond bond : originalmol.bonds())
		{ 
			bond.setID(Integer.toString(id));
			id++;
		}

		copymol.add(originalmol);

		int originalnumatoms=copymol.getAtomCount();
		int originalnumbonds=copymol.getBondCount();
		int cauchy =0;
		int numringsfound=0;
		int numcomponents = 0;

		stripTerminalAtoms(copymol);
		removeArtefacts(copymol);
	            	
		// need to partition into components first to calculate cauchy 
		if(ConnectivityChecker.isConnected(copymol))
		{
			if(copymol.getAtomCount()>0)
			{
				numcomponents = 1;
				cauchy = copymol.getBondCount() - copymol.getAtomCount() + numcomponents;                	
			}
			else
				cauchy=0;
		}
		else 
		{
			IMoleculeSet components = ConnectivityChecker.partitionIntoMolecules(copymol);
			numcomponents = components.getMoleculeCount();
			cauchy = copymol.getBondCount() - copymol.getAtomCount() + numcomponents;
		}

		if(cauchy>0)  
		{
			numringsfound=findSSSR
			(numcomponents,numringsfound,copymol,
					cauchy,originalnumatoms,cycles,originalmol,originalnumbonds);
		}

		if(cauchy==cycles.getAtomContainerCount())
		{
			return true; 		
		}
		else 
		{
			return false;
		}	
	}	
	
	private static void removeArtefacts(IAtomContainer newmol) {
		// deletes all atoms with no neighbours - these caused problems so needed resolving		
		
		int numbonds = 0;

       	for(IBond bond : newmol.bonds())
        {             	
       		IAtom a1 = bond.getAtom(0);
       		IAtom a2 = bond.getAtom(1);
       		
       		if((a1==null)||(a2==null))
       		{
   				newmol.removeBond(bond);       			
       		} 
        }      
       	
       	for(IAtom iatom1 : newmol.atoms())
        {             	
           	numbonds = newmol.getConnectedAtomsCount(iatom1);
   			if(numbonds==0)
   			{
   				newmol.removeAtom(iatom1);
   			}
        }      	
	}
	
	private static void stripTerminalAtoms (IAtomContainer mol)
	{
    	// find all one-connected atoms
    	// then remove neighbours iteratively until reach non-terminal atoms
		
		int numbonds = 0;
		IAtomContainer terminalatoms = 	DefaultChemObjectBuilder.getInstance().newAtomContainer();	
		
       	for(IAtom iatom1 : mol.atoms())
        {             	
           	numbonds = mol.getConnectedAtomsCount(iatom1);
   			if(numbonds==1)
   				terminalatoms.addAtom(iatom1);
        }

       	for(IAtom deleteatom : terminalatoms.atoms())
        {             			
   			numbonds=mol.getConnectedAtomsCount(deleteatom);
    		
      		IAtom iatom1 = deleteatom;
      		
   		   	while(numbonds==1)
           	{
           		List <IAtom> neighs = mol.getConnectedAtomsList(iatom1);
           		Iterator <IAtom> nextatom = neighs.iterator();
           		IAtom iatom2 = DefaultChemObjectBuilder.getInstance().newAtom();
           		while(nextatom.hasNext())
           		{
           			iatom2 = nextatom.next();
           		}

           		mol.removeBond(iatom1,iatom2);   
           		mol.removeAtom(iatom1);                     			

           		iatom1=iatom2;
       			numbonds=mol.getConnectedAtomsCount(iatom1);                   			           	
           	} 	
        }
	}
	
	private static int findSSSR(
			int numcomponents, int numringsfound, IAtomContainer mol, 
			int cauchy, 
			int originalnumatoms, 
			IMoleculeSet cycles,
			IAtomContainer originalmol, int originalnumbonds) {
			
		// Main function to find SSSR
		// fragments the molecule then finds simple rings
		// if unfound rings exist then do more complex ring finding

		int numfrags = 0;
	 	
		IAtomContainer cycle = DefaultChemObjectBuilder.getInstance().newMolecule();
			
		IMoleculeSet fragments = DefaultChemObjectBuilder.getInstance().newMoleculeSet();
		IAtomContainer fraggraph = DefaultChemObjectBuilder.getInstance().newMolecule();
		IAtomContainer origmol = DefaultChemObjectBuilder.getInstance().newMolecule();		

		origmol.add(mol);

		int [][] connectivity;
		connectivity = new int [originalnumatoms][2];

		int [][] cyclearray;
		cyclearray = new int [originalnumatoms*2][originalnumatoms*2]; 
		// can have > cauchy found and also > atoms! hence the size

		int [][] terminalatoms;
		terminalatoms = new int [originalnumatoms*2][3]; 
			
		int[][] atomfragmatrix; 
		atomfragmatrix = new int[originalnumatoms][originalnumatoms*2]; 
		// it is possible that the number frags can be more than the num original atoms
		// this happens in very dense graphs

		// optimisation - for checksplitting
		int [][] bondcyclearray;
		bondcyclearray=new int [originalnumbonds][originalnumbonds];
	    	
	    findConnectivities(connectivity,mol);
	    	
		while(mol.getAtomCount()>0)
		{
			// step one is to identify least connected atom
			// then make a cut to one of its bonds and chew away the open chains
			// this makes up the fragment and continues until nothing left in the molecule
		
			IAtom seedatom = findLeastConnectedAtom(connectivity,mol,originalnumatoms);
							
			IMolecule subgraph = DefaultChemObjectBuilder.getInstance().newMolecule();
				
			// now delete a bond - arbitrarily 
			List <IAtom> neighs = mol.getConnectedAtomsList(seedatom);
			Iterator <IAtom> nextatom = neighs.iterator();
			if(nextatom.hasNext())
			{
				IAtom iatom2 = nextatom.next();
				// need to do this first - just in case only one bond deleted
				subgraph.addAtom(seedatom);
				subgraph.addAtom(iatom2);
				subgraph.addBond(mol.getBond(seedatom,iatom2));

				// set up array for later
				atomfragmatrix[Integer.parseInt(seedatom.getID())][numfrags]=1;
				atomfragmatrix[Integer.parseInt(iatom2.getID())][numfrags]=1;
	         		
				mol.removeBond(seedatom,iatom2);   
				updateConnectivity(connectivity,seedatom,iatom2,mol);

				// now strip molecule of termini
				// and build a new molecule from what is removed - remembering to add deleted bond
				// and add a node to reduced edge molecule representing this subgraph
				stripMolecule(mol,subgraph,seedatom,iatom2,atomfragmatrix,numfrags,connectivity,terminalatoms);

				// now add the subgraph to our list of fragments
				fragments.addMolecule(subgraph);
	           		
				// now add the atom to a new fragment graph - edges added at last stage if it gets that far
				IAtom fragnode = DefaultChemObjectBuilder.getInstance().newAtom();
				fragnode.setID(Integer.toString(numfrags));

				// set atom type to C for default 
				fragnode.setAtomicNumber(6); 

				// check if subgraph contains a cycle
				if(checkCycle(subgraph)>0)
				{
					if(terminalatoms[numfrags][0]==1)
					{
						//Set node to type Ge = ring found in Type III fragment     
						fragnode.setAtomicNumber(32); 
	           		}
	           				
					else
					{
						// cyclic - so remove terminals - or rather number
						terminalatoms[numfrags][0]=0;
						
						 //Set node to Si = ring found - Type II fragment
						fragnode.setAtomicNumber(14);
	           		}
					cycle = addCycle(subgraph);
					cycles.addAtomContainer(cycle);
					addCycleArray(cycle,numringsfound,cyclearray,bondcyclearray);
					numringsfound++;
					if(numringsfound==cauchy)
					{   						
						return numringsfound;            				
	           		}
				}

				fraggraph.addAtom(fragnode);
				numfrags++;
	       	}
			else
			{
				mol.removeAtom(seedatom);
	       	}
		}

		// see if fragment has termini in the same fragment - if so resolve easy cycles
		int incommon=0;
		for(int y1=0;y1<numfrags;y1++)
		{
			for(int y2=y1+1;y2<numfrags;y2++)
			{
				IAtom frag1=fraggraph.getAtom(y1);
				IAtom frag2=fraggraph.getAtom(y2);        			
	      
				if((frag1.getAtomicNumber()==6)||(frag2.getAtomicNumber()==6))
				{
					incommon=countAtomsinCommon(atomfragmatrix,originalnumatoms,y1,y2,fragments.getAtomContainer(y1));
					if(incommon==2)
					{
						numringsfound+=resolveEasyCycle
						(fragments,y1,y2,cycles,fraggraph,originalmol,
								terminalatoms,cyclearray,numringsfound,bondcyclearray);
					}
					if((frag1.getAtomicNumber()==82)) // sets node as cycle solved = Pb
					{
						y2=numfrags;
	           		}
	       		}
	       	}
		}   

		if(numringsfound==cauchy)
			return numringsfound;        

		// if using to pre-process the cdk's SSSRFinder then the algorithm should end here
		// *******************************************************************************
		// return numringsfound;
		
		// next step is to find those unsolved linear fragments and 
		// find the shortest path in the molecule between the terminal atoms when fragment is removed

		removeArtefacts(origmol);
		
		numringsfound=resolveComplexCycles
		(numcomponents,terminalatoms,fraggraph,fragments,cycles,
				cauchy,origmol,originalmol,numringsfound,cyclearray,bondcyclearray,originalnumbonds);
		
		if(numringsfound==cauchy)
		{
			return numringsfound;            				
		}
			
		// now if it gets here it's because there the paths required to solve the unfound rings
		// are larger than the shortest paths which rarely happens or there are multiple shortest paths
		
		numringsfound=resolveRareCycles
		(numcomponents,terminalatoms,fraggraph,fragments,cycles,
				cauchy,origmol,originalmol,numringsfound,cyclearray,bondcyclearray,
				atomfragmatrix,originalnumatoms);
						
		return numringsfound;	    	
	}

	private static void findConnectivities(int[][] connectivity,
			IAtomContainer mol) {
		// generates connectivities for all atoms in molecule
		int connectivity_index;
		
		// first build up the list of connectivities
		for(IAtom atomit : mol.atoms())
		{
			connectivity_index=0;
			connectivity_index=findConnectivityIndex(mol.getConnectedAtomsCount(atomit));
			connectivity[Integer.parseInt(atomit.getID())][0]=connectivity_index; 

       		// do this otherwise unconnected or deleted atoms will be selected which is no good
       		if(connectivity_index==0)
       			connectivity_index=99999;
       		
			List <IAtom> neighs = mol.getConnectedAtomsList(atomit);
       		Iterator <IAtom> nextatom = neighs.iterator();
       		while(nextatom.hasNext())
       		{
       			IAtom iatom2 = nextatom.next();
       			connectivity_index+=findConnectivityIndex(mol.getConnectedAtomsCount(iatom2));
       		}
       		
       		connectivity[Integer.parseInt(atomit.getID())][1]=connectivity_index; 
       		
		}
	}	
		
	private static IAtom findLeastConnectedAtom(int[][] connectivity, IAtomContainer mol, int originalnumatoms) {
		// finds the atom least connected adapted from zamora paper JCICS 16 p40 1975
		// in zamora there is a multiplier of 64 but that is to find most connected atoms
		// removing this increases performance as it favours the least connected atom
		
		int leastconnectivity=9999999;
		IAtom leastatom=null;
		int leastconnected=mol.getAtomNumber(mol.getFirstAtom());
		for(IAtom atom : mol.atoms())
		{
			if(connectivity[Integer.parseInt(atom.getID())][1]<leastconnectivity)
			{
				leastconnectivity=connectivity[Integer.parseInt(atom.getID())][1];
				leastconnected=mol.getAtomNumber(atom);
			}
		}
		leastatom=mol.getAtom(leastconnected);
		return leastatom;
	}

	private static int findConnectivityIndex(int neighs) {
		if(neighs==2)
			return 1;
		if(neighs==3)
			return 8;
		if(neighs>=4)
			return 64;
		
		return 0;
	}

	private static void updateConnectivity(int[][] connectivity,
			IAtom seedatom, IAtom iatom2, IAtomContainer mol) {
		// when bond removed the neighbours' connectivity is updated

		int connectivity_index;
		int c1=0;
		int c2=0;
		int L1=0;
		int L2=0;
		
		L1=connectivity[Integer.parseInt(seedatom.getID())][1]-(connectivity[Integer.parseInt(seedatom.getID())][0]);
		connectivity_index=0;
		c1=connectivity[Integer.parseInt(seedatom.getID())][0];
		connectivity_index+=findConnectivityIndex(mol.getConnectedAtomsCount(seedatom));
		connectivity[Integer.parseInt(seedatom.getID())][0]=connectivity_index; // This is Ki		
		
		L2=connectivity[Integer.parseInt(iatom2.getID())][1]-(connectivity[Integer.parseInt(iatom2.getID())][0]);
		c2=connectivity[Integer.parseInt(iatom2.getID())][0];
		connectivity_index+=findConnectivityIndex(mol.getConnectedAtomsCount(iatom2));
		connectivity[Integer.parseInt(iatom2.getID())][0]=connectivity_index; // This is Ki
		
		connectivity[Integer.parseInt(seedatom.getID())][1]=connectivity[Integer.parseInt(seedatom.getID())][0]+L1; // Li
		connectivity[Integer.parseInt(iatom2.getID())][1]=connectivity[Integer.parseInt(iatom2.getID())][0]+L2; // Li

		List <IAtom> neighs = mol.getConnectedAtomsList(seedatom);
       	Iterator <IAtom> nextatom = neighs.iterator();
       	while(nextatom.hasNext())
       	{
       		IAtom iatom3 = nextatom.next();
       		connectivity[Integer.parseInt(iatom3.getID())][1]-=c1;
       		connectivity[Integer.parseInt(iatom3.getID())][1]+=connectivity[Integer.parseInt(seedatom.getID())][0];       		                                                 
       	}
       		
		List <IAtom> neighs2 = mol.getConnectedAtomsList(iatom2);
       	Iterator <IAtom> nextatom2 = neighs2.iterator();
       	while(nextatom.hasNext())
       	{
       		IAtom iatom4 = nextatom2.next();
       		connectivity[Integer.parseInt(iatom4.getID())][1]-=c2;
       		connectivity[Integer.parseInt(iatom4.getID())][1]+=connectivity[Integer.parseInt(iatom2.getID())][0];        		                                                 
       	}		
	}
	
	private static void stripMolecule(
			IAtomContainer mol, 
			IAtomContainer subgraph, 
			IAtom a, IAtom b, 
			int[][] atomfragmatrix, 
			int numfrags,
			int[][] connectivity,
			int[][] terminalatoms) {
		// same as strip terminals - except keep record of what is deleted
	
		int numbonds = 0;
		IAtomContainer deletionatoms = 	DefaultChemObjectBuilder.getInstance().newAtomContainer();	
		IAtom lastatom1 = 	DefaultChemObjectBuilder.getInstance().newAtom();	
		IAtom lastatom2 = 	DefaultChemObjectBuilder.getInstance().newAtom();	

		deletionatoms.addAtom(a);
		deletionatoms.addAtom(b);
		int numterminals=0;
		
		for(IAtom deleteatom : deletionatoms.atoms())
        {             	
   			numbonds=mol.getConnectedAtomsCount(deleteatom);
		
      		IAtom iatom1 = deleteatom;
      		
   		   	while(numbonds==1)
           	{
           		List <IAtom> neighs = mol.getConnectedAtomsList(iatom1);
           		Iterator <IAtom> nextatom = neighs.iterator();
           		IAtom iatom2 = DefaultChemObjectBuilder.getInstance().newAtom();
           		while(nextatom.hasNext())
           		{
           			iatom2 = nextatom.next();
           		}

       			if((iatom1!=a)&&(iatom1!=b))
       			{
       				subgraph.addAtom(iatom1);
       				if(iatom1.getID()!=null)
       					atomfragmatrix[Integer.parseInt(iatom1.getID())][numfrags]=1;
       			}
       			if((iatom2!=a)&&(iatom2!=b))
       			{
       				subgraph.addAtom(iatom2);
       				if(iatom2.getID()!=null)
       					atomfragmatrix[Integer.parseInt(iatom2.getID())][numfrags]=1;
       			}
       	
       			subgraph.addBond(mol.getBond(iatom1,iatom2));   			

           		mol.removeBond(iatom1,iatom2);   

           		if(iatom1.getID()!=null)
           		{
           			if(terminalatoms[numfrags][1]==Integer.parseInt(iatom1.getID()))
           			{
           				// then it has just seen an atom that is a terminal - ie a branching atom
           				// so reset in order to delete this and only have one terminal in these fragments
           				numterminals=0;
           			}
           		}
           		mol.removeAtom(iatom1);
           		
           		lastatom1=iatom1;
           		lastatom2=iatom2;
       
           		connectivity[Integer.parseInt(iatom1.getID())][1]=99999;
	
           		iatom1=iatom2;
       			numbonds=mol.getConnectedAtomsCount(iatom1);                   			           	
           	}
   		   	numterminals++;
   		   	terminalatoms[numfrags][0]=numterminals;
   		   	if(iatom1.getID()!=null)
   		   		terminalatoms[numfrags][numterminals]=Integer.parseInt(iatom1.getID());
        }
			
       	// efficient - update only for the last bond removed 
       	// and the first bond removed (done before come here)
       	// because everything else will have been deleted from the molecule and so can be ignored

       	if((lastatom1.getID()!=null)&&(lastatom2.getID()!=null))
      		updateConnectivity(connectivity,lastatom1,lastatom2,mol);	        
	}

	private static int checkCycle(IAtomContainer subgraph) {
		// simple cauchy to see if subgraph contains a cycle
		int subcauchy = 0;
		subcauchy = subgraph.getBondCount() - subgraph.getAtomCount() + 1; // only ever 1 component
		return subcauchy;
	}

	private static void addCycleArray(IAtomContainer cycle, int numringsfound,
			int[][] cyclearray, int[][] bondcyclearray) {
		// builds the array of cycles to make lookup a little faster
		// adds one when finds a cycle

		for(IAtom atom : cycle.atoms())
		{
			cyclearray[numringsfound][Integer.parseInt(atom.getID())]=1;
		}
		
		int count=0;
		for(IBond bond : cycle.bonds())
		{
			bondcyclearray[numringsfound][Integer.parseInt(bond.getID())]=1;
			count++;
		}
	}

	private static IAtomContainer addCycle(IAtomContainer subgraph) {
		// strip away any terminals and then copy what's left into the list of found cycles
		IAtomContainer newcycle = DefaultChemObjectBuilder.getInstance().newMolecule();
		newcycle.add(subgraph);
		
		stripTerminalAtoms(newcycle);
		
		return newcycle;
	}

	private static int countAtomsinCommon(int[][] atomfragmatrix, 
			int originalnumatoms, int i, int j, 
			IAtomContainer frag1) {
		// simply counts number of atoms shared in two fragments		
		
		int incommon=0;
		for(IAtom atom : frag1.atoms())
		{
			if((atomfragmatrix[Integer.parseInt(atom.getID())][j]==1))
				incommon++;
		}

		return incommon; 
	}

	private static int resolveEasyCycle(IMoleculeSet fragments, int i, int j,
			IMoleculeSet cycles,
			IAtomContainer fraggraph,
			IAtomContainer origmol, int[][] terminalatoms,
			int[][] cyclearray, int numringsfound, int[][] bondcyclearray) {

		// Finds simple 2-connected fragment rings 
		// - and also checks for splitting of rings already found 
					
		IAtomContainer origA = fragments.getMolecule(i);
		IAtomContainer origB = fragments.getMolecule(j);       						
		IAtomContainer cycle = DefaultChemObjectBuilder.getInstance().newMolecule();
		IAtomContainer fragA = DefaultChemObjectBuilder.getInstance().newMolecule();
		IAtomContainer fragB = DefaultChemObjectBuilder.getInstance().newMolecule();
		
		fragA.add(origA);
		fragB.add(origB);
		int splitflag=0;
		int foundring=0;
		
		if((checkCycle(fragA)==0)&&(checkCycle(fragB)==0))
       	{
       		// if both linear then this is easy to resolve
			cycle.add(fragA);
       		cycle.add(fragB); 
       		stripTerminalAtoms(cycle);
       		splitflag=checkSplitting(cycle,cycles,cyclearray,bondcyclearray);
       		
       		if(splitflag!=99)
       		{
       			numringsfound-=splitflag;
       		
       			if(checkUniqueCycle(cycle,cycles,cyclearray)==1)
       			{
       				cycles.addAtomContainer(cycle);	
       				addCycleArray(cycle,numringsfound,cyclearray,bondcyclearray);
       				foundring=1;

       				int pbfrag = findTerminals(i,j,fragA,fragB,terminalatoms,origmol);
       				fraggraph.getAtom(pbfrag).setAtomicNumber(82); //Pb - ie resolved linear
       			}
       		}
       		return foundring;
       	}
       	else // one must be a cycle - so need shortest path in this fragment
       	{
			if(checkCycle(fragA)==0)
			{
				IAtom linker1 = DefaultChemObjectBuilder.getInstance().newAtom();
				IAtom linker2 = DefaultChemObjectBuilder.getInstance().newAtom();

				if(terminalatoms[i][0]==2)
				{
					linker1=origmol.getAtom(terminalatoms[i][1]);
					linker2=origmol.getAtom(terminalatoms[i][2]);					
					
				}
				else
				{
					// never happens
					System.err.println("\nERROR - linear chain mismatch number = " + terminalatoms[i][0]);
				}

				List <IAtom> pathatoms=PathTools.getShortestPath(fragB,linker1,linker2);
				IAtomContainer shortestpath = DefaultChemObjectBuilder.getInstance().newMolecule();				

				// now build shortest path from fragB
				for(int h=0; h<pathatoms.size(); h++)
				{
					IAtom atom1=pathatoms.get(h);
								
					shortestpath.addAtom(atom1);
				}
				for(IBond bond1 : fragB.bonds())
				{
					if(shortestpath.contains(bond1.getAtom(0))&&(shortestpath.contains(bond1.getAtom(1))))
					{
						shortestpath.addBond(bond1);									
					}
				}
				shortestpath.add(shortestpath);
				cycle.add(shortestpath);
    			cycle.add(fragA);
    			//trim cycle
    			stripTerminalAtoms(cycle);
    			// check for splitting 
    			splitflag=checkSplitting(cycle,cycles,cyclearray,bondcyclearray);
           		
           		if(splitflag!=99)
           		{
           			numringsfound-=splitflag;
         			
           			if(checkUniqueCycle(cycle,cycles,cyclearray)==1)
           			{
           				cycles.addAtomContainer(cycle);
           				addCycleArray(cycle,numringsfound,cyclearray,bondcyclearray);
           				foundring=1;

           				fraggraph.getAtom(i).setAtomicNumber(82); //Pb - ie resolved linear
           			}
           		}
       			return foundring;          								
       		}
       		else
       		{
    			IAtom linker1 = DefaultChemObjectBuilder.getInstance().newAtom();
       			IAtom linker2 = DefaultChemObjectBuilder.getInstance().newAtom();
				if(terminalatoms[j][0]==2)
				{
					linker1=origmol.getAtom(terminalatoms[j][1]);
					linker2=origmol.getAtom(terminalatoms[j][2]);					
					
				}
				else
				{
					// never comes here
					System.err.println("\nERROR - linear chain mismatch number = " + terminalatoms[j][0]);
				}

       			// now find shortest path in other fragment
       			List <IAtom> pathatoms=PathTools.getShortestPath(fragA,linker1,linker2);
				IAtomContainer shortestpath = DefaultChemObjectBuilder.getInstance().newMolecule();
				
				for(int h=0; h<pathatoms.size(); h++)
				{
					IAtom atom1=pathatoms.get(h);	
					shortestpath.addAtom(atom1);
				}
				for(IBond bond1 : fragA.bonds())
				{
					if(shortestpath.contains(bond1.getAtom(0))&&(shortestpath.contains(bond1.getAtom(1))))
					{
						shortestpath.addBond(bond1);
					}					
				}
				shortestpath.add(shortestpath);
       			cycle.add(fragB);
       			cycle.add(shortestpath);
    			stripTerminalAtoms(cycle);
    			splitflag=checkSplitting(cycle,cycles,cyclearray,bondcyclearray);
           		
           		if(splitflag!=99)
           		{
           			numringsfound-=splitflag;

           			if(checkUniqueCycle(cycle,cycles,cyclearray)==1)
           			{
           				cycles.addAtomContainer(cycle);
           				addCycleArray(cycle,numringsfound,cyclearray,bondcyclearray);
           				foundring=1;

           				fraggraph.getAtom(j).setAtomicNumber(82); //Pb - ie resolved linear
           			}
           		}
         	    return foundring;
			}			
		}
	}

	private static int checkSplitting(IAtomContainer newcycle,
			IMoleculeSet cycles, int[][] cyclearray, int[][] bondcyclearray) {
		// checks that a cycle is not part of a larger ring found previously
		// if so the larger ring is divided into two rings and the new ring tested for uniqueness
		
		int deleted=0;
		int index=0;
		
		if(checkUniqueCycle(newcycle,cycles,cyclearray)==0)
		{		
			return 99;
		}	

		for(IAtomContainer cycle : cycles.molecules())
		{
			if(newcycle.getAtomCount()<cycle.getAtomCount())
			{
				// find bonds in common
				IAtomContainer commonbonds =  DefaultChemObjectBuilder.getInstance().newMolecule();
				IAtomContainer exclusivebonds = DefaultChemObjectBuilder.getInstance().newMolecule();
				IAtomContainer oldcycle = DefaultChemObjectBuilder.getInstance().newMolecule();				
				for(IBond newbond : newcycle.bonds())
				{
					if(bondcyclearray[index][Integer.parseInt(newbond.getID())]==1)
					{
						commonbonds.addBond(newbond);
					}
					else
					{
						exclusivebonds.addBond(newbond);
					}
				}
				
				// now check if the number of bonds is greater than number exclusive bonds
				if(commonbonds.getBondCount()>exclusivebonds.getBondCount())
				{
					oldcycle.add(cycle);
					
					for(IBond deletebond : commonbonds.bonds())
					{
						oldcycle.removeAtom(deletebond.getAtom(0));
						oldcycle.removeAtom(deletebond.getAtom(1));						
						oldcycle.removeBond(deletebond);						
					}

					for(IBond addbond : exclusivebonds.bonds())
					{
						oldcycle.addAtom(addbond.getAtom(0));
						oldcycle.addAtom(addbond.getAtom(1));						
						oldcycle.addBond(addbond);						
					}
										
					if(checkUniqueCycle(oldcycle,cycles,cyclearray)==0)
					{
						cycles.removeAtomContainer(cycle);
						removeCycleArray(cycle,index,cyclearray,bondcyclearray);					
						deleted++;
					}
					else
					{
						removeCycleArray(cycle,index,cyclearray,bondcyclearray);
						addCycleArray(oldcycle,index,cyclearray,bondcyclearray);
						replaceCycle(cycle,oldcycle);
					}
				}
			}
			
			index++;
		}	
		
	return deleted;
	
	}

	private static int checkUniqueCycle(IAtomContainer newcycle,
			IMoleculeSet cycles, int[][] cyclearray) {
		// simply check that this is not a cycle already found
	
		IAtomContainer cycle = null;

		for(int i=0;i<cycles.getAtomContainerCount();i++)
		{
			int incommon=0;
			cycle=cycles.getAtomContainer(i);
			if(cycle.getAtomCount()==newcycle.getAtomCount())
			{
				for(IAtom newatom : newcycle.atoms())
				{
					if(cyclearray[i][Integer.parseInt(newatom.getID())]==1)
						incommon++;
					else
						break;
				}
				if(incommon==newcycle.getAtomCount())
				{
					return 0;
				}
			}
		}
		
		return 1;
	}
	
	private static void replaceCycle(IAtomContainer cycle,
			IAtomContainer newcycle) {
		// need to replace the old split cycle with the new one

		IAtomContainer copy=DefaultChemObjectBuilder.getInstance().newMolecule();
		copy.add(cycle);
		cycle.remove(copy);
		cycle.add(newcycle);
		
	}

	private static void removeCycleArray(IAtomContainer cycle, int index,
			int[][] cyclearray, int[][] bondcyclearray) {
		// clear array out for replacing

		for(IAtom atom : cycle.atoms())
		{
			cyclearray[index][Integer.parseInt(atom.getID())]=0;
		}
		int count=0;
		for(IBond bond : cycle.bonds())
		{
			bondcyclearray[index][Integer.parseInt(bond.getID())]=0;
			count++;
		}
	}

	
	private static int findTerminals(int i, int j, IAtomContainer fragA,
			IAtomContainer fragB, int[][] terminalatoms, IAtomContainer origmol) {
		// finds which of the fragments has its terminal atoms in the other
		// if both the same - send the first

		IAtom linker1=origmol.getAtom(terminalatoms[i][1]);
		IAtom linker2=origmol.getAtom(terminalatoms[i][2]);					

		if((fragB.contains(linker1))&&(fragB.contains(linker2)))
		{
			return i;
		}
		else
		{
			return j;
		}		
	}
	
	private static int resolveComplexCycles(int numcomponents, int[][] terminalatoms, IAtomContainer fraggraph,
			IMoleculeSet fragments, IMoleculeSet cycles, int cauchy,
			IAtomContainer copymol, IAtomContainer originalmol, int numringsfound, int[][] cyclearray,
			int[][] bondcyclearray, int originalnumbonds) {
			
		// simply remove the fragment which is linear from mol
		// find shortest path to linkers
		// if closes ring then check ring formed

		int splitflag=0;
		int fragment=0;

		for(IAtom node : fraggraph.atoms())
		{
			if(node.getAtomicNumber()==6)
			{
				IAtomContainer cycle = DefaultChemObjectBuilder.getInstance().newMolecule();		
				IAtomContainer newmol = DefaultChemObjectBuilder.getInstance().newMolecule();		
				newmol.add(copymol);
					
				IAtomContainer deletefrag = DefaultChemObjectBuilder.getInstance().newMolecule();					
				deletefrag.add(fragments.getAtomContainer(fragment));
					
				IAtom linker1 = DefaultChemObjectBuilder.getInstance().newAtom();
				IAtom linker2 = DefaultChemObjectBuilder.getInstance().newAtom();
	      							
				// find linker atoms in fragment
				// use terminalatoms array
				// as must be linear must be two terminal atoms
				if(terminalatoms[fragment][0]==2)
				{
					linker1=originalmol.getAtom(terminalatoms[fragment][1]);
					linker2=originalmol.getAtom(terminalatoms[fragment][2]);						
				}
				//need to remove everything except the linker atoms
				for(IAtom fragatom : deletefrag.atoms())
				{
					// need to remove the bonds from the fragment only
					for(IBond fragbond : deletefrag.getConnectedBondsList(fragatom))
					{
						newmol.removeBond(fragbond);	
					}
				}

				removeArtefacts(newmol);
			
				// now find shortest path and add this to fragment
				// fails if linear creates disconnected fragments so test for this
				IMoleculeSet components = ConnectivityChecker.partitionIntoMolecules(newmol);
				for(IAtomContainer comp : components.molecules())
				{
					removeArtefacts(comp);
					if((comp.getAtomCount()>0)&&(comp.contains(linker1))&&(comp.contains(linker2)))
					{
						// now find shortest path in other fragment
						List <IAtom> pathatoms=PathTools.getShortestPath(newmol,linker1,linker2);
						IAtomContainer shortestpath = DefaultChemObjectBuilder.getInstance().newMolecule();					    				
						// now build shortest path from fragB
	    				
						for(int h=0; h<pathatoms.size(); h++)
						{
							IAtom atom1=pathatoms.get(h);
							shortestpath.addAtom(atom1);
	                	}
						// now connect the atoms
						for(IBond bond1 : newmol.bonds())
						{
							if(shortestpath.contains(bond1.getAtom(0))&&(shortestpath.contains(bond1.getAtom(1))))
							{
								shortestpath.addBond(bond1);									
	                		}
	                	}
						shortestpath.add(shortestpath);

						cycle.add(shortestpath);
						cycle.add(deletefrag);
						stripTerminalAtoms(cycle);
						// check for splitting here
						splitflag=0; 	      		

						splitflag=checkSplittingComplex
						(cycle,cycles,cyclearray,numringsfound,bondcyclearray,originalnumbonds);
	                				
	                	if(splitflag!=99)
	                	{
	                		if(splitflag==999)
	                		{
	                			numringsfound++;
	                		}
	                		else
	                		{
	                			numringsfound-=splitflag;
	       	
	                			if(checkUniqueCycle(cycle,cycles,cyclearray)==1)
	                			{
	                				cycles.addAtomContainer(cycle);
	                				addCycleArray(cycle,numringsfound,cyclearray,bondcyclearray);
	                				numringsfound++;
	                			}
	                		}
	                		if(numringsfound==cauchy)
	                			return numringsfound;
	                	
	                			// this should flag where found ring
	                			node.setAtomicNumber(82);               			
                		}	   					
	       			}
	    		}        								
			}
			fragment++;
		}
					
		return numringsfound;
	}		
	
	private static int resolveRareCycles(int numcomponents,
			int[][] terminalatoms, IAtomContainer fraggraph,
			IMoleculeSet fragments, IMoleculeSet cycles, int cauchy,
			IAtomContainer copymol, IAtomContainer originalmol,
			int numringsfound, int[][] cyclearray, int[][] bondcyclearray,
			int[][] atomfragmatrix, int originalnumatoms) {
		// to finish then we need to find sssr in fraggraph using cdk's sssrfinder itself
		// build and test the cycle formed
		// to build this we need to find shortest paths to linkers for each added fragment

		// we haven't built the fragraph connections yet so do this
		buildfraggraph(fraggraph,terminalatoms,atomfragmatrix,fragments,originalnumatoms);
			
		IRingSet fragRings = new SSSRFinder(fraggraph).findSSSR();     		
		List<IAtomContainer> ringlist = RingSetManipulator.getAllAtomContainers(fragRings);

		IAtomContainer sssrringatoms = DefaultChemObjectBuilder.getInstance().newMolecule();		
		for (int j = 0; j < fragRings.getAtomContainerCount(); j++) 
		{
			sssrringatoms = ringlist.get(j);	  			
			// build up the ring from the fragments
					
			int [] linkingatoms;
			linkingatoms = new int [sssrringatoms.getAtomCount()];
			  				
			int [] fraglist;
			fraglist = new int [sssrringatoms.getAtomCount()];		  				
			  				
			// first find the linkers between nodes in cycle
			int count=0;
			for(IAtom fragnode : sssrringatoms.atoms())
			{
				fraglist[count]=Integer.parseInt(fragnode.getID());
				count++;		  							  						
			}
			count=0;
			int y1,y2;
			for(int k=0;k<sssrringatoms.getAtomCount();k++)
			{		  							  			
				y1=fraglist[k];
						
				if(k==(sssrringatoms.getAtomCount()-1))
					y2=fraglist[0];
			  										
				else
					y2=fraglist[k+1];
			  					
				linkingatoms[count]=findAtomsinCommon(atomfragmatrix,originalnumatoms,y1,y2,fragments.getAtomContainer(y1));
			  	 		      		
				count++;
			}
			  				
			// now trace the shortest paths between these atoms in the fragments
			// really should ignore k1,3 interchange but as this doesn't make a ring forget it
			IAtomContainer cycle =  DefaultChemObjectBuilder.getInstance().newMolecule();	
			for(int i=0;i<count;i++)
			{
				IAtomContainer fragment =  DefaultChemObjectBuilder.getInstance().newMolecule();	
					
				y1=linkingatoms[i];		
					
				if(i==(count-1))
				{
					y2=linkingatoms[0];
					fragment=fragments.getAtomContainer(fraglist[0]);
				}
			  										
				else
				{
					y2=linkingatoms[i+1];				
					fragment=fragments.getAtomContainer(fraglist[i+1]);
				}		

				IAtom linker1 = originalmol.getAtom(y1);
				IAtom linker2 = originalmol.getAtom(y2);				      	      
				// now find shortest path in other fragment
					
				if(linker1!=linker2)					
				{
					List <IAtom> pathatoms=PathTools.getShortestPath(fragment,linker1,linker2);
								
					IAtomContainer shortestpath = DefaultChemObjectBuilder.getInstance().newMolecule();					
					// now build shortest path from fragB
					for(int h=0; h<pathatoms.size(); h++)
					{
						IAtom atom1=pathatoms.get(h);		
						shortestpath.addAtom(atom1);
					}
					// now connect the atoms
					for(IBond bond1 : fragment.bonds())
					{
						if(shortestpath.contains(bond1.getAtom(0))&&(shortestpath.contains(bond1.getAtom(1))))
						{
							shortestpath.addBond(bond1);									
						}
					}
					shortestpath.add(shortestpath);
					cycle.add(shortestpath);
				
				}
			}
			//trim cycle
			stripTerminalAtoms(cycle);
	    		
			int splitflag;
			
			if(cycle.getAtomCount()>0)
			{
				splitflag=checkSplitting(cycle,cycles,cyclearray,bondcyclearray);
	    	}
	    		
			else
				splitflag=99;
	    		
	    	if(splitflag!=99)
	    	{    		
	    		if(splitflag==999)
	    		{	
	    			numringsfound++;
	    		}
	    		else
	    		{
	    			numringsfound-=splitflag;
				
	    			if(checkUniqueCycle(cycle,cycles,cyclearray)==1)
	    			{
	    				cycles.addAtomContainer(cycle);
	    				addCycleArray(cycle,numringsfound,cyclearray,bondcyclearray);
	    				numringsfound++;
	    			}
	    		}
	    		if(numringsfound==cauchy)
	    			return numringsfound;				
	    	}
		}
			  					
		return numringsfound;
	}
		
	private static void buildfraggraph(IAtomContainer fraggraph,
			int[][] terminalatoms, int[][] atomfragmatrix,
			IMoleculeSet fragments, int originalnumatoms) {
		// makes the connections between the fraggraph nodes as we will need these now
			
		for(int i=0;i<(fragments.getAtomContainerCount()-1);i++)
		{	
			for(int j=i+1;j<fragments.getAtomContainerCount();j++)
			{	
				for(int k=0;k<originalnumatoms;k++)
				{
					if((atomfragmatrix[k][i]==1)&&(atomfragmatrix[k][j]==1))
					{
						IBond bond = DefaultChemObjectBuilder.getInstance().newBond();	
							
						IAtom node1 = DefaultChemObjectBuilder.getInstance().newAtom();	
						IAtom node2 = DefaultChemObjectBuilder.getInstance().newAtom();	
							
						node1=fraggraph.getAtom(i);
						node2=fraggraph.getAtom(j);
						
						bond.setAtom(node1,0);
						bond.setAtom(node2,1);
													
						fraggraph.addBond(bond);
							
						k=originalnumatoms;
					}			
				}		
			}					
		}			
	}

	private static int findAtomsinCommon(int[][] atomfragmatrix, 
			int originalnumatoms, int i, int j, 
			IAtomContainer frag1) {
		// simply counts number of atoms shared in two fragments

		for(IAtom atom : frag1.atoms())
		{
			if((atomfragmatrix[Integer.parseInt(atom.getID())][j]==1))
				return(Integer.parseInt(atom.getID()));
		}

		return -1; //should not come here
	}		
		
	private static int checkSplittingComplex(IAtomContainer newcycle,
			IMoleculeSet cycles, int[][] cyclearray, int numringsfound,
			int[][] bondcyclearray, int originalnumbonds) {
		// does check splitting - but this is where the query cycle can be larger than the test cycle
		
		int deleted=0;
		int index=0;
			
		if(checkUniqueCycle(newcycle,cycles,cyclearray)==0)
		{		
			return 99;
		}	
			
		for(IAtomContainer cycle : cycles.molecules())
		{
			if(newcycle.getAtomCount()<cycle.getAtomCount())
			{
				// now find bonds in common
				IAtomContainer commonbonds =  DefaultChemObjectBuilder.getInstance().newMolecule();
				IAtomContainer exclusivebonds = DefaultChemObjectBuilder.getInstance().newMolecule();
				IAtomContainer oldcycle = DefaultChemObjectBuilder.getInstance().newMolecule();				
					
				for(IBond newbond : newcycle.bonds())
				{
					if(bondcyclearray[index][Integer.parseInt(newbond.getID())]==1)
					{
						commonbonds.addBond(newbond);
					}
					else
					{
						exclusivebonds.addBond(newbond);
					}
				}

				// now check if the number of bonds is greater than number exclusive bonds					
				if(commonbonds.getBondCount()>exclusivebonds.getBondCount())
				{
					// then update the old cycle by removing commonbonds
					oldcycle.add(cycle);

					for(IBond deletebond : commonbonds.bonds())
					{
						oldcycle.removeAtom(deletebond.getAtom(0));
						oldcycle.removeAtom(deletebond.getAtom(1));						
						oldcycle.removeBond(deletebond);						
					}
					for(IBond addbond : exclusivebonds.bonds())
					{
						oldcycle.addAtom(addbond.getAtom(0));
						oldcycle.addAtom(addbond.getAtom(1));						
						oldcycle.addBond(addbond);						
					}			
					if(checkUniqueCycle(oldcycle,cycles,cyclearray)==0)
					{
						cycles.removeAtomContainer(cycle);
						removeCycleArray(cycle,index,cyclearray,bondcyclearray);					
						deleted++;
					}
					else
					{
						removeCycleArray(cycle,index,cyclearray,bondcyclearray);
						addCycleArray(oldcycle,index,cyclearray,bondcyclearray);
						replaceCycle(cycle,oldcycle);
					}						
				}	
			}
				
			else //the ring is bigger
			{			
				IAtomContainer commonbonds =  DefaultChemObjectBuilder.getInstance().newMolecule();
				IAtomContainer exclusivebonds = DefaultChemObjectBuilder.getInstance().newMolecule();
				IAtomContainer oldcycle = DefaultChemObjectBuilder.getInstance().newMolecule();				

				int [] commonbondarray;
				commonbondarray=new int [originalnumbonds];
				int id;
				for(IBond newbond : newcycle.bonds())
				{
					id = Integer.parseInt(newbond.getID());
					if(bondcyclearray[index][id]==1)
					{
						commonbonds.addBond(newbond);
						commonbondarray[id]=1;
					}					
				}
	
				for(IBond bond : cycle.bonds())
				{
					if(commonbondarray[Integer.parseInt(bond.getID())]==0)
					{
						exclusivebonds.addBond(bond);							
					}
				}				

				if(commonbonds.getBondCount()>exclusivebonds.getBondCount())
				{
					oldcycle.add(newcycle);
					for(IBond deletebond : commonbonds.bonds())
					{
						oldcycle.removeAtom(deletebond.getAtom(0));
						oldcycle.removeAtom(deletebond.getAtom(1));						
						oldcycle.removeBond(deletebond);						
					}
					for(IBond addbond : exclusivebonds.bonds())
					{
						oldcycle.addAtom(addbond.getAtom(0));
						oldcycle.addAtom(addbond.getAtom(1));						
						oldcycle.addBond(addbond);						
					}		
					if(checkUniqueCycle(oldcycle,cycles,cyclearray)==0)
					{
						return 99;
					}
					else
					{					
						newcycle=oldcycle;
						cycles.addAtomContainer(newcycle);
						addCycleArray(newcycle,numringsfound,cyclearray,bondcyclearray); 				
						return 999;
					}
				}						
			}
			index++;
		}
		
		return deleted;
	}
	
}

	
	
	

