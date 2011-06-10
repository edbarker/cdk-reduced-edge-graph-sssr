package BarkerSSSR;

import java.io.*;
import java.util.Iterator;
import java.util.List;

import org.openscience.cdk.*;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.io.iterator.IteratingSMILESReader;
import org.openscience.cdk.ringsearch.SSSRFinder;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.tools.manipulator.RingSetManipulator;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IMoleculeSet;
import org.openscience.cdk.graph.PathTools;

/**
 * Pre-processes a molecule to find the Smallest Set of Smallest Rings.
 * If it cannot find this set at the end of the process then the molecule
 * is deemed too complex and so passes to the full SSSR algorithm.  
 * This helps to speed up the full SSSR algorithm.
 * This is an implementation of an algorithm 
 * by Edward Barker
 * 
 *
 */

public class barkerSSSR {

	public static boolean findRings(IAtomContainer originalmol, IMoleculeSet cycles) throws Exception 
    	{
        	IAtomContainer copymol = DefaultChemObjectBuilder.getInstance().newMolecule();
		
        	// set id for each atom so can keep track of them when atoms are deleted
        	int id = 0;
        	IAtom atom = null;
        	for(int k=0; k<originalmol.getAtomCount(); k++)
        	{ 
        		atom = originalmol.getAtom(k);
        		atom.setID(Integer.toString(id));
        	    	id++;
        	}
               
		int numcomponents = 0;       	 	
        	copymol.add(originalmol);
            		
        	int originalnumatoms=copymol.getAtomCount();
       		int cauchy=0;
       		int numringsfound=0;

        	// need to find number of components to calculate cauchy 
        	if(ConnectivityChecker.isConnected(copymol))
        	{
        		numcomponents = 1;
	        }
	      	else 
	       	{
	       		IMoleculeSet components = ConnectivityChecker.partitionIntoMolecules(copymol);
	       		numcomponents = components.getMoleculeCount();
	       	}
        
	       	cauchy = copymol.getBondCount() - copymol.getAtomCount() + numcomponents;
  
	       	if(cauchy>0)
	       	{      	
	          	stripTerminalAtoms(copymol);
	      		numringsfound=reduceMolecule(copymol,cauchy,originalnumatoms,cycles,originalmol);
	       	}
       	   
	        if(cauchy==numringsfound)
	        {
	        	return true;
	        }
	        else
	        {
	        	return false;
	        }        
    	} 
    
	private static int reduceMolecule(
			IAtomContainer mol, 
			int cauchy, 
			int originalnumatoms, 
			IMoleculeSet cycles,
			IAtomContainer originalmol) {
		
		// Main function to reduce the molecule to the reduced edge graph
		// finds the least connected atom and removes a bond from this atom
		// then removes all acyclic branches and reconnects what's removed into a fragment
		// then this fragment is represented as a node in a reduced edge graph
		// then solves simple cycles - Type II and Type III
		// if cauchy not found - then solve all linear fragments
		// if still not found cauchy - then stage one ends and give up on these molecules

		int numfrags = 0;
	 	int numringsfound=0;
	 	
		IMoleculeSet fragments = DefaultChemObjectBuilder.getInstance().newMoleculeSet();
		IAtomContainer fraggraph = DefaultChemObjectBuilder.getInstance().newMolecule();
		IAtomContainer origmol = DefaultChemObjectBuilder.getInstance().newMolecule();		
	
		origmol.add(mol);
	
		int [][] cyclearray;
		cyclearray = new int [originalnumatoms][originalnumatoms]; // can have > cauchy found
		
		// where terminalatoms[fragnode][0]=num terminal atoms
		//  terminalatoms[fragnode][1]=terminal atom 1 etc. 
		int [][] terminalatoms;
		terminalatoms = new int [originalnumatoms*2][3]; 
		
		int[][] atomfragmatrix; 
    		atomfragmatrix = new int[originalnumatoms][originalnumatoms*2]; 
    		// it is possible that the number frags can be more than the num original atoms
    	
    		numringsfound=findFragments
			(numfrags,mol,originalnumatoms,fraggraph,
			cycles,fragments,atomfragmatrix,terminalatoms,cyclearray);

    		if(numringsfound==cauchy)
    			return numringsfound;
    	
    		numringsfound=findEasyRings
    			(numfrags,fraggraph,atomfragmatrix,originalnumatoms,
    			fragments,cyclearray,numringsfound,terminalatoms,originalmol,cycles);
 
		return numringsfound;	    	
	}

	private static int findEasyRings(int numfrags, IAtomContainer fraggraph,
			int[][] atomfragmatrix, int originalnumatoms,
			IMoleculeSet fragments, int[][] cyclearray, int numringsfound,
			int[][] terminalatoms, IAtomContainer originalmol,
			IMoleculeSet cycles) {
			   	
		// finds where a linear fragment has two terminal atoms connecting the same fragment
		// then resolves these cycles
		
	    	int incommon=0;
	    	for(int y1=0;y1<numfrags;y1++)
	    	{
	        	for(int y2=y1+1;y2<numfrags;y2++)
	        	{
	       			IAtom frag1=fraggraph.getAtom(y1);
	       			IAtom frag2=fraggraph.getAtom(y2);        			
      
	       			if((frag1.getAtomicNumber()==6)||(frag2.getAtomicNumber()==6))
	       			{
	       				// then flagged as linear
	       				incommon=countAtomsinCommon(atomfragmatrix,originalnumatoms,
       						y1,y2,fragments.getAtomContainer(y1));
	       				if(incommon==2)
	       				{
  						numringsfound+=resolveEasyCycle
   							(fragments,y1,y2,cycles,fraggraph,originalmol,
							terminalatoms,cyclearray,numringsfound);

	      				}
	       			}
	       		}
	       	}   
	    	    	
		return numringsfound;
	}

	private static int resolveEasyCycle(IMoleculeSet fragments, int i, int j,
			IMoleculeSet cycles,
			IAtomContainer fraggraph,
			IAtomContainer origmol, int[][] terminalatoms,
			int[][] cyclearray, int numringsfound) {
		// Finds cycle made from two fragments
		// also checks for splitting of rings already found and if cycle found is unique 

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
	       		splitflag=checkSplittingComplex(cycle,cycles,cyclearray);
       		
	       		if(splitflag!=99)
	       		{
	       			numringsfound-=splitflag;
       		
	       			if(checkUniqueCycle(cycle,cycles,cyclearray)==1)
	       			{
	       				cycles.addAtomContainer(cycle);	
	       				addCycleArray(cycle,numringsfound,cyclearray);
	       				foundring=1;
	       			}
	       		}
	       		return numringsfound;
	       	}
	       	else // one must be a cycle - so need shortest path
	       	{
	       		// find linear fragment
	       		// find linkers
	       		// find shortest path in other fragment - add this path to cycle
	       		// add linear frag to cycle
			if(checkCycle(fragA)==0)
			{
				numringsfound+=resolveCyclicFragment(i,j,numringsfound,fragA,fragB,
						terminalatoms,origmol,cycles,cyclearray,fraggraph);
      								
       			}
       			else
       			{
				numringsfound+=resolveCyclicFragment(j,i,numringsfound,fragB,fragA,
						terminalatoms,origmol,cycles,cyclearray,fraggraph);
	
       			}
			return numringsfound;
		}
	}
	
	private static int resolveCyclicFragment(
			int linearfragmentindex, 
			int cyclefragmentindex, int numringsfound, IAtomContainer linearfrag,
			IAtomContainer cyclicfrag, int[][] terminalatoms,
			IAtomContainer origmol, IMoleculeSet cycles, int[][] cyclearray,
			IAtomContainer fraggraph) {
		// finds the shortest path between linkers and adds this to create a ring
		
		IAtom linker1 = DefaultChemObjectBuilder.getInstance().newAtom();
		IAtom linker2 = DefaultChemObjectBuilder.getInstance().newAtom();
		IAtomContainer cycle = DefaultChemObjectBuilder.getInstance().newMolecule();
	
		int foundring=0;
		int splitflag=0;
		
		if(terminalatoms[linearfragmentindex][0]==2)
		{
			linker1=origmol.getAtom(terminalatoms[linearfragmentindex][1]);
			linker2=origmol.getAtom(terminalatoms[linearfragmentindex][2]);
						
		}
		else
		{
			//never happens
			//System.err.println("\nERROR - linear chain mismatch number = " 
			//		+ terminalatoms[linearfragmentindex][0]);
		}

		// now find shortest path in cyclic fragment
		List <IAtom> pathatoms=PathTools.getShortestPath(cyclicfrag,linker1,linker2);
		IAtomContainer shortestpath = DefaultChemObjectBuilder.getInstance().newMolecule();
				
		for(int h=0; h<pathatoms.size(); h++)
		{
			IAtom atom1=pathatoms.get(h);
			shortestpath.addAtom(atom1);
		}
		// now connect the atoms
		for(IBond bond1 : cyclicfrag.bonds())
		{
			if(shortestpath.contains(bond1.getAtom(0))&&(shortestpath.contains(bond1.getAtom(1))))
			{
				shortestpath.addBond(bond1);
									
			}
		}
		shortestpath.add(shortestpath);
		cycle.add(shortestpath);
		cycle.add(linearfrag);
		//trim cycle
		stripTerminalAtoms(cycle);
		// check for splitting here
		
		// not sure about this - should be checksplitting?**************************************
		if(linearfrag.getAtomCount()<pathatoms.size())
		{
			if(checkBranching(cyclefragmentindex,fraggraph)==1)
			{
				// if a branching cycle then need to check that linkers are in the cycle
				if(checkLinkers(cyclicfrag,linker1,linker2)==1)
					splitCycle(cycles,cyclicfrag,linearfrag,shortestpath);	    						
			}
			else
				splitCycle(cycles,cyclicfrag,linearfrag,shortestpath);	
		}

		splitflag=checkSplittingComplex(cycle,cycles,cyclearray);
   		
		// ****************************************************************************************
		
   		if(splitflag!=99)
   		{
   			numringsfound-=splitflag;
 			
   			if(checkUniqueCycle(cycle,cycles,cyclearray)==1)
   			{
   				cycles.addAtomContainer(cycle);
   				addCycleArray(cycle,numringsfound,cyclearray);
   				foundring=1;
   			}
   		}
		return foundring;    
		
	}

	private static int checkSplittingComplex(IAtomContainer newcycle,
			IMoleculeSet cycles, int[][] cyclearray) {
		// does check splitting - but this is where only have the cycle info
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
				
				int common;
				for(IBond newbond : newcycle.bonds())
				{
					common=0;
					for(IBond bond : cycle.bonds())
					{
						if(bond==newbond)
						{
							commonbonds.addBond(bond);
							common=1;
						}
					}			
					// if the newbond not in other cycle then add to exclusivebonds
					if(common==0)
					{
						exclusivebonds.addBond(newbond);
					}
				}
			
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
									
						deleted++;
					}
					else
					{
						removeCycleArray(cycle,index,cyclearray);
						addCycleArray(oldcycle,index,cyclearray);
						replaceCycle(cycle,oldcycle);
					}
					
				}	
				
			}
			
			index++;
		}		
		return deleted;
	}

	private static int countAtomsinCommon(int[][] atomfragmatrix, 
			int originalnumatoms, int i, int j, 
			IAtomContainer frag1) {
		// counts number of atoms shared in two fragments
		//optimised version - use fragatommatrix
		
		int incommon=0;
		
		for(IAtom atom : frag1.atoms())
		{
			if((atomfragmatrix[Integer.parseInt(atom.getID())][j]==1))
				incommon++;
		}

		return incommon; 
	}

	private static void removeCycleArray(IAtomContainer cycle, int index,
			int[][] cyclearray) {
		// clear array out for relacing
		for(IAtom atom : cycle.atoms())
		{
			cyclearray[index][Integer.parseInt(atom.getID())]=0;
		}
	}
	
	private static void replaceCycle(IAtomContainer cycle,
			IAtomContainer newcycle) {
		// need to replace the old split cycle with the new one
		IAtomContainer copy=DefaultChemObjectBuilder.getInstance().newMolecule();
		copy.add(cycle);
		cycle.remove(copy);
		cycle.add(newcycle);		
	}

	private static void splitCycle(
			IMoleculeSet cycles,
			IAtomContainer cyclefragment, 
			IAtomContainer smallpath, 
			IAtomContainer largepath) {
		// split the cycle previously found by removing largepath and replace with smallpath
		IAtomContainer newcycle = DefaultChemObjectBuilder.getInstance().newMolecule();
		stripTerminalAtoms(newcycle); //justincase branching type
		for(IAtomContainer cycle : cycles.molecules())
		{
			if(newcycle == cycle) 
			{
				// now remove largepath and replace with small path
				cycle.remove(largepath);
				cycle.add(smallpath);			
			}
		}
	}
	
	private static int checkBranching(int j, IAtomContainer fraggraph) {
		// simple check to see if has one connected atom
		// optimised - just label fragnode as Ge if identified as TypeIII during fragmentation
		IAtom fragnode=fraggraph.getAtom(j);
		if(fragnode.getAtomicNumber()==32)
		{
			return 1;			
		}

		return 0;
	}
	
	private static int checkLinkers(IAtomContainer fragment, IAtom linker1,
			IAtom linker2) {
		// checks branching cycles to see if the linkers are both still there when trim linear portion
		IAtomContainer newfrag = fragment;
		stripTerminalAtoms(newfrag);
		int hit=0;
		for(IAtom atom : newfrag.atoms())
		{
			if(atom==linker1)
				hit++;
			
			else if(atom==linker2)
				hit++;
		}
		
		if(hit==2)
			return 1;
		
		else
			return 0;
	}

	private static int findFragments(int numfrags, IAtomContainer mol,
			int originalnumatoms, IAtomContainer fraggraph,
			IMoleculeSet cycles, IMoleculeSet fragments, 
			int[][] atomfragmatrix, int[][] terminalatoms, int[][] cyclearray) {
	 
		int [][] connectivity;
		connectivity = new int [originalnumatoms][2];
		IAtomContainer cycle = DefaultChemObjectBuilder.getInstance().newMolecule();	
		int numringsfound=0;
    
		findConnectivities(connectivity,mol);
 			
		while(mol.getAtomCount()>0)
		{	
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

	           		stripMolecule		
				(mol,subgraph,seedatom,iatom2,atomfragmatrix,numfrags,connectivity,terminalatoms);

	           		// now add the subgraph to our list of fragments
	           		fragments.addMolecule(subgraph);
           		
	           		// now add the atom to a new fragment graph - we'll add the edges later
	           		IAtom fragnode = DefaultChemObjectBuilder.getInstance().newAtom();
	           		fragnode.setID(Integer.toString(numfrags));
	           		// set atom type to C for default - could also set for different classes
	           		// for example linear, cycle, branching
	           		fragnode.setAtomicNumber(6); // default is 6 - carbon

		               	// check if subgraph contains a cycle
	           		if(checkCycle(subgraph)>0)
	           		{
	           			if(terminalatoms[numfrags][0]==1)
	           			{
	           				fragnode.setAtomicNumber(32); //Ge = branching ring found      
	           			}
           				
	           			else
	           			{
		           			// cyclic - so remove terminals - or rather number
           					terminalatoms[numfrags][0]=0;
           					fragnode.setAtomicNumber(14); //Si = ring found
           				}
          	 			cycle = addCycle(subgraph);
           				cycles.addAtomContainer(cycle);

           				addCycleArray(cycle,numringsfound,cyclearray);
           				numringsfound++;
           			}

           			fraggraph.addAtom(fragnode);
           			numfrags++;
           	    
       			}
       			else{
       				mol.removeAtom(seedatom);
       			}
		}

		return numringsfound;
	}

	private static void addCycleArray(IAtomContainer cycle, int numringsfound,
			int[][] cyclearray) {
		// builds the array of cycles to make lookup a little faster

		for(IAtom atom : cycle.atoms())
		{
			cyclearray[numringsfound][Integer.parseInt(atom.getID())]=1;
		}
			
	}

	private static IAtomContainer addCycle(IAtomContainer subgraph) {
		// strip away any terminals and then copy what's left into the list of found cycles
		IAtomContainer newcycle = DefaultChemObjectBuilder.getInstance().newMolecule();
		newcycle.add(subgraph);
		
		stripTerminalAtoms(newcycle);
		
		return newcycle;
	}
	
	private static int checkCycle(IAtomContainer subgraph) {
		// simple cauchy to see if subgraph contains a cycle
		int subcauchy = 0;
		subcauchy = subgraph.getBondCount() - subgraph.getAtomCount() + 1; // only ever 1 component
		return subcauchy;
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
           			// now find the neighbouring atom
           			List <IAtom> neighs = mol.getConnectedAtomsList(iatom1);
           			Iterator <IAtom> nextatom = neighs.iterator();
           			IAtom iatom2 = DefaultChemObjectBuilder.getInstance().newAtom();
           			while(nextatom.hasNext())
           			{
           				iatom2 = nextatom.next();
           			}
           			
           			// now remove atom1
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
           					// so reset in order that only have one terminal in these fragments
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

       		// now update connectivity now deleted a bond
       		if((lastatom1.getID()!=null)&&(lastatom2.getID()!=null))
      			updateConnectivity(connectivity,lastatom1,lastatom2,mol);	       
	}
	
	private static void updateConnectivity(int[][] connectivity,
			IAtom seedatom, IAtom iatom2, IAtomContainer mol) {
		// when bond removed the neighbours' connectivity index is updated
	
		int connectivity_index;
		int c1=0;
		int c2=0;
		int L1=0;
		int L2=0;
		
		L1=connectivity[Integer.parseInt(seedatom.getID())][1]-(connectivity[Integer.parseInt(seedatom.getID())][0]);

		connectivity_index=0;
		c1=connectivity[Integer.parseInt(seedatom.getID())][0];
		connectivity_index+=findConnectivityIndex(mol.getConnectedAtomsCount(seedatom));
		connectivity[Integer.parseInt(seedatom.getID())][0]=connectivity_index;	
		
		L2=connectivity[Integer.parseInt(iatom2.getID())][1]-(connectivity[Integer.parseInt(iatom2.getID())][0]);
		c2=connectivity[Integer.parseInt(iatom2.getID())][0];
		connectivity_index+=findConnectivityIndex(mol.getConnectedAtomsCount(iatom2));
		connectivity[Integer.parseInt(iatom2.getID())][0]=connectivity_index;
		
		connectivity[Integer.parseInt(seedatom.getID())][1]=connectivity[Integer.parseInt(seedatom.getID())][0]+L1;
		connectivity[Integer.parseInt(iatom2.getID())][1]=connectivity[Integer.parseInt(iatom2.getID())][0]+L2; 

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
	
	private static void findConnectivities(int[][] connectivity,
			IAtomContainer mol) {
		// optimised - finds connectivities once for all atoms then updates this list as atoms removed
		// connectivity[atom][0]=atom's connectivity index.
		// connectivity[atom][1]=atom's total connectivity index - summed over neighbours
		
		int connectivity_index;
		
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
		// adapted from zamora paper
		
		if(neighs==2)
			return 1;
		if(neighs==3)
			return 8;
		if(neighs>=4)
			return 64;
		
		return 0;
	}
	

	private static void stripTerminalAtoms (IAtomContainer mol)
	{
	    	// finds all one-connected atoms to use as seed
    		// removes them and checks if neighbours become terminal 
		// traces along branch until reach atom with >1 neighbour
		
		int numbonds = 0;
		IAtomContainer terminalatoms = 	DefaultChemObjectBuilder.getInstance().newAtomContainer();	
		
		// find terminal atoms
	       	for(IAtom iatom1 : mol.atoms())
	        {             	
	           	numbonds = mol.getConnectedAtomsCount(iatom1);
   			if(numbonds==1)
   				terminalatoms.addAtom(iatom1);
	        }

	       	// use each terminal as a seed to chew away linear chain
	       	for(IAtom deleteatom : terminalatoms.atoms())
	        {             	
   			numbonds=mol.getConnectedAtomsCount(deleteatom);
	      		IAtom iatom1 = deleteatom;
      		
   		   	while(numbonds==1)
	           	{
	           		// now find the neighbouring atom
	           		List <IAtom> neighs = mol.getConnectedAtomsList(iatom1);
	           		Iterator <IAtom> nextatom = neighs.iterator();
	           		IAtom iatom2 = DefaultChemObjectBuilder.getInstance().newAtom();
	           		while(nextatom.hasNext())
	           		{
	           			iatom2 = nextatom.next();
	           		}

	           		mol.removeBond(iatom1,iatom2);   
	           		mol.removeAtom(iatom1);                     // don't do this?			

	           		iatom1=iatom2;
	       			numbonds=mol.getConnectedAtomsCount(iatom1);                   			
	           	
	           	} 	        
	        }
	}	
}

