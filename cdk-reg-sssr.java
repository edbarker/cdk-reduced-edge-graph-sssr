package SmilesRead;

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

/*
 * Algorithm to generate a reduced edge graph and solve sssr
 */

public class SMILESRead {
	 	
	public static void main(String[] args) throws Exception 
    	{
     	int numcomponents = 0;
     	int solved =0;
     	int mismatch=0;
     	int DIAGNOSTICS=0;
     	
        String filename = "/home/ed/chem/eclipse/stuff/test.smi";
        FileInputStream ins = new FileInputStream(filename);
        BufferedWriter out = new BufferedWriter(new FileWriter("/home/ed/chem/eclipse/stuff/output.ring"));
        try {
            IteratingSMILESReader reader = new IteratingSMILESReader(ins);
            int num = 0;

            while (reader.hasNext())
            {
                IAtomContainer originalmol = DefaultChemObjectBuilder.getInstance().newMolecule();
                IAtomContainer copymol = DefaultChemObjectBuilder.getInstance().newMolecule();
                
            	originalmol = (IMolecule)reader.next();
                // set id for each atom
            	// is this necessary? - try to see if can avoid this
                int id = 0;
               
                IAtom atom = null;
                // i think this is necessary - otherwise when atoms deleted hard to keep track of them
                for(int k=0; k<originalmol.getAtomCount(); k++)
                { 
       	        	atom = originalmol.getAtom(k);
                    atom.setID(Integer.toString(id));
                    id++;
                }
            	
                // Note error - for some reason cdk not recognising certain element types
                // not sure how to switch this off as element type is not important in finding rings
                // to solve this I replaced certain very heavy elements (Th and Sm eg) by hand in the 
                // NCI dataset - otherwise it crashes
                // but this is a bodge and not a fix
                
                // output mol
            	out.write("MOL ");
            	out.write(Integer.toString(num));
            	out.write("\n");             	
            	
            	copymol.add(originalmol);
            	
        		IMoleculeSet cycles = DefaultChemObjectBuilder.getInstance().newMoleculeSet();		
            	int originalnumatoms=copymol.getAtomCount();
            	if(DIAGNOSTICS==1)
            		showMolecule(originalmol);
            	int cauchy =0; // for the ring nullity calculation
            	int numringsfound=0;

            	// need to partition into components first
            	// I think this is necessary - certainly to calculate cauchy need to know num components 
            	if(ConnectivityChecker.isConnected(copymol))
            	{
            		numcomponents = 1;
                	cauchy = copymol.getBondCount() - copymol.getAtomCount() + numcomponents;        		
                	if(cauchy>0)
                	{
                    	// do initial strip terminal branches now - saves time later
                    	stripTerminalAtoms(copymol);
                    	
                		// first fragment the molecule
                		numringsfound+=reduceMolecule(copymol,cauchy,originalnumatoms,cycles);
                		
                	}
            	}
            	else 
            	{
            		IMoleculeSet components = ConnectivityChecker.partitionIntoMolecules(copymol);
            		numcomponents = components.getMoleculeCount();
                	cauchy = copymol.getBondCount() - copymol.getAtomCount() + numcomponents;
            
        			//IAtomContainer componentit = components.molecules();
                	if(cauchy>0)
                	{      	
                		for(int f=0; f<numcomponents; f++)
                		{
                			IAtomContainer component = components.getMolecule(f);
                			int compcauchy = component.getBondCount() - component.getAtomCount() + 1;
                			if(compcauchy>0)
                			{
                            	// do initial strip terminal branches now - saves time later
                            	stripTerminalAtoms(component);
                            	
                				numringsfound+=reduceMolecule(component,compcauchy,originalnumatoms,cycles);
                			}
                		}
                	}
            	}
            	if(DIAGNOSTICS==1)
            	{
            		showCycles(cycles);
            	}
            	if(cauchy==cycles.getAtomContainerCount())
            	{
               		// output num rings
               		out.write("NUMRINGS ");
               		out.write(Integer.toString(cauchy));
               		out.write("\n");              		
            		solved++;
            		for (int k = 0; k < cauchy; k++) {
                 		out.write("RING ");
                   		out.write(Integer.toString(k));
                   		out.write("\t"); 
                   		
                       	IAtomContainer cycle = cycles.getAtomContainer(k);
                       	                       	
                       	out.write(Integer.toString(cycle.getAtomCount()));
                   		out.write("   ");
                    		
                       	for(int m=0;m<cycle.getAtomCount();m++) {
                        		
                       		IAtom atom2 = cycle.getAtom(m);
                       		out.write(atom2.getID());  
                       		out.write(" ");                 		
                       	}
                       	out.write("\n");   
            	
            		}
            	}
            	// then just use sssr - 
            	// in case of using just phase I - happens in 7% of cpds in nci db
            	else 
            	{
                    // now find rings using sssr
                    IAtomContainer ringatoms = null;
                    IRingSet sssrRings = new SSSRFinder(originalmol).findSSSR();

                    // output num rings
              		out.write("NUMRINGS ");
               		out.write(Integer.toString(sssrRings.getAtomContainerCount()));
               		out.write("\n");  
                		
                    List<IAtomContainer> ring = RingSetManipulator.getAllAtomContainers(sssrRings);
                    for (int j = 0; j < sssrRings.getAtomContainerCount(); j++) {
                 		out.write("RING ");
                   		out.write(Integer.toString(j));
                   		out.write("\t"); 
                   		
                       	ringatoms = ring.get(j);
                     		
                       	out.write(Integer.toString(ringatoms.getAtomCount()));
                   		out.write("   ");
                    		
                       	for(int m=0;m<ringatoms.getAtomCount();m++) {
                        		
                       		atom = ringatoms.getAtom(m);
                       		out.write(atom.getID());  
                       		out.write(" ");                 		
                       	}
                        	
                       	out.write("\n");   
         	           	    	
                    }
            	} 
                num++;
           		      	
            	if(DIAGNOSTICS==1)
            	{
            		System.err.println();                   		
            		//System.err.println("Found rings = " + Integer.toString(numringsfound));
            		//System.err.println("Cauchy = " + Integer.toString(cauchy));    		
            		System.err.println();
            	}

                out.write("|\n");            
                System.err.println(Integer.toString(num));
                //if(num%10000==0)
                	//System.err.println(Integer.toString(num));
            }

            ins.close();
            out.close();
            System.err.print("\nNumber of molecules read in = ");
            System.err.print(Integer.toString(num));
            System.err.println();
            System.err.println("Number of compounds solved = " + solved);
            System.err.println();
        } catch (Exception e) {
            e.printStackTrace();     
            //System.out.println(e.toString());
        }
    }

	private static int reduceMolecule(
			IAtomContainer mol, 
			int cauchy, 
			int originalnumatoms, 
			IMoleculeSet cycles) {
		// Main function to reduce the molecule to the reduced edge graph
		// finds the least connected atom and removes a bond from this atom
		// then removes all acyclic branches and reconnects what's removed into a fragment
		// then this fragment is represented as a node in a reduced edge graph
		// then solves simple cycles - Type II and Type III
		// if cauchy not found - then solve all linear fragments
		// if still not found cauchy - find cycles in reduced edge graph
		// solve these cycles
		// if still not found cycles - go to sssrfinder
		int numringsfound=0;
		int numfrags = 0;
	 	int DIAGNOSTICS=0;
	 	
		IAtomContainer cycle = DefaultChemObjectBuilder.getInstance().newMolecule();
		
		IMoleculeSet fragments = DefaultChemObjectBuilder.getInstance().newMoleculeSet();
		IAtomContainer fraggraph = DefaultChemObjectBuilder.getInstance().newMolecule();
		IAtomContainer origmol = DefaultChemObjectBuilder.getInstance().newMolecule();		
	
		origmol.add(mol);
		
		int[][] atomfragmatrix; 
    	atomfragmatrix = new int[originalnumatoms][originalnumatoms*2]; 
    	//overkill - but ok - 
    	// actually no - it is possible that the number frags can be more than the num original atoms
    	// although these molecules really are too complex to resolve - but it does happen
    	// these molecules are densely connected - strange Boron containing compounds in NCI
    	// checked - java will initialise arrays to zero
    	
		while(mol.getAtomCount()>0)
		{
			// step one is to identify least connected atom
			// Q - if remove this - how much time is saved and how many fewer are solved?
			// but this is done each time - tested this and it is worth the time taken to calculate
			// however I believe that it could be possible to speed this up simply by
			// creating an array or linked list of the connectivity indices
			// which is updated whenever an atom is deleted
			// and least connected could be easily looked up
			
			//IAtom seedatom = findLeastConnectedAtom(mol);
			IAtom seedatom = mol.getFirstAtom();
						
			IMolecule subgraph = DefaultChemObjectBuilder.getInstance().newMolecule();
			
			// now delete a bond - arbitrarily (maybe least connected neigh is better?)
			List <IAtom> neighs = mol.getConnectedAtomsList(seedatom);
       		Iterator <IAtom> nextatom = neighs.iterator();
       		if(nextatom.hasNext())
       		{
       			IAtom iatom2 = nextatom.next();
       			
       			// need to do this first - just in case only one bond deleted
         		subgraph.addAtom(seedatom);
         		subgraph.addAtom(iatom2);
         		subgraph.addBond(mol.getBond(seedatom,iatom2));
         		
         		if(DIAGNOSTICS==1)
         		{
         			System.err.println("\nAtom : " + Integer.parseInt(seedatom.getID()) + " & " + Integer.parseInt(iatom2.getID()) + " frag : " + Integer.toString(numfrags));
         		}
         		// set up array for later
         		atomfragmatrix[Integer.parseInt(seedatom.getID())][numfrags]=1;
         		atomfragmatrix[Integer.parseInt(iatom2.getID())][numfrags]=1;
         		
       			mol.removeBond(seedatom,iatom2);   
           		// now strip molecule of termini
           		// and build a new molecule from what is removed - remembering to add deleted bond
           		// and add a node to reduced edge molecule representing this subgraph
           		stripMolecule(mol,subgraph,seedatom,iatom2,atomfragmatrix,numfrags);
          
            	if(DIAGNOSTICS==1)
            	{
            		System.err.println();
            		System.err.println("Fragment " + numfrags);
            		showMolecule(subgraph);
            		System.err.println();
            	}
           		// now add the subgraph to our list of fragments
           		fragments.addMolecule(subgraph);
           		
           		// now add the atom to a new fragment graph - we'll add the edges later
           		IAtom fragnode = DefaultChemObjectBuilder.getInstance().newAtom();
           		fragnode.setID(Integer.toString(numfrags));
           		// set atom type to C for default - could also set for different classes
           		// for example linear, cycle, branching
           		// but really need for unfound cycle
           		fragnode.setAtomicNumber(6); // default is 6 - carbon
           		// now add the atom to a new fragment graph - we'll add the edges later
         		
               	// check if subgraph contains a cycle
           		if(checkCycle(subgraph)>0)
           		{
           			fragnode.setAtomicNumber(14); //Si = ring found
           			cycle = addCycle(subgraph);
           			cycles.addAtomContainer(cycle);
           			numringsfound++;
           			if(numringsfound==cauchy)
           			{   						
           				if(DIAGNOSTICS==1)
							showCycles(cycles);
           				return numringsfound;            				
           			}

           			// problem if components as not know what found - fix this
           			// *** ERROR - can't end now as the ring may be split later
           			// - not true as if found cauchy there can be no further splitting! eejut
           		}

           		fraggraph.addAtom(fragnode);
           		numfrags++;
           	    
       		}
       		else{
       			mol.removeAtom(seedatom);
       		}
		}
		
    	if(DIAGNOSTICS==1)
    	{
    		System.err.println();
    		System.err.println("Generated all fragments*********************************");
    		System.err.println("Generated Fragments: " + fragments.getMoleculeCount());		
    	}
    	
    	// now build fragment graph
       	// first go through each fragment and mark atoms as being members, 
       	// also count numfrags atom is in

    	if(DIAGNOSTICS==1)
    	{
    		for(int i=0;i<origmol.getAtomCount();i++)
    		{
    			System.err.print("\n Atom " + Integer.toString(i) + " : "); 
    			for(int j=0;j<numfrags;j++)
    			{
    				System.err.print(" " + Integer.toString(atomfragmatrix[i][j]));
    			}
    		}
    	}

    	// optimise - instead of this maybe better to just look for linear fragments which 2-con 
    	// to single frag
    	// then if not all found then can build fragment graph afterwards 
    	// - this avoids wasting time building fraggraph when don't actually use it
    	    	
    	// find easy cycles - optimised for efficiency
    	// find linear fragment
    	// see if frag has both termini in the same fragment - if so resolve
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
       				incommon=countAtomsinCommon(atomfragmatrix,originalnumatoms,y1,y2);
       				if(incommon==2)
       				{
   						numringsfound+=resolveEasyCycle(fragments,y1,y2,cycles);
   						// there are times when can find a ring needing splitting that
   						// contains two previously found rings
   						// so don't return when find all cycles until tried all possibilities
   						
   						/*if(numringsfound==cauchy)
   	           			{
   							if(DIAGNOSTICS==1)
   								showCycles(cycles);
   	           				return numringsfound;            				
  	           			}*/  	
      				}
       			}
       		}
       	}   
    	
    	// possible optimisation -
    	// repeat stage one by recreating original molecule
    	// fragment as before - but somehow find a way to direct the fragmentation
    	// so the last fragment found is a ring we haven't yet seen
    	// tried this in my C version and it works - but not on all
    	// again may improve speed by doing this.
    	
    	if(numringsfound==cauchy)
    		return numringsfound;

    	int[][] matrix; 
    	matrix = new int[numfrags][numfrags];
       	
    	int fragcount=0;
    	for(int k=0; k< originalnumatoms; k++)
       	{    
       		//System.err.println("\nFragment " + frag1);
       		//IAtomContainer fragment = fragments.getMolecule(k);
       		//showMolecule(fragment);
       		fragcount=0;
       		for(int z=0; z<numfrags; z++)
       		{
       			// a fragment is linked by one-connected atoms so find these first
       			// this is wrong as sometimes two cycles fused - ring fusion - spiro
       			// so need to go through all atoms! Need to optimise this
       			//optimised version using atom array
       			if(atomfragmatrix[k][z]==1)
       			{
       				fragcount++;
       			}
       		}
       		if(fragcount>1)
       		{
       			// then we go through the array again and make link to neigbouring fragments
       			for(int m=0; m<numfrags; m++)
       	       	{          		
       				if(atomfragmatrix[k][m]==1)
       				{
       					for(int n=m+1; n<numfrags; n++)
       					{ 
       						if(atomfragmatrix[k][n]==1)
       						{
   	        					matrix[n][m]=1;
   	        					matrix[m][n]=1;       						
       						}
       					}

       				}       				
       			}
       		}       			
       	}
 
       	// now connect fragments together using matrix

       	for(int i=0;i<numfrags;i++)
       	{
       		IAtom fragtoadd1=fraggraph.getAtom(i);
       		//fragtoadd1.setID(Integer.toString(i)); // why did i do this?
       		
       		int numshared=0;
       		for(int j=i;j<numfrags;j++)
       		{
           		IAtom fragtoadd2=fraggraph.getAtom(j);
       			if(matrix[i][j]>0)
       			{
       				IBond fragbond = DefaultChemObjectBuilder.getInstance().newBond();

   					fragbond.setAtom(fragtoadd1,0);
   					fragbond.setAtom(fragtoadd2,1);
   						
   					numshared=countAtomsinCommon(fragments,i,j);
   					
   					if(numshared==1)
   						fragbond.setOrder(IBond.Order.SINGLE);
   					if(numshared>1)
   					{
   						fragbond.setOrder(IBond.Order.DOUBLE);   
   						matrix[i][j]=2;
   						matrix[j][i]=2;   				
   						
   						// now resolve this fragment - no done above now
   				      	if(DIAGNOSTICS==1)
   				    	{
   				    		int num=0;
   				    		for(IAtomContainer fragx : fragments.molecules())
   				    		{
   				    			System.err.println("\nCheck Fragment " + num);
   				    			showMolecule(fragx);
   				    			num++;
   				    		}
   				    	}/*
   						numringsfound+=resolveEasyCycle(fragments,i,j,cycles);
   				      	if(DIAGNOSTICS==1)
   				    	{
   				    		int num=0;
   				    		for(IAtomContainer fragx : fragments.molecules())
   				    		{
   				    			System.err.println("\nCheck Fragment 2" + num);
   				    			showMolecule(fragx);
   				    			num++;
   				    		}
   				    	}
   						if(numringsfound==cauchy)
   	           			{
   	           				return numringsfound;            				
   	           			}  */						
   					}
					       			
   					fraggraph.addBond(fragbond);	
       			}
       		}
       	}
       	
    	// to begin - simplify fraggraph
       	stripTerminalAtoms(fraggraph);
       	if(DIAGNOSTICS==1)
       		showMolecule(fraggraph);

		if(DIAGNOSTICS==1)
		{
			int fragcauchy = fraggraph.getBondCount()-fraggraph.getAtomCount()+1;
			System.err.println("\nFragcycles: cauchy = " + fragcauchy);
			
		}	
       	//now find cycles in this fraggraph
		
		// possible optimisation -
		// if only one cycle (which is majority of cases)
		// don't bother doing sssr finding as what's left is the only ring
		// this may possibly save time - not sure
		// Also possibly could do triangle finding first as this may be faster than sssrfinder
		// not sure - again majority of cycles are 3-membered rings so this would speed it up
		
        //IAtomContainer ringatoms = null;
        IRingSet sssrRings = new SSSRFinder(fraggraph).findSSSR();
		//AllRingsFinder allRingsFinder = new AllRingsFinder();
		//IRingSet sssrRings = allRingsFinder.findAllRings(atomContainer);
		List<IAtomContainer> fragcycle = RingSetManipulator.getAllAtomContainers(sssrRings);
        //now for each node in the cycle need to resolve ring
        numringsfound=resolveComplexCycles(fragcycle,fragments,cycles,cauchy,numringsfound);
		if(numringsfound==cauchy)
		{
			return numringsfound;            				
		}
				
        // if it still has not found all the cycles then we can try the next idea
        // - find the largest cycle - which is assumed to be one which needs splitting
        // - add on each fragment in turn and see if the cycles created are new
        // if not then remove the fragment and add the next fragment
        // hopefully the cycles will all be found - desperate I know!
        // not sure we actually need this now - was an error with split cycles
        //numringsfound=resolveUnfoundCycles(fragments,cycles,cauchy,numringsfound);
		
		else
		{
			// throw in the towel
			// if it gets here then weird stuff is happening 
			// sssr seems not to find cycles in complex graphs
			// so just do sssr on original molecule
			
			// reason not found at this stage I believe is due to problem where sssrfinder
			// does not see a ring that this algorithm does. Way around this may be to do 
			// allringsfinder - but this is going to be too slow so better to just send 
			// to sssrfinder
			
			return numringsfound;			
		}

	}

	private static int countAtomsinCommon(IMoleculeSet fragments, int i, int j) {
		// simply counts number of atoms shared in two fragments
		
		IAtomContainer frag1 = fragments.getMolecule(i);
		IAtomContainer frag2 = fragments.getMolecule(j);		
		
		// perhaps could optimise this with some sort of lookup table
		
		int incommon=0;
		for(IAtom frag1atom : frag1.atoms())
		{
			for(IAtom frag2atom : frag2.atoms())
			{			
				if(frag1atom == frag2atom)
				{
					incommon++;
				}
			}
		}
	
		return incommon; //must be at least 1 or 2 - doubt it can be more and never 0 as not come here otherwise
	}
	
	private static int resolveComplexCycles(
			List<IAtomContainer> fragcycles,
			IMoleculeSet fragments, 
			IMoleculeSet cycles, 
			int cauchy,
			int ringsfound) {
		// so now we need to find the fragment, find the linkers, build the cycle
		int DIAGNOSTICS=0;
		//int ringsfound=0;
		int count=0;
		Iterator <IAtomContainer> fragcycleit = fragcycles.iterator();
		while(fragcycleit.hasNext())
		{
			IAtomContainer fragcycle=fragcycleit.next();
			if(DIAGNOSTICS==1)
			{
				System.err.println("\nFragcycle: " + count);
				showMolecule(fragcycle);
			}
			count++;
			// now find the fragments
			int [] fragnode;
			fragnode = new int [fragcycle.getAtomCount()];
			int numfrags=0;
			for(IAtom node : fragcycle.atoms())
			{
				fragnode[numfrags]=Integer.parseInt(node.getID());
				numfrags++;
			}
			
			// now just add fragments together and trim it!!!
			// not that easy - what if one of frags is a cycle
			IAtomContainer newcycle = DefaultChemObjectBuilder.getInstance().newMolecule();
			for(int i=0;i<numfrags;i++)
			{
				// hmmmm complicated
				// if the fragment is linear then just add it
				// if cyclic then add shortest path through linkers
				// but then you need to find linkers
				
				// beware - k1,3 interchange! - deal with this before do this or will not work
				// - no just add as if linear - strip will delete all atoms so check before adding cycle

				// hmmmm - there must be a linear by def (can't be all 2connected or cycles
				// so could use this as a seed - first one come to
				// then find linker to next node and that node's neighbour
				
				// cheat method - use sssr to find cycles and see which one is unique - add this one
				// until find way to do it myself - which will be possible if can sort out structures in java
				// in my C version I had a quick way to do all this but for some reason
				// can't find out how to do data structures in java that could do this
				// so have to use sssr finder instead
				
				// add all frags together
				if(DIAGNOSTICS==1)
				{
					System.err.println("\nAdding fragment " + fragnode[i]);
					showMolecule(fragments.getAtomContainer(fragnode[i]));
				}
				newcycle.add(fragments.getAtomContainer(fragnode[i]));				
			}
			if(DIAGNOSTICS==1)
			{
				System.err.print("\nCycle made up of fragments: ");
				for(int p=0;p<numfrags;p++)
					System.err.print(" " + Integer.toString(fragnode[p]));
				System.err.println();
				System.err.println();					
				showMolecule(newcycle);
				System.err.println();
				System.err.println();				
			}
			// now strip terminal branches
			stripTerminalAtoms(newcycle);
			
			//check if actually anything left - k1,3 should be linear so all deleted
			// could do cycle check here but not sure if empty molecule would come back as ring=1
			if(newcycle.getAtomCount()>0)
			{		
				int numrings = checkCycle(newcycle);
				if(numrings==1)
				{		
					// then easy - just check it's unique and add
					if(checkUniqueCycle(newcycle,cycles)==1)
					{
						// check splitting too
			        	checkSplittingComplex(newcycle,cycles);
			        	ringsfound++;
			        	if(DIAGNOSTICS==1)
			        		System.err.println("FOUND A RING " + Integer.toString(ringsfound) + " of " +Integer.toString(cauchy));
			        	cycles.addAtomContainer(newcycle);		
   	           			if(ringsfound==cauchy)
   	           				return ringsfound; 
					}			
				}
				else if (numrings>1)
				{
			       	//now find cycles in the molecule
			        IRingSet fragrings = new SSSRFinder(newcycle).findSSSR();
			        List<IAtomContainer> fragringset = RingSetManipulator.getAllAtomContainers(fragrings);
					Iterator <IAtomContainer> fragringit = fragringset.iterator();
					
					while(fragringit.hasNext())
					{
						IAtomContainer fragcyc = fragringit.next();
						IAtomContainer copycyc = DefaultChemObjectBuilder.getInstance().newMolecule();
						copycyc.add(fragcyc);
				        // now check if ring is unique
				        if(checkUniqueCycle(copycyc,cycles)==1)
				        {
				        	
				        	if(DIAGNOSTICS==1)
				        		showMolecule(copycyc);
				        	
				        	checkSplittingComplex(copycyc,cycles);
				        	
				        	if(DIAGNOSTICS==1)
				        		showMolecule(copycyc);
				    
				        	cycles.addAtomContainer(copycyc);
				        	ringsfound++;
				        	if(DIAGNOSTICS==1)
				         		System.err.println("FOUND A RING " + Integer.toString(ringsfound) + " of " +Integer.toString(cauchy));
							   
				        	if(ringsfound==cauchy)
	   	           				return ringsfound; 
				        }							
					}		
				}
			}								
		}	
		
		return ringsfound;
	}

	private static void checkSplittingComplex(IAtomContainer newcycle,
			IMoleculeSet cycles) {
		// does check splitting - but this is where only have the cycle info
		int deleted=0;
		int DIAGNOSTICS=0;
		if(DIAGNOSTICS==1)
		{
			System.err.println("\nSplitting ring check");
		}

		for(IAtomContainer cycle : cycles.molecules())
		{
			// first check if the cycle is larger than newcycle - otherwise can't be split
			// (am I sure about this?)
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

				// now check if the number of bonds is greater than number exclusive bonds
				if(DIAGNOSTICS==1)
				{
					System.err.println("\nCommon bonds = " + commonbonds.getBondCount());
					System.err.println("\nExclusive bonds = " + exclusivebonds.getBondCount());
			
				}
				if(commonbonds.getBondCount()>exclusivebonds.getBondCount())
				{
					// then update the old cycle by removing commonbonds
					// does this work? - ie does altering cycle alter cycles?
					oldcycle=cycle;
					for(IBond deletebond : commonbonds.bonds())
					{
						oldcycle.removeAtom(deletebond.getAtom(0));
						oldcycle.removeAtom(deletebond.getAtom(1));						
						oldcycle.removeBond(deletebond);						
					}

					//showMolecule(oldcycle);
					for(IBond addbond : exclusivebonds.bonds())
					{
						oldcycle.addAtom(addbond.getAtom(0));
						oldcycle.addAtom(addbond.getAtom(1));						
						oldcycle.addBond(addbond);						
					}

					//showMolecule(oldcycle);					
					//cycles.removeAtomContainer(cycle);
					//cycles.addAtomContainer(oldcycle);
					cycle=oldcycle;
					if(checkUniqueCycle(cycle,cycles)==0)
					{
						cycles.removeAtomContainer(cycle);
						deleted++;
					}
				}		
			}
		}		
	
	}

	private static int checkUniqueCycle(IAtomContainer newcycle,
			IMoleculeSet cycles) {
		// simply check that this is not a cycle already found
		for(IAtomContainer cycle : cycles.molecules())
		{
			int incommon=0;
			for(IAtom newatom : newcycle.atoms())
			{
				for(IAtom cycleatom : cycle.atoms())
				{
					if(newatom==cycleatom)
						incommon++;
				}
			}
			if(incommon==cycle.getAtomCount())
				return 0;
			
			
			//if(newcycle == cycle) // seems not to work
				//return 0;
		}
		
		return 1;
	}

	private static int resolveEasyCycle(IMoleculeSet fragments, int i, int j,
			IMoleculeSet cycles) {
		// Finds simple 2-connected fragment rings 
		// - and also checks for splitting of rings already found 
					
		int DIAGNOSTICS=0;
		if(DIAGNOSTICS==1)
		{
			System.err.println();
			System.err.println("Solving easy cycle - fragments : " + Integer.toString(i) + " and " + Integer.toString(j));
			System.err.println();
		}
		IAtomContainer origA = fragments.getMolecule(i);
		IAtomContainer origB = fragments.getMolecule(j);       						
		IAtomContainer cycle = DefaultChemObjectBuilder.getInstance().newMolecule();
		IAtomContainer fragA = DefaultChemObjectBuilder.getInstance().newMolecule();
		IAtomContainer fragB = DefaultChemObjectBuilder.getInstance().newMolecule();
		
		fragA.add(origA);
		fragB.add(origB);
		
		int foundring=0;
		
		if((checkCycle(fragA)==0)&&(checkCycle(fragB)==0))
       	{
       		// if both linear then this is easy to resolve

			cycle.add(fragA);
       		cycle.add(fragB); //does this work? does it fuse shared atoms?
       		stripTerminalAtoms(cycle);
       		// check for splitting here - no need as can't split a linear chain
       		// you're wrong - it can and does
       		checkSplittingComplex(cycle,cycles);
 
       		if(checkUniqueCycle(cycle,cycles)==1)
       			cycles.addAtomContainer(cycle);	
       		foundring=1;
			if(DIAGNOSTICS==1)
			{
				System.err.println();
				System.err.println("Solving easy cycle - linear fragments: ");
				showMolecule(fragA);					
				System.err.println();
				showMolecule(fragB);					
				System.err.println();
				showMolecule(cycle);
			}
       		
       		return foundring;
       	}
       	else // one must be a cycle - so need shortest path
       	{
       		// find linear fragment
       		// find linkers
       		// find sp in other fragment - add this path to cycle
       		// add linear frag to cycle
			if(checkCycle(fragA)==0)
			{
				IAtom linker1 = DefaultChemObjectBuilder.getInstance().newAtom();
				IAtom linker2 = DefaultChemObjectBuilder.getInstance().newAtom();
      							
				// find linker atoms in fragment 
				for(IAtom linker : fragA.atoms())
				{
					if(fragA.getConnectedAtomsCount(linker)==1)
					{
						if(linker1.getID()!=null)
						{
							linker2 = linker;
						}
						else
						{
							linker1=linker;
						}
					}
				}

				if(DIAGNOSTICS==1)
				{
					System.err.println();
					System.err.println("Solving easy cycle - fragments: linkers: " + linker1.getID() + " & " + linker2.getID());
					showMolecule(fragA);					
					System.err.println();
					showMolecule(fragB);					
					System.err.println();
				}
				
				// now find shortest path in other fragment
				List <IAtom> pathatoms=PathTools.getShortestPath(fragB,linker1,linker2);
				IAtomContainer shortestpath = DefaultChemObjectBuilder.getInstance().newMolecule();
				if(DIAGNOSTICS==1)
				{
					System.err.println("\nShortest path checking frag: " + fragB.getID());
					showMolecule(shortestpath);
				}						
				// now build shortest path from fragB
				for(int h=0; h<pathatoms.size(); h++)
				{
					IAtom atom1=pathatoms.get(h);
					
					if(DIAGNOSTICS==1)
					{
						System.err.println("\nShortest path checking atom: " + atom1.getID());
					}			
					shortestpath.addAtom(atom1);
				}
				// now connect the atoms
				for(IBond bond1 : fragB.bonds())
				{
					if(shortestpath.contains(bond1.getAtom(0))&&(shortestpath.contains(bond1.getAtom(1))))
					{
						shortestpath.addBond(bond1);									
					}
				}
				shortestpath.add(shortestpath);
				if(DIAGNOSTICS==1)
				{
					System.err.println("\nShortest path after frag: " + fragB.getID());
					showMolecule(shortestpath);
				}		
				if(DIAGNOSTICS==1)
				{
					System.err.println();
					System.err.println("Solving easy cycle - shortestpath: ");
					showMolecule(shortestpath);					
					System.err.println();
				}
				cycle.add(shortestpath);
    			cycle.add(fragA);
    			//trim cycle
    			stripTerminalAtoms(cycle);
    			// check for splitting here
    			if(fragA.getAtomCount()<pathatoms.size())
    			{
    				if(checkBranching(fragB)==1)
    				{
    					// if a branching cycle then need to check that linkers are in the cycle
    					if(checkLinkers(fragB,linker1,linker2)==1)
        					splitCycle(cycles,fragB,fragA,shortestpath);	    						
    				}
    				else
    					splitCycle(cycles,fragB,fragA,shortestpath);	
    			}
    			checkSplittingComplex(cycle,cycles);
          		if(checkUniqueCycle(cycle,cycles)==1)
          			cycles.addAtomContainer(cycle);
        		foundring=1;
        		if(DIAGNOSTICS==1)
        			showMolecule(cycle);
       			return foundring;          								
       		}
       		else
       		{
    			IAtom linker1 = DefaultChemObjectBuilder.getInstance().newAtom();
       			IAtom linker2 = DefaultChemObjectBuilder.getInstance().newAtom();
       			
       			// find linker atoms in fragment 
       			for(IAtom linker : fragB.atoms())
       			{
       				if(fragB.getConnectedAtomsCount(linker)==1)
       				{
       					if(linker1.getID()!=null)
       					{
       						linker2 = linker;
       					}
       					else
       					{
       						linker1=linker;
       					}
       				}
       			}

				if(DIAGNOSTICS==1)
				{
					System.err.println();
					System.err.println("Solving easy cycle - fragments: linkers: " + linker1.getID() + " & " + linker2.getID());
					showMolecule(fragA);					
					System.err.println();
					showMolecule(fragB);					
					System.err.println();
				}
				
       			// now find shortest path in other fragment
       			List <IAtom> pathatoms=PathTools.getShortestPath(fragA,linker1,linker2);
				IAtomContainer shortestpath = DefaultChemObjectBuilder.getInstance().newMolecule();

				// now build shortest path from fragB
				if(DIAGNOSTICS==1)
				{
					System.err.println("\nShortest path checking frag: " + fragA.getID());
					showMolecule(shortestpath);
				}						
				for(int h=0; h<pathatoms.size(); h++)
				{
					IAtom atom1=pathatoms.get(h);
					
					if(DIAGNOSTICS==1)
					{
						System.err.println("\nShortest path checking atom: " + atom1.getID());
					}			
					shortestpath.addAtom(atom1);
				}
				// now connect the atoms
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
       			// check for splitting here
    			if(fragB.getAtomCount()<pathatoms.size())
    			{
    				if(checkBranching(fragA)==1)
    				{
    					// if a branching cycle then need to check that linkers are in the cycle
    					if(checkLinkers(fragA,linker1,linker2)==1)
        					splitCycle(cycles,fragA,fragB,shortestpath);	    						
    				}
    				else
    					splitCycle(cycles,fragA,fragB,shortestpath);	
    			}
    			checkSplittingComplex(cycle,cycles);
          		if(checkUniqueCycle(cycle,cycles)==1)
          			cycles.addAtomContainer(cycle);
                foundring=1;
                if(DIAGNOSTICS==1)
                	showCycles(cycles);
         	    return foundring;
			}			
		}
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

	private static int checkBranching(IAtomContainer fragment) {
		// simple check to see if has one connected atom
		// needed to avoid problem with splitting when two linear portions connect
		
		for(IAtom atom : fragment.atoms())
		{
			if(fragment.getConnectedAtomsCount(atom)==1)
				return 1;
		}
		return 0;
	}

	private static void splitCycle(
			IMoleculeSet cycles,
			IAtomContainer cyclefragment, 
			IAtomContainer smallpath, 
			IAtomContainer largepath) {
		// split the cycle previously found by removing largepath and replace with smallpath
		IAtomContainer newcycle = DefaultChemObjectBuilder.getInstance().newMolecule();
		stripTerminalAtoms(newcycle); //justincase branching type
		
		// first find cycle to replace
		// may not be one in case where branching cycle fragment has cycle outside its cycle
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

	private static int countAtomsinCommon(int[][] atomfragmatrix, int originalnumatoms, int i, int j) {
		// simply counts number of atoms shared in two fragments
		
		//optimised version - use fragatommatrix
		int incommon=0;
		
		for(int z=0;z<originalnumatoms;z++)
		{
			if((atomfragmatrix[z][i]==1)&&(atomfragmatrix[z][j]==1))
					incommon++;
		}
		
		/*
		IAtomContainer frag1 = fragments.getMolecule(i);
		IAtomContainer frag2 = fragments.getMolecule(j);		
		
		int incommon=0;
		for(IAtom frag1atom : frag1.atoms())
		{
			for(IAtom frag2atom : frag2.atoms())
			{			
				if(frag1atom == frag2atom)
				{
					incommon++;
				}
			}
		}
		 */
		return incommon; //must be at least 1 or 2 - doubt it can be more and never 0 as not come here otherwise
	}


	private static int checkCycle(IAtomContainer subgraph) {
		// simple cauchy to see if subgraph contains a cycle
		int subcauchy = 0;
		subcauchy = subgraph.getBondCount() - subgraph.getAtomCount() + 1; // only ever 1 component
		return subcauchy;

	}

	private static IAtomContainer addCycle(IAtomContainer subgraph) {
		// strip away any terminals and then copy what's left into the list of found cycles
		IAtomContainer newcycle = DefaultChemObjectBuilder.getInstance().newMolecule();
		newcycle.add(subgraph);
		
		stripTerminalAtoms(newcycle);
		
		return newcycle;
	}

	private static void stripMolecule(
			IAtomContainer mol, 
			IAtomContainer subgraph, 
			IAtom a, IAtom b, 
			int[][] atomfragmatrix, int numfrags) {
		// same as strip terminals - except keep record of what is deleted
	
	   	// find all one-connected atoms
    	// then remove its neighbour if it will become terminal - rpt till reach dead end
		// repeat for all one-connected atoms
		// problem - as deleting messes up position numbers - so delete bonds instead
		
		int numbonds = 0;
		IAtomContainer terminalatoms = 	DefaultChemObjectBuilder.getInstance().newAtomContainer();	
		
       	for(IAtom iatom1 : mol.atoms())
        {             	
           	//iatom1 = mol.getAtom(k);
           	numbonds = mol.getConnectedAtomsCount(iatom1);
   			if(numbonds==1)
   				terminalatoms.addAtom(iatom1);   				
        }

       	for(IAtom deleteatom : terminalatoms.atoms())
        {             	
      		//System.err.print("\nDeleting Atom from list - seeding ");
   			//System.err.print(deleteatom.getID());
  		
   			numbonds=mol.getConnectedAtomsCount(deleteatom);
      		//System.err.print("\nNumneighs....." + numbonds);
      		
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
       				atomfragmatrix[Integer.parseInt(iatom1.getID())][numfrags]=1;
       			}
       			if((iatom2!=a)&&(iatom2!=b))
       			{
       				subgraph.addAtom(iatom2);
       				atomfragmatrix[Integer.parseInt(iatom2.getID())][numfrags]=1;
       			}
       	
       			subgraph.addBond(mol.getBond(iatom1,iatom2));   			
       			
           		//System.err.print("\nDeleting Bond to Atom - ");
       			//System.err.print(iatom1.getID());

           		mol.removeBond(iatom1,iatom2);   
           		mol.removeAtom(iatom1);                     			

           		iatom1=iatom2;
       			numbonds=mol.getConnectedAtomsCount(iatom1);                   			
           	
           	} 	
        }
	}
		
	private static IAtom findLeastConnectedAtom(IAtomContainer mol) {
		// finds the atom least connected according to zamora
		// although this works i may have got this wrong 
		// i do Ki + Li but should be 64*Ki + Li - confused
		
		// step one is to identify least connected atom
		// Q - if remove this - how much time is saved and how many fewer are solved?
		// but this is done each time - tested this and it is worth the time taken to calculate
		// however I believe that it could be possible to speed this up simply by
		// creating an array or linked list of the connectivity indices
		// which is updated whenever an atom is deleted
		// and least connected could be easily looked up
		
		int connectivity;
		int leastconnectivity=9999999;
		IAtom leastatom=null;
		// first build up the list of connectivities
		for(IAtom atomit : mol.atoms())
		{
			connectivity=0;
			connectivity+=findConnectivityIndex(mol.getConnectedAtomsCount(atomit));
			List <IAtom> neighs = mol.getConnectedAtomsList(atomit);
       		Iterator <IAtom> nextatom = neighs.iterator();
       		while(nextatom.hasNext())
       		{
       			IAtom iatom2 = nextatom.next();
       			connectivity+=findConnectivityIndex(mol.getConnectedAtomsCount(iatom2));
       		}				
       		if(connectivity<leastconnectivity)
       		{
       			leastatom=atomit;
       			leastconnectivity=connectivity;
       		}
		}
		//System.err.print("\nLeast connected atom - " + leastatom.getID() + " " + Integer.toString(leastconnectivity));
		return leastatom;
	}

	private static int findConnectivityIndex(int neighs) {
		// from zamora paper
		
		if(neighs==2)
			return 1;
		if(neighs==3)
			return 8;
		if(neighs>=4)
			return 64;
		
		return 0;
	}

	private static void showMolecule(IAtomContainer mol) {
		// show all the atoms and neighbours
		//IBond bond;
		for(int x=0;x<mol.getAtomCount();x++)
		{
			IAtom iatom1 = mol.getAtom(x);
			System.err.print("\nAtom: ");
			System.err.print(iatom1.getID());
			System.err.print("\tNeighs: ");
			List <IAtom> neighs = mol.getConnectedAtomsList(iatom1);
       		Iterator <IAtom> nextatom = neighs.iterator();
       		while(nextatom.hasNext())
       		{
       			IAtom iatom2 = nextatom.next();
       			System.err.print("\t");
       			System.err.print(iatom2.getID());           			
       		}				

		}
	}
	
	private static void showCycles(IMoleculeSet cycles) {
		// show what cycles have been found
		int numcycles=0;
		for(IAtomContainer cycle : cycles.molecules())
		{
			System.err.print("\nCycle: ");
			System.err.print(Integer.toString(numcycles));

			for(IAtom atom : cycle.atoms())
			{
				System.err.print("\nAtom: ");
				System.err.print(atom.getID());
				System.err.print("\tNeighs: ");
				List <IAtom> neighs = cycle.getConnectedAtomsList(atom);
	       		Iterator <IAtom> nextatom = neighs.iterator();
	       		while(nextatom.hasNext())
	       		{
	       			IAtom iatom2 = nextatom.next();
	       			System.err.print("\t");
	       			System.err.print(iatom2.getID());           			
	       		}					
			}
			numcycles++;
		}
	}
	
	public static void stripTerminalAtoms (IAtomContainer mol)
	{
    	// find all one-connected atoms
    	// then remove its neighbour if it will become terminal - rpt till reach dead end
		// repeat for all one-connected atoms
		// problem - as deleting messes up position numbers - so delete bonds instead
		// I think I have optimised this by tracing back along branches
		// - quicker than just going through repeatedly deleting terminals 
		// as have to go back to start after each deletion if do that 
		
		int numbonds = 0;
		IAtomContainer terminalatoms = 	DefaultChemObjectBuilder.getInstance().newAtomContainer();	
		
       	for(IAtom iatom1 : mol.atoms())
        {             	
           	//iatom1 = mol.getAtom(k);
           	numbonds = mol.getConnectedAtomsCount(iatom1);
      		//System.err.print("\nChecking Atom - ");
   			//System.err.print(iatom1.getID() + " has neighs =" + numbonds);
   			if(numbonds==1)
   				terminalatoms.addAtom(iatom1);
        }

       	for(IAtom deleteatom : terminalatoms.atoms())
        {             	
      		//System.err.print("\nDeleting Atom from list - seeding ");
   			//System.err.print(deleteatom.getID());
  		
   			numbonds=mol.getConnectedAtomsCount(deleteatom);
      		//System.err.print("\nNumneighs....." + numbonds);
      		
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
           			
           		//System.err.print("\nDeleting Bond to Atom - ");
       			//System.err.print(iatom1.getID());

           		mol.removeBond(iatom1,iatom2);   
           		mol.removeAtom(iatom1);                     			

           		iatom1=iatom2;
       			numbonds=mol.getConnectedAtomsCount(iatom1);                   			
           	
           	} 	
        
        }
       	
	}
}

