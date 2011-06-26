package barkerSSSRFinder;

import java.io.IOException;
import java.util.List;

import org.openscience.cdk.*;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.ringsearch.AllRingsFinder;
import org.openscience.cdk.ringsearch.SSSRFinder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.RingSetManipulator;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IMoleculeSet;


/*
 * tests that the rings found by barkerSSSRFinder are correct
 * it does this by finding the correct sssr and seeing if the rings match
 * if they do not match this is not necessarily because it has found a wrong ring
 * so all possible rings are then found and checked to see if the unfound ring
 * matches any of the unallocated all rings. The sample molecule has just such a mismatching ring.
 * Rings mismatch because there are alternative shortest paths and the two methods choose 
 * a different one to create an equally valid ring.
 *   
 */

public class tester {

		public static void main(String[] args) throws Exception 
	    {
			boolean test;
			int num=0;
			
			try 
			{
				SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
				IMolecule mol = sp.parseSmiles("CCCCCCCC[N+]12CN3CN(CN(C3)C1)C2");
				IMoleculeSet cycles = DefaultChemObjectBuilder.getInstance().newMoleculeSet();
			   
				test=barkerSSSRFinder.findRings(mol,cycles);	   
				
				if(test)
				{
					System.err.println("\nMolecule " + num + " - Cycles found: " + cycles.getAtomContainerCount());
										           		
					int cycletest=testCycles(cycles,mol,cycles.getAtomContainerCount());

               		if(cycletest==1)
               		{
               			System.err.println("MISMATCH - found incorrect rings ");
               		}								
               		else
               		{
               			System.err.println("All cycles are correct ");               			
               		}
				}
				else
				{
					System.err.println("\nError - molecule failed in barkerSSSRFinder: " + num);				
				}
				num++;
			
			} 
			catch (InvalidSmilesException ise) 
			{
			
			}
			 
	    }
				
		private static int checkCycles(IMoleculeSet cycles, IRingSet sssrRings1, 
				int[][] cycle_size, int[][] sssr_cycle_size) {
			// simply tests to see if all cycles are found and identical -if not throws up error
			int match=0;
			int ringsfound=0;
			
			int count=0;
			
			for(IAtomContainer cycle2 : sssrRings1.atomContainers())
			{
				sssr_cycle_size[count][0]=cycle2.getAtomCount();
				sssr_cycle_size[count][1]=0;			
			}
			
			for(IAtomContainer cycle1 : cycles.atomContainers())
			{
				cycle_size[count][0]=cycle1.getAtomCount();
				cycle_size[count][1]=0;
			
				int sssr_count=0;
				for(IAtomContainer cycle2 : sssrRings1.atomContainers())
				{
					match=0;
					
					for(IAtom atom1 : cycle1.atoms())
					{

						for(IAtom atom2 : cycle2.atoms())
						{
							if(atom1.getID()==atom2.getID())
							{
								match++;
								
							}
							
							
						}					
					}
					if(match==cycle1.getAtomCount())
					{
						ringsfound++;
						cycle_size[count][1]=1;
						sssr_cycle_size[sssr_count][1]=1;
					}
					sssr_count++;			
				}
				count++;
			}
			
			if(ringsfound!=cycles.getAtomContainerCount())
			{
				return 1;
			}
			return 0;
		}
		
		private static int testCycles(IMoleculeSet cycles,
				IAtomContainer originalmol, int cauchy) throws IOException, CDKException {
			// tests that the cycles found are the correct cycles
			// first by comparing with the sssr
			// if mismatch then check all rings and see if it's in there
			   		
	  		System.err.print("\nNUMRINGS ");
	  		System.err.print(Integer.toString(cauchy));
	  		System.err.print("\n");              		
			for (int k = 0; k < cauchy; k++) 
			{
	    		System.err.print("RING ");
	      		System.err.print(Integer.toString(k));
	      		System.err.print("\t"); 
	      		
	          	IAtomContainer cycle = cycles.getAtomContainer(k);
	          	                       	
	          	System.err.print(Integer.toString(cycle.getAtomCount()));
	      		System.err.print("   ");
	       		
	          	for(int m=0;m<cycle.getAtomCount();m++) 
	          	{
	          		IAtom atom2 = cycle.getAtom(m);
	          		System.err.print(atom2.getID());  
	          		System.err.print(" ");                 		
	          	}
	           	
	          	System.err.print("\n");   
			}
		
			// check finds same cycles
		    IAtomContainer sssrringatoms = null;
	        IRingSet sssrRings1 = new SSSRFinder(originalmol).findSSSR();
			RingSetManipulator.sort(sssrRings1);
			
	        int mismatch=0;
	        
	        int [][] cycle_size;
	        cycle_size = new int [cycles.getAtomContainerCount()][2];
	        
	        int [][] sssr_cycle_size;
	        sssr_cycle_size = new int [cycles.getAtomContainerCount()][2];   

	        // output num rings
	  	
	        System.err.print("\nSSSR NUMRINGS ");
	        System.err.print(Integer.toString(sssrRings1.getAtomContainerCount()));
	  		System.err.print("\n");  
	     		
	  		List<IAtomContainer> ring1 = RingSetManipulator.getAllAtomContainers(sssrRings1);
	  		for (int j = 0; j < sssrRings1.getAtomContainerCount(); j++) 
	  		{
	  			System.err.print("SSSR RING ");
	  			System.err.print(Integer.toString(j));
	  			System.err.print("\t"); 
	        		
	  			sssrringatoms = ring1.get(j);
	          		
	  			System.err.print(Integer.toString(sssrringatoms.getAtomCount()));
	  			System.err.print("   ");
	         		
	  			for(int m=0;m<sssrringatoms.getAtomCount();m++) 
	  			{	
	  				IAtom atom1 = sssrringatoms.getAtom(m);
	  				System.err.print(atom1.getID());  
	  				System.err.print(" ");                 				
	  			}	             	
            	System.err.print("\n");            	           	    	
	  		}   	

			int test=0;
			test=checkCycles(cycles,sssrRings1,cycle_size,sssr_cycle_size);
			
			if(test==1)
			{
				// find max size
				int size=0;
				int maxsize=0;
				for (int j = 0; j < cycles.getAtomContainerCount(); j++) {
		        		IAtomContainer ring = cycles.getAtomContainer(j);
		        		size = ring.getAtomCount();
		          		if(size>maxsize)
		          			maxsize=size;
				}
				
				AllRingsFinder allRingsFinder = new AllRingsFinder();
				IRingSet allRings1;
				allRingsFinder.setTimeout(30000000);
				allRings1 = allRingsFinder.findAllRings(originalmol,maxsize);
				int found=0;
				int count=0;
				// check if any of rings are not in allrings

				int [][] all_cycle_size;
				all_cycle_size = new int [allRings1.getAtomContainerCount()][2];  
	           
				for(IAtomContainer cycle1 : allRings1.atomContainers())
				{
					all_cycle_size[count][0]=cycle1.getAtomCount();
					all_cycle_size[count][1]=0;
			
					int sssr_count=0;
					for(IAtomContainer cycle2 : sssrRings1.atomContainers())
					{
						int match=0;
					
						for(IAtom atom1 : cycle1.atoms())
						{

							for(IAtom atom2 : cycle2.atoms())
							{
								if(atom1.getID()==atom2.getID())
								{
									match++;
								}	
							}					
						}
						
						if(match==cycle1.getAtomCount())
						{
							all_cycle_size[count][1]=1;
						}
						sssr_count++;

					
					}
					count++;
				}
			   		
				count=0;
	           
				for(IAtomContainer cyclea : cycles.atomContainers())
				{                        	
					found=0;
	           		int all_count=0;
	           		if(cycle_size[count][1]==0)
	           		{           
	           			for(int k=0;k<allRings1.getAtomContainerCount();k++)
	               		{
	           				IAtomContainer cycleb = allRings1.getAtomContainer(k);            		
	           				int incommon=0;
	           				for(IAtom newatom : cyclea.atoms())
	               			{
	               				for(IAtom cycleatom : cycleb.atoms())
	               				{
	               					if(newatom==cycleatom)
	               						incommon++;
	               				}
	               			}
	               			if((all_cycle_size[k][1]==0)&&
	               					(incommon==cyclea.getAtomCount())&&
	               					(incommon==cycleb.getAtomCount()))
	               			{
	               				for(int i=0;i<cycles.getAtomContainerCount();i++)
	               				{
	               					if((cycle_size[i][1]==0)&&(cycle_size[i][0]==all_cycle_size[k][0]))
	               					{
	               						cycle_size[i][1]=1;
	               						all_cycle_size[k][1]=1;
	               						found=1;
	               					}
	               				}
	               			}
	               			all_count++;
	               		}
	             	                    
	                    if(found==0)
	                    {
	                    	mismatch++;
	                    	System.err.println("\nMISMATCH cycle not found molecule ");
	                    }           			
	           		}
	           		count++;
	    
	           }
	           
				// output num rings
	     		System.err.print("\nALLRINGS NUMRINGS ");
	     		System.err.print(Integer.toString(allRings1.getAtomContainerCount()));
	      		System.err.print("\n");  
	         		
	      		List<IAtomContainer> ringa = RingSetManipulator.getAllAtomContainers(allRings1);
	      		for (int j = 0; j < allRings1.getAtomContainerCount(); j++) 
	      		{
	      			System.err.print("ALL RING ");
	      			System.err.print(Integer.toString(j) + "\t");
	        		
	      			IAtomContainer ringatoms = ringa.get(j);
	              		
	      			System.err.print(Integer.toString(ringatoms.getAtomCount()));
	      			System.err.print("   ");
	             		
	      			for(int m=0;m<ringatoms.getAtomCount();m++) 
	      			{
	      				
	      				IAtom atom1 = ringatoms.getAtom(m);
	      				System.err.print(atom1.getID());  
	      				System.err.print(" ");                 		
	      			}
	  				System.err.print("\n");                         			
	      		}
			}
			
	        return mismatch;
		
		}
			
}
