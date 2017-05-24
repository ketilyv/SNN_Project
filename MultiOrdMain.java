package ordkilde;

import java.io.*;

public class MultiOrdMain {
	//These are the early and late target spike times, indicating respectively good and bad hyphen positions.
	static int early = 22;
	static int late = 30;
	//The number of input neurons, indicating also maximal number of characters in a word 
	//before it has to be "divided" by the window procedure (see thesis chap. 6).
	static int window = 9;
	static int arrSize = 0;
	//The wordlist on which the network should be trained.
	static String filnavn = "C:/Users/rulle/workspace/Orddeling/src/ordkilde/Orddel_rand.txt";
	static LetterArrays la;
	static String[] stringArrayNoHyph;
	static int[][][] desArray;;
	static int[][][][] wordlist;
	
	//Denne settes "true" dersom man ønsker å teste med et fritt valgt ord (f.eks. et som ikke er fra treningsmengden),
	//og "false" dersom man vil teste med et ord fra denne mengden.
	static boolean unknownWord = true;
	
	public static void main(String[] args) throws IOException {
		//The word on which the network should be tested, in order to find hyphen postions.
		//In this version: max. 8 letters.
		String testStr = "BEFRIR";
		int[] hid = {5};
		MultilayerSNN snn = new MultilayerSNN(window, 1, 35);
		snn.setminWeight(0.0);  
		snn.setmaxWeight(1.5); 
		
		snn.setTimeInterval(36);
		snn.setThreshold(1.0);
		snn.setLearnRate(0.01);
		snn.setIterations(800);
		snn.setBiasEnabled(true);
		snn.setBiasFactor(1.0);
		snn.setPrintErrorArray(false); 
		//snn.setRandomWeights(); //Set in the MultiTrainPatterns method
		wordlist = makeWordlistSpikes();
		
		long start = System.currentTimeMillis();
		snn.multiTrainPatternsLog(wordlist, desArray);
		long end = System.currentTimeMillis(); 
		double sek = (double)(end - start) /1000.0; 
		//Her velger man hvilket ord fra trenings-ordlisten man vil bruke til testing. 
		int wordnr = 9;
		
		if (unknownWord){
			//Her velger man fritt et test-ord (men foreløpig begrenset til maks. 8 bokstaver).
			int[][][] ia3 = makeSpikesFromUnknownWord(testStr);	
			for (int i = 0; i < testStr.length()-1; i++){
				int[] svar = snn.runForward(ia3[i]);
				System.out.println("Plass nr " + (i+1) + ": " + svar[0]);
			}
		}else{
			for (int i = 0; i < desArray[wordnr-1].length; i++){
				int[] svar = snn.runForward(wordlist[wordnr][i]);
				System.out.println("Plass nr " + (i+1) + ": " + svar[0]); 
			}
		}
		//System.out.println("Mean squared error start: " + snn.earlyMSError);
		//System.out.println("Mean squared error end: " + snn.lateMSError);
		System.out.println("Runtime for pattern training: " + sek + " seconds.");
	}
	
	/**
	 * Turns the words into simulated spike times, being the only acceptable input for the system.
	 * @return  A 4-dimensional array of spike times.
	 * @throws IOException
	 */
	public static int[][][][] makeWordlistSpikes() throws IOException{
		la = new LetterArrays(window);
		File f = new File(filnavn);
		FileReader fr = new FileReader(f);
		BufferedReader br = new BufferedReader(fr);
		int lines = countLinesAndPlus8(filnavn, window);
		stringArrayNoHyph = new String[lines]; 
		int[][][][] wordArray = new int[lines][window][][];
		desArray = new int[lines][][];
		int i = 0;
		String s;
		while ((s = br.readLine()) != null){
			if(s.length() > window){
				for (int d = 0; d <= (s.length() - window); d++){
					String sub = s.substring(d, window +d);
					if (sub.charAt(0) != '-' && sub.charAt(window -1) != '-'){
						int[] hyphens = la.makeHyphenArray(sub);
						System.out.print("Entry no: " + (i+1) + "  ");
						System.out.println(sub);
						System.out.print("Hyphen array: ");
						la.print1DArray(hyphens);
						if (sub.length() <= window){
							sub = la.takeAwayHyphens(sub);
							stringArrayNoHyph[i] = sub;
							int[][] ia = la.makeSpikeTrainsMaj(sub, window); 
							//System.out.println("2DArray length: " + ia.length);
							wordArray[i] = la.makeHyphPatterns(ia);
							desArray[i] = new int[sub.length()-1][1]; 
							for (int j = 0; j < desArray[i].length; j++){
								desArray[i][j][0] = late;
							}
							for (int j = 0; j < hyphens.length; j++){
								desArray[i][hyphens[j]-(j+1)][0] = early; 
							}
							i++;
						}else{
							System.out.println("Substring too long!");
						}
					}
				}
				System.out.println("Longer word treated.");
			}else{
				System.out.print("Entry no: " + (i+1) + "  ");
				System.out.println(s);
				int[] hyphens = la.makeHyphenArray(s);
				//System.out.print("Hyphen array: ");
				//la.print1DArray(hyphens);
				/*if (s.length() < window-1){
					for (int j:hyphens){
						hyphens[j] = hyphens[j] +1;
					}
				}*/
				s = la.takeAwayHyphens(s);
				stringArrayNoHyph[i] = s;
				int[][] ia = la.makeSpikeTrainsMaj(s, window); 				
				wordArray[i] = la.makeHyphPatterns(ia);
				
				desArray[i] = new int[window -1][1];
				for (int j = 0; j < desArray[i].length; j++){
					desArray[i][j][0] = late;
				}
				for (int j = 0; j < hyphens.length; j++){
					desArray[i][hyphens[j]][0] = early; 
				}
				i++;
			}
		}
		br.close();
		return wordArray;
	}
	/**
	 * Counting the lines in the input textfile, but counting also the extra characters if a word
	 * has more than 8 (or "window" size -1) characters.  
	 * @param filename
	 * @param window If an input string exceeds this size, we need to count the extra characters.
	 * @return The total number of lines and extra characters.
	 * @throws IOException
	 */
	public static int countLinesAndPlus8(String filename, int window) throws IOException{
		File f = new File(filename);
		FileReader fr = new FileReader(f);
		BufferedReader br = new BufferedReader(fr);
		String s;
		int i = 0;
		while ((s = br.readLine()) != null){
			if (s.length() > window){
				i += s.length() - window;
			}
			i++;
		}
		br.close();
		System.out.println("Word array size: " + i);
		return i;
	}  
	
	/**
	 * Makes a pattern of simulated spikes from a word (String)
	 * that is previouslu unknown for the system.
	 * @param uw  A String input
	 * @return  A three-dimensional array of spikes. 
	 */
	public static int[][][] makeSpikesFromUnknownWord(String uw){
		LetterArrays la = new LetterArrays(window);
		int[][] ia2 = la.makeSpikeTrainsMaj(uw, window);
		int[][][] ia3 = la.makeHyphPatterns(ia2);
		return ia3;
	}
}
