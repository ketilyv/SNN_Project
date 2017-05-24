package ordkilde;
/**
 * @author rulle
 * Det sentrale i denne klassen er en to-dimensjonal tabell (kalt 'abc')
 * der hver bokstav tildeles en sekvens av "spikes" (altså et "spike-tog").
 * Hver bokstav tildeles altså en unik sekvens bestående av min. 1 og max. 5 spikes. 
 */
public class LetterArrays {
	//The following integers are neuron firing times (called 'spikes').
	private int s1 = 1;
	private int s2 = 5;
	private int s3 = 9;
	private int s4 = 13;
	private int s5 = 17;
	
	private int[][] abc;
	private int window;
	
	public LetterArrays(int win){
		 abc = new int[31][];
		 abc[0] = new int[1];//a
		 abc[0][0] = s1;
		 abc[1] = new int[1];//b
		 abc[1][0] = s2;
		 abc[2] = new int[1];//c
		 abc[2][0] = s3;
		 abc[3] = new int[1];//d
		 abc[3][0] = s4;
		 abc[4] = new int[1];//e
		 abc[4][0] = s5;
		 abc[5] = new int[2];//f
		 abc[5][0] = s1;
		 abc[5][1] = s2;
		 abc[6] = new int[2];//g
		 abc[6][0] = s1;
		 abc[6][1] = s3;
		 abc[7] = new int[2];//h
		 abc[7][0] = s1;
		 abc[7][1] = s4;
		 abc[8] = new int[2];//i
		 abc[8][0] = s1;
		 abc[8][1] = s5;
		 abc[9] = new int[2];//j
		 abc[9][0] = s2;
		 abc[9][1] = s3;
		 abc[10] = new int[2];//k
		 abc[10][0] = s2;
		 abc[10][1] = s4;
		 abc[11] = new int[2];//l
		 abc[11][0] = s2;
		 abc[11][1] = s5;
		 abc[12] = new int[2];//m
		 abc[12][0] = s3;
		 abc[12][1] = s4;
		 abc[13] = new int[2];//n
		 abc[13][0] = s3;
		 abc[13][1] = s5;
		 abc[14] = new int[2];//o
		 abc[14][0] = s4;
		 abc[14][1] = s5;
		 abc[15] = new int[3];//p
		 abc[15][0] = s1;
		 abc[15][1] = s2;
		 abc[15][2] = s3;
		 abc[16] = new int[3];//q
		 abc[16][0] = s1;
		 abc[16][1] = s2;
		 abc[16][2] = s4;
		 abc[17] = new int[3];//r
		 abc[17][0] = s1;
		 abc[17][1] = s2;
		 abc[17][2] = s5;
		 abc[18] = new int[4];//s
		 abc[18][0] = s1;
		 abc[18][1] = s2;
		 abc[18][2] = s3;
		 abc[18][3] = s5;
		 abc[19] = new int[4];//t
		 abc[19][0] = s1;
		 abc[19][1] = s2;
		 abc[19][2] = s3;
		 abc[19][3] = s4;
		 abc[20] = new int[4];//u
		 abc[20][0] = s1;
		 abc[20][1] = s2;
		 abc[20][2] = s4;
		 abc[20][3] = s5;
		 abc[21] = new int[5];//v
		 abc[21][0] = s1;
		 abc[21][1] = s2;
		 abc[21][2] = s3;
		 abc[21][3] = s4;
		 abc[21][4] = s5;
		 abc[22] = new int[3];//w
		 abc[22][0] = s1;
		 abc[22][1] = s3;
		 abc[22][2] = s4;
		 abc[23] = new int[3];//x
		 abc[23][0] = s1;
		 abc[23][1] = s3;
		 abc[23][2] = s5;
		 abc[24] = new int[3];//y
		 abc[24][0] = s1;
		 abc[24][1] = s4;
		 abc[24][2] = s5;
		 abc[25] = new int[2];//z
		 abc[25][0] = s1;
		 abc[25][1] = s5;
		 abc[26] = new int[3];//æ
		 abc[26][0] = s2;
		 abc[26][1] = s3;
		 abc[26][2] = s4;
		 abc[27] = new int[3];//ø
		 abc[27][0] = s2;
		 abc[27][1] = s3;
		 abc[27][2] = s5;
		 abc[28] = new int[3];//å
		 abc[28][0] = s2;
		 abc[28][1] = s4;
		 abc[28][2] = s5;
		 abc[29] = new int[3];//-
		 abc[29][0] = s3;
		 abc[29][1] = s4;
		 abc[29][2] = s5;
		 //The empty space corresponds to a neuron that does not fire at all.
		 abc[30] = new int[1];//Empty space
		 abc[30][0] = -1;
	}
	
	public int[][] getABC(){
		return abc;
	}
	
	/**
	 * This method converts strings of minuscules to integer arrays.
	 * @param str The string that should be converted
	 * @return Two-dimensional array that is a representation of the given string
	 */
	public int[][] makeSpikeTrainsMin(String str){
		int[][] spikes = new int[str.length()][];
		for (int i = 0; i < str.length(); i++){
			Character c = (Character)str.charAt(i);
			int ci = (int)c.charValue();
			if (ci > 96 && ci < 123){
				spikes[i] = abc[ci -97];
			}
			else if (ci == 45){
				spikes[i] = abc[29];
			}
			else if (ci == 230){
				spikes[i] = abc[26];
			}
			else if (ci == 248){
				spikes[i] = abc[27];
			} 
			else if (ci == 229){
				spikes[i] = abc[28];
			}else{
				spikes[i] = new int[1];
				spikes[i][0] = -1;
			}
			
		}
		return spikes;
	}
	/**
	 * This method converts strings of majuscules to integer arrays.
	 * @param str The string that should be converted
	 * @return Two-dimensional array that is a representation of the given string
	 */
	public int[][] makeSpikeTrainsMaj(String str, int win){
		int[][] spikes = new int[str.length()][];
		for (int i = 0; i < str.length(); i++){
			Character c = (Character)str.charAt(i);
			int ci = (int)c.charValue();
			if (ci > 9 && ci < 36){
				spikes[i] = abc[ci -10];
			}
			else if (ci == 45){
				spikes[i] = abc[29];
			}
			else if (ci == 198){
				spikes[i] = abc[26];
			}
			else if (ci == 216){
				spikes[i] = abc[27];
			} 
			else if (ci == 197){
				spikes[i] = abc[28];
			}else{
				spikes[i] = new int[1];
				spikes[i][0] = -1;
			}
		}
		return spikes;
	}
	/**
	 * This method makes every 'word' into an array of 'words', in which a hyphen is placed in all possible positions.
	 * @param word A two-dimensional array produced by any of the 'makeSpikeTrains()' methods.
	 * @return A three-dimnesional integer array representing the different hyphenizations.
	 */
	public int[][][] makeHyphPatterns(int[][] word){
		int[][][] mlr = new int[word.length -1][][];
		for (int i = 0; i < mlr.length;i++){
			mlr[i] = new int[word.length +1][];
			mlr[i][0] = word[0];
			mlr[i][word.length] = word[word.length -1];
			for (int j = 1; j < word.length; j++){
				if (j == i+1){
					mlr[i][j] = abc[29];
				}else if (j > i+1){
					mlr[i][j] = word[j-1];
				}else{
					mlr[i][j] = word[j];
				}
			}
		}
		if (mlr.length < (window -2)){
			int[][][] ml2 = new int[mlr.length][window][];
			for (int i = 0; i < ml2.length; i++){
				ml2[i][0] = new int[1];
				ml2[i][0][0] = -1; 
				for (int j = 1; j <= mlr[0].length; j++){
					ml2[i][j] = mlr[i][j-1];
				}
				for (int j = mlr[0].length +1 ; j < window; j++){
					ml2[i][j] = new int[1];
					ml2[i][j][0] = -1;
				}
			}
			//System.out.println("Returning ml2");
			return ml2;
		} else{
		//System.out.println("returning mlr");
		return mlr;
		}
	}
	
	//Checking integer value of characters	
	public void printcharValue(char c){
		int ci = c; //Character.getNumericValue(c);
		System.out.println(ci);
	}
	//Also for checking purposes
	public void printSpikeTrains(int[][] sTrains){
		for (int i = 0; i < sTrains.length; i++){
			for (int j = 0; j < sTrains[i].length; j++){
				System.out.print(sTrains[i][j] + "  ");
			}
			System.out.println();
		}
	}
	
	public int[] makeHyphenArray(String s){
		int[] ha = new int[1];
		int hyphens = 0;
		for (int i = 1; i < s.length() -1; i++){
			if (s.charAt(i) == '-'){
				if (hyphens+1 > ha.length){
					int[] temp = ha;
					ha = new int[temp.length +1];
					for (int j = 0; j < temp.length; j++)	{
						ha[j] = temp[j];
					}
				}
				ha[hyphens] = i;
				hyphens++;
			}
		}
		return ha;
	}
	
	public String takeAwayHyphens(String s){
		int[] ha = makeHyphenArray(s);
		if (ha[0] < 1){
			return s;
		} else {
			int[] ha2 = new int[ha.length+1];
			ha2[ha.length] = s.length();
			for (int i = 0; i < ha.length; i++){
				ha2[i] = ha[i];
			}
			String sub1 = s.substring(0, ha[0]);
			if (ha.length == 1){
				sub1+= s.substring(ha[0] +1); 
			} else { 
				for (int j = 0; j < ha.length; j++){
					String sub2 = s.substring(ha2[j]+1, ha2[j+1]);
					sub1 += sub2;
				}
			}
			return sub1;
		}
	}
	
	public void print1DArray(int[] arr){
		for (int i = 0; i < arr.length; i++){
			System.out.print(arr[i] + "  ");
		}
		System.out.println();
	}
}
