 package ordkilde;

//import java.util.Arrays;

public class MultilayerSNN {

		/*
		 * The neurons variable is three-dimensional. Example: neurons[0][2].length = 4 means that 
		 * the third neuron of the input layer has a spike train of length 4, and neurons[0][2][1] = 7 means
		 * that the second spike of this train was fired at time t = 7.  
		 */
		private int[][][] neurons;
		
		private double[][] potential;
		/*
		 * The weights and delays variables are four-dimensional. In addition to the three dimensions
		 * that are known from traditional neural networks, there is also the fourth one reflecting the fact
		 * that each connection consists of multiple synapses.
		 */
		private double[][][][] weights;
		
		private double[][][][] deltaWeights; 
		 
		private double[][][] biasWeights;
		
		private double[][][] deltaBiasWeights;
		 
		//private int[][] desiredOutputs;
		
		private int[] calculatedOutputs;
		
		private double minWeight;
		
		private double maxWeight;
		
		//private double lowerWeightBoundary = -2.0;
		
		//The bias neurons have a tendency to generate too many spikes, thus this biasfactor to reduce the bias impact. 
		private double biasFactor = 1.0;
		
		private double learnRate = 0.005;
		
		private int iterations = 100;
		
		private int hiddenSpikes;
		private int outputSpikes;
		
		private int biasSpikes;
		
		private int iterationCount = 0;
		
		//The number of connections (or sub-synapses) between any pair of (connected) neurons.
		private int sub_k = 20;
		
		//The time window being used.
		private int timeInterval = 16;
		
		//The threshold potential value for all neurons.
		private double TH = 1.0;
		
		//The following three are time constants used in the calculation of neuron potential.
		private final double TM = 4.0;
		private final double TS = 2.0;
		private final double TR = 20.0;
		
		private final double MINGRAD = 0.2;
		
		private boolean biasEnabled = false;
		
		public int callsForNP = 0;
		public int callsForCDW = 0;
		public double earlyMSError = 0.0;
		public double lateMSError = 0.0;
		private int testIter = 20;
		private boolean printErrorArray = false;

		/**
		 * Constructor for a SNN with no hidden layers.
		 * @param inputs	The number of neurons in the input layer
		 * @param outputs	The number of neurons in the output layer
		 */
		public MultilayerSNN(int inputs, int outputs, int subConnect){
			sub_k = subConnect;
			neurons = new int[2][][];
			neurons[0] = new int[inputs][];
			neurons[1] = new int[outputs][1];
			weights = new double[1][inputs][outputs][sub_k];
			deltaWeights = new double[1][inputs][outputs][sub_k];
			potential = new double[2][outputs];
			for (int i = 0; i < outputs; i++){
				potential[1][i] = 0.0;
			}
			biasWeights = new double[1][outputs][sub_k];
			deltaBiasWeights = new double[1][outputs][sub_k];
			hiddenSpikes = 0; 
			biasSpikes = 0;
		}

		/**
		 * Constructor for a MultiLayerSNN with one or more hidden layer(s).
		 * @param inputs	The number of neurons in the input layer
		 * @param hidden	An array of integers, in which the length og the array indicates the number of hidden layers,
		 * and each element indicates the number of neurons in each layer
		 * @param outputs	The number of neurons in the output layer
		 */
		public MultilayerSNN(int inputs, int[] hidden, int outputs, int subConnect){
			sub_k = subConnect;
			neurons = new int[2 + hidden.length][][];
			neurons[0] = new int[inputs][];
			potential = new double[2 + hidden.length][];
			for (int i = 0; i < hidden.length; i++){
				neurons[i+1] = new int[hidden[i]][1];
				potential[i+1] = new double[hidden[i]];
				for (int j = 0; j < hidden[i]; j++){
					potential[i+1][j] = 0.0;
				}
			} 
			potential[neurons.length -1] = new double[outputs];
			for (int i = 0; i < outputs; i++){
				potential[neurons.length -1][i] = 0.0;
			}
			neurons[neurons.length -1] = new int[outputs][1];//The outputneurons
			weights = new double[1 + hidden.length][][][];
			deltaWeights = new double[1 + hidden.length][][][];
			biasWeights = new double[1 + hidden.length][][];
			deltaBiasWeights = new double[1 + hidden.length][][];
			//delays = new int[1 + hidden.length][][][];
			for (int i = 0; i <= hidden.length; i++){
				weights[i] = new double[neurons[i].length][neurons[i+1].length][sub_k];
				deltaWeights[i] = new double[neurons[i].length][neurons[i+1].length][sub_k];
				biasWeights[i] = new double[neurons[i+1].length][sub_k];
				deltaBiasWeights[i] = new double[neurons[i+1].length][sub_k];
			}
			hiddenSpikes = 0;
			outputSpikes = 0;
			biasSpikes = 0;
		}
		
		/*
		 * Initiates each weight with at random value in the interval between minWeight and maxWeight
		 */
		public void setRandomWeights(){
			for ( int i = 0; i < neurons.length -1; i++){
				for (int j = 0; j < neurons[i].length; j++){
					for (int m = 0; m < neurons[i+1].length; m++){
						for (int n = 0; n < sub_k; n++){
							weights[i][j][m][n] = makeRandomValue(minWeight, maxWeight);
						}
					}
				}
			}
			if (biasEnabled){
				for ( int i = 0; i < neurons.length -1; i++){
					for (int j = 0; j < neurons[i+1].length; j++){
						for (int n = 0; n < sub_k; n++){
							biasWeights[i][j][n] = makeRandomValue(minWeight, maxWeight);
						}
					}
				}
			}
		}
		
		/**
		 * One training cycle, in which the input and the desired output is used to change all the weights 
		 * of the network (each weight once) in the wanted direction.
		 * The errors of the n first and n last (n = testIter) runs are also added up here.
		 * @param inputs	The 2-dimensional array of firing times (spikes) of all the input neurons.
		 * @param desiredOut	The array of firing times of the output neurons, 
		 * counting only the first spike (therefore 1-dimensional only). 
		 */
		public void trainingCycle(int[][] inputs, int[] desiredOut){
			iterationCount++;
			calculatedOutputs = runForward(inputs);
			double squaredError = computeDeltaWeights(desiredOut);
			//normalizeWeights();
			if (iterationCount <= testIter) earlyMSError += squaredError;
			else if (iterationCount > iterations - testIter) lateMSError += squaredError;
		}
		public void trainingCycleLog(int[][] inputs, int[] desiredOut){
			calculatedOutputs = runForward(inputs);
			double squaredError = computeDeltaWeights(desiredOut);
		}
		
		/**
		 * This method trains a complete set of patterns the wanted number of iterations. At each ieration, a pattern is picked at random 
		 * from the total set of patterns.
		 * The early and late error measure is divided by testIter, which is the number of iterations over which we average. 
		 * @param inputs	3-dimensional array of input firing times (An array of the 2-D array described in the above method).
		 * @param desiredOut 	2-dimensional array of output firing times (An array of the 1-D array described in the above method).
		 */
		public void trainPatterns(int[][][] inputs, int[][] desiredOut){
			for (int i = 0; i < iterations; i++){
				int j = (int) Math.floor(Math.random() * inputs.length);
				trainingCycle(inputs[j], desiredOut[j]);
			}
			earlyMSError = earlyMSError / (double)testIter;
			lateMSError = lateMSError / (double)testIter;
		}
		
		//See "multiTrainPatternsLog" method (below) for explanation
		public void trainPatternsLog(int[][][] inputs, int[][] desiredOut, int iters){
			for (int i = 0; i < iters; i++){
				int j = (int) Math.floor(Math.random() * inputs.length);
				trainingCycleLog(inputs[j], desiredOut[j]);
			}
		}
		/**
		 * A wordlist is trained, one entry at a time.
		 * @param wordlist
		 * @param desired  An array depicting desired hyphen positions
		 */
		public void multiTrainPatterns(int[][][][] wordlist, int[][][] desired){
			int i = 0;
			setRandomWeights();
			while(wordlist[i][0] != null){
				trainPatterns(wordlist[i], desired[i]);
				System.out.println("Gjennomført trening for ord nr " + (i+1));
				i++;
			}
		}
		
		//Same as above, but this time we take first 1/2 of the iterations, then 1/4, then 1/8 and so on
		//in order to not let the latter part of the list influence training more than the first part.
		public void multiTrainPatternsLog(int[][][][] wordlist, int[][][] desired){
			int i = 0;
			System.out.println("Runde " + i);
			int iters = iterations/2;
			setRandomWeights();
			while (iters > 10){
				for (int j = 0; j < wordlist.length -1; j++){
					if (wordlist[j][0] != null){
						trainPatternsLog(wordlist[j], desired[j], iters);
						System.out.println("Gjennomført trening for ord nr " + (j+1));
					}
				}
				iters /= 2;
			}
			i++;	
		}
		/**
		 * On basis of the inputs and synaptic weights, computing the potentials for the entire network
		 * @param inputs	The two-dimensional array of input spikes for the network
		 * @return		The computed one-dimensional array of output spikes from the network
		 */
		public int[] runForward(int[][] inputs){
			for (int i = 0; i < inputs.length; i++){
				neurons[0][i] = new int[inputs[i].length];
				for ( int j = 0; j < inputs[i].length; j++){
					neurons[0][i][j] = inputs[i][j]; 
				}
			}
			resetNeurons();//Deleting spikes from previous runs, and setting all potentials to zero.
			for (int t = 1; t <= timeInterval; t++){//For every discrete time step of the interval
				for (int i = 1; i < neurons.length; i++){//For every layer of neurons (except the input layer).
					for ( int j = 0; j < neurons[i].length; j++){// For every neuron of the current layer
						neuronPotential(i, j, t);
					}
				}
			}
			int outputneurons = neurons[neurons.length -1].length; 
			int[] outputArray = new int[outputneurons];
			for (int i = 0; i < outputneurons; i++){
				outputArray[i] = neurons[neurons.length -1][i][0];
				if (outputArray[i] == -1) outputArray[i] = timeInterval+1; 
			}
			return outputArray;
		} 
		
		/**
		 * Computing the "potential" of each neuron (except in the input layer) in the network 
		 * @param i		The layer to which the neuron belongs
		 * @param j		The number of the neuron in layer i, thus identifying the neuron
		 */ 
		private void neuronPotential(int i, int j, int t){
			callsForNP++;
			potential[i][j] = 0.0;
		    if (i == neurons.length -1 && neurons[i][j][0] > 0){
				 return;  //Then this output neuron has already fired a spike (and we count only the first one).
			}
			int spikes = neurons[i][j].length -1; //The number of spikes that the current neuron has fired.
			if (spikes > 0 && t > neurons[i][j][spikes -1]){ //Subtracting the refractory term.
				 potential[i][j] -= TH * Math.exp(-((double)(t- neurons[i][j][spikes -1]))/TR);
			 }
			if (biasEnabled){
				for (int k = 0; k < sub_k; k++){
					double s = (double)(t - k);//A bias neuron always fires at t=0, so we don't need a third term to compute 's'.
					if (s > 0.0){
						potential[i][j] += biasFactor * biasWeights[i-1][j][k] * (Math.exp(-s/TM) - Math.exp(-s/TS));
	  					if (potential[i][j] >= TH){ //if the potential reaches the threshold
	  						if (spikes > 0){
	  							int[] arr = new int[spikes]; 
	  							for (int x = 0; x < spikes; x++){ 
	  								arr[x] = neurons[i][j][x];
	  							}
	  							neurons[i][j] = new int[spikes +1];
	  							for (int x = 0; x < spikes; x++){
	  								neurons[i][j][x] = arr[x];
	  							}
	  							neurons[i][j][spikes] = t;
	  						}else{
	  							neurons[i][j][0] = t;
	  						}
							spikes++;
							if (i == neurons.length -1)
								outputSpikes++;
							else
								biasSpikes++;
							return;
						} 
		 			}
				}
			}
			for (int a = 0; a < neurons[i-1].length; a++){//For every neuron of the preceding layer
				if (neurons[i-1][a][0] > -1){ //if the spike train of the neuron is not empty
					for (int b = 0; b < neurons[i-1][a].length; b++){//For every spike of the above neuron
						for (int k = 0; k < sub_k; k++){//For every sub-connection
							double s = (double) (t - neurons[i-1][a][b] - k);//The variable of the epsilon-function.
							if (s > 0.0){
								if (i > 1 && a == 0)//making the first neuron of a hidden layer inhibitory
									 potential[i][j] -= weights[i-1][a][j][k] * (Math.exp(-s/TM) - Math.exp(-s/TS));
								else // The other hidden layer neurons are all excitatory.
									 potential[i][j] += weights[i-1][a][j][k] * (Math.exp(-s/TM) - Math.exp(-s/TS)); 
			  					if (potential[i][j] >= TH){ //if the potential reaches the threshold
			  						if (spikes > 0){//The size of the array of spikes needs to be incremented.
			  							int[] arr = new int[spikes]; 
			  							for (int x = 0; x < spikes; x++){
			  								arr[x] = neurons[i][j][x];
			  							}
			  							neurons[i][j] = new int[spikes +1];
			  							for (int x = 0; x < spikes; x++){
			  								neurons[i][j][x] = arr[x];
			  				 			}
			  							neurons[i][j][spikes] = t;
			  						}else{//Then this is the first spike of this neuron, and the spike array already has size of 1.
			  							neurons[i][j][0] = t;
			  						}
									spikes++;
									if (i == neurons.length -1)
										outputSpikes++;
									else
										hiddenSpikes++;
									return;
								}
				 			}
						}   
					}
				}
			}
		}

		/**
		 * This is the learning algorithm for a multi-layered SNN.
		 * @param desiredOut The desired results for the output layer
		 * @return The mean squared error of the network  
		 */		
		private double computeDeltaWeights(int[] desiredOut){
			callsForCDW++;
			double[] errorArray = new double[desiredOut.length];
			double meanSqErr, mse = 0.0; 
			if (iterationCount == testIter) System.out.println();
			for ( int n = 0; n < desiredOut.length; n++){
				errorArray[n] = (double)(calculatedOutputs[n] - desiredOut[n]);	
				mse += (errorArray[n] * errorArray[n]);
			}
			meanSqErr = mse / (double)desiredOut.length;
			if (printErrorArray &&(iterationCount < testIter || iterationCount >= iterations -testIter)){
				System.out.print("Iteration " +  iterationCount +" Error array: ");
				for ( int n = 0; n < desiredOut.length; n++){
					System.out.print(errorArray[n] + "  "); 
				}
				System.out.println();
			}
			//First, we will compute the weight change from the output neuron(s) to their preceding layer (Ca 45 following code lines)
			int i = neurons.length -1;
			for (int j = 0; j < neurons[i].length; j++){//For every neuron of the output layer
				double gradient = 0.0;// The denominator of eq. 3.13
				double s = 0.0; //This is the variable for the epsilon-function.
				for (int a = 0; a < neurons[i-1].length; a++){//For every neuron of the preceding layer
					for (int b = 0; b < neurons[i-1][a].length; b++){//For every spike of the neuron
						for (int k = 0; k < sub_k; k++){
							 s = (double)(neurons[i][j][0] - neurons[i-1][a][b] - k);
							 if (s > 0.0)
							 gradient += weights[i-1][a][j][k] * ((-1.0/TM)* Math.exp(-s/TM) + (1.0/TS) * Math.exp(-s/TS));
						}
					}
				} 
				if (biasEnabled){
					for (int k = 0; k < sub_k; k++){
						 s = (double)(neurons[i][j][0] - k);
						 if (s > 0.0)
						 gradient += biasWeights[i-1][j][k] * ((-1.0/TM)* Math.exp(-s/TM) + (1.0/TS) * Math.exp(-s/TS));
					}
				}
				if (Math.abs(gradient) < MINGRAD){ //To avoid the gradient getting too close to 0.
					if (gradient >= 0.0) gradient = MINGRAD;
					else gradient = -MINGRAD;
				} 
				for (int a = 0; a < neurons[i-1].length; a++){ //Running through the neurons of the preceding layer.
					for (int k = 0; k < sub_k; k++){//for the delay of every sub-connection
						double numer = 0.0; //The numerator of equation 3.13
						for (int b = 0; b < neurons[i][j].length; b++){//For every spike of the neuron
							s = (double)(neurons[i][j][0] - neurons[i-1][a][b] - k);
							if (s > 0.0)
								numer -= (Math.exp(-s/TM) - Math.exp(-s/TS));
						}
						deltaWeights[i-1][a][j][k] = -errorArray[j] * learnRate * (numer/gradient);
					}
				}
				if (biasEnabled){
					for (int k = 0; k < sub_k; k++){
						s = (double)(neurons[i][j][0] - k);
						if (s > 0.0){
							double numer = -(Math.exp(-s/TM) - Math.exp(-s/TS));
							deltaBiasWeights[i-1][j][k] = -errorArray[j] * learnRate * (numer/gradient);
						}
					}
				}
			} 
			//Then, we will compute the weight changes between hidden layers (if more than one) and between a hidden layer and the input layer.
			if (neurons.length > 2){//If there's (one or more) hidden layer(s)
				for (int m = neurons.length -2; m > 0; m--){ //for each hidden layer
					for (int n = 0; n < neurons[m].length; n++){ //for every neuron of the hidden layer
						double previousSW = 0.0; //Saving the spike-weight partial derivative from the previous spike
						for (int p = 0; p < neurons[m][n].length; p++){ //for every spike of the neuron 
							double pdErrorSpike = 0.0; //Partial derivative of net error with regard to spike (Eq. 4.11)
							for (int q = 0; q < neurons[m+1].length; q++){//For every neuron of the subsequent layer
				 				double num = 0.0, denom = 0.0; // Numerator and denominator of eq. 4.11
								for (int k = 0; k < sub_k; k++){
						 			double s = (double)(calculatedOutputs[q] - neurons[m][n][p] - k); //Kernel of epsilon function
						 			if (s > 0.0)
									num += weights[m][n][q][k] * ((-1.0/TM)* Math.exp(-s/TM) + (1.0/TS) * Math.exp(-s/TS));
								} 
								for (int nn = 0; nn < neurons[m].length; nn++){
									for (int r = neurons[m][n].length -1; r >= 0; r--){
										for (int k = 0; k < sub_k; k++){
											double s = (double)(calculatedOutputs[q] - neurons[m][n][r] - k);
											if (s > 0.0)
											denom += weights[m][n][q][k] * ((-1.0/TM)* Math.exp(-s/TM) + (1.0/TS) * Math.exp(-s/TS));
				 						}
									}
								 }
								if (Math.abs(denom) < MINGRAD){//Avoiding that the denominator gets too close to zero.
									if (denom >= 0.0) denom = MINGRAD;
									else denom = -MINGRAD;
								} 
								pdErrorSpike += (num / denom) * errorArray[q]; 
							}
							double pdSpikeWeight = 0.0; //Derivative of spike with regard to weight (eq. 4.16)
							double num = 0.0, denom = 0.0; //Numerator and denominator of eq. 4.16
							for (int q =  0; q < neurons[m-1].length; q++){ //for every neuron of the previous layer
								for (int k = 0; k < sub_k; k++){ //for every weight of the sub-connections
									for (int r = 0; r < neurons[m-1][q].length; r++){ //for every spike of neuron (of previous layer)
										double s = (double)(neurons[m][n][p] - neurons[m-1][q][r] -k);
				 						num -= (Math.exp(-s/TM) - Math.exp(-s/TS));
										for (int x = 0; x < neurons[m-1].length; x++){//for each neuron h of the previous layer
											for (int y = 0; y < sub_k; y++){//for each sub-connection
												s = (double)(neurons[m][n][p] - neurons[m-1][q][r] -y);
												if (s > 0.0)
												denom += weights[m-1][x][n][y] * ((-1.0/TM)* Math.exp(-s/TM) + (1.0/TS) * Math.exp(-s/TS));
											}
										}  
									}
									if (p > 0){ //We also need the refractory terms if this is not the first spike
										double sr = (double)(neurons[m][n][p] - neurons[m][n][p-1]);//The kernel of the eta-function
										num += (TH /TR) * Math.exp(-sr/TR) * previousSW; //Mult. w/ previous "pdSpikeWeight" since eq. 4.16 is recursive
										denom += (TH /TR) * Math.exp(-sr/TR);
									}
									if (Math.abs(denom) < MINGRAD){
										if (denom >= 0.0) denom = MINGRAD;
										else denom = -MINGRAD;
									}
									pdSpikeWeight = num / denom;//Partial derivative of spiketime with regard to weight
									/*if (iterationCount == iterations -1){
										System.out.print("pdSpikeWeight: ");
										System.out.format("%.5f%n", pdSpikeWeight);
									}*/ 
									deltaWeights[m-1][q][n][k] -= learnRate* pdSpikeWeight * pdErrorSpike;
									previousSW = pdSpikeWeight;
								}
							}
							if (biasEnabled){
								for (int k = 0; k < sub_k; k++){
									double s = (double)(neurons[m][n][p] -k); //Each biasneuron fires only once, and always at t = 0. 
									num -= (Math.exp(-s/TM) - Math.exp(-s/TS));
									for (int y = 0; y < sub_k; y++){//for each sub-connection
										s = (double)(neurons[m][n][p] -y);
										if (s > 0.0)
										denom += biasWeights[m-1][n][y] * ((-1.0/TM)* Math.exp(-s/TM) + (1.0/TS) * Math.exp(-s/TS));
									}
									if (Math.abs(denom) < MINGRAD){
										if (denom >= 0.0) denom = MINGRAD;
										else denom = -MINGRAD;
									}
									pdSpikeWeight = num / denom;
									deltaBiasWeights[m-1][n][k] -= learnRate* pdSpikeWeight * pdErrorSpike;
								}
							}
						}  
					}
				}
			}
			//if (iterationCount %10 == 0) printDeltaWeights();
			updateWeights();
			return meanSqErr;
		}
		
		public int getTestIter() {
			return testIter;
		}

		public void setTestIter(int testIter) {
			this.testIter = testIter;
		}

		private void updateWeights(){
			for (int m = 0; m < neurons.length -1; m++){
				for (int n = 0; n < neurons[m].length; n++){
					for (int p = 0; p < neurons[m+1].length; p++){
						for (int k = 0; k < sub_k; k++){
							//System.out.format("%.3f%n", deltaWeights[m][n][p][k]);
							weights[m][n][p][k] += deltaWeights[m][n][p][k];
							//Weights should not be negative (acc. to Bohte)
							//if (weights[m][n][p][k] < 0.0)
								//weights[m][n][p][k] = - weights[m][n][p][k];
							deltaWeights[m][n][p][k] = 0.0;
						}
					}
				}
			}
			if (biasEnabled){
				for (int m = 0; m < neurons.length -1; m++){
					for (int p = 0; p < neurons[m+1].length; p++){
						for (int k = 0; k < sub_k; k++){
							biasWeights[m][p][k] += deltaBiasWeights[m][p][k];
							if (biasWeights[m][p][k] < 0.0)
								biasWeights[m][p][k] = - biasWeights[m][p][k];
							deltaBiasWeights[m][p][k] = 0.0;
						}
					}
				}
			}
		}
		
		private void normalizeWeights(){
			double maxSum = 30.0;
			for (int m = 0; m < neurons.length -1; m++){
				for (int n = 0; n < neurons[m].length; n++){
					for (int p = 0; p < neurons[m+1].length; p++){
						double subweightsSum = 0.0;
						for (int k = 0; k < sub_k; k++){
							subweightsSum += weights[m][n][p][k];
						}
						if (subweightsSum > maxSum || subweightsSum < -maxSum){
							for (int k = 0; k < sub_k; k++){
								weights[m][n][p][k] = (weights[m][n][p][k] / subweightsSum) * sub_k; 
							}
						}
					}
				}
			}
			if (biasEnabled){
				for (int m = 0; m < neurons.length -1; m++){
					for (int p = 0; p < neurons[m+1].length; p++){
						double subweightsSum = 0.0;
						for (int k = 0; k < sub_k; k++){
							subweightsSum += biasWeights[m][p][k];
						}
						if (subweightsSum > maxSum || subweightsSum < -maxSum){
							for (int k = 0; k < sub_k; k++){
								biasWeights[m][p][k] = (biasWeights[m][p][k] / subweightsSum) * sub_k; 
							}
						}
					}
				}
			}
		}
		/**
		 * The method deletes all previous spikes from the neurons, giving each neuron a spike array of size 1, 
		 * and the value -1, showing that no spike has been fired yet.
		 * The potential of each neuron is likewise nullified.
		 */
		private void resetNeurons(){
			for (int i = 1; i < neurons.length; i++){
				for (int j = 0; j < neurons[i].length; j++){
					neurons[i][j] = new int[1];
					neurons[i][j][0] = -1;
					potential[i][j] = 0.0;	
				}
			}
		}
		/**
		 *Just for testing the weights. Prints out each weight. 
		 */
		public void printWeights(){
			for (int i = 0; i < weights.length; i++){
				for (int j = 0; j < weights[i].length; j++){
					for (int m = 0; m < weights[i][j].length; m++){
						System.out.println();
						for (int n = 0; n < sub_k; n++){
							System.out.format("%.3f%n", weights[i][j][m][n]);
						}
					}
				}
			}
			//System.out.println("Number of weights: " + printedWeights);
		}
		
		/**
		 * @return	The number of weights
		 */
		public int getNumberOfWeights(){
			int weightsNum = 0;
			for (int i = 0; i < weights.length; i++){
				for (int j = 0; j < weights[i].length; j++){
					for (int m = 0; m < weights[i][j].length; m++){
						weightsNum += sub_k;
					}
				}
			}
			return weightsNum;
		}
		
		/**
		 * @return	The number of biasWeights
		 */
		public int getNumberOfBiasWeights(){
			int biasWeightsNum = 0;
			for (int i = 0; i < biasWeights.length; i++){
				for (int j = 0; j < biasWeights[i].length; j++){
					biasWeightsNum += sub_k;
				}
			}
			return biasWeightsNum;
		}
		
		//Also just for testing purposes
		public void printPotential(){
			for (int i = 0; i < potential.length; i++){
				for (int j = 0; j < potential[i+1].length; j++){
					System.out.format("%.3f%n", potential[i+1][j]);
				}
			}
		}
		/**
		 * @return The weights arranged in a one-dimensional array
		 */
		public double[] get1DWeightsArray(){
			double[] arr = new double[getNumberOfWeights()];
			int a = 0;
			for (int i = 0; i < weights.length; i++){
				for (int j = 0; j < weights[i].length; j++){
					for (int m = 0; m < weights[i][j].length; m++){
						for (int n = 0; n < sub_k; n++){
							arr[a] = weights[i][j][m][n];
							a++;
						}
					}
				}
			}
			return arr;
		}
		
		//Same as above for the biasWeights.
		public double[] get1DBiasWeightsArray(){
			double[] arr = new double[getNumberOfBiasWeights()];
			int a = 0;
			for (int i = 0; i < biasWeights.length; i++){
				for (int j = 0; j < biasWeights[i].length; j++){
					for (int m = 0; m < sub_k; m++){
						arr[a] = biasWeights[i][j][m];
						a++;
					}
				}
			}
			return arr;
		}
	
		//For testing purposes
		public void printDeltaWeights(){
			System.out.println();
			System.out.println("Deltaweights for iteration " + iterationCount + ": ");
			double[] arr = new double[getNumberOfWeights()];
			int a = 0;
			for (int i = 0; i < deltaWeights.length; i++){
				for (int j = 0; j < deltaWeights[i].length; j++){
					for (int m = 0; m < deltaWeights[i][j].length; m++){
						for (int n = 0; n < sub_k; n++){
							arr[a] = deltaWeights[i][j][m][n];
							if (Math.abs(arr[a]) > 0.0){
								System.out.print("Deltaweights " + (a+1) + ": ");
								System.out.format("%.7f%n", arr[a]);
							}
							a++;
						}
						System.out.println();
					}
				}
			}			
		}
		
		//public int[][] makeSpikeTrain()
		
		//Setting the potential threshold for spikes
		public void setThreshold(double threshold){
			TH = threshold;
		}
		
		public double getminWeight() {
			return minWeight;
		}

		public void setminWeight(double minWeight) {
			this.minWeight = minWeight;
		}

		public double getmaxWeight() {
			return maxWeight;
		}

		public void setmaxWeight(double maxWeight) {
			this.maxWeight = maxWeight;
		}

		public void setBiasEnabled(boolean enabled){
			biasEnabled = enabled;
		}
		
		public boolean getBiasEnabled(){
			return biasEnabled;
		}
		
		public int getIterations() {
			return iterations;
		}

		public void setIterations(int iterations) {
			this.iterations = iterations;
		}
		
		public int getIterationCount() {
			return iterationCount;
		}

		public void setIterationCount(int iterationCount) {
			this.iterationCount = iterationCount;
		} 

		public double getLearnRate() {
			return learnRate;
		}

		public void setLearnRate(double learnRate) {
			this.learnRate = learnRate;
		}

		private double makeRandomValue(double lower, double upper){
			return Math.random() * (upper - lower) + lower;
		}
		
		public void setBiasFactor(double biasFactor){
			this.biasFactor = biasFactor;
		}
		
		public double getBiasFactor(){
			return biasFactor; 
		}
		
		public int getHiddenSpikes() {
			return hiddenSpikes;
		}
		public int getOutputSpikes() {
			return outputSpikes;
		}
		public void setHiddenSpikes(int hiddenSpikes) {
			this.hiddenSpikes = hiddenSpikes;
		}

		public int getBiasSpikes() {
			return biasSpikes;
		}

		public void setBiasSpikes(int biasSpikes) {
			this.biasSpikes = biasSpikes;
		}
		
		public int getSub_k(){
			return sub_k;
		}
		
		public void setPrintErrorArray(boolean printErrorArray) {
			this.printErrorArray = printErrorArray;
		}
		public void setTimeInterval(int t){
			timeInterval = t;
		}
	}


