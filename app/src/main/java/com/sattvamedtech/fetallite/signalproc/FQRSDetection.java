package com.sattvamedtech.fetallite.signalproc;

public class FQRSDetection {
	private static int mLength;
	MatrixFunctions mMatrixFunctions = new MatrixFunctions();
	static boolean i1Flag = false;
	static boolean i2Flag = false;
	static boolean i3Flag = false;
	static boolean i4Flag = false;

	double[] aChannel1;
	double[] aChannel2;
	double[] aChannel3;
	double[] aChannel4;

	double[][] mInput;

	int[] mQRS1, mQRS2, mQRS3, mQRS4;
	/**
	 * 
	 * @param iInput
	 * @param iQrsM
	 * @param iInterpolatedLength
	 * @param iQRSLast
	 * @param iRRMeanLast
	 * @param iNoDetectionFlag
	 * @return Object[] { aQRSFinal, aInterpolatedLength, aNoDetectionFLag };
	 * @throws Exception
	 */
	public Object[] fQRS(double[][] iInput, int[] iQrsM, int iInterpolatedLength, int iQRSLast, double iRRMeanLast,
			int iNoDetectionFlag) throws Exception {

		mLength = iInput.length;

		mInput = iInput;
		aChannel1 = new double[mLength];
		aChannel2 = new double[mLength];
		aChannel3 = new double[mLength];
		aChannel4 = new double[mLength];

//		for (int i = 0; i < mLength; i++) {
//			aChannel1[i] = iInput[i][0];
//			aChannel2[i] = iInput[i][1];
//			aChannel3[i] = iInput[i][2];
//			aChannel4[i] = iInput[i][3];
//		}

//		int qrs1[] = fetalQRS(channel1);
//		int qrs2[] = fetalQRS(channel2);
//		int qrs3[] = fetalQRS(channel3);
//		int qrs4[] = fetalQRS(channel4);
		
		Thread aQRSDet1 = new Thread(new Runnable() {
			
			@Override
			public void run() {
				for (int i = 0; i < mLength; i++) {
					aChannel1[i] = mInput[i][0];
				}
				try {
					mQRS1 = fetalQRS(aChannel1);
				} catch (Exception e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				i1Flag = true;
				
			}
		});
		
		Thread aQRSDet2 = new Thread(new Runnable() {
			
			@Override
			public void run() {
				for (int i = 0; i < mLength; i++) {
					aChannel2[i] = mInput[i][1];
				}
				try {
					mQRS2 = fetalQRS(aChannel2);
				} catch (Exception e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				i2Flag = true;
				
			}
		});
		Thread aQRSDet3 = new Thread(new Runnable() {
			
			@Override
			public void run() {
				for (int i = 0; i < mLength; i++) {
					aChannel3[i] = mInput[i][2];
				}
				try {
					mQRS3 = fetalQRS(aChannel3);
				} catch (Exception e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				i3Flag = true;
				
			}
		});
		Thread aQRSDet4 = new Thread(new Runnable() {
			
			@Override
			public void run() {
				for (int i = 0; i < mLength; i++) {
					aChannel4[i] = mInput[i][3];
				}
				try {
					mQRS4 = fetalQRS(aChannel4);
				} catch (Exception e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				i4Flag = true;
				
			}
		});
		aQRSDet1.start();
		aQRSDet2.start();
		aQRSDet3.start();
		aQRSDet4.start();
		
		while(true){
			Thread.sleep(1);
			if (i1Flag && i2Flag && i3Flag && i4Flag) {
				i1Flag = false;
				i2Flag = false;
				i3Flag = false;
				i4Flag = false;
				break;
			}
		}

		Object[] qrsSelectionInputs = mMatrixFunctions.channelSelection_Feb17(mQRS1, mQRS2, mQRS3, mQRS4,
				Constants.FQRS_VARIANCE_THRESHOLD, Constants.FQRS_RR_LOW_TH, Constants.FQRS_RR_HIGH_TH);

		// Object[] qrsSelectionInputs2 = Matrix.channelSelection(qrs1, qrs2,
		// qrs3, qrs4, Constants.FQRS_VARIANCE_THRESHOLD);

		QRSSelection aQrsSelect = new QRSSelection();
		Object[] aQrsSelected = aQrsSelect.qrsSelection((int[]) qrsSelectionInputs[0], (int) qrsSelectionInputs[1],
				iQrsM, iInterpolatedLength, iQRSLast, iRRMeanLast, iNoDetectionFlag);

		// qrsFSelectionQueue qrsFetal = new qrsFSelectionQueue();
		// int[] qrsF = qrsFetal.qrs((int[]) qrsSelectionInputs2[0], (int)
		// qrsSelectionInputs[1], iQrsM);

		return aQrsSelected;
	}

	/**
	 * Change the filter [a,b] values for 2 filters
	 * 
	 */

	public int[] fetalQRS(double[] channel) throws Exception{

		// differentiate and square
		mMatrixFunctions.convolutionQRSDetection(channel, Constants.QRS_DERIVATIVE);

		/**
		 * FIltering 0.8- 3Hz
		 */

		double bhigh[] = new double[2];
		for (int i0 = 0; i0 < 2; i0++) {
			bhigh[i0] = Constants.FQRS_BHIGH0 + Constants.FQRS_BHIGH_SUM * (double) i0;
		}

		mMatrixFunctions.filterLoHi(channel, Constants.FQRS_AHIGH, bhigh, Constants.FQRS_ZHIGH);

		// Have to add 6th order filter

		mMatrixFunctions.filterLoHi(channel, Constants.FQRS_ALOW, Constants.FQRS_BLOW, Constants.FQRS_ZLOW);

		/**
		 * Integrator
		 */

		double[] integrator = new double[mLength];

		double sum = 0;

		for (int j = 0; j < Constants.FQRS_WINDOW; j++) {
			sum = sum + channel[Constants.FQRS_WINDOW - j - 1];
		}
		integrator[Constants.FQRS_WINDOW - 1] = sum / Constants.FQRS_WINDOW;

		for (int i = Constants.FQRS_WINDOW; i < mLength; i++) {
			integrator[i] = integrator[i - 1]
					+ (-channel[i - Constants.FQRS_WINDOW] + channel[i]) / Constants.FQRS_WINDOW;
		}
		/**
		 * Find the 90% and 10% value to find the threshold
		 */

		double threshold = mMatrixFunctions.setIntegratorThreshold(integrator, Constants.FQRS_INTEGRATOR_THRESHOLD_SCALE);

		/**
		 * Peak Detection , not sure about return type have to change it Just
		 * return the first column of the Maxtab. No need the magnitudes.
		 */
		int peakLoc[] = mMatrixFunctions.peakDetection(integrator, threshold);

		int delay = Constants.FQRS_WINDOW / 2;
		int peakLength = peakLoc.length;
		// Check the starting peak is greater than delay/2 or remove nIt
		int count = 0;
		for (int i = 0; i < peakLength; i++) {
			if (peakLoc[i] < delay + 2) {
				count = count + 1;
			}
		}

		int lenQrs = peakLength - count;
		int[] qrs = new int[lenQrs];
		for (int i = 0; i < lenQrs; i++) {
			qrs[i] = peakLoc[i + count] - delay;
		}

		return qrs;
	}

}