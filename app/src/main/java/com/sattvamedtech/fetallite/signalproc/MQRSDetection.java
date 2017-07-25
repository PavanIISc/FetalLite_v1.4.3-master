package com.sattvamedtech.fetallite.signalproc;

public class MQRSDetection {
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
	public int[] mQRS(double[][] iInput) throws Exception {

		mLength = iInput.length;
		mInput = iInput;

		aChannel1 = new double[mLength];
		aChannel2 = new double[mLength];
		aChannel3 = new double[mLength];
		aChannel4 = new double[mLength];

//		for (int i = 0; i < length; i++) {
//			aChannel1[i] = iInput[i][0];
//			aChannel2[i] = iInput[i][1];
//			aChannel3[i] = iInput[i][2];
//			aChannel4[i] = iInput[i][3];
//		}

//		int aQRS1[] = maternalQRS(aChannel1);
//		int aQRS2[] = maternalQRS(aChannel2);
//		int aQRS3[] = maternalQRS(aChannel3);
//		int aQRS4[] = maternalQRS(aChannel4);
//		
		Thread aQRSDet1 = new Thread(new Runnable() {
			
			@Override
			public void run() {
				for (int i = 0; i < mLength; i++) {
					aChannel1[i] = mInput[i][0];
				}
				try {
					mQRS1 = maternalQRS(aChannel1);
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
					mQRS2 = maternalQRS(aChannel2);
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
					mQRS3 = maternalQRS(aChannel3);
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
					mQRS4 = maternalQRS(aChannel4);
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
				Constants.MQRS_VARIANCE_THRESHOLD, Constants.MQRS_RR_LOW_TH, Constants.MQRS_RR_HIGH_TH);

		// qrsmSelection qrsMaternal = new qrsmSelection();
		// int[] qrsM = qrsMaternal.qrsSelection((int[])qrsSelectionInputs[0],
		// (int)qrsSelectionInputs[1]);

		QRSSelection aQrsSelect = new QRSSelection();
		Object[] aQrsSelected = aQrsSelect.qrsSelection((int[]) qrsSelectionInputs[0], (int) qrsSelectionInputs[1],
				new int[] {}, 0, 0, 0, -1);

		return (int[]) aQrsSelected[0];
	}

	public int[] maternalQRS(double[] iChannel) throws Exception{
		int aLength = iChannel.length;
		// differentiate and square
		mMatrixFunctions.convolutionQRSDetection(iChannel, Constants.QRS_DERIVATIVE);

		/**
		 * FIltering 0.8- 3Hz
		 */

		double aBhigh[] = new double[2];
		double aAhigh[] = new double[2];
		for (int i0 = 0; i0 < 2; i0++) {
			aBhigh[i0] = Constants.MQRS_BHIGH0 + Constants.MQRS_BHIGH_SUM * (double) i0;
			aAhigh[i0] = 1.0 + Constants.MQRS_AHIGH_SUM * (double) i0;
		}

		mMatrixFunctions.filterLoHi(iChannel, aAhigh, aBhigh, Constants.MQRS_ZHIGH);

		double aAlow[] = new double[2];
		for (int i0 = 0; i0 < 2; i0++) {
			aAlow[i0] = 1.0 + Constants.MQRS_ALOW_SUM * (double) i0;
		}

		mMatrixFunctions.filterLoHi(iChannel, aAlow, Constants.MQRS_BLOW, Constants.MQRS_ZLOW);

		/**
		 * Integrator
		 */

		double[] aIntegrator = new double[aLength];

		double aSum = 0;

		for (int j = 0; j < Constants.MQRS_WINDOW; j++) {
			aSum = aSum + iChannel[Constants.MQRS_WINDOW - j - 1];
		}
		aIntegrator[Constants.MQRS_WINDOW - 1] = aSum / Constants.MQRS_WINDOW;

		for (int i = Constants.MQRS_WINDOW; i < aLength; i++) {
			aIntegrator[i] = aIntegrator[i - 1]
					+ (-iChannel[i - Constants.MQRS_WINDOW] + iChannel[i]) / Constants.MQRS_WINDOW;
		}

		/**
		 * Find the 90% and 10% value to find the threshold
		 */

		double aThreshold = mMatrixFunctions.setIntegratorThreshold(aIntegrator,
				Constants.MQRS_INTEGRATOR_THRESHOLD_SCALE);

		/**
		 * Peak Detection , Determines the peaks of the signal Just return the
		 * first column of the Maxtab. No need the magnitudes.
		 */
		int aPeakLoc[] = mMatrixFunctions.peakDetection(aIntegrator, aThreshold);

		int aDelay = Constants.MQRS_WINDOW / 2;
		int aPeakLength = aPeakLoc.length;
		// Check the starting peak is greater than delay/2 or remove nIt
		int aCount = 0;
		for (int i = 0; i < aPeakLength; i++) {
			if (aPeakLoc[i] < aDelay + 2) {
				aCount = aCount + 1;
			}
		}

		int aLenQrs = aPeakLength - aCount;
		int[] aQrs = new int[aLenQrs];
		for (int i = 0; i < aLenQrs; i++) {
			aQrs[i] = aPeakLoc[i + aCount] - aDelay;
		}

		return aQrs;
	}

}