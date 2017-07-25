package com.sattvamedtech.fetallite.signalproc;

public class MQRSCancelI {
	MatrixFunctions mMatrixFunctions = new MatrixFunctions();

	public double[][] cancel(double[][] iInput, int[] iQRSm) throws Exception {
		// input = Nx4

		int aNoQRSm = iQRSm.length;

		double[] aRRms = new double[aNoQRSm - 1]; // RR in milli seconds
		for (int i = 1; i < aNoQRSm; i++) {
			aRRms[i - 1] = (iQRSm[i] - iQRSm[i - 1]) / (double) Constants.FS;
		}

		double aRRmean = mMatrixFunctions.findMeanBetweenDistributionTails(aRRms, 4, 4);

		/**
		 * Initialize the no of points before and after QRS
		 */
		int aNoSamplesBeforeQRS = 0;
		int aNoSamplesAfterQRS = 0;
		if (Math.ceil(Constants.CANCEL_QRS_BEFORE_PERC * Constants.FS) == Constants.CANCEL_QRS_BEFORE_PERC
				* Constants.FS) {
			aNoSamplesBeforeQRS = (int) Math.ceil(Constants.CANCEL_QRS_BEFORE_PERC * Constants.FS);
		} else {
			aNoSamplesBeforeQRS = (int) Math.ceil(Constants.CANCEL_QRS_BEFORE_PERC * Constants.FS) - 1;
		}
		double aNoSamplesTemp = Constants.CANCEL_QRS_AFTER_PERC * (aRRmean - 0.1);

		if (aNoSamplesTemp > Constants.CANCEL_QRS_AFTER_TH) {
			aNoSamplesTemp = Constants.CANCEL_QRS_AFTER_TH;
		}
		if (Math.ceil(aNoSamplesTemp * Constants.FS) == aNoSamplesTemp * Constants.FS) {
			aNoSamplesAfterQRS = (int) Math.ceil(aNoSamplesTemp * Constants.FS);
		} else {
			aNoSamplesAfterQRS = (int) Math.ceil(aNoSamplesTemp * Constants.FS) - 1;
		}

		int aNoSamplesQRS = 1 + aNoSamplesBeforeQRS + aNoSamplesAfterQRS;

		/**
		 * Extend signals to manage first QRS
		 */
		int aInitialQrsIndexArr = 1;

		/**
		 * Extend Signals to manage first QRS
		 */
		// Always the first qrsM > 120, so take only else condition
		int aNoSamplesToLeft = 0;
		if (aNoSamplesBeforeQRS + 1 - iQRSm[0] > 0) {
			aNoSamplesToLeft = aNoSamplesBeforeQRS + 1 - iQRSm[0];
		}

		double[][] aRowExtract = mMatrixFunctions.subMatrix(iInput, 0, 0, 0, Constants.NO_OF_CHANNELS - 1);

		double[][] aRowExtension = mMatrixFunctions.repmat(aRowExtract, aNoSamplesToLeft);
		double[][] aInputExtension = mMatrixFunctions.verticalConcat(aRowExtension, iInput);

		/**
		 * Extend signals to manage last QRS
		 */
		// Always the last qrsM < len - 140, so take only else condition

		int aFinalQrsIndexArr = aNoQRSm;
		double[] aRRMedTemp = new double[Constants.CANCEL_NO_SAMPLES_END];
		double aRRMeanSum = 0;
		for (int i = 0; i < Constants.CANCEL_NO_SAMPLES_END; i++) {
			aRRMedTemp[i] = aRRms[aRRms.length - 1 - i];
			aRRMeanSum = aRRMedTemp[i] + aRRMeanSum;
		}
		double aRRmsMean = aRRMeanSum / Constants.CANCEL_NO_SAMPLES_END;
		double aRRmsMedian = mMatrixFunctions.findMedian(aRRMedTemp);

		double aTempD = (1 - Constants.CANCEL_LASTQRS_TH_HIGH_PERC) * Constants.FS * aRRmsMedian;
		int aTempI = 0;
		if (Math.ceil(aTempD) == aTempD) {
			aTempI = (int) Math.ceil(aTempD);
		} else {
			aTempI = (int) Math.ceil(aTempD) - 1;
		}
		double aNoSamplesAddedEndTemp = 0;
		int aNoSamplesAddedEnd = 0;

		if (iQRSm[aFinalQrsIndexArr - 1] + aTempI < Constants.NO_OF_SAMPLES) {
			// find max
			if (Constants.CANCEL_LASTQRS_TH_LOW_PERC * Constants.FS > Constants.CANCEL_LASTQRS_TH_HIGH_PERC
					* Constants.FS * aRRmsMean) {
				aNoSamplesAddedEndTemp = Constants.CANCEL_LASTQRS_TH_LOW_PERC * Constants.FS;
			} else {
				aNoSamplesAddedEndTemp = Constants.CANCEL_LASTQRS_TH_HIGH_PERC * Constants.FS * aRRmsMean;
			}

			if (Math.ceil(aNoSamplesAddedEndTemp) == aNoSamplesAddedEndTemp) {
				aNoSamplesAddedEnd = (int) Math.ceil(aNoSamplesAddedEndTemp);
			} else {
				aNoSamplesAddedEnd = (int) Math.ceil(aNoSamplesAddedEndTemp) - 1;
			}

			int aRowToExtend = aInputExtension.length - 1 - aNoSamplesAddedEnd - 1;

			// Do replicate the row and add it to the input extension
			double[][] aRowExtended = mMatrixFunctions.subMatrix(aInputExtension, aRowToExtend, aRowToExtend, 0,
					aInputExtension[0].length - 1);

			double[][] aInputExtendFinal = mMatrixFunctions.repmat(aRowExtended, aNoSamplesAddedEnd);
			for (int i = 0; i < aNoSamplesAddedEnd; i++) {
				for (int j = 0; j < aInputExtension[0].length; j++) {
					aInputExtension[i + aRowToExtend + 2][j] = aInputExtendFinal[i][j];
				}
			}

		} // end if qrsm[qf]
		/**
		 * no of samples to add to right of the signal
		 */
		int aNoSamplesToRight = 0;

		if (iQRSm[aFinalQrsIndexArr - 1] + aNoSamplesAfterQRS - Constants.NO_OF_SAMPLES > -1) {
			aNoSamplesToRight = iQRSm[aFinalQrsIndexArr - 1] + aNoSamplesAfterQRS - Constants.NO_OF_SAMPLES + 1;
		}
		double[][] aInputSVD = new double[aInputExtension.length + aNoSamplesToRight][aInputExtension[0].length];
		for (int i = 0; i < aInputExtension.length; i++) {
			for (int j = 0; j < aInputExtension[0].length; j++) {
				aInputSVD[i][j] = aInputExtension[i][j];
			}
		}
		// Do extension if required.
		if (aNoSamplesToRight > 0) {
			double[][] aRowExtractRight = mMatrixFunctions.subMatrix(iInput, Constants.NO_OF_SAMPLES - 1,
					Constants.NO_OF_SAMPLES - 1, 0, Constants.NO_OF_CHANNELS - 1);
			double[][] aReplicateSamplesRight = mMatrixFunctions.repmat(aRowExtractRight, aNoSamplesToRight);

			for (int i = 0; i < aNoSamplesToRight; i++) {
				for (int j = 0; j < Constants.NO_OF_CHANNELS; j++) {
					aInputSVD[i + aInputExtension.length][j] = aReplicateSamplesRight[i][j];
				}
			}
		}

		/**
		 * Added Samples to right :: ALL extensions of signal is done.
		 */

		int aNoSamplesExtend = aInputSVD.length;
		int aNoQrs = aFinalQrsIndexArr - aInitialQrsIndexArr + 1;

		int[] aInitQrsLocations = new int[aNoQrs];
		for (int i = aInitialQrsIndexArr; i <= aFinalQrsIndexArr; i++) {
			aInitQrsLocations[i - aInitialQrsIndexArr] = iQRSm[i - 1];
		}

		/**
		 * Start and end of QRS window
		 */
		int[] aQrsIndexArrI = new int[aNoQrs];
		int[] aQrsIndexArrF = new int[aNoQrs];

		for (int i = 0; i < aNoQrs; i++) {
			aQrsIndexArrI[i] = aInitQrsLocations[i] + aNoSamplesToLeft - aNoSamplesBeforeQRS;
			aQrsIndexArrF[i] = aInitQrsLocations[i] + aNoSamplesToLeft + aNoSamplesAfterQRS;
		}

		double aSvdExtract[][] = new double[aNoSamplesQRS][aNoQrs];

		// add weight function
		double[][] aWeightWindowFunction = mMatrixFunctions.weightFunction(aNoSamplesBeforeQRS, aNoSamplesAfterQRS,
				Constants.FS);
		// double[][] wwg = Matrix.transpose(wwgT);
		double[][] aApproxmSignal = new double[aNoSamplesExtend][Constants.NO_OF_CHANNELS];

		/**
		 * Start loop for doing SVD and substraction
		 */
		double[][] aRowextract = new double[1][aNoSamplesQRS];
		double[][] aRowWeighted = new double[1][aNoSamplesQRS];
		for (int is = 0; is < Constants.NO_OF_CHANNELS; is++) {

			
			for (int iq = 0; iq < aNoQrs; iq++) {
				int iQrsIndex = aQrsIndexArrI[iq];
				int fQrsIndex = aQrsIndexArrF[iq];

				aRowextract = mMatrixFunctions.subMatrix(aInputSVD, iQrsIndex, fQrsIndex, is, is);
				aRowWeighted = mMatrixFunctions.elementWiseMultiply(aRowextract, aWeightWindowFunction);
				for (int j = 0; j < aNoSamplesQRS; j++) {
					aSvdExtract[j][iq] = aRowWeighted[j][0];
				}
			} // end extracting the matrix A

			
			ImplicitSVDusingEVD aISVD = new ImplicitSVDusingEVD();

			double[][] aApproxSignal = aISVD.impcitsvd(aSvdExtract);
			/**
			 * We did svd for A^T. A = ui * vi^T * sigma.
			 */

			// putting back the approximation into a single channel
			double[][] aApproxSignalTemp = new double[aNoSamplesQRS][1];
			double[][] aApproxSignalTemp1 = new double[aNoSamplesQRS][1];
			for (int iq = 0; iq < aNoQrs; iq++) {
				int aIwq = aQrsIndexArrI[iq];
				int aFwq = aQrsIndexArrF[iq];
				aApproxSignalTemp = mMatrixFunctions.subMatrix(aApproxSignal, 0, aNoSamplesQRS - 1, iq, iq);
				aApproxSignalTemp1 = mMatrixFunctions.elementWiseDivide(aApproxSignalTemp, aWeightWindowFunction);

				for (int i = aIwq; i <= aFwq; i++) {
					aApproxmSignal[i][is] = aApproxSignalTemp1[i - aIwq][0];
				}
			} // end approx single channel

			// smoothening connectioons btw sucessive channels
			double aDifferenceValue = 0;
			double aPercentValue = 0;
			for (int iq = 1; iq < aNoQrs; iq++) {
				int aQrsIndexF = aQrsIndexArrF[iq - 1];
				int aQrsIndexI = aQrsIndexArrI[iq];
				if (aQrsIndexI > aQrsIndexF) {
					aDifferenceValue = aApproxmSignal[aQrsIndexI][is] - aApproxmSignal[aQrsIndexF][is];
					aPercentValue = aDifferenceValue / (aQrsIndexI - aQrsIndexF);
					for (int it = aQrsIndexF + 1; it < aQrsIndexI; it++) {
						aApproxmSignal[it][is] = aApproxmSignal[aQrsIndexF][is] + aPercentValue * (it - aQrsIndexF);
					}
				}
			} // end smoothening

		} // end SVD substraction loop

		double[][] aResidueOutput = new double[Constants.NO_OF_SAMPLES][Constants.NO_OF_CHANNELS];
		for (int i = 0; i < Constants.NO_OF_SAMPLES; i++) {
			for (int j = 0; j < Constants.NO_OF_CHANNELS; j++) {
				aResidueOutput[i][j] = aInputSVD[i + aNoSamplesToLeft][j] - aApproxmSignal[i + aNoSamplesToLeft][j];
			}
		}

		return aResidueOutput;
	} // end main

} // end class