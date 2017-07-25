package com.sattvamedtech.fetallite.signalproc;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;

public class TextFileReader {
	private static double mCheck = Math.pow(2, 23);
	private static double mVref = 4.5;
	private static double mGain = 24;
	private static double mCheckDivide = 2 * mCheck;
	private Queue<String> mStringArray = new LinkedList();
	private double mInput[][];
	private List<String> mStringCheck;
	public double[][] readFileNumber(){
		try {
			BufferedReader br = new BufferedReader(new FileReader("android_in.txt"));
//			StringBuilder sb = new StringBuilder();
			String line;
			ArrayList<String> aStrings = new ArrayList<>();

//			int counter = 0;
			while ((line = br.readLine()) != null) {
				aStrings.add(line);
			}
			
//			double[][] aDoubleValues = new double[aStrings.size()][4];
			double[][] aDoubleValues = new double[15000][4];
			
//			for (int i = 0; i < aStrings.size(); i++) {
			for (int i = 1; i<=15000; i++){
				String[] aSplitLine = aStrings.get(i).split(",");
				for (int j=0; j< 4; j++) {
					aDoubleValues[i][j] = Double.parseDouble(aSplitLine[j]);
				}
			}

			br.close();
			return aDoubleValues;
		} catch (Exception e) {
			e.printStackTrace();
		}
		return new double[0][0];
	}

	public double[][] readFile() {
		

		try {
			BufferedReader br = new BufferedReader(new FileReader("/Users/kishoresubramanian/Downloads/android_in.txt"));
			StringBuilder sb = new StringBuilder();
			String line;

			int counter = 0;
			while ((line = br.readLine()) != null) {
				String[] aSplitLine = line.split("\\+");
				mStringArray.addAll(Arrays.asList(aSplitLine).subList(1, aSplitLine.length));
				mStringCheck = Arrays.asList(aSplitLine).subList(1, aSplitLine.length);
				if (mStringCheck == null)
				{
					System.out.println("count = " + counter);
					break;
				}
				counter++;
			}
			System.out.println("count = " + counter);

			br.close();
		} catch (Exception e) {
			e.printStackTrace();
		}

		mInput = new double[mStringArray.size()][4];
//		for (int i = 0; i < mStringArray.size(); i++) {
//			for (int aChannleCount = 0; aChannleCount < 4; aChannleCount++) {
//				mInput[i][aChannleCount] = stringToDouble(
//						mStringArray.get(i).substring(6 * aChannleCount + 1, 6 * aChannleCount + 6));
//			}
//		}
		
		populateInputArray();
		/**
		 * To write into a file
		 */
//		try {
//			BufferedWriter br = new BufferedWriter(new FileWriter("uc_file_sample_Java.csv"));
//			StringBuilder sb = new StringBuilder();
//			for (int i = 0; i < mStringArray.size(); i++) {
//			for (int j = 0; j < 4; j++) {
//			 sb.append(mInput[i][j]);
//			 sb.append(",");
//			}
//			sb.append("\n");
//			}
//
//			br.write(sb.toString());
//			br.close();
//		} catch (Exception e) {
//			e.printStackTrace();
//		}

		return mInput;

	}
	
	private void feedInputArray(String iInputString, int iInputIndex) {
	    for (int aInputChannelIndex = 0; aInputChannelIndex < 4; aInputChannelIndex++) {
	        mInput[iInputIndex][aInputChannelIndex] = stringToDouble(iInputString.substring(6 * aInputChannelIndex + 1, 6 * aInputChannelIndex + 7));
//	        if (aInputChannelIndex == 0) {
//	            ApplicationUtils.mInputArrayUc[iInputIndex] = ApplicationUtils.mInputArray[iInputIndex][aInputChannelIndex];
//	        }
	    }
	}

	private void populateInputArray() {
		int aInputArrayCounter = 0;
		String aSample = getNextValidSample();

		int aLastIndex = Character.getNumericValue(aSample.charAt(0));
		feedInputArray(aSample, aInputArrayCounter);
		for (aInputArrayCounter++; aInputArrayCounter < mStringArray.size(); aInputArrayCounter++) {
			aSample = getNextValidSample();
			int aCurrentIndex = Character.getNumericValue(aSample.charAt(0));
			int aIndexDiff = aCurrentIndex - aLastIndex;
			if (aIndexDiff <= 0)
				aIndexDiff += 10;
			// DO CHANGE IN THIS LINE
			if (aIndexDiff == 1 || aIndexDiff == -9) {
				feedInputArray(aSample, aInputArrayCounter);
			} else {
				aInputArrayCounter += aIndexDiff-1;
				feedInputArray(aSample, aInputArrayCounter);
				interpolate(aInputArrayCounter, aInputArrayCounter - aIndexDiff, aLastIndex);
			}
			aLastIndex = Character.getNumericValue(aSample.charAt(0));
		}

	}

	private String getNextValidSample() {
		String aValidSample = "";
		if (mStringArray.size() > 0) {
			do { 
				aValidSample = mStringArray.remove();
			} while (aValidSample.length() != 25 && mStringArray.size() > 0);
		}
		return aValidSample;
	}

	private void interpolate(int iCurrentInputIndex, int iStartIndex, int iEndIndex) {
		System.out.println("iCurrentInputIndex: " + iCurrentInputIndex);
		System.out.println("iStartIndex: " + iStartIndex);
		System.out.println("iEndIndex: " + iEndIndex);
		for (int k = iStartIndex + 1; k < iCurrentInputIndex; k++) {
			for (int aInputChannelIndex = 0; aInputChannelIndex < 4; aInputChannelIndex++) {
				// DO CHANGE IN THIS LINE
				mInput[k][aInputChannelIndex] = mInput[iStartIndex][aInputChannelIndex] + (mInput[iCurrentInputIndex][aInputChannelIndex] - mInput[iStartIndex][aInputChannelIndex]) / (iCurrentInputIndex - iStartIndex) * (k-iStartIndex);
			}
		}
	}

	private double stringToDouble(String iChannelInput) {
		return doubleConv(new BigInteger(iChannelInput, 16).doubleValue());
	}

	private double doubleConv(double iDoubleValue) {
		double aOut;
		if (iDoubleValue >= mCheck) {
			aOut = (iDoubleValue - mCheckDivide) * mVref / (mCheck - 1) / mGain;
		} else {
			aOut = iDoubleValue / (mCheck - 1) / mGain * mVref;
		}
		return aOut;
	}
}

