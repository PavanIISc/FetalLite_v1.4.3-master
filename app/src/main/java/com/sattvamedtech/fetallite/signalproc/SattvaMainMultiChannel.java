package com.sattvamedtech.fetallite.signalproc;

import java.io.BufferedWriter;
import java.io.FileWriter;


public class SattvaMainMultiChannel 
{
	public static void main(String[]args)
	{
		System.out.println("***Started At : "+(new java.text.SimpleDateFormat("H:mm:ss:SSS")).format(java.util.Calendar.getInstance().getTime()));

		TextFileReader readTxt = new TextFileReader();
		double[][] input = readTxt.readFile();
//		double[][] inputAndroid = readTxt.readFileNumber();

		AlgorithmMain Algo = new AlgorithmMain();

		
		
		for (int it = 0; it <28; it++)
		{
			double[][] input1 = new double[15000][4];
	
				for (int i = 0; i<15000; i++)
				{
					for (int j = 0; j<4; j++)
					{
						input1[i][j] = input[i+Constants.QRS_SHIFT*it][j];
					}
				}
				
				Object[] Final;
				try {
					long T1 = System.currentTimeMillis();
					
					Final = Algo.algoStart(input1, it);
					
					long T2 = System.currentTimeMillis();
					System.out.println("Time for Algo to complete: " + (T2 - T1) + " ms");					
					
					
				} catch (Exception e) {
					e.printStackTrace();
				}
				System.out.println("***Completed At : "+(new java.text.SimpleDateFormat("H:mm:ss:SSS")).format(java.util.Calendar.getInstance().getTime()));
				
				System.gc();
		}
		
		int lengthQRSF = Constants.QRS_FETAL_LOCATION.size();
		int lengthQRSM = Constants.QRS_MATERNAL_LOCATION.size();
		double[][] QRSF = new double[lengthQRSF][2];
		for (int i =0; i<lengthQRSF; i++)
		{
			QRSF[i][0] = Constants.QRS_FETAL_LOCATION.get(i);
			QRSF[i][1] = Constants.HR_FETAL.get(i);
			
		}
		double[][] QRSM = new double[lengthQRSM][2];
		for (int i =0; i<lengthQRSM; i++)
		{
			QRSM[i][0] = Constants.QRS_MATERNAL_LOCATION.get(i);
			QRSM[i][1] = Constants.HR_MATERNAL.get(i);
		}

		
		
		
		try {
			BufferedWriter br = new BufferedWriter(new FileWriter("qrs.csv"));
			StringBuilder sb = new StringBuilder();
			for (int i = 0; i < lengthQRSF; i++) {
			 sb.append(QRSF[i][0]);
			 sb.append(",");
			 sb.append(QRSF[i][1]);
			 sb.append("\n");
			}
//			sb.append("\n");
			

			br.write(sb.toString());
			br.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		
//		try {
//			BufferedWriter br = new BufferedWriter(new FileWriter("qrsfRandSVD.csv"));
//			StringBuilder sb = new StringBuilder();
//			for (int i = 0; i < lengthQRSF1; i++) {
//			 sb.append(qrsF1.get(i));
////			 sb.append(",");
////			 sb.append(QRSM[i][1]);
//			 sb.append("\n");
//			}
////			sb.append("\n");
//			
//
//			br.write(sb.toString());
//			br.close();
//		} catch (Exception e) {
//			e.printStackTrace();
//		}
		
		System.gc();
	} // end main

} // end class