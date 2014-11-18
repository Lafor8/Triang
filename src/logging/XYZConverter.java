package logging;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.text.DecimalFormat;

import org.opencv.core.Mat;
import org.opencv.highgui.Highgui;

import android.os.Environment;
import android.util.Log;

public class XYZConverter {
	private static String ROOT = Environment.getExternalStorageDirectory().toString();
	private static File dir;
	private static File xyzFile;

	private static String dirPath = "/_thesis/";
	private static String fileName = "try.xyz";

	public XYZConverter() {
	}
	
	public static void writeAllToXYZFile(Mat mat, Mat img1, Mat img2, Mat points4D, int divr,int divc, Mat rgb, String fileName) {
		StringBuilder sb = new StringBuilder();
		DecimalFormat df = new DecimalFormat("0.000");

		sb.append("0 0 0 255 0 0");
		sb.append('\n');
		sb.append("1 0 0 130 0 0");
		sb.append('\n');

		sb.append("0 0 1 0 130 0");
		sb.append('\n');
		sb.append("0 0 -1 0 130 0");
		sb.append('\n');

		sb.append("0 1 0 0 0 130");
		sb.append('\n');
		sb.append("0 -1 0 0 0 130");
		sb.append('\n');
		
		double max = Double.MIN_VALUE;
		double min = Double.MAX_VALUE;
		double off = 0;
		double num;
		for (int i = 0; i < mat.rows(); ++i) {
			for (int j = 0; j < mat.cols(); ++j) {
				num = mat.get(i,j)[0];
				if(num > max)
					max = num;
				if(num < min)
					min = num;
			}
		}
		
		if(min < 0){
			off = min;
		}
		else{
			off = -min;
		}

		min += off;
		max += off;
		
		for (int i = 0; i < mat.rows(); ++i) {
			for (int j = 0; j < mat.cols(); ++j) {
				num = (int)(((mat.get(i, j)[0] + off)/max )*255);
				
				// Previous Image
				
//				sb.append(i-mat.rows());
//				sb.append(' ');
//				sb.append('0');
//				sb.append(' ');
//				sb.append(j+mat.cols());
//				sb.append(' ');
//				sb.append((int)img1.get(i, j)[0]);
//				sb.append(' ');
//				sb.append((int)img1.get(i, j)[0]);
//				sb.append(' ');
//				sb.append((int)img1.get(i, j)[0]);
//				sb.append('\n');
//				
//				// Current Image
//				
//				sb.append(i-mat.rows());
//				sb.append(' ');
//				sb.append('0');
//				sb.append(' ');
//				sb.append(j-mat.cols());
//				sb.append(' ');
//				sb.append((int)img2.get(i, j)[0]);
//				sb.append(' ');
//				sb.append((int)img2.get(i, j)[0]);
//				sb.append(' ');
//				sb.append((int)img2.get(i, j)[0]);
//				sb.append('\n');
//				
//				// Optical Flow Map
//				sb.append(i-mat.rows());
//				sb.append(' ');
//				sb.append('0');
//				sb.append(' ');
//				sb.append(j);
//				sb.append(' ');
//				sb.append((int)num);
//				sb.append(" 0 0");
//				sb.append('\n');
//				
				
			}
		}
		
		mul = 10;
		int tra = 0;
		double x,y,z;
		//int skip=0;
		for (int i = 0, k = 0; i < rgb.rows()/divr; ++i) {
			for (int j = 0; j < rgb.cols()/divc; ++j, ++k) {
				
				x =points4D.get(0, k)[0] / points4D.get(3, k)[0];
				y =points4D.get(2, k)[0] / points4D.get(3, k)[0];
				z =points4D.get(1, k)[0] / points4D.get(3, k)[0];
//				if(x +y+z > 500 ){
//					skip++;
//					continue;
//				}
				sb.append(df.format(x*mul+tra));
				sb.append(' ');
				sb.append(df.format(y*mul+tra));
				sb.append(' ');
				sb.append(df.format(z*mul+tra));
				sb.append(' ');
				sb.append((int) (rgb.get(i, j)[0]));
				sb.append(' ');
				sb.append((int) (rgb.get(i, j)[1]));
				sb.append(' ');
				sb.append((int) (rgb.get(i, j)[2]));
				sb.append('\n');
			}
		}
	//Log.i("Triangule", skip + " points were left out.");

		writeToFile(sb.toString(),fileName);
	}
	
	public static void writeAllToXYZFile(Mat mat, Mat img1, Mat img2, Mat points4D, int div, Mat rgb) {
		StringBuilder sb = new StringBuilder();
		DecimalFormat df = new DecimalFormat("0.000");

		sb.append("0 0 0 255 0 0");
		sb.append('\n');
		sb.append("1 0 0 130 0 0");
		sb.append('\n');

		sb.append("0 0 1 0 130 0");
		sb.append('\n');
		sb.append("0 0 -1 0 130 0");
		sb.append('\n');

		sb.append("0 1 0 0 0 130");
		sb.append('\n');
		sb.append("0 -1 0 0 0 130");
		sb.append('\n');
		
		double max = Double.MIN_VALUE;
		double min = Double.MAX_VALUE;
		double off = 0;
		double num;
		for (int i = 0; i < mat.rows(); ++i) {
			for (int j = 0; j < mat.cols(); ++j) {
				num = mat.get(i,j)[0];
				if(num > max)
					max = num;
				if(num < min)
					min = num;
			}
		}
		
		if(min < 0){
			off = min;
		}
		else{
			off = -min;
		}

		min += off;
		max += off;
		
		for (int i = 0; i < mat.rows(); ++i) {
			for (int j = 0; j < mat.cols(); ++j) {
				num = (int)(((mat.get(i, j)[0] + off)/max )*255);
				
				// Previous Image
				
//				sb.append(i-mat.rows());
//				sb.append(' ');
//				sb.append('0');
//				sb.append(' ');
//				sb.append(j+mat.cols());
//				sb.append(' ');
//				sb.append((int)img1.get(i, j)[0]);
//				sb.append(' ');
//				sb.append((int)img1.get(i, j)[0]);
//				sb.append(' ');
//				sb.append((int)img1.get(i, j)[0]);
//				sb.append('\n');
//				
//				// Current Image
//				
//				sb.append(i-mat.rows());
//				sb.append(' ');
//				sb.append('0');
//				sb.append(' ');
//				sb.append(j-mat.cols());
//				sb.append(' ');
//				sb.append((int)img2.get(i, j)[0]);
//				sb.append(' ');
//				sb.append((int)img2.get(i, j)[0]);
//				sb.append(' ');
//				sb.append((int)img2.get(i, j)[0]);
//				sb.append('\n');
//				
//				// Optical Flow Map
//				sb.append(i-mat.rows());
//				sb.append(' ');
//				sb.append('0');
//				sb.append(' ');
//				sb.append(j);
//				sb.append(' ');
//				sb.append((int)num);
//				sb.append(" 0 0");
//				sb.append('\n');
//				
				
			}
		}
		
		mul = 1;
		int tra = 0;
		double x,y,z;
		//int skip=0;
		for (int i = 0, k = 0; i < rgb.rows()/div; ++i) {
			for (int j = 0; j < rgb.cols()/div; ++j, ++k) {
				
				x =points4D.get(0, k)[0] / points4D.get(3, k)[0];
				y =points4D.get(2, k)[0] / points4D.get(3, k)[0];
				z =points4D.get(1, k)[0] / points4D.get(3, k)[0];
//				if(x +y+z > 500 ){
//					skip++;
//					continue;
//				}
				sb.append(df.format(x*mul+tra));
				sb.append(' ');
				sb.append(df.format(y*mul+tra));
				sb.append(' ');
				sb.append(df.format(z*mul+tra));
				sb.append(' ');
				sb.append((int) (rgb.get(i, j)[0]));
				sb.append(' ');
				sb.append((int) (rgb.get(i, j)[1]));
				sb.append(' ');
				sb.append((int) (rgb.get(i, j)[2]));
				sb.append('\n');
			}
		}
	//Log.i("Triangule", skip + " points were left out.");

		writeToFile(sb.toString());
	}


	public static void writePointsToXYZFile(Mat points4D) {
		StringBuilder sb = new StringBuilder();
		DecimalFormat df = new DecimalFormat("############0000.0000############");

		sb.append("0 0 0 255 0 0");
		sb.append('\n');
		sb.append("1 0 0 130 0 0");
		sb.append('\n');

		sb.append("0 0 1 0 130 0");
		sb.append('\n');
		sb.append("0 0 -1 0 130 0");
		sb.append('\n');

		sb.append("0 1 0 0 0 130");
		sb.append('\n');
		sb.append("0 -1 0 0 0 130");
		sb.append('\n');

		for (int i = 0; i < points4D.cols(); ++i) {
			sb.append(df.format(-points4D.get(2, i)[0] / points4D.get(3, i)[0]));
			sb.append(' ');
			sb.append(df.format(-points4D.get(0, i)[0] / points4D.get(3, i)[0]));
			sb.append(' ');
			sb.append(df.format(-points4D.get(1, i)[0] / points4D.get(3, i)[0]));
			sb.append(' ');
			sb.append("100 100 100");
			sb.append('\n');
		}

		writeToFile(sb.toString());

	}
	static int mul;
	public static void writeRGBPointsToXYZFile(Mat points4D, Mat mat, int div) {
		StringBuilder sb = new StringBuilder();
		DecimalFormat df = new DecimalFormat("############0.000############");

		
		mul = 20;
		sb.append("0 0 0 255 0 0");
		sb.append('\n');
		sb.append("1 0 0 130 0 0");
		sb.append('\n');

		sb.append("0 0 1 0 130 0");
		sb.append('\n');
		sb.append("0 0 -1 0 130 0");
		sb.append('\n');

		sb.append("0 1 0 0 0 130");
		sb.append('\n');
		sb.append("0 -1 0 0 0 130");
		sb.append('\n');

		for (int i = 0, k = 0; i < mat.rows()/div; ++i) {
			for (int j = 0; j < mat.cols()/div; ++j, ++k) {

				sb.append(df.format((points4D.get(0, k)[0] / points4D.get(3, k)[0])*mul));
				sb.append(' ');
				sb.append(df.format((points4D.get(2, k)[0] / points4D.get(3, k)[0])*mul));
				sb.append(' ');
				sb.append(df.format((points4D.get(1, k)[0] / points4D.get(3, k)[0])*mul));
				sb.append(' ');
				sb.append((int) (mat.get(i, j)[0]));
				sb.append(' ');
				sb.append((int) (mat.get(i, j)[1]));
				sb.append(' ');
				sb.append((int) (mat.get(i, j)[2]));
				sb.append('\n');
			}
		}

		writeToFile(sb.toString());

	}

	// Given Points with 32 decimal places and 32 real no. places, the String
	// can hold N lines
	// N = ((2^31)-1)/(3*(32+1+32)+3*(3)+5) = 1,027,5041 lines/points only
	// Image resolution now limited to either 1280 �~ 720 (16:9) or 1152 x 86
	// (4:3)
	private static void writeToFile(String s) {

		dir = new File(ROOT + dirPath);
		xyzFile = new File(ROOT + dirPath + fileName);

		FileOutputStream outputStream = null;

		try {
			if (!dir.exists())
				if (dir.mkdirs())
					throw new IOException();

			if (!xyzFile.exists()) {
				xyzFile.createNewFile();
				outputStream = new FileOutputStream(xyzFile);
			} else
				outputStream = new FileOutputStream(xyzFile, false);

			outputStream.write(s.getBytes());
		} catch (Exception e) {

		} finally {
			if (outputStream != null)
				try {
					outputStream.close();
				} catch (IOException e) {
					Log.e("XYZConverter", "Problem closing stream", e);
				}
		}
	}
	
	
	private static void writeToFile(String s, String cFileName) {

		dir = new File(ROOT + dirPath);
		xyzFile = new File(ROOT + dirPath + cFileName);

		FileOutputStream outputStream = null;

		try {
			if (!dir.exists())
				if (dir.mkdirs())
					throw new IOException();

			if (!xyzFile.exists()) {
				xyzFile.createNewFile();
				outputStream = new FileOutputStream(xyzFile);
			} else
				outputStream = new FileOutputStream(xyzFile, false);

			outputStream.write(s.getBytes());
		} catch (Exception e) {

		} finally {
			if (outputStream != null)
				try {
					outputStream.close();
				} catch (IOException e) {
					Log.e("XYZConverter", "Problem closing stream", e);
				}
		}
	}

	public static void writeGrayMatToXYZFile(Mat mat) {
		StringBuilder sb = new StringBuilder();
		DecimalFormat df = new DecimalFormat("############0000.0000############");

		sb.append("0 0 0 255 0 0");
		sb.append('\n');
		sb.append("1 0 0 130 0 0");
		sb.append('\n');

		sb.append("0 0 1 0 130 0");
		sb.append('\n');
		sb.append("0 0 -1 0 130 0");
		sb.append('\n');

		sb.append("0 1 0 0 0 130");
		sb.append('\n');
		sb.append("0 -1 0 0 0 130");
		sb.append('\n');

		for (int i = 0; i < mat.rows(); ++i) {
			for (int j = 0; j < mat.cols(); ++j) {
				sb.append(i);
				sb.append(' ');
				sb.append('0');
				sb.append(' ');
				sb.append(j);
				sb.append(' ');
				sb.append((int) (mat.get(i, j)[0]));
				sb.append(' ');
				sb.append((int) (mat.get(i, j)[0]));
				sb.append(' ');
				sb.append((int) (mat.get(i, j)[0]));
				sb.append('\n');
			}
		}

		writeToFile(sb.toString());
	}
	
	public static void writeOFMatToXYZFile(Mat mat) {
		StringBuilder sb = new StringBuilder();
		DecimalFormat df = new DecimalFormat("############0.000############");

		sb.append("0 0 0 255 0 0");
		sb.append('\n');
		sb.append("1 0 0 130 0 0");
		sb.append('\n');

		sb.append("0 0 1 0 130 0");
		sb.append('\n');
		sb.append("0 0 -1 0 130 0");
		sb.append('\n');

		sb.append("0 1 0 0 0 130");
		sb.append('\n');
		sb.append("0 -1 0 0 0 130");
		sb.append('\n');

		double max = Double.MIN_VALUE;
		double min = Double.MAX_VALUE;
		double off = 0;
		double num;
		for (int i = 0; i < mat.rows(); ++i) {
			for (int j = 0; j < mat.cols(); ++j) {
				num = mat.get(i,j)[0];
				if(num > max)
					max = num;
				if(num < min)
					min = num;
			}
		}
		
		if(min < 0){
			off = min;
		}
		else{
			off = -min;
		}

		min += off;
		max += off;
		
		for (int i = 0; i < mat.rows(); ++i) {
			for (int j = 0; j < mat.cols(); ++j) {
				num = (int)(((mat.get(i, j)[0] + off)/max )*255);
				
				sb.append(i);
				sb.append(' ');
				sb.append('0');
				sb.append(' ');
				sb.append(j);
				sb.append(' ');
				sb.append(num);
				sb.append(" 0 0");
				sb.append('\n');
				
//				sb.append(i);
//				sb.append(' ');
//				sb.append('0');
//				sb.append(' ');
//				sb.append(j);
//				sb.append(' ');
//				sb.append("100 0 0");
//				sb.append('\n');
//				
//				sb.append(i + mat.get(i, j)[0]);
//				sb.append(' ');
//				sb.append('0');
//				sb.append(' ');
//				sb.append(j + mat.get(i, j)[1]);
//				sb.append(' ');
//				sb.append("0 0 100");
//				sb.append('\n');
			}
		}
		

	

		writeToFile(sb.toString());
	}
	
	
	public static void writeGrayOFMatToXYZFile(Mat mat,Mat img1, Mat img2) {
		StringBuilder sb = new StringBuilder();
		DecimalFormat df = new DecimalFormat("############0.000############");

		sb.append("0 0 0 255 0 0");
		sb.append('\n');
		sb.append("1 0 0 130 0 0");
		sb.append('\n');

		sb.append("0 0 1 0 130 0");
		sb.append('\n');
		sb.append("0 0 -1 0 130 0");
		sb.append('\n');

		sb.append("0 1 0 0 0 130");
		sb.append('\n');
		sb.append("0 -1 0 0 0 130");
		sb.append('\n');

		double max = Double.MIN_VALUE;
		double min = Double.MAX_VALUE;
		double off = 0;
		double num;
		for (int i = 0; i < mat.rows(); ++i) {
			for (int j = 0; j < mat.cols(); ++j) {
				num = mat.get(i,j)[0];
				if(num > max)
					max = num;
				if(num < min)
					min = num;
			}
		}
		
		if(min < 0){
			off = min;
		}
		else{
			off = -min;
		}

		min += off;
		max += off;
		
		for (int i = 0; i < mat.rows(); ++i) {
			for (int j = 0; j < mat.cols(); ++j) {
				num = (int)(((mat.get(i, j)[0] + off)/max )*255);
				
				sb.append(i);
				sb.append(' ');
				sb.append('0');
				sb.append(' ');
				sb.append(j);
				sb.append(' ');
				sb.append((int)num);
				sb.append(" 0 0");
				sb.append('\n');
				
				sb.append(i);
				sb.append(' ');
				sb.append('0');
				sb.append(' ');
				sb.append(j-mat.cols());
				sb.append(' ');
				sb.append((int)img1.get(i, j)[0]);
				sb.append(' ');
				sb.append((int)img1.get(i, j)[0]);
				sb.append(' ');
				sb.append((int)img1.get(i, j)[0]);
				sb.append('\n');
				
				sb.append(i);
				sb.append(' ');
				sb.append('0');
				sb.append(' ');
				sb.append(j+mat.cols());
				sb.append(' ');
				sb.append((int)img2.get(i, j)[0]);
				sb.append(' ');
				sb.append((int)img2.get(i, j)[0]);
				sb.append(' ');
				sb.append((int)img2.get(i, j)[0]);
				sb.append('\n');
			}
		}
		

	

		writeToFile(sb.toString());
	}
	
	public static void writeRGBMatToXYZFile(Mat mat) {
		StringBuilder sb = new StringBuilder();
		DecimalFormat df = new DecimalFormat("############0000.0000############");

		sb.append("0 0 0 255 0 0");
		sb.append('\n');
		sb.append("1 0 0 130 0 0");
		sb.append('\n');

		sb.append("0 0 1 0 130 0");
		sb.append('\n');
		sb.append("0 0 -1 0 130 0");
		sb.append('\n');

		sb.append("0 1 0 0 0 130");
		sb.append('\n');
		sb.append("0 -1 0 0 0 130");
		sb.append('\n');

		for (int i = 0; i < mat.rows(); ++i) {
			for (int j = 0; j < mat.cols(); ++j) {
				sb.append(i);
				sb.append(' ');
				sb.append('0');
				sb.append(' ');
				sb.append(j);
				sb.append(' ');
				sb.append((int) (mat.get(i, j)[0]));
				sb.append(' ');
				sb.append((int) (mat.get(i, j)[1]));
				sb.append(' ');
				sb.append((int) (mat.get(i, j)[2]));
				sb.append('\n');
			}
		}

		writeToFile(sb.toString());
	}

	public static void writeRgbOFMatToXYZFile(Mat mat, Mat img1, Mat img2) {
		StringBuilder sb = new StringBuilder();
		DecimalFormat df = new DecimalFormat("############0.000############");

		sb.append("0 0 0 255 0 0");
		sb.append('\n');
		sb.append("1 0 0 130 0 0");
		sb.append('\n');

		sb.append("0 0 1 0 130 0");
		sb.append('\n');
		sb.append("0 0 -1 0 130 0");
		sb.append('\n');

		sb.append("0 1 0 0 0 130");
		sb.append('\n');
		sb.append("0 -1 0 0 0 130");
		sb.append('\n');
		
		for (int i = 0; i < mat.rows(); ++i) {
			for (int j = 0; j < mat.cols(); ++j) {
				
				sb.append(i+mat.get(i, j)[0]);
				sb.append(' ');
				sb.append("30");
				sb.append(' ');
				sb.append(j+mat.get(i, j)[1]-mat.cols());
				sb.append(' ');
				sb.append((int)img1.get(i, j)[0]);
				sb.append(' ');
				sb.append((int)img1.get(i, j)[1]);
				sb.append(' ');
				sb.append((int)img1.get(i, j)[2]);
				sb.append('\n');
				
				sb.append(i);
				sb.append(' ');
				sb.append('0');
				sb.append(' ');
				sb.append(j+mat.cols());
				sb.append(' ');
				sb.append((int)img1.get(i, j)[0]);
				sb.append(' ');
				sb.append((int)img1.get(i, j)[1]);
				sb.append(' ');
				sb.append((int)img1.get(i, j)[2]);
				sb.append('\n');
				
				sb.append(i);
				sb.append(' ');
				sb.append('0');
				sb.append(' ');
				sb.append(j-mat.cols());
				sb.append(' ');
				sb.append((int)img2.get(i, j)[0]);
				sb.append(' ');
				sb.append((int)img2.get(i, j)[1]);
				sb.append(' ');
				sb.append((int)img2.get(i, j)[2]);
				sb.append('\n');
			}
		}
		

	

		writeToFile(sb.toString());
	}
}
