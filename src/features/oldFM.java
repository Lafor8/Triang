package features;

public class oldFM {

}
/*
 * package features;

import java.util.ArrayList;
import java.util.List;

import logging.XYZConverter;
import motionestimation.DevicePose;

import org.opencv.android.BaseLoaderCallback;
import org.opencv.android.CameraBridgeViewBase;
import org.opencv.android.CameraBridgeViewBase.CvCameraViewFrame;
import org.opencv.android.CameraBridgeViewBase.CvCameraViewListener2;
import org.opencv.android.LoaderCallbackInterface;
import org.opencv.android.OpenCVLoader;
import org.opencv.calib3d.Calib3d;
import org.opencv.core.Core;
import org.opencv.core.CvType;
import org.opencv.core.Mat;
import org.opencv.core.MatOfByte;
import org.opencv.core.MatOfFloat;
import org.opencv.core.MatOfKeyPoint;
import org.opencv.core.MatOfPoint2f;
import org.opencv.core.Point;
import org.opencv.core.Scalar;
import org.opencv.core.Size;
import org.opencv.features2d.FeatureDetector;
import org.opencv.features2d.KeyPoint;
import org.opencv.video.Video;

import android.app.Activity;
import android.util.Log;
import android.view.SurfaceView;
import dlsu.vins.R;
import ekf.PointDouble;

public class FeatureManager implements CvCameraViewListener2 {
	private static final String TAG = "Feature Manager";
	private final Scalar BLACK = new Scalar(0);
	private final Scalar WHITE = new Scalar(255);

	private BaseLoaderCallback loaderCallback;

	private FeatureDetector detector;

	private int frames = 0; // TODO: ivan sir what is this for, i dunno
	private CameraBridgeViewBase cameraView;

	// Optical flow fields
	private MatOfPoint2f prevCurrent;
	private MatOfPoint2f prevNew;
	private Mat prevImage;
	private Mat currentImage;

	// Triangulation fields
	private Size imageSize;
	private Mat cameraMatrix, distCoeffs, Rot2, Rot1, T2, T1;
	private Mat R1, R2, P1, P2, Q;
	private Mat points4D;
	private Mat F, E, W;
	private Mat u, w, vt;
	private Mat nullMatF, tempMat, RotW, RotFinal;

	static boolean get = false;

	private FeatureManagerListener listener;

	private Activity caller;

	public FeatureManager(Activity caller, FeatureManagerListener listener) {
		Log.i(TAG, "constructed");

		this.caller = caller;
		Log.i(TAG, "Trying to load OpenCV library");
		this.listener = listener;
		initLoader(caller);

		cameraView = (CameraBridgeViewBase) caller.findViewById(R.id.surface_view);
		// http://stackoverflow.com/a/17872107
		// cameraView.setMaxFrameSize(720, 1280); // sets to 720 x 480
		cameraView.setMaxFrameSize(400, 1280); // sets to 320 x 240
		cameraView.setVisibility(SurfaceView.VISIBLE);

		cameraView.setCvCameraViewListener(this);

		if (!OpenCVLoader.initAsync(OpenCVLoader.OPENCV_VERSION_2_4_9, caller, loaderCallback)) {
			Log.e(TAG, "Cannot connect to OpenCV Manager");
		}
	}

	private void initLoader(Activity caller) {
		loaderCallback = new BaseLoaderCallback(caller) {
			@Override
			public void onManagerConnected(int status) {
				switch (status) {
				case LoaderCallbackInterface.SUCCESS: {
					Log.i(TAG, "OpenCV loaded successfully");
					cameraView.enableView();
					prevCurrent = new MatOfPoint2f();
					prevNew = new MatOfPoint2f();
					prevImage = new Mat();
					detector = FeatureDetector.create(FeatureDetector.FAST);
					listener.initDone();
				}
					break;
				default: {
					super.onManagerConnected(status);
				}
					break;
				}
			}
		};
	}

	private void initRectifyVariables() {
		// INITIALIZATION FOR STEREORECTIFY()

		// INPUT VARIABLES

		cameraMatrix = Mat.zeros(3, 3, CvType.CV_64F);
		distCoeffs = Mat.zeros(5, 1, CvType.CV_64F);
		imageSize = new Size(240, 320);

		Rot1 = Mat.eye(3, 3, CvType.CV_64F);
		Rot2 = Mat.eye(3, 3, CvType.CV_64F);
		T1 = Mat.ones(3, 1, CvType.CV_64F);
		T2 = Mat.ones(3, 1, CvType.CV_64F);

		// CALIBRATION RESULTS FOR 320 x 240
		cameraMatrix.put(0, 0, 287.484405747163);
		cameraMatrix.put(0, 2, 119.5);
		cameraMatrix.put(1, 1, 287.484405747163);
		cameraMatrix.put(1, 2, 159.5);
		cameraMatrix.put(2, 2, 1);

		distCoeffs.put(0, 0, 0.1831508618865668);
		distCoeffs.put(1, 0, -0.8391135375141514);
		distCoeffs.put(2, 0, 0);
		distCoeffs.put(3, 0, 0);
		distCoeffs.put(4, 0, 1.067914298622483);

		// OUTPUT VARIABLES

		R1 = Mat.zeros(3, 3, CvType.CV_64F);
		R2 = Mat.zeros(3, 3, CvType.CV_64F);
		P1 = Mat.zeros(3, 4, CvType.CV_64F);
		P2 = Mat.zeros(3, 4, CvType.CV_64F);
		Q = Mat.zeros(4, 4, CvType.CV_64F);

		// INITIALIZATION END

		// CALL STEREORECTIFY EACH FRAME AFTER THE FIRST
		// JUST PASS A NEW ROTATION AND TRANSLATION MATRIX

		// CALIBRATION RESULTS FOR 1920 x 1080
		// cameraMatrix.put(0, 0, 1768.104971372035, 0, 959.5);
		// cameraMatrix.put(1, 0, 0, 1768.104971372035, 539.5);
		// cameraMatrix.put(2, 0, 0, 0, 1);
		//
		// distCoeffs.put(0, 0, 0.1880897270445046);
		// distCoeffs.put(1, 0, -0.7348187497379466);
		// distCoeffs.put(2, 0, 0);
		// distCoeffs.put(3, 0, 0);
		// distCoeffs.put(4, 0, 0.6936210153459164);
	}

	public void onCameraViewStarted(int width, int height) {
		nullMatF = Mat.zeros(0, 0, CvType.CV_64F);
		initRectifyVariables();
	}

	public void onCameraViewStopped() {
	}

	static boolean firstImg = true;
	static Mat rgbMap;

	public Mat onCameraFrame(CvCameraViewFrame inputFrame) {
		Log.d("VINS", "onCameraFrame");

		// Mat mat = Mat.zeros(0, 0,CvType.CV_32F), temp = inputFrame.gray().t();
		// Core.flip(temp, mat, 1);

		Mat mat, temp;

		mat = Mat.zeros(0, 0, CvType.CV_32F);
		temp = inputFrame.rgba().t();
		Core.flip(temp, mat, 1);
		rgbMap = mat;

		mat = Mat.zeros(0, 0, CvType.CV_32F);
		temp = inputFrame.gray().t();
		Core.flip(temp, mat, 1);
		currentImage = mat;
		

		if (firstImg) {// there is a case that the whole image is filed with zeros
			if (inputFrame.gray().get(0, 0)[0] != 0) {
				// XYZConverter.writeRGBMatToXYZFile(mat);
				firstImg = false;
			} else {
				currentImage = null;
			}
		}

		frames++;
		return inputFrame.gray();
	}

	public FeatureUpdate getFeatureUpdate(DevicePose devicePose) {
		Log.d(TAG, "Getting Feature Update");

		if (currentImage == null)
			return null;

		Mat detectMask = currentImage.clone();
		detectMask.setTo(WHITE);

		FeatureUpdate update = new FeatureUpdate();
		Log.i(TAG, "prevCurrent: " + prevCurrent.size() + "\nprevNew: " + prevNew.size());

		MatOfPoint2f oldPt, newPt;
		Log.i("Image", "Image Captured");
		int div = 2;
		if (!prevImage.empty()) {// prevCurrent.size().height +
									// prevNew.size().height > 0) {
			Log.i("ENTERED", "here");
			// // Optical Flow

			MatOfByte status = new MatOfByte();
			MatOfFloat err = new MatOfFloat();
			MatOfPoint2f nextFeatures = new MatOfPoint2f();

			int prevCurrentSize = (int) prevCurrent.size().height; // whut
			if (prevNew != null && prevNew.size().height > 0)
				prevCurrent.push_back(prevNew); // combined

			Log.i("ENTERED", "prefarne");
			Mat flow = Mat.zeros(0, 0, CvType.CV_32FC2);
			Video.calcOpticalFlowFarneback(prevImage.submat(0, prevImage.rows() / div, 0, prevImage.cols() / div),
					currentImage.submat(0, currentImage.rows() / div, 0, currentImage.cols() / div), flow, 0.5, 5, 150, 60, 7, 1.5, Video.OPTFLOW_FARNEBACK_GAUSSIAN);

//			 if (!get) {
//			 //Log.i("FILE WRITE", "written " + points4D.width() + " points");
//			
//			 XYZConverter.writeGrayOFMatToXYZFile(flow,prevImage.submat(0, prevImage.rows() / div, 0, prevImage.cols() / div), currentImage.submat(0, currentImage.rows() / div, 0, currentImage.cols() / div));//.writeRGBPointsToXYZFile(points4D, rgbMap,div);
//			 get = true;
//			 }

			Log.i("ENTERED", "postfarne");
			// 60 80
			Log.i("Farne", flow.rows() + " " + flow.cols());

			if (!get) {
			//	Log.i("DERP", flow.dump());
			}

			double off[];
			oldPt = new MatOfPoint2f();
			newPt = new MatOfPoint2f();
			Point oldPoint = new Point();
			Point newPoint = new Point();
			for (int i = 0; i < flow.rows(); ++i)
				for (int j = 0; j < flow.cols(); ++j) {
					off = flow.get(i, j);
					oldPoint.x = i;
					oldPoint.y = j;
					newPoint.x = i + off[0];
					newPoint.y = j + off[1];
					oldPt.push_back(new MatOfPoint2f(oldPoint));
					newPt.push_back(new MatOfPoint2f(newPoint));
				}

			// Video.calcOpticalFlowPyrLK(prevImage, currentImage, prevCurrent,
			// nextFeatures, status, err);

			// Use status to filter out good points from bad

			// List<Point> oldPoints = prevCurrent.toList();
			// List<Point> newPoints = nextFeatures.toList();
			List<Point> oldPoints = oldPt.toList();
			List<Point> newPoints = newPt.toList();

			// List<Point> goodOldList = new ArrayList<>();
			// List<Point> goodNewList = new ArrayList<>();
			List<Point> goodOldList = oldPoints;
			List<Point> goodNewList = newPoints;
			List<Integer> badPointsIndex = new ArrayList<>();

			int index = 0;
			int currentSize = 0;
			// for (Byte item : status.toList()) {
			// if (item.intValue() == 1) {
			// if (index < prevCurrentSize)
			// currentSize++;
			// goodOldList.add(oldPoints.get(index));
			// goodNewList.add(newPoints.get(index));
			// Core.circle(detectMask, newPoints.get(index), 10, BLACK, -1); //
			// mask
			// // out
			// // during
			// // detection
			// } else if (index < prevCurrentSize) { // TODO: double check sir
			// badPointsIndex.add(Integer.valueOf(index));
			// }
			// index++;
			// }

			// Convert from List to OpenCV matrix for triangulation

			MatOfPoint2f goodOld = new MatOfPoint2f();
			MatOfPoint2f goodNew = new MatOfPoint2f();
			goodOld.fromList(goodOldList);
			goodNew.fromList(goodNewList);

			Log.i("pointsassign", goodOld.size() + " " + goodNew.size());

			// Triangulation

			// TODO: might want to initialize points4D with a large Nx4 Array
			// so that both memory and time will be saved (instead of
			// reallocation each time)
			// TODO: consider separating triangulation into different class

			if (!goodOld.empty() && !goodNew.empty()) {
				// SOLVING FOR THE ROTATION AND TRANSLATION MATRICES

				// GETTING THE FUNDAMENTAL MATRIX
				Log.i("fun","1");
				do {
					F = Calib3d.findFundamentalMat(goodOld, goodNew);

					cameraMatrix = cameraMatrix.clone();

					tempMat = nullMatF.clone();
					E = nullMatF.clone();

					// GETTING THE ESSENTIAL MATRIX

					Core.gemm(cameraMatrix.t(), F, 1, nullMatF, 0, tempMat);
					Core.gemm(tempMat, cameraMatrix, 1, nullMatF, 0, E);
					
					Log.i("fun","2");
					
					if (Math.abs(Core.determinant(E)) > 1e-07)
						continue;

					W = Mat.zeros(3, 3, CvType.CV_64F);
					W.put(0, 0, 0, -1, 0);
					W.put(1, 0, 1, 0, 0);
					W.put(2, 0, 0, 0, 1);
					u = nullMatF.clone();
					w = nullMatF.clone();
					vt = nullMatF.clone();

					// DECOMPOSING ESSENTIAL MATRIX TO GET THE ROTATION AND
					// TRANSLATION MATRICES

					// Decomposing Essential Matrix to obtain Rotation and
					// Translation Matrices

					Core.SVDecomp(E, w, u, vt);

					Log.i("fun","3");
					
					// check if first and second singular values are the same (as they should be)
					double singular_values_ratio = Math.abs(w.get(0, 0)[0]) / Math.abs(w.get(1, 0)[0]);
					if (singular_values_ratio > 1.0)
						singular_values_ratio = 1.0 / singular_values_ratio; // flip ratio to keep it [0,1]
					if (singular_values_ratio < 0.7) {
						Log.i("Der", "Singular values too far apart");
						continue;
					}
					
					Log.i("fun","4");

					Core.gemm(u, W, 1, nullMatF, 0, tempMat);
					Core.gemm(tempMat, vt, 1, nullMatF, 0, Rot1);
					T1 = u.col(2);

					Core.gemm(u, W.inv(), 1, nullMatF, 0, tempMat);
					Core.gemm(tempMat, vt, 1, nullMatF, 0, Rot2);
					T2 = u.col(2).mul(Mat.ones(u.col(2).size(), u.col(2).type()), -1);

					Log.i("fun","5");
					
					if (Core.determinant(R1) + 1.0 < 1e-09) {
						// according to http://en.wikipedia.org/wiki/Essential_matrix#Showing_that_it_is_valid
						E.mul(Mat.ones(E.size(), E.type()), -1);

						Core.SVDecomp(E, w, u, vt);

						Log.i("fun","6");
						
						// check if first and second singular values are the same (as they should be)
						singular_values_ratio = Math.abs(w.get(0, 0)[0]) / Math.abs(w.get(1, 0)[0]);
						if (singular_values_ratio > 1.0)
							singular_values_ratio = 1.0 / singular_values_ratio; // flip ratio to keep it [0,1]
						if (singular_values_ratio < 0.7) {
							Log.i("Der", "Singular values too far apart");
							continue;
						}
						
						Log.i("fun","7");

						Core.gemm(u, W, 1, nullMatF, 0, tempMat);
						Core.gemm(tempMat, vt, 1, nullMatF, 0, Rot1);
						T1 = u.col(2);

						Core.gemm(u, W.inv(), 1, nullMatF, 0, tempMat);
						Core.gemm(tempMat, vt, 1, nullMatF, 0, Rot2);
						T2 = u.col(2).mul(Mat.ones(u.col(2).size(), u.col(2).type()), -1);
					}
					
					Log.i("fun","8");
					
					if (Math.abs(Core.determinant(R1)) - 1.0 > 1e-07) {
						Log.i("Der", "This is not a rotation Matrix");
						continue;
					}

					P1.put(0, 0, 0, 0, 0);
					P1.put(1, 0, 0, 0, 0);
					P1.put(2, 0, 0, 0, 0);

					P2.put(0, 0, Rot2.get(0, 0)[0], Rot2.get(0, 1)[0], Rot2.get(0, 2)[0], T2.get(0, 0)[0]);
					P2.put(1, 0, Rot2.get(1, 0)[0], Rot2.get(1, 1)[0], Rot2.get(1, 2)[0], T2.get(1, 0)[0]);
					P2.put(2, 0, Rot2.get(2, 0)[0], Rot2.get(2, 1)[0], Rot2.get(2, 2)[0], T2.get(2, 0)[0]);
				
					Log.i("fun","9");
					break;
					
					
				} while (true);
				// (DEBUG) LOGGING THE VARIOUS MATRICES

				// Log.i("E", E.dump());
				// Log.i("K", cameraMatrix.dump());
				// Log.i("F", F.dump());
				// Log.i("K", cameraMatrix.dump());
				// Log.i("E", E.dump());
				// Log.i("w", w.dump());
				// Log.i("u", u.dump());
				// Log.i("vt", vt.dump());
				// Log.i("r", Rot.dump());
				Log.i("t", T2.dump());
				// Log.i("nullMatF", nullMatF.dump());

				Log.i("fun","10");
				// RotW = Mat.zeros(3, 3, CvType.CV_64F);
				// RotW.put(0, 0, devicePose.getRotWorld().getArray()[0]);
				// RotW.put(1, 0, devicePose.getRotWorld().getArray()[1]);
				// RotW.put(2, 0, devicePose.getRotWorld().getArray()[2]);
				// RotFinal = Mat.zeros(3, 3, CvType.CV_64F);
				// Core.gemm(RotW, Rot, 1, Mat.zeros(0, 0, CvType.CV_64F), 0, RotFinal);

				points4D = Mat.zeros(0, 4, CvType.CV_64F);

				// RotFinal = Mat.zeros(3, 3, CvType.CV_64F);
				//
				// RotFinal.put(0, 0, Math.cos(Math.PI / 4));
				// RotFinal.put(0, 2, Math.sin(Math.PI / 4));
				// RotFinal.put(1, 1, 1);
				// RotFinal.put(2, 0, -Math.sin(Math.PI / 4));
				// RotFinal.put(2, 2, Math.cos(Math.PI / 4));
				// T.put(0, 0, 500);
				// T.put(1, 0, 0);
				// T.put(2, 0, 500);

				// Point a = new Point(119, 159);
				// goodOld = new MatOfPoint2f();
				// goodNew = new MatOfPoint2f();
				// goodOld.fromArray(a);
				// goodNew.fromArray(a);

				// Calib3d.stereoRectify(cameraMatrix, distCoeffs, cameraMatrix.clone(), distCoeffs.clone(), imageSize, RotFinal, T, R1, R2, P1, P2, Q);
				// Calib3d.triangulatePoints(P1, P2, goodOld, goodNew,
				// points4D);
				points4D = triangulatePoints(goodOld, goodNew, cameraMatrix, P1, P2);

				// double transPixel[] = T.t().get(0, 0);
				// double transMetric[] = { devicePose.get_xPos(),
				// devicePose.get_yPos(), devicePose.get_zPos() };
				// double metricScale = 0;
				//
				// for (int i = 0; i < transPixel.length; ++i)
				// metricScale += transMetric[i] / transPixel[i];
				// metricScale /= 3;

				// TODO: maybe this method is more optimized??
				// Mat points3D = new Mat();
				// Calib3d.convertPointsFromHomogeneous(points4D, points3D);

				// becomes 2n: n (with method above) + n (iterating to split
				// into current and new
				// I mean, sure, 2n

				// Split points to current and new PointDouble
				// TODO verify this shit
				// TODO yass corrected

				Log.i(TAG, "points4D size: " + points4D.size().height + " " + points4D.size().width + " " + points4D.channels());
				Log.i(TAG, "points4D size: " + points4D.rows() + " " + points4D.cols());
				// Log.i(TAG, points4D.dump());
				// Log.i(TAG, T.dump());
				// Log.i(TAG, devicePose.toString());

				PointDouble point = null;
				List<PointDouble> current2d = new ArrayList<>();
				List<PointDouble> new2d = new ArrayList<>();

				if (!get) {
					Log.i("FILE WRITE", "written " + points4D.width() + " points");

					XYZConverter.writeRGBPointsToXYZFile(points4D, rgbMap, div);
					get = true;
				}
//				for (int i = 0; i < goodOld.height(); i++) {
//					// Log.i("trinaufng", goodOld.height() + "");
//					double x = points4D.get(0, i)[0] / points4D.get(3, i)[0];
//					double y = points4D.get(2, i)[0] / points4D.get(3, i)[0];
//
//					point = new PointDouble(x, y);
//					if (i < currentSize) {
//						current2d.add(point);
//					} else {
//						new2d.add(point);
//					}
//
//					final String str;
//					String strTemp = "";
//					Mat MatTemp = points4D.t();
//					// write 3D homogenous point
//					strTemp += (i + 1) + ": [" + MatTemp.get(0, 0)[0] + ", " + MatTemp.get(0, 1)[0] + ", " + MatTemp.get(0, 2)[0] + ", " + MatTemp.get(0, 3)[0] + "]\n";
//
//					// write 3D euclidean point
//					strTemp += point + "\n";
//					strTemp += cameraMatrix.dump() + "\n";
//					strTemp += P1.dump() + "\n";
//					strTemp += P2.dump();
//					str = strTemp;
//
//					// caller.runOnUiThread(new Runnable() {
//					// @Override
//					// public void run() {
//					// TextView tv = (TextView)
//					// caller.findViewById(R.id.debugTextView);
//					//
//					// tv.setText(str + "\n" + tv.getText());
//					// }
//					// });
//				}
//
//				Log.i("FM Feats", point.toString() + "");
//
//				update.setCurrentPoints(current2d);
//				update.setNewPoints(new2d);
			}
//			update.setBadPointsIndex(badPointsIndex);
			goodNew.copyTo(prevCurrent);
		}

		// Detect new points based on optical flow mask
		MatOfKeyPoint newFeatures = new MatOfKeyPoint();
		detector.detect(currentImage, newFeatures, detectMask);
		if (newFeatures.size().height > 0) {
			prevNew = convert(newFeatures);
		}

		currentImage.copyTo(prevImage);
		return update;
	}

	private Mat triangulatePoints(MatOfPoint2f goodOld, MatOfPoint2f goodNew, Mat k, Mat p1, Mat p2) {
		k = k.inv();

		Mat points4D = Mat.zeros(goodNew.height(), 0, CvType.CV_64F);

		Mat u1 = Mat.zeros(1, 3, CvType.CV_64F);
		Mat u2 = Mat.zeros(1, 3, CvType.CV_64F);
		Mat um1 = Mat.zeros(3, 1, CvType.CV_64F);
		Mat um2 = Mat.zeros(3, 1, CvType.CV_64F);
		Mat x;

		for (int i = 0; i < goodOld.height(); i++) {
			// Log.i("triangU",
			// goodOld.get(i + 1, 0)[0] + " " + goodOld.get(i + 1, 0)[1] + " " +
			// goodNew.get(i + 1, 0)[0] + " "
			// + goodNew.get(i + 1, 0)[1]);
			u1.put(0, 0, goodOld.get(i, 0)[0], goodOld.get(i, 0)[1], 1);

			// Log.i("triangU " + i, k.rows() + "x" + k.cols() + " " +
			// u1.t().rows() + "x" + u1.t().cols());
			Core.gemm(k, u1.t(), 1, nullMatF, 0, um1);

			// Log.i("triangU " + i, k.rows() + "x" + k.cols() + " " +
			// u1.t().rows() + "x" + u1.t().cols());

			u2.put(0, 0, goodNew.get(i, 0)[0], goodNew.get(i, 0)[1], 1);
			Core.gemm(k, u2.t(), 1, nullMatF, 0, um2);

			x = linearLSTriangulation(um1.t(), p1, um2.t(), p2);
			points4D.push_back(x.t());
		}

		// 1x3 3 1 3
		// Log.i("Trinaug", points4D.size() + " " + points4D.rows() + " " +
		// points4D.width() + " " + points4D.height());

		return points4D.t();
	}

	private Mat linearLSTriangulation(Mat u, Mat p1, Mat u1, Mat p2) {
		// build matrix A for homogenous equation system Ax = 0
		// assume X = (x,y,z,1), for Linear-LS method
		// which turns it into a AX = B system, where A is 4x3, X is 3x1 and B
		// is 4x1

		// Matx43d A( u.x*P(2,0)-P(0,0), u.x*P(2,1)-P(0,1), u.x*P(2,2)-P(0,2),
		// u.y*P(2,0)-P(1,0), u.y*P(2,1)-P(1,1), u.y*P(2,2)-P(1,2),
		// u1.x*P1(2,0)-P1(0,0), u1.x*P1(2,1)-P1(0,1), u1.x*P1(2,2)-P1(0,2),
		// u1.y*P1(2,0)-P1(1,0), u1.y*P1(2,1)-P1(1,1), u1.y*P1(2,2)-P1(1,2)
		// );
		// Mat_ B = (Mat_(4,1) << -(u.x*P(2,3) -P(0,3)),
		// -(u.y*P(2,3) -P(1,3)),
		// -(u1.x*P1(2,3) -P1(0,3)),
		// -(u1.y*P1(2,3) -P1(1,3)));
		//
		// Mat_ X;
		// solve(A,B,X,DECOMP_SVD);
		//
		// return X;

		double arrMatA[][] = {
				{ u.get(0, 0)[0] * p1.get(2, 0)[0] - p1.get(0, 0)[0], u.get(0, 0)[0] * p1.get(2, 1)[0] - p1.get(0, 1)[0], u.get(0, 0)[0] * p1.get(2, 2)[0] - p1.get(0, 2)[0] },

				{ u.get(0, 1)[0] * p1.get(2, 0)[0] - p1.get(1, 0)[0], u.get(0, 1)[0] * p1.get(2, 1)[0] - p1.get(1, 1)[0], u.get(0, 1)[0] * p1.get(2, 2)[0] - p1.get(1, 2)[0] },

				{ u1.get(0, 0)[0] * p2.get(2, 0)[0] - p2.get(0, 0)[0], u1.get(0, 0)[0] * p2.get(2, 1)[0] - p2.get(0, 1)[0], u1.get(0, 0)[0] * p2.get(2, 2)[0] - p2.get(0, 2)[0] },

				{ u1.get(0, 1)[0] * p2.get(2, 0)[0] - p2.get(1, 0)[0], u1.get(0, 1)[0] * p2.get(2, 1)[0] - p2.get(1, 1)[0], u1.get(0, 1)[0] * p2.get(2, 2)[0] - p2.get(1, 2)[0] } };

		Mat A = Mat.zeros(4, 3, CvType.CV_64F);
		A.put(0, 0, arrMatA[0]);
		A.put(1, 0, arrMatA[1]);
		A.put(2, 0, arrMatA[2]);
		A.put(3, 0, arrMatA[3]);

		double arrMatB[][] = { { -(u.get(0, 0)[0] * p1.get(2, 3)[0] - p1.get(0, 3)[0]) }, { -(u.get(0, 1)[0] * p1.get(2, 3)[0] - p1.get(1, 3)[0]) },
				{ -(u1.get(0, 0)[0] * p2.get(2, 3)[0] - p2.get(0, 3)[0]) }, { -(u1.get(0, 1)[0] * p2.get(2, 3)[0] - p2.get(1, 3)[0]) }, };

		Mat B = Mat.zeros(4, 1, CvType.CV_64F);
		B.put(0, 0, arrMatB[0]);
		B.put(1, 0, arrMatB[1]);
		B.put(2, 0, arrMatB[2]);
		B.put(3, 0, arrMatB[3]);

		Mat x = Mat.zeros(0, 0, CvType.CV_64F);
		Core.solve(A, B, x, Core.DECOMP_SVD);
		// Log.i("tiranugn", x.dump());
		x.push_back(Mat.ones(1, 1, CvType.CV_64F));
		// Log.i("tiranugn", x.dump());
		return x;
	}

	private MatOfPoint2f convert(MatOfKeyPoint keyPoints) {
		KeyPoint[] keyPointsArray = keyPoints.toArray();
		Point[] pointsArray = new Point[keyPointsArray.length];

		for (int i = 0; i < keyPointsArray.length; i++) {
			pointsArray[i] = (Point) keyPointsArray[i].pt;
		}

		return new MatOfPoint2f(pointsArray);
	}

}

 */