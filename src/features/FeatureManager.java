package features;

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

	private BaseLoaderCallback loaderCallback;
	private CameraBridgeViewBase cameraView;

	// Optical flow fields
	private Mat prevImage;
	private Mat currentImage;

	// Triangulation fields
	private Size imageSize;
	private Mat cameraMatrix, distCoeffs, Rot2, Rot1, T2, T1;
	private Mat R1, R2, P1, P2, Q;
	private Mat points4D;
	private Mat F, E, W;
	private Mat u, w, vt;
	private Mat nullMatF, tempMat;

	private FeatureManagerListener listener;

	private Activity caller;

	public FeatureManager(Activity caller, FeatureManagerListener listener) {
		Log.i(TAG, "constructed");

		this.caller = caller;
		Log.i(TAG, "Trying to load OpenCV library");
		this.listener = listener;
		initLoader(caller);

		cameraView = (CameraBridgeViewBase) caller.findViewById(R.id.surface_view);
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
					prevImage = null;
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
	}

	public void onCameraViewStarted(int width, int height) {
		nullMatF = Mat.zeros(0, 0, CvType.CV_64F);
		initRectifyVariables();
	}

	public void onCameraViewStopped() {
	}

	static boolean firstImg = true;
	static Mat rgbMap;
	static Mat prevRgbMap;
	static Mat rgbFrame;
	static Mat currentFrame;

	public Mat onCameraFrame(CvCameraViewFrame inputFrame) {
		Log.d("VINS", "onCameraFrame");

		Mat mat, temp;

		mat = Mat.zeros(0, 0, CvType.CV_32F);
		temp = inputFrame.rgba().t();
		Core.flip(temp, mat, 1);
		rgbFrame = mat;

		mat = Mat.zeros(0, 0, CvType.CV_32F);
		temp = inputFrame.gray().t();
		Core.flip(temp, mat, 1);
		currentFrame = mat;

		// if (firstImg) {// there is a case that the whole image is filed with zeros
		// if (currentImage.get(0, 0)[0] != currentImage.get(currentImage.rows()-1, 0)[0]) {
		// // XYZConverter.writeRGBMatToXYZFile(mat);
		// Log.i("Camera Frame", "First valid image");
		// firstImg = false;
		// } else {
		// currentImage = null;
		// }
		// }
		// else{
		// if(currentImage.get(0, 0)[0] == currentImage.get(currentImage.rows()-1, 0)[0]){
		// firstImg = true;
		// currentImage = null;
		// }
		// }

		if (currentFrame.get(0, 0)[0] == currentFrame.get(currentFrame.rows() - 1, 0)[0]) {
			currentFrame = null;
			Log.i("Image Capture", "Invalid Image");
		}

		return inputFrame.gray();
	}

	public static boolean get = false;
	int divr = 2, divc = 1;
	Mat flow;

	public FeatureUpdate getFeatureUpdate(DevicePose devicePose) {
		if (currentFrame == null) {
			prevImage = null;
			return null;
		} else {
			currentImage = currentFrame;
			rgbMap = rgbFrame;
		}
		Log.i("Image", "Valid Image");

		MatOfPoint2f oldPt, newPt;

		if (prevImage != null) {

			// Optical Flow

			Log.i("ENTERED", "pre Optical Flow");

			flow = Mat.zeros(prevImage.size(), CvType.CV_32FC2);

			Mat prevSnippet = prevImage.submat(0, prevImage.rows() / divr, 0, prevImage.cols() / divc);
			Mat currentSnippet = currentImage.submat(0, currentImage.rows() / divr, 0, currentImage.cols() / divc);

			//Video.calcOpticalFlowFarneback(prevSnippet, currentSnippet, flow, 0.5, 3, 15, 3, 5, 1.2, Video.OPTFLOW_FARNEBACK_GAUSSIAN);

			// TODO: do the thigns
			//Video.calcOpticalFlowPyrLK(prevSnippet, currentSnippet, prevPts, nextPts, status, err);
			
			Mat RigidT = Video.estimateRigidTransform(prevSnippet, currentSnippet, false);
			Log.i("ARr",RigidT.dump()+"\n"+RigidT.type());
			Mat moved = Mat.zeros(1, 3, CvType.CV_64F);
			moved = moved.t();
			Mat movedPt = Mat.zeros(0, 0,CvType.CV_64F);
			for(int i = 0; i < flow.rows();++i){
				Log.i("Arr", i + "/" + flow.rows());
				for(int j=0;j<flow.cols();++i){
					moved = moved.t();
					moved.put(0, 0, i,j,1);
					moved = moved.t();
					//Log.i("ARr",RigidT.dump()+"\n"+RigidT.type());
					Core.gemm(RigidT, moved, 1, nullMatF, 0,movedPt);
					double flowPt[] = {movedPt.get(0, 0)[0]-i,movedPt.get(1, 0)[0]-j};
					flow.put(i, j, flowPt);
				}
		}
			/*
			for (int x=0; x<img_1.cols; x++) {
				for (int y=0; y<img_1.rows; y++) {
	//				Mat_<double> moved = H * (Mat_<double>(3,1) << x , y , 1);
					Mat_<double> moved = T * (Mat_<double>(3,1) << x , y , 1);
					Point2f movedpt(moved(0),moved(1));
					flow_from_features(y,x) = Point2f(movedpt.x-x,movedpt.y-y);
			*/
			 //Video.calcOpticalFlowFarneback(prevSnippet, currentSnippet, flow, 0.5, 5, 150, 60, 7, 1.5, Video.OPTFLOW_FARNEBACK_GAUSSIAN);
			 Video.calcOpticalFlowFarneback(prevSnippet, currentSnippet, flow, 0.5, 2, 40, 40, 5, 0.5, Video.OPTFLOW_USE_INITIAL_FLOW);
			 Video.calcOpticalFlowFarneback(prevSnippet, currentSnippet, flow, 0.5, 0, 25, 40, 3, 0.25, Video.OPTFLOW_USE_INITIAL_FLOW);

			if (!get) {
				 //Log.i("FILE WRITE", "written " + points4D.width() + " points");
				Log.i("FILE WRITE", "written " );
				XYZConverter.writeRgbOFMatToXYZFile(flow, this.prevRgbMap.submat(0, prevImage.rows() / divr, 0, prevImage.cols() / divc),
						this.rgbMap.submat(0, currentImage.rows() / divr, 0, currentImage.cols() / divc));// .writeRGBPointsToXYZFile(points4D, rgbMap,div);
				
//				XYZConverter.writeGrayOFMatToXYZFile(flow, prevImage.submat(0, prevImage.rows() / divr, 0, prevImage.cols() / divc),
//						currentImage.submat(0, currentImage.rows() / divr, 0, currentImage.cols() / divc));// .writeRGBPointsToXYZFile(points4D, rgbMap,div);
				get = true;
			}

			Log.i("ENTERED", "post Optical Flow");

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

			// Convert from List to OpenCV matrix for triangulation

			MatOfPoint2f goodOld = new MatOfPoint2f();
			MatOfPoint2f goodNew = new MatOfPoint2f();
			goodOld.fromList(oldPt.toList());
			goodNew.fromList(newPt.toList());

			Log.i("fun", "pre-triangulation");
			// Triangulation

			if (!goodOld.empty() && !goodNew.empty()) {
				// SOLVING FOR THE ROTATION AND TRANSLATION MATRICES

				// GETTING THE FUNDAMENTAL MATRIX

				// There is a case that fundamental matrix is not found
				Log.i("fun", "1");

				do {
					F = Calib3d.findFundamentalMat(goodOld, goodNew);

					cameraMatrix = cameraMatrix.clone();

					tempMat = nullMatF.clone();
					E = nullMatF.clone();

					// GETTING THE ESSENTIAL MATRIX

					Core.gemm(cameraMatrix.t(), F, 1, nullMatF, 0, tempMat);
					Core.gemm(tempMat, cameraMatrix, 1, nullMatF, 0, E);

					Log.i("fun", "2");

					if (Math.abs(Core.determinant(E)) > 1e-07)
						continue;

					W = Mat.zeros(3, 3, CvType.CV_64F);
					W.put(0, 0, 0, -1, 0);
					W.put(1, 0, 1, 0, 0);
					W.put(2, 0, 0, 0, 1);
					u = nullMatF.clone();
					w = nullMatF.clone();
					vt = nullMatF.clone();
					Core.SVDecomp(E, w, u, vt);

					Log.i("fun", "3");

					// check if first and second singular values are the same (as they should be)
					double singular_values_ratio = Math.abs(w.get(0, 0)[0]) / Math.abs(w.get(1, 0)[0]);
					if (singular_values_ratio > 1.0)
						singular_values_ratio = 1.0 / singular_values_ratio; // flip ratio to keep it [0,1]
					if (singular_values_ratio < 0.7) {
						Log.i("Der", "Singular values too far apart");
						continue;
					}

					Log.i("fun", "4");

					Core.gemm(u, W, 1, nullMatF, 0, tempMat);
					Core.gemm(tempMat, vt, 1, nullMatF, 0, Rot1);
					T1 = u.col(2);

					Core.gemm(u, W.inv(), 1, nullMatF, 0, tempMat);
					Core.gemm(tempMat, vt, 1, nullMatF, 0, Rot2);
					T2 = u.col(2).mul(Mat.ones(u.col(2).size(), u.col(2).type()), -1);

					Log.i("fun", "5");

					if (Core.determinant(R1) + 1.0 < 1e-09) {
						// according to http://en.wikipedia.org/wiki/Essential_matrix#Showing_that_it_is_valid
						E.mul(Mat.ones(E.size(), E.type()), -1);

						Core.SVDecomp(E, w, u, vt);

						Log.i("fun", "6");

						// check if first and second singular values are the same (as they should be)
						singular_values_ratio = Math.abs(w.get(0, 0)[0]) / Math.abs(w.get(1, 0)[0]);
						if (singular_values_ratio > 1.0)
							singular_values_ratio = 1.0 / singular_values_ratio; // flip ratio to keep it [0,1]
						if (singular_values_ratio < 0.7) {
							Log.i("Der", "Singular values too far apart");
							continue;
						}

						Log.i("fun", "7");

						Core.gemm(u, W, 1, nullMatF, 0, tempMat);
						Core.gemm(tempMat, vt, 1, nullMatF, 0, Rot1);
						T1 = u.col(2);

						Core.gemm(u, W.inv(), 1, nullMatF, 0, tempMat);
						Core.gemm(tempMat, vt, 1, nullMatF, 0, Rot2);
						T2 = u.col(2).mul(Mat.ones(u.col(2).size(), u.col(2).type()), -1);
					}

					Log.i("fun", "8");

					if (Math.abs(Core.determinant(R1)) - 1.0 > 1e-07) {
						Log.i("Der", "This is not a rotation Matrix");
						continue;
					}

					P1.put(0, 0, 1, 0, 0, 0);
					P1.put(1, 0, 0, 1, 0, 0);
					P1.put(2, 0, 0, 0, 1, 0);

					Log.i("fun", "9");
					break;

				} while (true);

				Log.i("t", T2.dump());

				Log.i("fun", "10");

				if (!get) {

					// Combination 1
					// P2.put(0, 0, Rot1.get(0, 0)[0], Rot1.get(0, 1)[0], Rot1.get(0, 2)[0], T1.get(0, 0)[0]);
					// P2.put(1, 0, Rot1.get(1, 0)[0], Rot1.get(1, 1)[0], Rot1.get(1, 2)[0], T1.get(1, 0)[0]);
					// P2.put(2, 0, Rot1.get(2, 0)[0], Rot1.get(2, 1)[0], Rot1.get(2, 2)[0], T1.get(2, 0)[0]);
					//
					// points4D = Mat.zeros(0, 4, CvType.CV_64F);
					// points4D = triangulatePoints(goodOld, goodNew, cameraMatrix, P1, P2);
					//
					// XYZConverter.writeAllToXYZFile(flow, prevImage, currentImage, points4D, divr,divc, rgbMap, "tryC1.xyz");
					//
					// // Combination 2
					// P2.put(0, 0, Rot1.get(0, 0)[0], Rot1.get(0, 1)[0], Rot1.get(0, 2)[0], T2.get(0, 0)[0]);
					// P2.put(1, 0, Rot1.get(1, 0)[0], Rot1.get(1, 1)[0], Rot1.get(1, 2)[0], T2.get(1, 0)[0]);
					// P2.put(2, 0, Rot1.get(2, 0)[0], Rot1.get(2, 1)[0], Rot1.get(2, 2)[0], T2.get(2, 0)[0]);
					//
					// points4D = Mat.zeros(0, 4, CvType.CV_64F);
					// points4D = triangulatePoints(goodOld, goodNew, cameraMatrix, P1, P2);
					// XYZConverter.writeAllToXYZFile(flow, prevImage, currentImage, points4D, divr, divc, rgbMap, "tryC2.xyz");
					//
					// // Combination 3
					// P2.put(0, 0, Rot2.get(0, 0)[0], Rot2.get(0, 1)[0], Rot2.get(0, 2)[0], T1.get(0, 0)[0]);
					// P2.put(1, 0, Rot2.get(1, 0)[0], Rot2.get(1, 1)[0], Rot2.get(1, 2)[0], T1.get(1, 0)[0]);
					// P2.put(2, 0, Rot2.get(2, 0)[0], Rot2.get(2, 1)[0], Rot2.get(2, 2)[0], T1.get(2, 0)[0]);
					//
					// points4D = Mat.zeros(0, 4, CvType.CV_64F);
					// points4D = triangulatePoints(goodOld, goodNew, cameraMatrix, P1, P2);
					// XYZConverter.writeAllToXYZFile(flow, prevImage, currentImage, points4D, divr,divc, rgbMap, "tryC3.xyz");
					//
					// // Combination 4
					// P2.put(0, 0, Rot2.get(0, 0)[0], Rot2.get(0, 1)[0], Rot2.get(0, 2)[0], T2.get(0, 0)[0]);
					// P2.put(1, 0, Rot2.get(1, 0)[0], Rot2.get(1, 1)[0], Rot2.get(1, 2)[0], T2.get(1, 0)[0]);
					// P2.put(2, 0, Rot2.get(2, 0)[0], Rot2.get(2, 1)[0], Rot2.get(2, 2)[0], T2.get(2, 0)[0]);
					//
					// points4D = Mat.zeros(0, 4, CvType.CV_64F);
					// points4D = triangulatePoints(goodOld, goodNew, cameraMatrix, P1, P2);
					// XYZConverter.writeAllToXYZFile(flow, prevImage, currentImage, points4D, divr,divc, prevRgbMap, "tryC4.xyz");
					//
					Log.i("FILE WRITE", "written " + points4D.width() + " points");
					get = true;
				}

				// check RoT and T

				// Calib3d.triangulatePoints(P1,P2, goodOld, goodNew, points4D);
				// if (!get) {
				// Log.i("FILE WRITE", "written " + points4D.width() + " points");
				//
				// XYZConverter.writeRGBPointsToXYZFile(points4D, rgbMap, div);
				// get = true;
				// }

			}
		}

		prevImage = currentImage.clone();
		prevRgbMap = rgbMap.clone();
		return null;
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
			// Log.i("Triangulation", "Triangulating, attempt "+ i);
			u1.put(0, 0, goodOld.get(i, 0)[0], goodOld.get(i, 0)[1], 1);

			Core.gemm(k, u1.t(), 1, nullMatF, 0, um1);

			u2.put(0, 0, goodNew.get(i, 0)[0], goodNew.get(i, 0)[1], 1);
			Core.gemm(k, u2.t(), 1, nullMatF, 0, um2);

			// x = iterativeLinearLSTriangulation(um1.t(), p1, um2.t(), p2);
			x = linearLSTriangulation(um1.t(), p1, um2.t(), p2);
			points4D.push_back(x.t());
		}
		return points4D.t();
	}

	private Mat linearLSTriangulation(Mat u, Mat p1, Mat u1, Mat p2) {

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
				{ -(u1.get(0, 0)[0] * p2.get(2, 3)[0] - p2.get(0, 3)[0]) }, { -(u1.get(0, 1)[0] * p2.get(2, 3)[0] - p2.get(1, 3)[0]) } };

		Mat B = Mat.zeros(4, 1, CvType.CV_64F);
		B.put(0, 0, arrMatB[0]);
		B.put(1, 0, arrMatB[1]);
		B.put(2, 0, arrMatB[2]);
		B.put(3, 0, arrMatB[3]);

		Mat x = Mat.zeros(0, 0, CvType.CV_64F);
		Core.solve(A, B, x, Core.DECOMP_SVD);
		x.push_back(Mat.ones(1, 1, CvType.CV_64F));

		return x;
	}

	double EPSILON = 0.0001;

	private Mat iterativeLinearLSTriangulation(Mat u, Mat p1, Mat u1, Mat p2) {

		double wi = 1, wi1 = 1;

		Mat x = null;
		Mat X = linearLSTriangulation(u, p1, u1, p2);

		Mat p2xm, temp;

		double p2x, p2x1;

		// Log.i("Tri", X.dump());

		for (int i = 0; i < 10; i++) { // Hartley suggests 10 iterations at most

			// recalculate weights
			p2xm = Mat.zeros(0, 0, CvType.CV_64F);
			temp = Mat.zeros(1, 4, CvType.CV_64F);
			temp.put(0, 0, p1.get(2, 0)[0], p1.get(2, 1)[0], p1.get(2, 2)[0], p1.get(2, 3)[0]);
			Core.gemm(temp, X, 1, nullMatF, 0, p2xm);

			// Log.i("Tri", "Temp "+temp.dump() + "\nP1 " + p1.dump() + "\nP2xm "+p2xm.dump());
			//
			p2x = p2xm.get(0, 0)[0];

			p2xm = Mat.zeros(0, 0, CvType.CV_64F);
			temp = Mat.zeros(1, 4, CvType.CV_64F);
			temp.put(0, 0, p2.get(2, 0)[0], p2.get(2, 1)[0], p2.get(2, 2)[0], p2.get(2, 3)[0]);
			Core.gemm(temp, X, 1, nullMatF, 0, p2xm);

			p2x1 = p2xm.get(0, 0)[0];

			// breaking point
			if (Math.abs(wi - p2x) <= EPSILON && Math.abs(wi1 - p2x1) <= EPSILON)
				break;

			wi = p2x;
			wi1 = p2x1;

			// Log.i("Tri", wi + " " + wi1 + " " + p2x + " " + p2x1);

			double arrMatA[][] = {
					{ (u.get(0, 0)[0] * p1.get(2, 0)[0] - p1.get(0, 0)[0]) / wi, (u.get(0, 0)[0] * p1.get(2, 1)[0] - p1.get(0, 1)[0]) / wi,
							(u.get(0, 0)[0] * p1.get(2, 2)[0] - p1.get(0, 2)[0]) / wi },

					{ (u.get(0, 1)[0] * p1.get(2, 0)[0] - p1.get(1, 0)[0]) / wi, (u.get(0, 1)[0] * p1.get(2, 1)[0] - p1.get(1, 1)[0]) / wi,
							(u.get(0, 1)[0] * p1.get(2, 2)[0] - p1.get(1, 2)[0]) / wi },

					{ (u1.get(0, 0)[0] * p2.get(2, 0)[0] - p2.get(0, 0)[0]) / wi1, (u1.get(0, 0)[0] * p2.get(2, 1)[0] - p2.get(0, 1)[0]) / wi1,
							(u1.get(0, 0)[0] * p2.get(2, 2)[0] - p2.get(0, 2)[0]) / wi1 },

					{ (u1.get(0, 1)[0] * p2.get(2, 0)[0] - p2.get(1, 0)[0]) / wi1, (u1.get(0, 1)[0] * p2.get(2, 1)[0] - p2.get(1, 1)[0]) / wi1,
							(u1.get(0, 1)[0] * p2.get(2, 2)[0] - p2.get(1, 2)[0]) / wi1 } };

			Mat A = Mat.zeros(4, 3, CvType.CV_64F);
			A.put(0, 0, arrMatA[0]);
			A.put(1, 0, arrMatA[1]);
			A.put(2, 0, arrMatA[2]);
			A.put(3, 0, arrMatA[3]);

			double arrMatB[][] = { { -((u.get(0, 0)[0] * p1.get(2, 3)[0] - p1.get(0, 3)[0])) / wi }, { -((u.get(0, 1)[0] * p1.get(2, 3)[0] - p1.get(1, 3)[0])) / wi },
					{ -((u1.get(0, 0)[0] * p2.get(2, 3)[0] - p2.get(0, 3)[0])) / wi1 }, { -((u1.get(0, 1)[0] * p2.get(2, 3)[0] - p2.get(1, 3)[0])) / wi1 } };

			Mat B = Mat.zeros(4, 1, CvType.CV_64F);
			B.put(0, 0, arrMatB[0]);
			B.put(1, 0, arrMatB[1]);
			B.put(2, 0, arrMatB[2]);
			B.put(3, 0, arrMatB[3]);

			x = Mat.zeros(0, 0, CvType.CV_64F);
			Core.solve(A, B, x, Core.DECOMP_SVD);
			x.push_back(Mat.ones(1, 1, CvType.CV_64F));

			// Log.i("Tri", x.dump());
			X = x.clone();
		}
		return X;

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
