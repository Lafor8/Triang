package dlsu.vins;

import java.util.ArrayList;
import java.util.List;

import org.opencv.calib3d.Calib3d;
import org.opencv.core.Core;
import org.opencv.core.CvType;
import org.opencv.core.Mat;
import org.opencv.core.MatOfPoint2f;
import org.opencv.core.Size;

import ekf.PointDouble;

public class TriTest {
	// Triangulation fields
	private static Size imageSize;
	private static Mat cameraMatrix, distCoeffs, Rot, T;
	private static Mat R1, R2, P1, P2, Q;
	private static Mat points4D;
	private static Mat F, E, W;
	private static Mat u, w, vt;
	private static Mat nullMatF, tempMat, RotW, RotFinal;

	public static void main() {
		TriTest t = new TriTest();
		t.initRectifyVariables();

		// Convert from List to OpenCV matrix for triangulation

		MatOfPoint2f goodOld = new MatOfPoint2f();
		MatOfPoint2f goodNew = new MatOfPoint2f();

		// Triangulation

		// TODO: might want to initialize points4D with a large Nx4 Array
		// so that both memory and time will be saved (instead of
		// reallocation each time)
		// TODO: consider separating triangulation into different class

		if (!goodOld.empty() && !goodNew.empty()) {
			// SOLVING FOR THE ROTATION AND TRANSLATION MATRICES

			// GETTING THE FUNDAMENTAL MATRIX

			F = Calib3d.findFundamentalMat(goodOld, goodNew);

			cameraMatrix = cameraMatrix.clone();

			tempMat = nullMatF.clone();
			E = nullMatF.clone();

			// GETTING THE ESSENTIAL MATRIX

			Core.gemm(cameraMatrix.t(), F, 1, nullMatF, 0, tempMat);
			Core.gemm(tempMat, cameraMatrix, 1, nullMatF, 0, E);

			W = Mat.zeros(3, 3, CvType.CV_64F);
			W.put(0, 1, -1);
			W.put(1, 0, 1);
			W.put(2, 2, 1);
			u = nullMatF.clone();
			w = nullMatF.clone();
			vt = nullMatF.clone();

			// DECOMPOSING ESSENTIAL MATRIX TO GET THE ROTATION AND
			// TRANSLATION MATRICES

			// Decomposing Essential Matrix to obtain Rotation and
			// Translation Matrices

			Core.SVDecomp(E, w, u, vt);

			Core.gemm(u, W, 1, nullMatF, 0, tempMat);
			Core.gemm(tempMat, vt, 1, nullMatF, 0, Rot);
			T = u.col(2);

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
			// Log.i("t", T.dump());
			// Log.i("nullMatF", nullMatF.dump());

			// WORLD COORDINATE SYSTEM
			RotW = Mat.eye(3, 3, CvType.CV_64F);
			// RotW = Mat.zeros(3, 3, CvType.CV_64F);
			// RotW.put(0, 0, devicePose.getRotWorld().getArray()[0]);
			// RotW.put(1, 0, devicePose.getRotWorld().getArray()[1]);
			// RotW.put(2, 0, devicePose.getRotWorld().getArray()[2]);
			RotFinal = Mat.zeros(3, 3, CvType.CV_64F);
			Core.gemm(RotW, Rot, 1, Mat.zeros(0, 0, CvType.CV_64F), 0, RotFinal);
			RotFinal = Mat.eye(3, 3, CvType.CV_64F);
			
			RotFinal.put(0, 0, Math.cos(Math.PI/4), 0, Math.sin(Math.PI/4));
			RotFinal.put(2, 0, -Math.sin(Math.PI/4), 0, Math.cos(Math.PI/4));
			T.put(0, 0, 1);
			T.put(1, 0, 0);
			T.put(2, 0, 1);
			//goodOld = Mat.zeros(rows, cols, type)
			
			points4D = Mat.zeros(1, 4, CvType.CV_64F);
			Calib3d.stereoRectify(cameraMatrix, distCoeffs, cameraMatrix.clone(), distCoeffs.clone(), imageSize, RotFinal, T, R1, R2, P1,
					P2, Q);
			Calib3d.triangulatePoints(P1, P2, goodOld, goodNew, points4D);

			// alrun code
			// double observedDistance =
			// deviceCoords.computeDistanceTo(observedFeatureCoords);
			// double observedHeading = Math.atan((observedFeatureCoords.getY()
			// - deviceCoords.getY())
			// / (observedFeatureCoords.getX() - deviceCoords.getX()))
			// - this.getHeadingRadians();

			// METRIC SCALE
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

			// Log.i(TAG, "points4D size: " + points4D.size().width);
			// Log.i(TAG, T.dump());
			// Log.i(TAG, devicePose.toString());

			PointDouble point = null;
			List<PointDouble> current2d = new ArrayList<>();
			List<PointDouble> new2d = new ArrayList<>();

			System.out.println(points4D.dump());
			for (int i = 0; i < goodOld.height(); i++) {
				// double x = points4D.get(0, i)[0] * metricScale /
				// points4D.get(3, i)[0];
				// double y = points4D.get(2, i)[0] * metricScale /
				// points4D.get(3, i)[0];
				double x = points4D.get(0, i)[0] / points4D.get(3, i)[0];
				double y = points4D.get(2, i)[0] / points4D.get(3, i)[0];

				point = new PointDouble(x, y);

				System.out.println(point);
				// if (i < currentSize) {
				// current2d.add(point);
				// } else {
				// new2d.add(point);
				// }
			}
			System.out.println();
		}

	}

	private void initRectifyVariables() {
		// INITIALIZATION FOR STEREORECTIFY()

		// INPUT VARIABLES

		cameraMatrix = Mat.zeros(3, 3, CvType.CV_64F);
		distCoeffs = Mat.zeros(5, 1, CvType.CV_64F);
		imageSize = new Size(1920, 1080);
		Rot = Mat.zeros(3, 3, CvType.CV_64F);
		T = Mat.ones(3, 1, CvType.CV_64F);

		// CALIBRATION RESULTS FOR 320 x 240
		cameraMatrix.put(0, 0, 287.484405747163, 0, 159.5);
		cameraMatrix.put(1, 0, 0, 287.484405747163, 119.5);
		cameraMatrix.put(2, 0, 0, 0, 1);

		distCoeffs.put(0, 0, 0.1831508618865668);
		distCoeffs.put(1, 0, -0.8391135375141514);
		distCoeffs.put(2, 0, 0);
		distCoeffs.put(3, 0, 0);
		distCoeffs.put(4, 0, 1.067914298622483);

		Rot.put(0, 0, 1, 0, 0);
		Rot.put(1, 0, 0, 1, 0);
		Rot.put(2, 0, 0, 0, 1);

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
}
