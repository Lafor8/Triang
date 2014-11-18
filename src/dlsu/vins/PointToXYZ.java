package dlsu.vins;

public class PointToXYZ {
	private String printPoint(){
		StringBuilder sb = new StringBuilder();
		sb.append(0 + " ");
		sb.append(0 + " ");
		sb.append(0 + " ");
		sb.append(255 + " ");
		sb.append(255 + " ");
		sb.append(255 + "\n");
		
		return sb.toString();
	} 
}

/* 
usage:
java -jar PotreeConverter.jar "path to xyz" "path to output directory" [max points]

ex.:
java -jar PotreeConverter.jar "C:/dev/pointclouds/skatepark/skate/skate.xyz" "D:/temp/potree/skatepark/converterTest/02" 40000000

*/
