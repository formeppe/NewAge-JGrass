/*
 * This file is part of JGrasstools (http://www.jgrasstools.org)
 * (C) HydroloGIS - www.hydrologis.com 
 * 
 * JGrasstools is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */
package org.jgrasstools.hortonmachine.modules.hydrogeomorphology.adige;

import org.jgrasstools.gears.io.rasterreader.OmsRasterReader;
import org.jgrasstools.gears.io.timedependent.OmsTimeSeriesIteratorReader;
import org.jgrasstools.gears.io.timedependent.OmsTimeSeriesIteratorWriter;
import org.jgrasstools.gears.libs.modules.JGTConstants;

import java.awt.image.RenderedImage;
import java.awt.image.WritableRaster;
import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Set;

import javax.media.jai.iterator.RandomIterFactory;
import javax.media.jai.iterator.WritableRandomIter;

import static org.jgrasstools.gears.libs.modules.JGTConstants.doubleNovalue;
import static org.jgrasstools.gears.libs.modules.JGTConstants.isNovalue;
import oms3.annotations.Author;
import oms3.annotations.Bibliography;
import oms3.annotations.Description;
import oms3.annotations.Documentation;
import oms3.annotations.Execute;
import oms3.annotations.In;
import oms3.annotations.Keywords;
import oms3.annotations.Label;
import oms3.annotations.License;
import oms3.annotations.Name;
import oms3.annotations.Out;
import oms3.annotations.Status;
import oms3.annotations.Unit;

import org.geotools.coverage.grid.GridCoverage2D;
import org.geotools.coverage.grid.GridGeometry2D;
import org.geotools.data.simple.SimpleFeatureCollection;
import org.geotools.feature.FeatureIterator;
import org.geotools.feature.SchemaException;
import org.geotools.geometry.DirectPosition2D;
import org.geotools.referencing.crs.DefaultGeographicCRS;
import org.hamcrest.core.IsNull;
import org.jgrasstools.gears.libs.modules.JGTModel;
import org.jgrasstools.gears.libs.monitor.IJGTProgressMonitor;
import org.jgrasstools.gears.libs.monitor.LogProgressMonitor;
import org.jgrasstools.gears.utils.CrsUtilities;
import org.jgrasstools.gears.utils.RegionMap;
import org.jgrasstools.gears.utils.coverage.CoverageUtilities;
import org.jgrasstools.gears.utils.geometry.GeometryUtilities;
import org.joda.time.DateTime;
import org.joda.time.DateTimeZone;
import org.joda.time.format.DateTimeFormat;
import org.joda.time.format.DateTimeFormatter;
import org.opengis.feature.simple.SimpleFeature;
import org.opengis.geometry.DirectPosition;
import org.opengis.geometry.MismatchedDimensionException;
import org.opengis.referencing.crs.CoordinateReferenceSystem;
import org.opengis.referencing.operation.MathTransform;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.Point;

@Description("Calculate the amount of power incident on a surface in a period of time.")
@Documentation("Insolation.html")
@Author(name = "Giuseppe Formetta", contact = "")
@Keywords("Hydrology, Snow Model")
@Label(JGTConstants.HYDROGEOMORPHOLOGY)
@Name("insolation")
@Status(Status.CERTIFIED)
@License("General Public License Version 3 (GPLv3)")
public class SnowMelting extends JGTModel {

	@Description("The map of the skyview factor.")
	@In
	public GridCoverage2D inSkyview = null;

	@Description("The map of the skyview factor.")
	@In
	public double pm1 = 1.0;

	@Description("Parameter of the snowmelting.")
	@In
	public int pDimMeas;

	@Description("The first day of the simulation.")
	@In
	public String tStartDate = null;
	@Description("The time step in minutes of the measurement.")
	@In
	public int inTimestep;
	@Description("The last day of the simulation.")
	@In
	public String tEndDate = null;

	@Description("The map of the dem.")
	@In
	public GridCoverage2D inDem = null;

	@Description("The current time.")
	@In
	public String time = null;

	@Description("The Tmelt.")
	@In
	@Unit("C")
	public double pTmelt;

	@Description("The pr of lmax.")
	@In
	@Unit("C")
	public double pR;

	@Description("The map of the temperature.")
	@In
	public GridCoverage2D inTempGrid = null;
	@In
	public GridCoverage2D incondIniIGrid = null;
	@In
	public GridCoverage2D incondIniLGrid = null;
	@In
	public GridCoverage2D inRainGrid = null;
	@In
	public GridCoverage2D inInsFeb = null;

	@In
	public GridCoverage2D inInsJan = null;

	@In
	public GridCoverage2D inInsMar = null;

	@In
	public GridCoverage2D inInsApr = null;

	@In
	public GridCoverage2D inInsMay = null;

	@In
	public GridCoverage2D inInsJun = null;

	@In
	public SimpleFeatureCollection inStations = null;

	@In
	public String fStationsid = null;

	@In
	public boolean doRaster = true;

	@In
	public HashMap<Integer, double[]> inTemp;

	@In
	public HashMap<Integer, double[]> inRainfall;

	@In
	public double pCmf;
	@Description("Parameter of the snowmelting.")
	@In
	public double pCrf;

	@In
	public double pCff;

	@In
	public boolean doDaily;

	@In
	public double pCr;
	@In
	public double pCs;

	@Description("Parameter of the snowmelting.")
	@In
	public String pathTemp;

	@Description("Parameter of the snowmelting.")
	@In
	public int numSitesWhereCalibrate = 1;

	@Description("Parameter of the snowmelting.")
	@In
	public String pathToSolarRad;

	@Description("Parameter of the snowmelting.")
	@In
	public String pathRainf;

	@Description("0=cazorzi; 1=Classical;2=Hoock")
	@In
	public int pMode;
	@Description("The progress monitor.")
	@In
	public boolean pDoReadMeas = false;

	@In
	public String pathRainfMaps;
	@In
	public String pathTempMaps;

	@Description("The progress monitor.")
	@In
	public IJGTProgressMonitor pm = new LogProgressMonitor();

	@Out
	public GridCoverage2D outGrid = null;

	@In
	public String pathToSwe = null;
	@In
	public String pathToFreezing = null;
	@In
	public String pathToSolidWater = null;
	@In
	public String pathToLiquidWater = null;

	@Description("The map of the insolation for March.")
	@In
	public String pPathtoMeas = null;
	@Description("The interpolated gridded data (for mode 2 and 3.")
	@Out
	public double[] outMeasured = null;

	@Description("The interpolated data (for mode 0 and 1).")
	@In
	public String pathToMelting = null;

	@Description("The interpolated data (for mode 0 and 1).")
	@Out
	public GridCoverage2D outMeltingDataGrid = null;

	@Description("The interpolated data (for mode 0 and 1).")
	@Out
	public GridCoverage2D outSweDataGrid = null;

	@Description("The interpolated data (for mode 0 and 1).")
	@Out
	public double[] outSwevector = null;

	private static final double pLapse = -.0065;

	private HashMap<String, Double> attribute;
	private int height = 0;
	private int width = 0;
	public HashMap<Integer, double[]> outMeltingData = null;
	public HashMap<Integer, double[]> outSWEData = null;
	public HashMap<Integer, double[]> outFreezingData = null;

	private WritableRaster insFebWR = null;
	private WritableRaster insJanWR = null;
	private WritableRaster intempWR = null;
	private WritableRaster inrainWR = null;

	private WritableRaster demwr = null;
	private WritableRaster insMarWR = null;
	private WritableRaster insAprWR = null;
	private WritableRaster insJunWR = null;
	private WritableRaster insMayWR = null;
	public double lambda;
	private double lambdaVect[];
	private double[] xStation;
	private double[] yStation;
	private int[] idStation;
	private int[] colnumvetVect;
	private int[] rownumvetVect;
	private boolean flag = true;
	private WritableRaster skyviewfactorWR = null;
	private HashMap<Integer, double[]> rain_values;
	private HashMap<Integer, double[]> solar_values;
	private HashMap<Integer, double[]> temp_values;
	private HashMap<Integer, double[]> condiniI_values;
	private HashMap<Integer, double[]> condiniL_values;
	boolean init2 = true;
	private double hoursToJan = 720;
	private double hoursToFeb = 1416;
	private double hoursToMar = 2160;
	private double hoursToApr = 2880;
	private double hoursToMay = 3631;
	private double hoursToJune = 4351;
	public HashMap<Integer, double[]> outMelting;
	public HashMap<Integer, double[]> outSWE;
	private WritableRaster outWRSWE = null;
	private WritableRaster outWRMEL = null;
	private double[] outMELTINGVECTOR = null;
	private double[] outSWEVECTOR = null;

	private double[] outMELTINGVECTOR_forMaps = null;
	private double[] outSWEVECTOR_forMaps = null;
	public HashMap<Integer, double[]> inCondIniI;
	public HashMap<Integer, double[]> inCondIniL;

	public File folder = null;
	public File[] P_listOfFiles = null;
	public File[] T_listOfFiles;

	public int conta = 0;
	private LinkedHashMap<Integer, Coordinate> pointsToComputeId2Coordinates = null;
	private int cols;
	private int rows;
	private double south;
	private double west;
	private double xres;
	private double yres;
	private boolean doOneTime = true;

	@Execute
	public void process() throws Exception { // transform the

		if (init2) {
			outMeasured = new double[pDimMeas];
			init2 = false;
			if (pDoReadMeas) {
				// outSimulated=new double[pDimMeas+1];
				// outMeasured = new double[pDimMeas];
				int dim = pDimMeas;
				double[] portate = new double[dim];
				int cont_portate = 0;
				// System.out.println("SONOENTRATO");
				// lettura portate//
				try {

					String str = new String();
					str = pPathtoMeas;

					FileInputStream fstream = new FileInputStream(str);
					DataInputStream in = new DataInputStream(fstream);
					BufferedReader br = new BufferedReader(
							new InputStreamReader(in));
					String strLine;

					double aa = 0;
					while ((strLine = br.readLine()) != null) {
						aa = Double.parseDouble(strLine);
						portate[cont_portate] = aa;
						// System.out.println(aa);
						cont_portate += 1;

					}
					in.close();
				} catch (Exception e) {
					// System.err.println("Errore: " + e.getMessage());
				}

				outMeasured = portate;
				pDoReadMeas = false;

			}
		}

		/*
		 * The models use only one value of the latitude. So I have decided to
		 * set it to the center of the raster. Extract the CRS of the
		 * GridCoverage and transform the value of a WGS84 latitude.
		 */
		CoordinateReferenceSystem sourceCRS = inDem
				.getCoordinateReferenceSystem2D();

		int numPointToCompute = 0;

		GridGeometry2D inRainGridGeo = null;

		if (doRaster == false) {
			pointsToComputeId2Coordinates = getCoordinate(numPointToCompute,
					inStations, fStationsid);
			numPointToCompute = pointsToComputeId2Coordinates.size();

		} else if (doRaster == true) {
			if (inRainGrid != null) {
				inRainGridGeo = inRainGrid.getGridGeometry();
			}
			pointsToComputeId2Coordinates = getCoordinate(inRainGridGeo);
			numPointToCompute = pointsToComputeId2Coordinates.size();
		}
		Set<Integer> pointsToInterpolateIdSet = pointsToComputeId2Coordinates
				.keySet();
		Iterator<Integer> idIterator = pointsToInterpolateIdSet.iterator();
		int j = 0;
		xStation = new double[numPointToCompute];
		yStation = new double[numPointToCompute];
		idStation = new int[numPointToCompute];
		colnumvetVect = new int[numPointToCompute];
		rownumvetVect = new int[numPointToCompute];
		lambdaVect = new double[numPointToCompute];
		CoordinateReferenceSystem targetCRS = DefaultGeographicCRS.WGS84;
		while (idIterator.hasNext()) {
			int id = idIterator.next();
			idStation[j] = id;
			Coordinate coordinate = (Coordinate) pointsToComputeId2Coordinates
					.get(id);
			xStation[j] = coordinate.x;
			yStation[j] = coordinate.y;

			double srcPts[] = new double[] { (xStation[j]), (yStation[j]) };
			Coordinate source = new Coordinate(srcPts[0], srcPts[1]);

			Point[] so = new Point[] { GeometryUtilities.gf().createPoint(
					source) };
			CrsUtilities.reproject(sourceCRS, targetCRS, so);
			// the latitude value
			lambdaVect[j] = Math.toRadians(so[0].getY());

			j += 1;

		}
		MathTransform transf = inDem.getGridGeometry().getCRSToGrid2D();
		for (int i = 0; i < xStation.length; i++) {

			DirectPosition point = new DirectPosition2D(sourceCRS, xStation[i],
					yStation[i]);
			DirectPosition gridPoint = transf.transform(point, null);

			colnumvetVect[i] = (int) gridPoint.getCoordinate()[0];
			rownumvetVect[i] = (int) gridPoint.getCoordinate()[1];
			// System.out.println(idStation[i] + " "
			// + gridPoint.getCoordinate()[0] + " "
			// + gridPoint.getCoordinate()[1]);
		}

		double minimofeb = 0;
		double minimomar = 0;
		double minimoapr = 0;
		double minimomagg = 0;
		double minimogiu = 0;
		double dx = 0;
		DateTimeFormatter formatter = DateTimeFormat.forPattern(
				"yyyy-MM-dd HH:mm").withZone(DateTimeZone.UTC);
		DateTime startcurrentDatetime = formatter.parseDateTime(tStartDate);

		DateTime endcurrentDatetime = formatter.parseDateTime(tEndDate);
		long diff = 0;
		if (doDaily == false) {
			diff = (endcurrentDatetime.getMillis() - startcurrentDatetime
					.getMillis()) / (1000 * 60 * 60);
		}
		if (doDaily) {
			diff = (endcurrentDatetime.getMillis() - startcurrentDatetime
					.getMillis()) / (1000 * 60 * 60 * 24);
		}
		DateTime array[] = new DateTime[(int) diff];
		if (doDaily == false) {
			for (int i = 0; i < array.length; i++) {
				array[i] = startcurrentDatetime.plusHours(i);
			}
		}
		if (doDaily) {
			for (int i = 0; i < array.length; i++) {
				array[i] = startcurrentDatetime.plusDays(i);
			}
		}
		if (doOneTime) {
			attribute = CoverageUtilities
					.getRegionParamsFromGridCoverage(inSkyview);
			dx = attribute.get(CoverageUtilities.XRES);

			double srcPts[] = new double[] {
					attribute.get(CoverageUtilities.EAST),
					attribute.get(CoverageUtilities.SOUTH) };
			// CoordinateReferenceSystem targetCRS = DefaultGeographicCRS.WGS84;
			Coordinate source = new Coordinate(srcPts[0], srcPts[1]);
			Point[] so = new Point[] { GeometryUtilities.gf().createPoint(
					source) };
			CrsUtilities.reproject(sourceCRS, targetCRS, so);
			// the latitude value
			lambda = Math.toRadians(so[0].getY());
			CoverageUtilities.getRegionParamsFromGridCoverage(inSkyview);

			if (pMode == 0) {

				RenderedImage insFebTmpRI = inInsFeb.getRenderedImage();
				width = insFebTmpRI.getWidth();
				height = insFebTmpRI.getHeight();
				insFebWR = CoverageUtilities.replaceNovalue(insFebTmpRI,
						-9999.0);

				RenderedImage insJanTmpRI = inInsJan.getRenderedImage();
				width = insJanTmpRI.getWidth();
				height = insJanTmpRI.getHeight();
				insJanWR = CoverageUtilities.replaceNovalue(insJanTmpRI,
						-9999.0);

				RenderedImage insMarTmpRI = inInsMar.getRenderedImage();
				insMarWR = CoverageUtilities.replaceNovalue(insMarTmpRI,
						-9999.0);
				insMarTmpRI = null;

				RenderedImage insAprTmpRI = inInsApr.getRenderedImage();
				insAprWR = CoverageUtilities.replaceNovalue(insAprTmpRI,
						-9999.0);
				insAprTmpRI = null;

				RenderedImage insMayTmpRI = inInsMay.getRenderedImage();
				insMayWR = CoverageUtilities.replaceNovalue(insMayTmpRI,
						-9999.0);
				insMayTmpRI = null;

				RenderedImage insJuneTmpRI = inInsJun.getRenderedImage();
				insJunWR = CoverageUtilities.replaceNovalue(insJuneTmpRI,
						-9999.0);
				insJuneTmpRI = null;
				minimofeb = findminumum(insFebWR, hoursToFeb);
				minimomar = findminumum(insMarWR, hoursToMar);
				minimoapr = findminumum(insAprWR, hoursToApr);
				minimomagg = findminumum(insMayWR, hoursToMay);
				minimogiu = findminumum(insJunWR, hoursToJune);

			}
			RenderedImage drmri = inDem.getRenderedImage();
			demwr = CoverageUtilities.replaceNovalue(drmri, -9999.0);
			drmri = null;
			RenderedImage SkyTmpRI = inSkyview.getRenderedImage();
			skyviewfactorWR = CoverageUtilities.renderedImage2WritableRaster(
					SkyTmpRI, true);

			inCondIniL = new HashMap<Integer, double[]>();
			inCondIniI = new HashMap<Integer, double[]>();
			outSWEData = new HashMap<Integer, double[]>();
			outFreezingData = new HashMap<Integer, double[]>();

			doOneTime = false;
		}
		folder = null;
		P_listOfFiles = null;
		T_listOfFiles = null;
		if (doRaster == true) {

			folder = new File(pathRainfMaps);
			P_listOfFiles = folder.listFiles();
			for (int i = 0; i < P_listOfFiles.length; i++) {
				System.out.println(P_listOfFiles[i]);
			}
			folder = new File(pathTempMaps);
			T_listOfFiles = folder.listFiles();
			for (int i = 0; i < T_listOfFiles.length; i++) {
				System.out.println(T_listOfFiles[i]);
			}

		}
		OmsTimeSeriesIteratorReader reader_temp = new OmsTimeSeriesIteratorReader();
		OmsTimeSeriesIteratorReader reader_solar = new OmsTimeSeriesIteratorReader();

		OmsTimeSeriesIteratorReader reader_rainf = new OmsTimeSeriesIteratorReader();
		OmsTimeSeriesIteratorWriter writer = new OmsTimeSeriesIteratorWriter();
		OmsTimeSeriesIteratorWriter writer2 = new OmsTimeSeriesIteratorWriter();
		OmsTimeSeriesIteratorWriter writer3 = new OmsTimeSeriesIteratorWriter();
		OmsTimeSeriesIteratorWriter writer4 = new OmsTimeSeriesIteratorWriter();
		OmsTimeSeriesIteratorWriter writer5 = new OmsTimeSeriesIteratorWriter();

		if (doRaster == false) {

			if (!((pathTemp == null))) {
				reader_temp.file = pathTemp;
				reader_temp.idfield = "ID";
				reader_temp.tStart = tStartDate;
				reader_temp.tEnd = tEndDate;
				reader_temp.fileNovalue = "-9999";
				reader_temp.tTimestep = inTimestep;
			}

			if (!((pathToSolarRad == null))) {
				reader_solar.file = pathToSolarRad;
				reader_solar.idfield = "ID";
				reader_solar.tStart = tStartDate;
				reader_solar.tEnd = tEndDate;
				reader_solar.fileNovalue = "-9999";
				reader_solar.tTimestep = inTimestep;
			}

			if (!(pathRainf == null)) {

				reader_rainf.file = pathRainf;
				reader_rainf.idfield = "ID";
				reader_rainf.tStart = tStartDate;
				reader_rainf.tEnd = tEndDate;
				reader_rainf.fileNovalue = "-9999";
				reader_rainf.tTimestep = inTimestep;

			}

			if (pathToMelting != null) {
				writer.file = pathToMelting;
				writer.tStart = tStartDate;
				writer.tTimestep = inTimestep;
			}
			if (pathToSwe != null) {
				writer2.file = pathToSwe;
				writer2.tStart = tStartDate;
				writer2.tTimestep = inTimestep;

			}
			if (pathToFreezing != null) {
				writer3.file = pathToFreezing;
				writer3.tStart = tStartDate;
				writer3.tTimestep = inTimestep;

			}
			if (pathToSolidWater != null) {
				writer4.file = pathToSolidWater;
				writer4.tStart = tStartDate;
				writer4.tTimestep = inTimestep;
			}
			if (pathToLiquidWater != null) {
				writer5.file = pathToLiquidWater;
				writer5.tStart = tStartDate;
				writer5.tTimestep = inTimestep;

			}

			if (numSitesWhereCalibrate != 1) {
				outSWEVECTOR = new double[array.length * numSitesWhereCalibrate];
				outMELTINGVECTOR = new double[array.length
						* numSitesWhereCalibrate];
			} else {
				outSWEVECTOR = new double[array.length];
				outMELTINGVECTOR = new double[array.length];

			}
			conta = 0;
		}

		outMELTINGVECTOR_forMaps = new double[numPointToCompute];
		outSWEVECTOR_forMaps = new double[numPointToCompute];

		for (int i = 0; i < array.length; i++) {
			outSWE = new HashMap<Integer, double[]>();
			outMeltingData = new HashMap<Integer, double[]>();
			// System.out.println(" data=" + array[i]);
			DateTime currentime = array[i];
			if (doRaster == false) {
				if (!(pathTemp == null)) {
					reader_temp.nextRecord();
					temp_values = reader_temp.outData;
				}

				if (!(pathRainf == null)) {
					reader_rainf.nextRecord();
					rain_values = reader_rainf.outData;
				}
				if (!(pathToSolarRad == null)) {
					reader_solar.nextRecord();
					solar_values = reader_solar.outData;
				}
			}
			int month = currentime.getMonthOfYear();
			switch (month) {

			case 1:
				calcMelting(insJanWR, demwr, dx, currentime, hoursToJan,
						minimofeb, i);
				break;
			case 10:
				calcMelting(insFebWR, demwr, dx, currentime, hoursToFeb,
						minimofeb, i);
				break;
			case 11:
				calcMelting(insFebWR, demwr, dx, currentime, hoursToFeb,
						minimofeb, i);
				break;
			case 12:
				calcMelting(insFebWR, demwr, dx, currentime, hoursToFeb,
						minimofeb, i);
				break;
			case 2:
				calcMelting(insFebWR, demwr, dx, currentime, hoursToFeb,
						minimofeb, i);
				break;
			case 3:
				calcMelting(insMarWR, demwr, dx, currentime, hoursToMar,
						minimomar, i);
				break;
			case 4:
				calcMelting(insAprWR, demwr, dx, currentime, hoursToApr,
						minimoapr, i);
				break;
			case 5:
				calcMelting(insMayWR, demwr, dx, currentime, hoursToMay,
						minimomagg, i);
				break;
			case 6:
				calcMelting(insJunWR, demwr, dx, currentime, hoursToJune,
						minimogiu, i);
				break;
			case 7:
				calcMelting(insJunWR, demwr, dx, currentime, hoursToJune,
						minimogiu, i);
				break;
			case 8:
				calcMelting(insJunWR, demwr, dx, currentime, hoursToJune,
						minimogiu, i);
				break;
			case 9:
				calcMelting(insJunWR, demwr, dx, currentime, hoursToJune,
						minimogiu, i);
				break;
			default:
				break;
			}

			//
			// calcInsolation(lambda, pitWR, gradientWR, insolationWR,
			// staoWR,
			// diffuseWR, dx, currentime);
			if (doRaster == false) {
				if (pathToMelting != null) {
					writer.inData = outMeltingData;

					writer.writeNextLine();
				}
				if (pathToSolidWater != null) {
					writer4.inData = inCondIniI;

					writer4.writeNextLine();
				}
				if (pathToLiquidWater != null) {
					writer5.inData = inCondIniL;

					writer5.writeNextLine();
				}
				if (pathToSwe != null) {
					writer2.inData = outSWEData;
					writer2.writeNextLine();
				}
				if (pathToFreezing != null) {
					writer3.inData = outFreezingData;
					writer3.writeNextLine();
				}
			} else if (doRaster) {

				storeResult(outSWEVECTOR_forMaps, outMELTINGVECTOR_forMaps,
						pointsToComputeId2Coordinates);
				// }

			}
			// pm.worked(i - startDay);
		}
		if (pathToMelting != null) {
			writer.close();
		}
		if (pathToSwe != null) {
			writer2.close();
		}
		if (pathToFreezing != null) {
			writer3.close();
		}
		if (pathToSolidWater != null) {
			writer4.close();
		}
		if (pathToLiquidWater != null) {
			writer5.close();
		}
		outSwevector = outSWEVECTOR;

	}

	// System.out.print(time + "   ");

	// DateTimeFormatter formatter = DateTimeFormat.forPattern(
	// "yyyy-MM-dd HH:mm").withZone(DateTimeZone.UTC);
	// DateTime currentDatetime = formatter.parseDateTime(time);
	// int month = currentDatetime.getMonthOfYear();
	// switch (month) {
	//
	// case 1:
	// calcInsolation(insFebWR, demwr, dx, currentDatetime, hoursToFeb,
	// minimofeb);
	// break;
	// case 10:
	// calcInsolation(insFebWR, demwr, dx, currentDatetime, hoursToFeb,
	// minimofeb);
	// break;
	// case 11:
	// calcInsolation(insFebWR, demwr, dx, currentDatetime, hoursToFeb,
	// minimofeb);
	// break;
	// case 12:
	// calcInsolation(insFebWR, demwr, dx, currentDatetime, hoursToFeb,
	// minimofeb);
	// break;
	// case 2:
	// calcInsolation(insFebWR, demwr, dx, currentDatetime, hoursToFeb,
	// minimofeb);
	// break;
	// case 3:
	// calcInsolation(insMarWR, demwr, dx, currentDatetime, hoursToMar,
	// minimomar);
	// break;
	// case 4:
	// calcInsolation(insAprWR, demwr, dx, currentDatetime, hoursToApr,
	// minimoapr);
	// break;
	// case 5:
	// calcInsolation(insMayWR, demwr, dx, currentDatetime, hoursToMay,
	// minimomagg);
	// break;
	// case 6:
	// calcInsolation(insJunWR, demwr, dx, currentDatetime, hoursToJune,
	// minimogiu);
	// break;
	// case 7:
	// calcInsolation(insJunWR, demwr, dx, currentDatetime, hoursToJune,
	// minimogiu);
	// break;
	// case 8:
	// calcInsolation(insJunWR, demwr, dx, currentDatetime, hoursToJune,
	// minimogiu);
	// break;
	// case 9:
	// calcInsolation(insJunWR, demwr, dx, currentDatetime, hoursToJune,
	// minimogiu);
	// break;
	// default:
	// break;
	// }
	// if (doRaster == false) {
	// storeResult(outSWEVECTOR, outMELTINGVECTOR, idStation);
	// } else {
	// storeResult(outSWEVECTOR, outMELTINGVECTOR,
	// pointsToComputeId2Coordinates);
	// }

	private void calcMelting(WritableRaster energyWR, WritableRaster demWR,
			double dx, DateTime time, double ore, double minimoEI, int iii)
			throws Exception {

		CoverageUtilities.getRegionParamsFromGridCoverage(inSkyview);

		double condiniI = doubleNovalue;
		double condiniL = doubleNovalue;

		for (int contastaz = 0; contastaz < xStation.length; contastaz++) {

			int colnuber = colnumvetVect[contastaz];
			int rownumber = rownumvetVect[contastaz];
			int i = colnuber;
			int j = rownumber;
			int id = idStation[contastaz];
			// System.out.println("STATIOn ID. ========> " + id);

			double temperatura = doubleNovalue;
			double pioggia = doubleNovalue;
			double solar = doubleNovalue;

			if (doRaster == false) {
				if (temp_values != null) {
					temperatura = temp_values.get(id)[0];
				}
				if (rain_values != null) {
					pioggia = rain_values.get(id)[0];
				}
				if (solar_values != null) {
					solar = solar_values.get(id)[0];
				}
				if (flag) {

					if (!(condiniI_values != null)) {
						condiniI = 0;
					}
					if (!(condiniL_values != null)) {
						condiniL = 0;
					}
					if (contastaz == xStation.length - 1) {
						flag = false;
					}
				} else {
					condiniI = inCondIniI.get(id)[0];
					condiniL = inCondIniL.get(id)[0];

				}
			} else if (doRaster == true) {

				// String fileT=pathTempMaps+iii+".asc";
				// String fileP=pathRainfMaps+iii+".asc";

				OmsRasterReader readersT = new OmsRasterReader();
				readersT.file = T_listOfFiles[iii].toString();
				readersT.fileNovalue = -9999.0;
				readersT.geodataNovalue = Double.NaN;
				readersT.process();
				GridCoverage2D t = readersT.outRaster;

				OmsRasterReader readersP = new OmsRasterReader();
				readersP.file = P_listOfFiles[iii].toString();
				readersP.fileNovalue = -9999.0;
				readersP.geodataNovalue = Double.NaN;
				readersP.process();
				GridCoverage2D p = readersP.outRaster;

				RenderedImage inRainRI = p.getRenderedImage();
				inrainWR = CoverageUtilities.replaceNovalue(inRainRI, -9999.0);
				inRainRI = null;

				RenderedImage inTempRI = t.getRenderedImage();
				intempWR = CoverageUtilities.replaceNovalue(inTempRI, -9999.0);
				inTempRI = null;

				if (!(temp_values != null)) {
					temperatura = intempWR.getSampleDouble(i, j, 0);
				}
				if (!(rain_values != null)) {
					pioggia = inrainWR.getSampleDouble(i, j, 0);
				}
				if (flag) {

					if (!(condiniI_values != null)) {
						condiniI = 0;
					}
					if (!(condiniL_values != null)) {
						condiniL = 0;
					}
					if (contastaz == xStation.length - 1) {
						flag = false;
					}
				} else {
					condiniI = inCondIniI.get(id)[0];
					condiniL = inCondIniL.get(id)[0];

				}
			}
			double[] aaa = calcMelting(i, j, energyWR, demWR, temperatura,
					pioggia, solar, ore, time, minimoEI, condiniI, condiniL,
					lambdaVect[contastaz]);

			outSWEVECTOR[conta] = aaa[0];
			outSWEVECTOR_forMaps[contastaz] = aaa[0];

			inCondIniI.put(id, new double[] { aaa[1] });
			inCondIniL.put(id, new double[] { aaa[2] });

			outMeltingData.put(id, new double[] { aaa[3] });
			outSWEData.put(id, new double[] { aaa[0] });
			outFreezingData.put(id, new double[] { aaa[4] });

			outMELTINGVECTOR[conta] = aaa[3];
			outMELTINGVECTOR_forMaps[contastaz] = aaa[3];

			conta++;
		}

	}

	private double[] calcMelting(int i, int j, WritableRaster enWR,
			WritableRaster demWR, double tem, double pio, double sol,
			double ore, DateTime time, double minimoEI, double condiniI,
			double condiniL, double lll) throws IOException {

		double[] risultato = new double[5];
		double melting = 0;
		double prain = 0;
		double psnow = 0;

		int day = time.getDayOfYear();

		double dayangb = (360 / 365.25) * (day - 79.436);

		dayangb = Math.toRadians(dayangb);

		// Evaluate the declination of the sun.
		double delta = getDeclination(dayangb);
		// Evaluate the radiation in this day.
		double ss = Math.acos(-Math.tan(delta) * Math.tan(lll));

		// double ssDaniele = Math.acos(-Math.tan(delta) * Math.tan(lambda));
		double sunrise = 12 * (1.0 - ss / Math.PI);
		double sunset = 12 * (1.0 + ss / Math.PI);
		double hhh = (double) time.getMillisOfDay() / (1000 * 60 * 60);
		// System.out.println(time+"  day=" + day + "  delta=  " +
		// delta+"  sunrise= "+sunrise+"  sunset="+sunset);
		if (doDaily == false) {

			if ((hhh >= (sunrise) && hhh <= (sunset))) {
				if (skyviewfactorWR.getSampleDouble(i, j, 0) != -9999.0) {

					double ei = 0;
					if (pMode == 0) {
						ei = enWR.getSampleDouble(i, j, 0) / (ore / 24);
					}

					double temp = tem;

					if (isNovalue(tem)) {
						double zzz = demwr.getSampleDouble(i, j, 0);
						temp = 273 + pLapse * (zzz - 4000);
					}
					if (isNovalue(pio) || pio < 0) {
						pio = 0;
					}
					if (pMode == 0) {

						if (temp > pTmelt) {

							melting = ei * pCmf * (temp - pTmelt)
									* skyviewfactorWR.getSampleDouble(i, j, 0);
							melting = Math.min(melting, (condiniI + condiniL));
						} else {
							melting = 0;
						}

					} else if (pMode == 1) {
						if (temp > pTmelt) {

							melting = pCmf * (temp - pTmelt)
									* skyviewfactorWR.getSampleDouble(i, j, 0);
						} else {
							melting = 0;
						}
						melting = Math.max(melting, 0);

					} else if (pMode == 2) {
						if (temp > pTmelt) {
							melting = (pCmf + sol * pCrf) * (temp - pTmelt)
									* skyviewfactorWR.getSampleDouble(i, j, 0);
						} else {
							melting = 0;
						}
						melting = Math.max(melting, 0);

					}

					else {
						System.out.println("ERROREEEE");
					}

					prain = (pio / Math.PI) * Math.atan((temp - pTmelt) / pm1)
							+ pio * 0.5;

					psnow = pio - prain;
					prain = prain * pCr;
					psnow = pCs * psnow;

					// freezing
					double freezing = 0;

					if (temp < pTmelt) {
						freezing = pCff * (pTmelt - temp);
					} else {
						freezing = 0;
					}

					double I = 0;
					double deltat = 1.0;
					// if (temp > pTmelt) {
					I = condiniI + deltat * (psnow + freezing - melting);
					// } else {
					// I = condiniI + deltat * (psnow + freezing);
					// }
					if (I < 0) {
						I = 0;
						melting = 0;
					}

					double L = 0;
					L = condiniL + deltat * (prain + melting - freezing);
					double L_max = pR * I;
					double melting_discharge = 0;
					if (L < 0)
						L = 0.0;
					if (L > L_max) {
						melting_discharge = L - L_max;
						L = L_max;
					}
					risultato[0] = L + I;
					risultato[1] = I;
					risultato[2] = L;
					risultato[3] = melting_discharge;
					risultato[4] = freezing;

				}

				else {

					risultato[0] = -9999;
					risultato[1] = condiniI;
					risultato[2] = condiniL;
					risultato[3] = -9999;
					risultato[4] = -9999;

				}
			} else {

				if (skyviewfactorWR.getSampleDouble(i, j, 0) != -9999.0) {

					double z = 0;
					if (pMode == 0) {
						z = minimoEI / ore;
					}

					double temp = tem;

					if (isNovalue(tem)) {
						double zzz = demWR.getSampleDouble(i, j, 0);
						temp = pLapse * (zzz - 4000);
					}

					if (pMode == 0) {
						if (temp > pTmelt) {
							melting = z * pCmf * (temp - pTmelt)
									* skyviewfactorWR.getSampleDouble(i, j, 0);
							melting = Math.min(melting, (condiniI + condiniL));

						} else {
							melting = 0;
						}
						// melting = Math.max(melting, 0);

					} else if (pMode == 1) {
						if (temp > pTmelt) {
							melting = pCmf * (temp - pTmelt)
									* skyviewfactorWR.getSampleDouble(i, j, 0);
						} else {
							melting = 0;
						}
						melting = Math.max(melting, 0);

					}

					else if (pMode == 2) {
						if (temp > pTmelt) {
							melting = (pCmf + sol * pCrf) * (temp - pTmelt)
									* skyviewfactorWR.getSampleDouble(i, j, 0);
						} else {
							melting = 0;
						}
						melting = Math.max(melting, 0);

					} else {
						System.out.println("ERROREEEE");
					}

					if (isNovalue(pio) || pio < 0) {
						pio = 0;
					}

					prain = (pio / Math.PI) * Math.atan((temp - pTmelt) / pm1)
							+ pio * 0.5;

					psnow = pio - prain;
					prain = prain * pCr;
					psnow = pCs * psnow;

					// freezing
					double freezing = 0;

					if (temp < pTmelt) {
						freezing = pCff * (pTmelt - temp);

					} else {
						freezing = 0;
					}

					double I = 0;
					double deltat = 1.0;
					// if (temp > pTmelt) {
					I = condiniI + deltat * (psnow + freezing - melting);
					// } else {
					// I = condiniI + deltat * (psnow + freezing);
					// }
					if (I < 0) {
						I = 0;
						melting = 0;
					}

					double L = 0;
					L = condiniL + deltat * (prain + melting - freezing);
					double L_max = pR * I;
					double melting_discharge = 0;
					if (L < 0)
						L = 0.0;
					if (L > L_max) {
						melting_discharge = L - L_max;
						L = L_max;
					}
					risultato[0] = L + I;
					risultato[1] = I;
					risultato[2] = L;
					risultato[3] = melting_discharge;
					risultato[4] = freezing;

				}

				else {

					risultato[0] = -9999;
					risultato[1] = condiniI;
					risultato[2] = condiniL;
					risultato[3] = -9999;
					risultato[4] = -9999;

				}

			}

		} else {

			if (skyviewfactorWR.getSampleDouble(i, j, 0) != -9999.0) {
				double ei = 0;
				if (pMode == 0) {
					ei = enWR.getSampleDouble(i, j, 0) / (ore / 24);
				}
				double temp = tem;

				if (isNovalue(tem)) {
					double zzz = demWR.getSampleDouble(i, j, 0);
					temp = 273 + pLapse * (zzz - 4000);
				}
				if (isNovalue(pio)) {
					pio = 0;
				}
				if (pMode == 0) {
					if ((temp > pTmelt)) {
						melting = ei * pCmf * (temp - pTmelt)
								* skyviewfactorWR.getSampleDouble(i, j, 0);
					} else {
						melting = 0;
					}
					melting = Math.max(melting, 0);

				} else if (pMode == 1) {
					if ((temp > pTmelt)) {
						melting = pCmf * (temp - pTmelt)
								* skyviewfactorWR.getSampleDouble(i, j, 0);
					} else {
						melting = 0;
					}
					melting = Math.max(melting, 0);

				}

				else if (pMode == 2) {
					if (temp > pTmelt) {
						melting = (pCmf + sol * pCrf) * (temp - pTmelt)
								* skyviewfactorWR.getSampleDouble(i, j, 0);
					} else {
						melting = 0;
					}
					melting = Math.max(melting, 0);

				}

				else {
					System.out.println("ERROREEEE");
				}

				prain = (pio / Math.PI) * Math.atan((temp - pTmelt) / pm1)
						+ pio * 0.5;

				psnow = pio - prain;
				prain = prain * pCr;
				psnow = pCs * psnow;

				// freezing
				double freezing = 0;

				if (temp < pTmelt) {
					freezing = pCff * (pTmelt - temp);
				} else {
					freezing = 0;
				}

				double I = 0;
				double deltat = 1.0;
				// if (temp > pTmelt) {
				I = condiniI + deltat * (psnow + freezing - melting);
				// } else {
				// I = condiniI + deltat * (psnow + freezing);
				// }
				if (I < 0) {
					I = 0;
					melting = 0;
				}

				double L = 0;
				L = condiniL + deltat * (prain + melting - freezing);
				double L_max = pR * I;
				double melting_discharge = 0;
				if (L < 0)
					L = 0.0;
				if (L > L_max) {
					melting_discharge = L - L_max;
					L = L_max;
				}
				risultato[0] = L + I;
				risultato[1] = I;
				risultato[2] = L;
				risultato[3] = melting_discharge;
				risultato[4] = freezing;

			}

			else {

				risultato[0] = -9999;
				risultato[1] = condiniI;
				risultato[2] = condiniL;
				risultato[3] = -9999;
				risultato[4] = -9999;

			}

		}

		return risultato;

	}

	private double findminumum(WritableRaster w, double num) {

		double minimo = 10000000;
		for (int j = 0; j < height; j++) {
			for (int i = 0; i < width; i++) {
				if (w.getSample(i, j, 0) != -9999) {
					if (w.getSample(i, j, 0) / num > 0.0
							&& w.getSample(i, j, 0) / num < minimo) {
						minimo = w.getSample(i, j, 0) / num;

					}
				}
			}
		}

		return minimo;

	}

	private LinkedHashMap<Integer, Coordinate> getCoordinate(GridGeometry2D grid) {
		LinkedHashMap<Integer, Coordinate> out = new LinkedHashMap<Integer, Coordinate>();
		int count = 0;
		RegionMap regionMap = CoverageUtilities
				.gridGeometry2RegionParamsMap(grid);
		cols = regionMap.getCols();
		rows = regionMap.getRows();
		south = regionMap.getSouth();
		west = regionMap.getWest();
		xres = regionMap.getXres();
		yres = regionMap.getYres();

		outWRSWE = CoverageUtilities.createDoubleWritableRaster(cols, rows,
				null, null, null);
		outWRMEL = CoverageUtilities.createDoubleWritableRaster(cols, rows,
				null, null, null);

		double northing = south;
		double easting = west;
		for (int i = 0; i < cols; i++) {
			easting = easting + xres;
			for (int j = 0; j < rows; j++) {
				northing = northing + yres;
				Coordinate coordinate = new Coordinate();
				coordinate.x = west + i * xres;
				coordinate.y = south + j * yres;
				out.put(count, coordinate);
				count++;
			}
		}

		return out;
	}

	/**
	 * Extract the coordinate of a FeatureCollection in a HashMap with an ID as
	 * a key.
	 * 
	 * @param nStaz
	 * @param collection
	 * @throws Exception
	 *             if a fiel of elevation isn't the same of the collection
	 */
	private LinkedHashMap<Integer, Coordinate> getCoordinate(int nStaz,
			SimpleFeatureCollection collection, String idField)
			throws Exception {
		LinkedHashMap<Integer, Coordinate> id2CoordinatesMap = new LinkedHashMap<Integer, Coordinate>();
		FeatureIterator<SimpleFeature> iterator = collection.features();
		Coordinate coordinate = null;
		try {
			while (iterator.hasNext()) {
				SimpleFeature feature = iterator.next();
				int name = ((Number) feature.getAttribute(idField)).intValue();
				coordinate = ((Geometry) feature.getDefaultGeometry())
						.getCentroid().getCoordinate();
				double z = 0;

				coordinate.z = z;
				id2CoordinatesMap.put(name, coordinate);
			}
		} finally {
			iterator.close();
		}

		return id2CoordinatesMap;
	}

	// private void storeResult(double[] resultSwe, double[] resultMel, int[]
	// id)
	// throws SchemaException {
	// outSWEData = new HashMap<Integer, double[]>();
	// outMeltingData = new HashMap<Integer, double[]>();
	// for (int i = 0; i < resultSwe.length; i++) {
	// outSWEData.put(id[i], new double[] { resultSwe[i] });
	// outMeltingData.put(id[i], new double[] { resultMel[i] });
	// }
	// }

	private void storeResult(double[] resultSwe, double[] resultMel,
			HashMap<Integer, Coordinate> interpolatedCoordinatesMap)
			throws MismatchedDimensionException, Exception {

		WritableRandomIter outIterSWE = RandomIterFactory.createWritable(
				outWRSWE, null);
		WritableRandomIter outIterMEL = RandomIterFactory.createWritable(
				outWRMEL, null);

		Set<Integer> pointsToInterpolateIdSett = interpolatedCoordinatesMap
				.keySet();
		Iterator<Integer> idIterator = pointsToInterpolateIdSett.iterator();
		int c = 0;
		MathTransform transf = inRainGrid.getGridGeometry().getCRSToGrid2D();

		final DirectPosition gridPoint = new DirectPosition2D();

		while (idIterator.hasNext()) {
			int id = idIterator.next();
			Coordinate coordinate = (Coordinate) interpolatedCoordinatesMap
					.get(id);

			DirectPosition point = new DirectPosition2D(
					inRainGrid.getCoordinateReferenceSystem(), coordinate.x,
					coordinate.y);
			transf.transform(point, gridPoint);

			double[] gridCoord = gridPoint.getCoordinate();
			int x = (int) gridCoord[0];
			int y = (int) gridCoord[1];

			outIterSWE.setSample(x, y, 0, resultSwe[c]);
			outIterMEL.setSample(x, y, 0, resultMel[c]);
			c++;

		}

		RegionMap regionMap = CoverageUtilities
				.gridGeometry2RegionParamsMap(inRainGrid.getGridGeometry());

		outMeltingDataGrid = CoverageUtilities.buildCoverage("gridded",
				outWRMEL, regionMap, inRainGrid.getGridGeometry()
						.getCoordinateReferenceSystem());
		outSweDataGrid = CoverageUtilities.buildCoverage("gridded", outWRSWE,
				regionMap, inRainGrid.getGridGeometry()
						.getCoordinateReferenceSystem());

	}

	private double getDeclination(double dayangb) {
		double delta = .3723 + 23.2567 * Math.sin(dayangb) - .758
				* Math.cos(dayangb) + .1149 * Math.sin(2 * dayangb) + .3656
				* Math.cos(2 * dayangb) - .1712 * Math.sin(3 * dayangb) + .0201
				* Math.cos(3 * dayangb);
		return Math.toRadians(delta);
	}

}
