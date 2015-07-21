/*
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package org.jgrasstools.hortonmachine.modules.hydrogeomorphology.lwrb;

import static org.jgrasstools.gears.libs.modules.ModelsEngine.calcInverseSunVector;
import static org.jgrasstools.gears.libs.modules.ModelsEngine.calcNormalSunVector;
import static org.jgrasstools.gears.libs.modules.ModelsEngine.calculateFactor;

import org.jgrasstools.gears.io.rasterwriter.OmsRasterWriter;
import org.jgrasstools.gears.io.timedependent.OmsTimeSeriesIteratorReader;
import org.jgrasstools.gears.io.timedependent.OmsTimeSeriesIteratorWriter;
import org.jgrasstools.gears.libs.modules.JGTConstants;

import java.awt.image.RenderedImage;
import java.awt.image.SampleModel;
import java.awt.image.WritableRaster;
import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import static org.jgrasstools.gears.libs.modules.JGTConstants.doubleNovalue;
import static org.jgrasstools.gears.libs.modules.JGTConstants.isNovalue;

import javax.media.jai.RasterFactory;
import javax.media.jai.iterator.RandomIter;
import javax.media.jai.iterator.RandomIterFactory;
import javax.media.jai.iterator.WritableRandomIter;

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
import org.geotools.data.simple.SimpleFeatureCollection;
import org.geotools.feature.FeatureIterator;
import org.geotools.geometry.DirectPosition2D;
import org.geotools.referencing.crs.DefaultGeographicCRS;
import org.jgrasstools.gears.libs.modules.JGTModel;
import org.jgrasstools.gears.libs.monitor.IJGTProgressMonitor;
import org.jgrasstools.gears.libs.monitor.LogProgressMonitor;
import org.jgrasstools.gears.utils.CrsUtilities;
import org.jgrasstools.gears.utils.coverage.CoverageUtilities;
import org.jgrasstools.gears.utils.geometry.GeometryUtilities;
import org.joda.time.DateTime;
import org.joda.time.DateTimeZone;
import org.joda.time.format.DateTimeFormat;
import org.joda.time.format.DateTimeFormatter;
import org.opengis.feature.simple.SimpleFeature;
import org.opengis.geometry.DirectPosition;
import org.opengis.referencing.crs.CoordinateReferenceSystem;
import org.opengis.referencing.operation.MathTransform;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.Point;

@Description("Calculate longwave upwelling and downwelling solar radiation.")
@Documentation("LWRB.html")
@Author(name = "Giuseppe Formetta and Riccardo Rigon", contact = "giuseppe.formetta@unical.it")
@Keywords("Hydrology, Radiation, SkyviewFactor, Hillshade")
@Label(JGTConstants.HYDROGEOMORPHOLOGY)
@Name("lwrb")
@Status(Status.CERTIFIED)
@License("General Public License Version 3 (GPLv3)")
public class OmsLongwaveRadiationBalance extends JGTModel {
	@Description("The map of the elevation.")
	@In
	public GridCoverage2D inElev = null;

	@Description("The map of the elevation.")
	@In
	public double pYs = -999;

	@Description("The map of the elevation.")
	@In
	public double pZs = -999;

	@Description("The map of the elevation.")
	@In
	public double pXs = -999;

	@Description("The map of the skyview factor.")
	@In
	public GridCoverage2D inskyview = null;

	@Description("The vector of the measurement points, containing the position of the stations.")
	@In
	public SimpleFeatureCollection inStations = null;

	@Description("The field of the vector of stations, defining the id.")
	@In
	public String fStationsid = null;

	@Description("The field of the vector of stations, defining the albedo.")
	@In
	public String fStationsAlb = null;

	@Description("The first day of the simulation.")
	@In
	public String tStartDate = null;

	@Description("The last day of the simulation.")
	@In
	public String tEndDate = null;
	@Description("The time step in minutes of the measurement.")
	@In
	public int inTimestep;

	@Description("The path to diffuse output.")
	@In
	public String pathToDiffuse;

	@Description("The path to atmospheric trasmittance pruduct output.")
	@In
	public String pathToTau;
	@Description("The path to atmospheric trasmittance pruduct output.")
	@In
	public double pEpsilonS=0.98;	
	
	

	@Description("The path to beam component output.")
	@In
	public String pathToBeam;

	@Description("The path to Longwave Downwelling output.")
	@In
	public String pathToLongwaveDownwelling;

	@Description("The path to Longwave Upwelling output.")
	@In
	public String pathToLongwaveUpwelling;

	@Description("The path to Longwave  output.")
	@In
	public String pathToLongwave;

	@Description("The path to top atm output.")
	@In
	public String pathToTopatm;

	@Description("Raster Mode=true; Vector Mode=F.")
	@In
	public boolean workWithRaster = true;

	@Description("Raster Mode=true; Vector Mode=F.")
	@In
	public boolean doLongwave = true;

	@Description("The calculation time.")
	@In
	public int cumulationTime = 0;

	@Description("Path to the temperature file in input.")
	@In
	public String inPathToTemp = null;

	@Description("Path to the temperature file in input.")
	@In
	public String inPathToClearness = null;

	@Description("Path to the humidity file in input.")
	@In
	public String inPathToHumidity = null;

	@Description("Path to the humidity file in input.")
	@In
	public String inPathToSoilTemp = null;

	@Description("parameters Cm o3")
	@In
	public double pCmO3 = -1;

	@Description("parameters Visibility [km]")
	@In
	public double pVisibility = -1;

	@Description("parameters FC")
	@In
	public double pFC = -1;

	@Description("The progress monitor.")
	@In
	public IJGTProgressMonitor pm = new LogProgressMonitor();

	@Description("The map of total insolation in case workWithRaster=T.")
	@Out
	public GridCoverage2D outIns;

	@Description("The soil albedo.")
	@In
	public double pAlphag = -1.0;

	@Description("The soil albedo.")
	@In
	public int pModeDownCLR = 0;

	@Description("The HashMap of direct radiation.")
	@Out
	HashMap<Integer, double[]> outHMdirect;

	@Description("The HashMap of direct radiation.")
	@In
	public double pC_up = 0;

	@Description("The HashMap of direct radiation.")
	@In
	public double pD_up = 0;

	@Description("The HashMap of diffuse radiation.")
	@Out
	public boolean doDaily = false;

	@Description("The HashMap of diffuse radiation.")
	@In
	public double pA_Cloud = 0;

	@Description("The HashMap of diffuse radiation.")
	@In
	public double pB_Cloud = 0;

	@Description("The HashMap of diffuse radiation.")
	@Out
	HashMap<Integer, double[]> outHMdiffuse;

	@Description("The matrix of top atmosphere radiation.")
	@Out
	HashMap<Integer, double[]> outHMtopatm;

	@Description("The matrix of top atmosphere radiation.")
	@Out
	HashMap<Integer, double[]> outHMlongwaveDownwelling;

	@Description("The matrix of top atmosphere radiation.")
	@Out
	HashMap<Integer, double[]> outHMlongwaveUpwelling;

	@Description("The measured file.")
	@In
	@Unit("")
	public String pPathtoMeas = null;

	@Description("Read the measured file.")
	@In
	@Unit("")
	public boolean pDoReadMeas = false;

	@Description("MeasuredData.")
	@Out
	public double[] outMeasured;

	@Description("Simulated.")
	@Out
	public double[] outSimulated;

	@In
	@Unit("-")
	public int pDimMeas = 0;

	@Description("The matrix of top atmosphere radiation.")
	@Out
	HashMap<Integer, double[]> outHMlongwave;

	private static final double pRH = 0.7;

	private static final double pLapse = -.0065;

	private static final double SOLARCTE = 1370.0;

	private static final double ATM = 1013.25;

	private double lambda;
	private HashMap<String, Double> attribute;
	private double delta;
	private int height = 0;
	private int width = 0;
	private double omega;
	private double vetgio[];
	private int contaore = 1;
	private WritableRaster resultstaoWR = null;
	private WritableRaster resultinsWR = null;
	private WritableRaster resultdiffWR = null;
	private WritableRaster pitWR = null;
	private int contastampe = 0;
	private double[] xStation;
	private double[] yStation;
	private double[] zStation;
	private double[] albedoStation;
	private static final double ConstBoltz = 5.670373 * Math.pow(10, -8);
	private int[] idStation;
	private int[] colnumvetVect;
	private int[] rownumvetVect;
	public int contaconta = 0;

	private HashMap<Integer, double[]> tempvalues;
	private HashMap<Integer, double[]> clearindexvalue;
	private HashMap<Integer, double[]> soiltempvalues;
	private HashMap<Integer, double[]> umivalues;
	private WritableRaster skyviewfactorWR = null;

	boolean init = true;

	boolean init2 = true;

	@Execute
	public void process() throws Exception { // transform the

		if (init) {
			outSimulated = new double[pDimMeas];
			outMeasured = new double[pDimMeas];

			init = false;

		}
		if (init2) {
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

		// extract some attributes of the map
		attribute = CoverageUtilities.getRegionParamsFromGridCoverage(inElev);
		double dx = attribute.get(CoverageUtilities.XRES);

		/*
		 * The models use only one value of the latitude. So I have decided to
		 * set it to the center of the raster. Extract the CRS of the
		 * GridCoverage and transform the value of a WGS84 latitude.
		 */
		CoordinateReferenceSystem sourceCRS = inElev
				.getCoordinateReferenceSystem2D();
		CoordinateReferenceSystem targetCRS = DefaultGeographicCRS.WGS84;

		DateTimeFormatter formatter = DateTimeFormat.forPattern(
				"yyyy-MM-dd HH:mm").withZone(DateTimeZone.UTC);
		DateTime startcurrentDatetime = formatter.parseDateTime(tStartDate);

		DateTime endcurrentDatetime = formatter.parseDateTime(tEndDate);

		long diff = (endcurrentDatetime.getMillis() - startcurrentDatetime
				.getMillis()) / (1000 * 60 * 60);
		DateTime array[] = new DateTime[(int) diff];
		if (doDaily == false) {
			for (int i = 0; i < array.length; i++) {
				array[i] = startcurrentDatetime.plusHours(i);
			}
		}
		if (doDaily == true) {
			for (int i = 0; i < array.length; i++) {
				array[i] = startcurrentDatetime.plusDays(i);
			}
		}
		if (workWithRaster == false) {
			List<Double> xStationList = new ArrayList<Double>();
			List<Double> yStationList = new ArrayList<Double>();
			List<Double> zStationList = new ArrayList<Double>();
			List<Double> iddList = new ArrayList<Double>();
			List<Double> albedoList = new ArrayList<Double>();

			/*
			 * counter for the number of station with measured value !=0.
			 */
			/*
			 * Store the station coordinates and measured data in the array.
			 */
			FeatureIterator<SimpleFeature> stationsIter = inStations.features();
			try {
				while (stationsIter.hasNext()) {
					SimpleFeature feature = stationsIter.next();
					double zzz = 0;
					int id = ((Number) feature.getAttribute(fStationsid))
							.intValue();
					if (fStationsAlb != null) {
						double albedo = ((Number) feature
								.getAttribute(fStationsAlb)).doubleValue();
						albedoList.add(albedo);
					}
					Coordinate coordinate = ((Geometry) feature
							.getDefaultGeometry()).getCentroid()
							.getCoordinate();
					xStationList.add(coordinate.x);
					yStationList.add(coordinate.y);
					zStationList.add(zzz);

					iddList.add((double) id);

				}
			} finally {
				stationsIter.close();
			}

			int nStaz = xStationList.size();
			/*
			 * The coordinates of the station points plus in last position a
			 * place for the coordinate of the point to interpolate.
			 */
			xStation = new double[nStaz];
			yStation = new double[nStaz];
			zStation = new double[nStaz];
			idStation = new int[nStaz];
			colnumvetVect = new int[nStaz];
			rownumvetVect = new int[nStaz];
			albedoStation = new double[nStaz];
			if (nStaz != 0) {
				for (int i = 0; i < nStaz; i++) {
					double xTmp = xStationList.get(i);
					double yTmp = yStationList.get(i);
					double zTmp = zStationList.get(i);
					int idTmp = iddList.get(i).intValue();
					xStation[i] = xTmp;
					yStation[i] = yTmp;
					zStation[i] = zTmp;
					idStation[i] = idTmp;
					if (albedoList.size() != 0l) {
						albedoStation[i] = albedoList.get(i);
					} else {
						if (!isNovalue(pAlphag)) {
							albedoStation[i] = pAlphag;
						} else {
							albedoStation[i] = 0.5;
						}
					}
				}
			}

			MathTransform transf = inElev.getGridGeometry().getCRSToGrid2D();
			for (int i = 0; i < xStation.length; i++) {

				DirectPosition point = new DirectPosition2D(sourceCRS,
						xStation[i], yStation[i]);
				DirectPosition gridPoint = transf.transform(point, null);

				colnumvetVect[i] = (int) gridPoint.getCoordinate()[0];
				rownumvetVect[i] = (int) gridPoint.getCoordinate()[1];
				System.out.println(idStation[i] + "  "
						+ gridPoint.getCoordinate()[0] + " "
						+ gridPoint.getCoordinate()[1]);
			}
		}

		double srcPts[] = new double[] { attribute.get(CoverageUtilities.EAST),
				attribute.get(CoverageUtilities.SOUTH) };

		Coordinate source = new Coordinate(srcPts[0], srcPts[1]);
		Point[] so = new Point[] { GeometryUtilities.gf().createPoint(source) };
		CrsUtilities.reproject(sourceCRS, targetCRS, so);
		// the latitude value
		lambda = Math.toRadians(so[0].getY());
		CoverageUtilities.getRegionParamsFromGridCoverage(inElev);
		RenderedImage pitTmpRI = inElev.getRenderedImage();
		width = pitTmpRI.getWidth();
		height = pitTmpRI.getHeight();
		pitWR = CoverageUtilities.replaceNovalue(pitTmpRI, -9999.0);
		pitTmpRI = null;

		WritableRaster insolationWR = CoverageUtilities
				.createDoubleWritableRaster(width, height, null,
						pitWR.getSampleModel(), 0.0);

		WritableRaster staoWR = CoverageUtilities.createDoubleWritableRaster(
				width, height, null, pitWR.getSampleModel(), 0.0);

		WritableRaster diffuseWR = CoverageUtilities
				.createDoubleWritableRaster(width, height, null,
						pitWR.getSampleModel(), 0.0);
		WritableRaster gradientWR = normalVector(pitWR, dx);

		//
		vetgio = new double[array.length];

		if (workWithRaster) {
			resultstaoWR = CoverageUtilities.createDoubleWritableRaster(width,
					height, null, pitWR.getSampleModel(), 0.0);
			resultinsWR = CoverageUtilities.createDoubleWritableRaster(width,
					height, null, pitWR.getSampleModel(), 0.0);
			resultdiffWR = CoverageUtilities.createDoubleWritableRaster(width,
					height, null, pitWR.getSampleModel(), 0.0);

			for (int i = 0; i < array.length; i++) {
				DateTime currentime = array[i];
				calcInsolation(lambda, pitWR, gradientWR, insolationWR, staoWR,
						diffuseWR, dx, currentime);

				// pm.worked(i - startDay);
			}
		}
		if (workWithRaster == false) {

			OmsTimeSeriesIteratorReader reader_temp = new OmsTimeSeriesIteratorReader();
			if (!((inPathToTemp == null))) {
				reader_temp.file = inPathToTemp;
				reader_temp.idfield = "ID";
				reader_temp.tStart = tStartDate;
				reader_temp.tEnd = tEndDate;
				reader_temp.fileNovalue = "-9999";
				reader_temp.tTimestep = inTimestep;
			}

			OmsTimeSeriesIteratorReader reader_umi = new OmsTimeSeriesIteratorReader();

			if (!(inPathToHumidity == null)) {

				reader_umi.file = inPathToHumidity;
				reader_umi.idfield = "ID";
				reader_umi.tStart = tStartDate;
				reader_umi.tEnd = tEndDate;
				reader_umi.fileNovalue = "-9999";
				reader_umi.tTimestep = inTimestep;

			}

			OmsTimeSeriesIteratorReader reader_soiltemp = new OmsTimeSeriesIteratorReader();

			if (!(inPathToSoilTemp == null)) {

				reader_soiltemp.file = inPathToSoilTemp;
				reader_soiltemp.idfield = "ID";
				reader_soiltemp.tStart = tStartDate;
				reader_soiltemp.tEnd = tEndDate;
				reader_soiltemp.fileNovalue = "-9999";
				reader_soiltemp.tTimestep = inTimestep;

			}

			OmsTimeSeriesIteratorReader reader_clear = new OmsTimeSeriesIteratorReader();

			if (!(inPathToClearness == null)) {

				reader_clear.file = inPathToClearness;
				reader_clear.idfield = "ID";
				reader_clear.tStart = tStartDate;
				reader_clear.tEnd = tEndDate;
				reader_clear.fileNovalue = "-9999";
				reader_clear.tTimestep = inTimestep;

			}

			OmsTimeSeriesIteratorWriter writer4 = new OmsTimeSeriesIteratorWriter();
			OmsTimeSeriesIteratorWriter writer5 = new OmsTimeSeriesIteratorWriter();
			OmsTimeSeriesIteratorWriter writer6 = new OmsTimeSeriesIteratorWriter();

			writer4.file = pathToLongwaveDownwelling;
			writer5.file = pathToLongwaveUpwelling;
			writer6.file = pathToLongwave;

			writer4.tStart = tStartDate;
			writer5.tStart = tStartDate;
			writer6.tStart = tStartDate;

			writer4.tTimestep = inTimestep;
			writer5.tTimestep = inTimestep;
			writer6.tTimestep = inTimestep;

			for (int i = 0; i < array.length; i++) {
				outHMlongwaveDownwelling = new HashMap<Integer, double[]>();
				outHMlongwaveUpwelling = new HashMap<Integer, double[]>();
				outHMlongwave = new HashMap<Integer, double[]>();

				System.out.println(" data=" + array[i]);
				DateTime currentime = array[i];
				if (!(inPathToTemp == null)) {
					reader_temp.nextRecord();
					tempvalues = reader_temp.outData;
				}

				if (!(inPathToClearness == null)) {
					reader_clear.nextRecord();
					clearindexvalue = reader_clear.outData;
				}

				if (!(inPathToHumidity == null)) {
					reader_umi.nextRecord();
					umivalues = reader_umi.outData;
				}

				if (!(inPathToSoilTemp == null)) {
					reader_soiltemp.nextRecord();
					soiltempvalues = reader_soiltemp.outData;
				}

				calcInsolation(lambda, pitWR, gradientWR, insolationWR, staoWR,
						diffuseWR, dx, currentime);

				writer4.inData = outHMlongwaveDownwelling;
				writer4.writeNextLine();
				writer5.inData = outHMlongwaveUpwelling;
				writer5.writeNextLine();
				writer6.inData = outHMlongwave;
				writer6.writeNextLine();

				// pm.worked(i - startDay);
			}
			if (pathToLongwaveDownwelling != null) {
				writer4.close();
			}
			if (pathToLongwaveUpwelling != null) {
				writer5.close();
			}
			if (pathToLongwave != null) {
				writer6.close();
			}
		}

	}

	/**
	 * Evaluate the radiation.
	 * 
	 * @param lambda
	 *            the latitude.
	 * @param demWR
	 *            the raster of elevation
	 * @param gradientWR
	 *            the raster of the gradient value of the dem.
	 * @param insolationWR
	 *            the wr where to store the result.
	 * @param the
	 *            day in the year.
	 * @throws Exception
	 * @paradx the resolutiono of the dem.
	 */
	private void calcInsolation(double lambda, WritableRaster demWR,
			WritableRaster gradientWR, WritableRaster insolationWR,
			WritableRaster staoWR, WritableRaster diffuseWR, double dx,
			DateTime time) throws Exception {

		int day = time.getDayOfYear();

		double dayangb = (360 / 365.25) * (day - 79.436);

		dayangb = Math.toRadians(dayangb);

		// Evaluate the declination of the sun.
		delta = getDeclination(dayangb);
		// Evaluate the radiation in this day.
		double ss = Math.acos(-Math.tan(delta) * Math.tan(lambda));
		double sunrise = 12 * (1.0 - ss / Math.PI);
		double sunset = 12 * (1.0 + ss / Math.PI);
		double hhh = (double) time.getMillisOfDay() / (1000 * 60 * 60);

		double tetaj = 2 * Math.PI * (day - 1.0) / 365.0;
		double ecccorr = 1.00011 + 0.034221 * Math.cos(tetaj) + 0.00128
				* Math.sin(tetaj) + 0.000719 * Math.cos(2 * tetaj) + 0.000077
				* Math.sin(2 * tetaj);

		double hourangle = (hhh / 12.0 - 1.0) * Math.PI;

		// CoverageUtilities.getRegionParamsFromGridCoverage(inskyview);
		// RenderedImage SkyTmpRI = inskyview.getRenderedImage();
		// int width = SkyTmpRI.getWidth();
		// int height = SkyTmpRI.getHeight();

		if (workWithRaster) {

			for (int j = 0; j < height; j++) {
				for (int i = 0; i < width; i++) {

					if (hhh <= (sunrise)) {

						insolationWR.setSample(i, j, 0, 0);
						staoWR.setSample(i, j, 0, 0);
						diffuseWR.setSample(i, j, 0, 0);

					}
					if (hhh >= (sunset)) {
						// System.out.println(0);
						vetgio[contaore] = 0;
						insolationWR.setSample(i, j, 0, 0);
						staoWR.setSample(i, j, 0, 0);
						diffuseWR.setSample(i, j, 0, 0);

					}
					if (hhh > (sunrise) && hhh < (sunset)) {
						if (demWR.getSampleDouble(i, j, 0) != -9999.0) {
							omega = hourangle;
							// calculating the vector related to the sun
							double sunVector[] = calcSunVector();
							double zenith = calcZenith(sunVector[2]);

							double[] inverseSunVector = calcInverseSunVector(sunVector);
							double[] normalSunVector = calcNormalSunVector(sunVector);

							height = demWR.getHeight();
							width = demWR.getWidth();
							WritableRaster sOmbraWR = calculateFactor(height,
									width, sunVector, inverseSunVector,
									normalSunVector, demWR, dx);

							double mr = 1 / (sunVector[2] + 0.15 * Math.pow(
									(93.885 - (zenith * 180 / 3.14)), (-1.253)));

							double temp = doubleNovalue;
							double umi = doubleNovalue;
							// evaluate the radiation.
							// double[] aaa = calcRadiation(i, j, demWR,
							// sOmbraWR,
							// insolationWR, sunVector, gradientWR, mr,
							// ecccorr, temp, umi);
							// vetgio[contaore] = aaa * 0.0864 / 24;
							// System.out.println(aaa);
							// staoWR.setSample(i, j, 0, aaa[0] * 0.0864 / 24);
							// insolationWR.setSample(i, j, 0,
							// aaa[1] * 0.0864 / 24);
							// diffuseWR.setSample(i, j, 0, aaa[2] * 0.0864 /
							// 24);

						}

						else {
							staoWR.setSample(i, j, 0, -9999);
							insolationWR.setSample(i, j, 0, -9999);
							diffuseWR.setSample(i, j, 0, -9999);

						}
					}
				}
			}
			//
		} else {

			for (int contastaz = 0; contastaz < xStation.length; contastaz++) {
				System.out.println("STAZIONE N. ========> " + contastaz);
				int colnuber = colnumvetVect[contastaz];
				int rownumber = rownumvetVect[contastaz];
				int i = colnuber;
				int j = rownumber;
				int id = idStation[contastaz];

				// evaluate the radiation.
				double temperatura = doubleNovalue;
				if (tempvalues != null) {
					temperatura = tempvalues.get(id)[0];
				}

				double temperaturaSuolo = doubleNovalue;
				if (soiltempvalues != null) {
					temperaturaSuolo = soiltempvalues.get(id)[0];
				} else {
					temperaturaSuolo = temperatura;
				}

				double clearness = doubleNovalue;
				if (clearindexvalue != null) {
					clearness = clearindexvalue.get(id)[0];
				}

				double umidita = doubleNovalue;
				if (umivalues != null) {
					umidita = umivalues.get(id)[0];
				}
				System.out.println("ID===" + id);
				double[] aaa = calcLongwaveRadiation(i, j, demWR, insolationWR,
						gradientWR, ecccorr, temperatura, umidita, clearness,
						temperaturaSuolo);
				outSimulated[contaconta] = aaa[1];
				contaconta++;
				outHMlongwaveUpwelling.put(id, new double[] { aaa[2] });
				outHMlongwaveDownwelling.put(id, new double[] { aaa[1] });
				outHMlongwave.put(id, new double[] { aaa[2] + aaa[1] });

			}
		}

		//
		System.out.println("========contaora========= " + contaore);
		contaore += 1;
		if (workWithRaster) {
			printmap(insolationWR, staoWR, diffuseWR, cumulationTime, contaore);
		}
		// System.out.println("___________________");
	}

	public void printmap(WritableRaster rasterIsolation,
			WritableRaster rasterstao, WritableRaster rasterdiffuse,
			int cumulata, int ore) throws Exception {

		if (ore <= cumulata) {
			// System.out.println(ore+"  "+contaore);
			for (int j = 0; j < height; j++) {
				for (int i = 0; i < width; i++) {
					// evaluate the radiation.
					double valueins = rasterIsolation.getSampleDouble(i, j, 0)
							+ resultinsWR.getSampleDouble(i, j, 0);
					double valuestau = rasterstao.getSampleDouble(i, j, 0)
							+ resultstaoWR.getSampleDouble(i, j, 0);
					double valuediff = rasterdiffuse.getSampleDouble(i, j, 0)
							+ resultdiffWR.getSampleDouble(i, j, 0);
					// System.out.println(i+"  "+j+"  "+value);
					// System.out.println("somma=" + value);
					resultstaoWR.setSample(i, j, 0, valuestau);
					resultdiffWR.setSample(i, j, 0, valuediff);
					resultinsWR.setSample(i, j, 0, valueins);

				}
			}
		} else {

			for (int y = 2; y < height - 2; y++) {
				for (int x = 2; x < width - 2; x++) {
					if (pitWR.getSampleDouble(x, y, 0) == -9999.0) {
						resultdiffWR.setSample(x, y, 0, Double.NaN);
						resultinsWR.setSample(x, y, 0, Double.NaN);
						resultstaoWR.setSample(x, y, 0, Double.NaN);
					}

				}
			}
			contaore = 0;

			contastampe += 1;
			String resstau = "map_STAU" + contastampe;
			String percorsostau = "/Users/giuseppeformetta/Desktop/e/"
					+ resstau + ".asc";
			String resdiff = "map_DIFF" + contastampe;
			String percorsodiff = "/Users/giuseppeformetta/Desktop/e/"
					+ resdiff + ".asc";
			String resins = "map_INS" + contastampe;
			String percorsoins = "/Users/giuseppeformetta/Desktop/e/" + resins
					+ ".asc";

			GridCoverage2D coveragestau = CoverageUtilities.buildCoverage(
					"resultstao", resultstaoWR, attribute,
					inElev.getCoordinateReferenceSystem());
			GridCoverage2D coverageins = CoverageUtilities.buildCoverage(
					"resultins", resultinsWR, attribute,
					inElev.getCoordinateReferenceSystem());
			GridCoverage2D coveragediff = CoverageUtilities.buildCoverage(
					"resultdiff", resultdiffWR, attribute,
					inElev.getCoordinateReferenceSystem());

			OmsRasterWriter writer1 = new OmsRasterWriter();
			writer1.file = percorsostau;
			writer1.inRaster = coveragestau;
			writer1.process();

			OmsRasterWriter writer2 = new OmsRasterWriter();
			writer2.file = percorsoins;
			writer2.inRaster = coverageins;
			writer2.process();

			OmsRasterWriter writer3 = new OmsRasterWriter();
			writer3.file = percorsodiff;
			writer3.inRaster = coveragediff;
			writer3.process();

			resultdiffWR = null;
			resultdiffWR = CoverageUtilities.createDoubleWritableRaster(width,
					height, null, pitWR.getSampleModel(), 0.0);

			resultinsWR = null;
			resultinsWR = CoverageUtilities.createDoubleWritableRaster(width,
					height, null, pitWR.getSampleModel(), 0.0);

			resultstaoWR = null;
			resultstaoWR = CoverageUtilities.createDoubleWritableRaster(width,
					height, null, pitWR.getSampleModel(), 0.0);

		}

	}

	/*
	 * Evaluate the declination.
	 */
	private double getDeclination(double dayangb) {
		double delta = .3723 + 23.2567 * Math.sin(dayangb) - .758
				* Math.cos(dayangb) + .1149 * Math.sin(2 * dayangb) + .3656
				* Math.cos(2 * dayangb) - .1712 * Math.sin(3 * dayangb) + .0201
				* Math.cos(3 * dayangb);
		return Math.toRadians(delta);
	}

	/*
	 * evaluate several component of the radiation and then multiply by the
	 * sOmbra factor.
	 */

	protected double[] calcSunVector() {
		double sunVector[] = new double[3];
		sunVector[0] = -Math.sin(omega) * Math.cos(delta);
		sunVector[1] = Math.sin(lambda) * Math.cos(omega) * Math.cos(delta)
				- Math.cos(lambda) * Math.sin(delta);
		sunVector[2] = Math.cos(lambda) * Math.cos(omega) * Math.cos(delta)
				+ Math.sin(lambda) * Math.sin(delta);

		return sunVector;

	}

	protected WritableRaster normalVector(WritableRaster pitWR, double res) {

		int minX = pitWR.getMinX();
		int minY = pitWR.getMinY();
		int rows = pitWR.getHeight();
		int cols = pitWR.getWidth();

		RandomIter pitIter = RandomIterFactory.create(pitWR, null);
		/*
		 * Initializa the Image of the normal vector in the central point of the
		 * cells, which have 3 components so the Image have 3 bands..
		 */
		SampleModel sm = RasterFactory
				.createBandedSampleModel(5, cols, rows, 3);
		WritableRaster tmpNormalVectorWR = CoverageUtilities
				.createDoubleWritableRaster(cols, rows, null, sm, 0.0);
		WritableRandomIter tmpNormalIter = RandomIterFactory.createWritable(
				tmpNormalVectorWR, null);
		/*
		 * apply the corripio's formula (is the formula (3) in the article)
		 */
		for (int j = minY; j < minX + rows - 1; j++) {
			for (int i = minX; i < minX + cols - 1; i++) {
				double zij = pitIter.getSampleDouble(i, j, 0);
				double zidxj = pitIter.getSampleDouble(i + 1, j, 0);
				double zijdy = pitIter.getSampleDouble(i, j + 1, 0);
				double zidxjdy = pitIter.getSampleDouble(i + 1, j + 1, 0);
				double firstComponent = res * (zij - zidxj + zijdy - zidxjdy);
				double secondComponent = res * (zij + zidxj - zijdy - zidxjdy);
				double thirthComponent = 2 * (res * res);
				double den = Math.sqrt(firstComponent * firstComponent
						+ secondComponent * secondComponent + thirthComponent
						* thirthComponent);
				tmpNormalIter.setPixel(i, j, new double[] {
						firstComponent / den, secondComponent / den,
						thirthComponent / den });

			}
		}
		pitIter.done();

		return tmpNormalVectorWR;

	}

	private double calcZenith(double sunVector2) {
		return Math.acos(sunVector2);
	}

	private double[] calcLongwaveRadiation(int i, int j, WritableRaster demWR,
			WritableRaster insolationWR, WritableRaster gradientWR,
			double eccentricita, double tem, double umi, double clearness,
			double tsuolo) throws IOException {

		double risultato[] = new double[3];
		if (demWR.getSampleDouble(i, j, 0) != -9999.0) {
			double z = demWR.getSampleDouble(i, j, 0);
			double pressure = ATM * Math.exp(-0.0001184 * z);

			double temp = tem;
			double pRHH = umi;

			if (isNovalue(umi)) {
				pRHH = pRH;
			}

			if (isNovalue(tem)) {
				temp = pLapse * (z - 4000);
			}
			if (isNovalue(tsuolo)) {
				tsuolo = tem;
			}

			if (isNovalue(clearness)) {
				clearness = 1.0;
			}

			// humidity should be between 0-1
			// double es = 6.11 * Math.exp((17.3 * temp) / (237.3 + temp));
			// double e = pRHH * es;// [kPa]
			// transform temperature in K
			// temp = temp + 273.15;
			//
			// double vap_psat = Math.exp(26.23 - 5416.0 / temp);
			// double wPrec = 0.493 * (pRHH / 100) * vap_psat / temp;

			double downwellingclearsky = computedownwellingCS(temp, pRHH / 100,
					pModeDownCLR);

			double cloudnessIndex = 1 + pA_Cloud
					* Math.pow(clearness, pB_Cloud);

			double downwellingallsky = downwellingclearsky * cloudnessIndex;
			double t=temp+ 273.15;
			System.out.println("t= "+t);

			t=pC_up+pD_up*t;
			System.out.println("tmodified= "+t);

			double upwelling = pEpsilonS * ConstBoltz
					* Math.pow(t, 4);
			System.out.println("upwelling= "+upwelling);

			risultato[0] = downwellingclearsky;
			risultato[1] = downwellingallsky;
			risultato[2] = upwelling;

		}

		else {

			risultato[0] = -9999.0;
			risultato[1] = -9999.0;
			risultato[2] = -9999.0;

		}
		return risultato;

	}

	public double computedownwellingCS(double temperature, double humidity,
			int mode) {

		double res = 0;
		// es and e should be in hPa
		// humidity should be between 0-1
		double es = 6.11 * Math.exp((17.3 * temperature)
				/ (237.3 + temperature)) / 10;

		es = 6.11 * Math.pow(10, (7.5 * temperature) / (237.3 + temperature)) / 10;

		double e = humidity * es;
		// e KPa
		// transform temperature in K
		temperature = temperature + 273.15;

		// ALLSKY
		if (mode == 1) {
			// MODE 1: Angstrom 1918;The coefficients of the Angstrom scheme
			// used in
			// this study were fitted using summer data from Sodankyla, Finland,
			// 1997:

			// FLERCHINGER ET AL.: ATMOSPHERIC LONG-WAVE RADIATION ALGORITHMS

			double uno = Math.pow(10, pZs * e);

			// res = (0.83 - 0.18 * uno) * ConstBoltz * Math.pow(temperature,
			// 4);

			res = (pXs - pYs * uno) * ConstBoltz * Math.pow(temperature, 4);

		}

		if (mode == 2) {

			// Brunt�s 1932 scheme
			// cite
			// Estimation of daytime downward longwave radiation under
			// clear and cloudy skies conditions over a sub-humid region Facundo
			// Carmona & Raúl Rivas & Vicente Caselles

			// FLERCHINGER ET AL.: ATMOSPHERIC LONG-WAVE RADIATION ALGORITHMS

			double uno = Math.pow(e, 0.5);
			res = (pXs + pYs * uno) * ConstBoltz * Math.pow(temperature, 4);

		}

		if (mode == 3) {

			// Swinbank 1963 and FLERCHINGER ET AL.: ATMOSPHERIC LONG-WAVE
			// RADIATION ALGORITHMS
			res = (pXs * Math.pow(10, -13.0)) * Math.pow(temperature, 6);
		}

		if (mode == 4) {

			// Idso and Jackson 1969
			res = (1.0 - pXs
					* Math.exp((pYs * Math.pow(10, -4) * (273 - temperature) * (273 - temperature))))
					* ConstBoltz * Math.pow(temperature, 4);
		}
		if (mode == 5) {

			// Brutsaert 1975

			double uno = Math.pow(e / temperature, (1.0 / 7.0));
			res = pXs * uno * ConstBoltz * Math.pow(temperature, 4);
		}

		if (mode == 6) {

			// Idso 1981
			// FLERCHINGER ET AL.: ATMOSPHERIC LONG-WAVE RADIATION ALGORITHMS
			double uno = Math.exp(1500 / temperature);
			res = (pXs + pYs * Math.pow(10, -4) * e * uno) * ConstBoltz
					* Math.pow(temperature, 4);
		}

		if (mode == 7) {

			// Monteith and Unsworth (1990)

			res = (1.0 / (ConstBoltz * Math.pow(temperature, 4)))
					* (pXs + pYs * ConstBoltz * Math.pow(temperature, 4))
					* ConstBoltz * Math.pow(temperature, 4);
		}

		if (mode == 8) {

			// Konzelmann et al 1994

			double uno = Math.pow((e / 10) / temperature, (1.0 / 8.0));
			res = (pXs + pYs * uno) * ConstBoltz * Math.pow(temperature, 4);
		}

		if (mode == 9) {

			double wp = 4650 * e / temperature;
			res = (1.0 - (pXs + wp)
					* Math.exp((-Math.pow((pYs + pZs * wp), 0.5))))
					* ConstBoltz * Math.pow(temperature, 4);

		}

		if (mode == 10) {

			// Dilley and OBrien 1998
			double w = 4650 * e / temperature;
			res = (pXs + pYs * Math.pow(temperature / 273.16, 6) + pZs
					* Math.pow((w / 25), 0.5));
		}
		if (mode == 11) {
			res = ConstBoltz * Math.pow(temperature, 4)
					* (1 - pXs * Math.exp(-pYs * e / temperature));
			System.out.println("ciao");
		}

		if (mode == 12) {
			res = ConstBoltz * Math.pow(temperature, 4)
					* (pXs - pYs * Math.pow(e, pZs));
		}

		return res;
	}

}