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
package org.jgrasstools.hortonmachine.models.hm;

import java.util.HashMap;

import org.geotools.coverage.grid.GridCoverage2D;
import org.geotools.data.simple.SimpleFeatureCollection;
import org.jgrasstools.gears.io.rasterreader.OmsRasterReader;
import org.jgrasstools.gears.io.shapefile.OmsShapefileFeatureReader;
import org.jgrasstools.gears.utils.HMTestCase;
import org.jgrasstools.gears.utils.coverage.CoverageUtilities;
import org.jgrasstools.hortonmachine.modules.hydrogeomorphology.lwrb.OmsLongwaveRadiationBalance;
import org.jgrasstools.hortonmachine.utils.HMTestMaps;
import org.opengis.referencing.crs.CoordinateReferenceSystem;

/**
 * Test ab.
 * 
 * @author Andrea Antonello (www.hydrologis.com)
 */
public class TestLWRBNODEM extends HMTestCase {

	public void testAb() throws Exception {

		OmsShapefileFeatureReader stationsReader = new OmsShapefileFeatureReader();
		stationsReader.file = "/Users/giuseppeformetta/Desktop/8/stazioni_tagliate.shp";
		stationsReader.readFeatureCollection();
		SimpleFeatureCollection hydrometersFC = stationsReader.geodata;

		OmsRasterReader reader = new OmsRasterReader();
		reader.file = "/Users/giuseppeformetta/Desktop/8/dtm1.asc";
		reader.fileNovalue = -9999.0;
		reader.geodataNovalue = Double.NaN;
		reader.process();
		GridCoverage2D readCoverage = reader.outRaster;

		OmsLongwaveRadiationBalance l = new OmsLongwaveRadiationBalance();
		// l.pPathtoMeas =
		// "/Users/giuseppeformetta/Desktop/MarialauraWork/LWRB/data/Downwellingm_CLEAR";
		// l.pDimMeas = 1138;
		// l.pDoReadMeas = true;
		l.inPathToHumidity = "/Users/giuseppeformetta/Desktop/8/H.csv";
		l.inPathToTemp = "/Users/giuseppeformetta/Desktop/8/Taria.csv";
		l.inPathToClearness = "/Users/giuseppeformetta/Desktop/8/ClearnessIndex.csv";
		//l.inPathToSoilTemp = "/Users/giuseppeformetta/Desktop/8/Tsuolo.csv";
		l.fStationsid = "cat";
		l.inStations = hydrometersFC;
		l.pathToLongwaveUpwelling = "/Users/giuseppeformetta/Desktop/8/UPWELLING.csv";
		l.inTimestep = 60;
		l.pModeDownCLR = 1;
		l.pA_Cloud = 0;
		l.pB_Cloud = 0;
		l.pXs = 0.963;
		l.pYs = 0.323;
		l.pZs = -0.247;
		l.pC_up = 1;
		l.pD_up = 1;
		// l.pYs = 0.91;
		l.inElev = readCoverage;
		l.workWithRaster = false;
		l.tStartDate = "2010-01-01 00:00";
		l.tEndDate = "2012-12-31 00:00";

		l.process();

	}

}
