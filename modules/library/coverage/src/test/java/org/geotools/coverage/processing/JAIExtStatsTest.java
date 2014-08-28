package org.geotools.coverage.processing;

import static java.lang.Math.round;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.awt.geom.AffineTransform;
import java.awt.geom.NoninvertibleTransformException;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;

import javax.imageio.ImageIO;
import javax.media.jai.PlanarImage;
import javax.media.jai.ROI;
import javax.media.jai.ROIShape;

import it.geosolutions.imageioimpl.plugins.tiff.TIFFImageReader;
import it.geosolutions.imageioimpl.plugins.tiff.TIFFImageReaderSpi;
import it.geosolutions.jaiext.range.Range;
import it.geosolutions.jaiext.range.RangeFactory;
import it.geosolutions.jaiext.stats.Extrema;
import it.geosolutions.jaiext.stats.Statistics;
import it.geosolutions.jaiext.stats.Statistics.StatsType;

import org.geotools.TestData;
import org.geotools.coverage.CoverageFactoryFinder;
import org.geotools.coverage.GridSampleDimension;
import org.geotools.coverage.grid.GridCoverage2D;
import org.geotools.coverage.grid.GridEnvelope2D;
import org.geotools.coverage.grid.GridGeometry2D;
import org.geotools.coverage.grid.ViewType;
import org.geotools.data.DataStore;
import org.geotools.data.FeatureSource;
import org.geotools.data.FileDataStoreFinder;
import org.geotools.data.WorldFileReader;
import org.geotools.feature.FeatureCollection;
import org.geotools.feature.FeatureIterator;
import org.geotools.geometry.GeneralEnvelope;
import org.geotools.referencing.crs.DefaultGeographicCRS;
import org.geotools.referencing.operation.transform.ProjectiveTransform;
import org.geotools.util.logging.Logging;
import org.junit.Before;
import org.junit.Test;
import org.opengis.feature.simple.SimpleFeature;
import org.opengis.feature.simple.SimpleFeatureType;
import org.opengis.geometry.BoundingBox;
import org.opengis.metadata.spatial.PixelOrientation;
import org.opengis.parameter.ParameterValueGroup;
import org.opengis.referencing.operation.MathTransform;
import org.opengis.referencing.operation.TransformException;

import com.vividsolutions.jts.geom.CoordinateSequence;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.LineString;
import com.vividsolutions.jts.geom.MultiPolygon;
import com.vividsolutions.jts.geom.Polygon;

public class JAIExtStatsTest extends GridProcessingTestBase {
    private final static double DELTA = 10E-4;

    private final static Logger LOGGER = Logging.getLogger(JAIExtStatsTest.class);

    /**
     * The grid coverage to test.
     */
    private GridCoverage2D coverage;

    /**
     * Set up common objects used for all tests.
     */
    @Before
    public void setUp() throws Exception {
        coverage = EXAMPLES.get(0);
        TestData.unzipFile(this, "test.zip");
    }

    @Test
    public void testHistogram() throws Exception {
        GridCoverage2D source = coverage.view(ViewType.NATIVE);
        CoverageProcessor processor = CoverageProcessor.getInstance();

        ParameterValueGroup param = processor.getOperation("Stats").getParameters();
        param.parameter("Source").setValue(source);
        param.parameter("stats").setValue(new StatsType[] { StatsType.HISTOGRAM });

        GridCoverage2D processed = (GridCoverage2D) processor.doOperation(param);
        Statistics[][] stats = (Statistics[][]) processed.getProperty(Statistics.STATS_PROPERTY);
        assertNotNull(stats);
        double[] counts = ((double[]) stats[0][0].getResult());
        assertEquals(256, counts.length);
        assertEquals(0.0, counts[12], 0d);
        assertEquals(92.0, counts[13], 0d);
        assertEquals(763.0, counts[14], 0d);
        assertEquals(2281.0, counts[15], 0d);
        assertEquals(5401.0, counts[16], 0d);
        assertEquals(8280.0, counts[17], 0d);
        assertEquals(10253.0, counts[18], 0d);
        assertEquals(8254.0, counts[19], 0d);
        assertEquals(5682.0, counts[20], 0d);
        assertEquals(9869.0, counts[21], 0d);
        assertEquals(11010.0, counts[22], 0d);
        assertEquals(12304.0, counts[23], 0d);
        assertEquals(17020.0, counts[24], 0d);
        assertEquals(17307.0, counts[25], 0d);
        assertEquals(28534.0, counts[26], 0d);
        assertEquals(14165.0, counts[27], 0d);
        assertEquals(13878.0, counts[28], 0d);
        assertEquals(3159.0, counts[29], 0d);
        assertEquals(102.0, counts[30], 0d);
        assertEquals(0.0, counts[31], 0d);
    }

    @Test
    public void testHistogramWithNumBins() throws Exception {
        GridCoverage2D source = coverage.view(ViewType.NATIVE);
        CoverageProcessor processor = CoverageProcessor.getInstance();

        ParameterValueGroup param = processor.getOperation("Stats").getParameters();
        param.parameter("Source").setValue(source);
        param.parameter("stats").setValue(new StatsType[] { StatsType.HISTOGRAM });
        param.parameter("numBins").setValue(new int[] { 10 });

        GridCoverage2D processed = (GridCoverage2D) processor.doOperation(param);

        Statistics[][] stats = (Statistics[][]) processed.getProperty(Statistics.STATS_PROPERTY);
        assertNotNull(stats);
        double[] counts = ((double[]) stats[0][0].getResult());

        assertEquals(10, counts.length);

        assertEquals(99645.0, counts[0], 0d);
        assertEquals(68709.0, counts[1], 0d);
        assertEquals(0.0, counts[2], 0d);
    }

    @Test
    public void testHistogramWithBounds() throws Exception {
        GridCoverage2D source = coverage.view(ViewType.NATIVE);
        CoverageProcessor processor = CoverageProcessor.getInstance();

        ParameterValueGroup param = processor.getOperation("Stats").getParameters();
        param.parameter("Source").setValue(source);
        param.parameter("stats").setValue(new StatsType[] { StatsType.HISTOGRAM });
        param.parameter("numBins").setValue(new int[] { 10 });
        param.parameter("minBounds").setValue(new double[] { 100 });
        param.parameter("maxBounds").setValue(new double[] { 200 });
        GridCoverage2D processed = (GridCoverage2D) processor.doOperation(param);

        Statistics[][] stats = (Statistics[][]) processed.getProperty(Statistics.STATS_PROPERTY);
        assertNotNull(stats);
        double[] counts = ((double[]) stats[0][0].getResult());
        assertEquals(10, counts.length);
        assertEquals(0.0, counts[0], 0d);
        assertEquals(0.0, counts[9], 0d);
    }

    private static class MultiRange extends Range {
        private Range[] ranges;

        private Number max;

        private Number min;

        private boolean isPoint;

        public MultiRange(Range[] ranges) {
            super();
            this.ranges = ranges;
            isPoint = true;
            for (Range r : ranges) {
                if (!r.isPoint()) {
                    isPoint = false;
                }
                if (max == null || r.getMax().doubleValue() > max.doubleValue()) {
                    max = r.getMax();
                }
                if (min == null || r.getMin().doubleValue() < min.doubleValue()) {
                    min = r.getMin();
                }
            }
        }

        @Override
        public DataType getDataType() {
            return ranges[0].getDataType();
        }

        @Override
        public boolean isPoint() {
            return isPoint;
        }

        @Override
        public Number getMax() {
            return max;
        }

        @Override
        public Number getMin() {
            return min;
        }

        @Override
        public boolean contains(byte value) {
            for (Range r : ranges) {
                if (r.contains(value)) {
                    return true;
                }
            }
            return false;
        }

        @Override
        public boolean contains(short value) {
            for (Range r : ranges) {
                if (r.contains(value)) {
                    return true;
                }
            }
            return false;
        }

        @Override
        public boolean contains(int value) {
            for (Range r : ranges) {
                if (r.contains(value)) {
                    return true;
                }
            }
            return false;
        }

        @Override
        public boolean contains(float value) {
            for (Range r : ranges) {
                if (r.contains(value)) {
                    return true;
                }
            }
            return false;
        }

        @Override
        public boolean contains(double value) {
            for (Range r : ranges) {
                if (r.contains(value)) {
                    return true;
                }
            }
            return false;
        }

        @Override
        public boolean contains(long value) {
            for (Range r : ranges) {
                if (r.contains(value)) {
                    return true;
                }
            }
            return false;
        }

        @Override
        public <T extends Number & Comparable<T>> boolean contains(T value) {
            for (Range r : ranges) {
                if (r.contains(value)) {
                    return true;
                }
            }
            return false;
        }

    }

    private static class StatisticsTool {

        /*
         * external user params
         */
        private Set<StatsType> statisticsSet;

        private int[] bands;

        private GridCoverage2D gridCoverage2D;

        private List<SimpleFeature> featureList;

        private Range noDataRange;

        /*
         * results
         */
        private Map<String, Map<StatsType, Statistics[]>> feature2StatisticsMap = new HashMap<String, Map<StatsType, Statistics[]>>();

        private StatisticsTool(Set<StatsType> statisticsSet, GridCoverage2D gridCoverage2D,
                int[] bands, List<SimpleFeature> polygonsList, Range noDataRange) {
            this.statisticsSet = statisticsSet;
            this.gridCoverage2D = gridCoverage2D;
            this.bands = bands;
            this.featureList = polygonsList;
            this.noDataRange = noDataRange;
        }

        /**
         * Run the requested analysis.
         * 
         * <p>
         * This is the moment in which the analysis takes place. This method is intended to give the user the possibility to choose the moment in
         * which the workload is done.
         * 
         * @throws Exception
         */
        public void run() throws Exception {
            processPolygonMode();
        }

        private void processPolygonMode() throws TransformException {
            final AffineTransform gridToWorldTransformCorrected = new AffineTransform(
                    (AffineTransform) ((GridGeometry2D) gridCoverage2D.getGridGeometry())
                            .getGridToCRS2D(PixelOrientation.UPPER_LEFT));
            final MathTransform worldToGridTransform;
            try {
                worldToGridTransform = ProjectiveTransform.create(gridToWorldTransformCorrected
                        .createInverse());
            } catch (NoninvertibleTransformException e) {
                throw new IllegalArgumentException(e.getLocalizedMessage());
            }

            for (SimpleFeature feature : featureList) {
                final String fid = feature.getID();
                final Geometry geometry = (Geometry) feature.getDefaultGeometry();
                if (geometry instanceof Polygon || geometry instanceof MultiPolygon) {
                    final BoundingBox bbox = feature.getBounds();

                    /*
                     * crop on region of interest
                     */
                    final CoverageProcessor processor = CoverageProcessor.getInstance();
                    final ParameterValueGroup param = processor.getOperation("CoverageCrop")
                            .getParameters();
                    param.parameter("Source").setValue(gridCoverage2D);
                    param.parameter("Envelope").setValue(new GeneralEnvelope(bbox));
                    final GridCoverage2D cropped = (GridCoverage2D) processor.doOperation(param);
                    ROI roi = null;
                    final int numGeometries = geometry.getNumGeometries();
                    for (int i = 0; i < numGeometries; i++) {
                        Geometry geometryN = geometry.getGeometryN(i);
                        java.awt.Polygon awtPolygon = toAWTPolygon((Polygon) geometryN,
                                worldToGridTransform);
                        if (roi == null) {
                            roi = new ROIShape(awtPolygon);
                        } else {
                            ROI newRoi = new ROIShape(awtPolygon);
                            roi.add(newRoi);
                        }
                    }

                    StatsType[] statsTypes = statisticsSet.toArray(new StatsType[statisticsSet
                            .size()]);
                    final OperationJAI op = new OperationJAI("Stats");
                    ParameterValueGroup params = op.getParameters();
                    params.parameter("Source").setValue(cropped);
                    params.parameter("stats").setValue(statsTypes);
                    params.parameter("bands").setValue(bands);
                    params.parameter("roi").setValue(roi);
                    params.parameter("noData").setValue(noDataRange);
                    final GridCoverage2D coverage = (GridCoverage2D) op.doOperation(params, null);
                    final Statistics[][] allStats = (Statistics[][]) coverage
                            .getProperty(Statistics.STATS_PROPERTY);
                    final Map<StatsType, Statistics[]> statsMap = new HashMap<StatsType, Statistics[]>();
                    for (int i = 0; i < statsTypes.length; i++) {
                        Statistics[] perBandStats = new Statistics[bands.length];
                        for (int b = 0; b < bands.length; b++) {
                            perBandStats[b] = allStats[b][i];
                        }
                        statsMap.put(statsTypes[i], perBandStats);
                    }
                    feature2StatisticsMap.put(fid, statsMap);
                }
            }

        }

        private java.awt.Polygon toAWTPolygon(final Polygon roiInput,
                MathTransform worldToGridTransform) throws TransformException {
            final boolean isIdentity = worldToGridTransform.isIdentity();
            final java.awt.Polygon retValue = new java.awt.Polygon();
            final double coords[] = new double[2];
            final LineString exteriorRing = roiInput.getExteriorRing();
            final CoordinateSequence exteriorRingCS = exteriorRing.getCoordinateSequence();
            final int numCoords = exteriorRingCS.size();
            for (int i = 0; i < numCoords; i++) {
                coords[0] = exteriorRingCS.getX(i);
                coords[1] = exteriorRingCS.getY(i);
                if (!isIdentity)
                    worldToGridTransform.transform(coords, 0, coords, 0, 1);
                retValue.addPoint((int) round(coords[0] + 0.5d), (int) round(coords[1] + 0.5d));
            }
            return retValue;
        }

        /**
         * Gets the performed statistics.
         * 
         * @param fId the id of the feature used as region for the analysis.
         * @return the results of the analysis for all the requested {@link Statistic} for the requested bands. Note that the result contains for
         *         every {@link Statistic} a result value for every band.
         */
        public Map<StatsType, Statistics[]> getStatistics(String fId) {
            return feature2StatisticsMap.get(fId);
        }
    }

    @Test
    public void testPolygonStats() throws Exception {
        // by combining the inclusion ranges into a multi- nodata range, this is comparable to the global ZonalStats test
        final File tiff = TestData.file(this, "test.tif");
        final File tfw = TestData.file(this, "test.tfw");

        final TIFFImageReader reader = (it.geosolutions.imageioimpl.plugins.tiff.TIFFImageReader) new TIFFImageReaderSpi()
                .createReaderInstance();
        reader.setInput(ImageIO.createImageInputStream(tiff));
        final BufferedImage image = reader.read(0);
        reader.dispose();

        final MathTransform transform = new WorldFileReader(tfw).getTransform();
        final GridCoverage2D coverage2D = CoverageFactoryFinder.getGridCoverageFactory(null)
                .create("coverage",
                        image,
                        new GridGeometry2D(new GridEnvelope2D(PlanarImage.wrapRenderedImage(image)
                                .getBounds()), transform, DefaultGeographicCRS.WGS84),
                        new GridSampleDimension[] { new GridSampleDimension("coverage") }, null,
                        null);

        final File fileshp = TestData.file(this, "testpolygon.shp");
        final DataStore store = FileDataStoreFinder.getDataStore(fileshp.toURI().toURL());
        FeatureSource<SimpleFeatureType, SimpleFeature> featureSource = store
                .getFeatureSource(store.getNames().get(0));
        FeatureCollection<SimpleFeatureType, SimpleFeature> featureCollection = featureSource
                .getFeatures();
        List<SimpleFeature> polygonList = new ArrayList<SimpleFeature>();
        FeatureIterator<SimpleFeature> featureIterator = featureCollection.features();
        if (featureIterator.hasNext()) {
            SimpleFeature feature = featureIterator.next();
            polygonList.add(feature);
        }
        featureIterator.close();

        // choose the stats
        Set<StatsType> statsSet = new LinkedHashSet<StatsType>();
        statsSet.add(StatsType.MIN);
        statsSet.add(StatsType.MAX);
        statsSet.add(StatsType.MEAN);
        statsSet.add(StatsType.VARIANCE);
        statsSet.add(StatsType.DEV_STD);
        statsSet.add(StatsType.EXTREMA);

        // select the bands to work on
        int[] bands = new int[] { 0 };
        // no data ranges
        Range noDataRange = new MultiRange(
                new Range[] {
                        RangeFactory.create(Float.valueOf(1300), false, Float.valueOf(1370), false,
                                true),
                        RangeFactory.create(Float.valueOf(1600), false, Float.POSITIVE_INFINITY,
                                true, true) });

        // create the proper instance
        StatisticsTool statisticsTool = new StatisticsTool(statsSet, coverage2D, bands,
                polygonList, noDataRange);

        // do analysis
        statisticsTool.run();

        // get the results
        String id = "testpolygon.1";
        Map<StatsType, Statistics[]> statistics = statisticsTool.getStatistics(id);
        Extrema extrema = (Extrema) statistics.get(StatsType.EXTREMA)[0];
        double[] extremaResult = (double[]) extrema.getResult();
        assertEquals(extremaResult[1] - extremaResult[0], 343.0, DELTA);
        assertEquals((double) statistics.get(StatsType.DEV_STD)[0].getResult(), 88.7358, DELTA);
        assertEquals((double) statistics.get(StatsType.MIN)[0].getResult(), 1255.0, DELTA);
        assertEquals((double) statistics.get(StatsType.MEAN)[0].getResult(), 1380.5423, DELTA);
        assertEquals((double) statistics.get(StatsType.VARIANCE)[0].getResult(), 7874.0598, DELTA);
        assertEquals((double) statistics.get(StatsType.MAX)[0].getResult(), 1598.0, DELTA);

        id = "testpolygon.2";
        statistics = statisticsTool.getStatistics(id);
        LOGGER.info(id + statistics.toString());

        extrema = (Extrema) statistics.get(StatsType.EXTREMA)[0];
        extremaResult = (double[]) extrema.getResult();

        assertEquals(extremaResult[1] - extremaResult[0], 216.0, DELTA);
        assertEquals((double) statistics.get(StatsType.DEV_STD)[0].getResult(), 36.7996, DELTA);
        assertEquals((double) statistics.get(StatsType.MIN)[0].getResult(), 1192.0, DELTA);
        assertEquals((double) statistics.get(StatsType.MEAN)[0].getResult(), 1248.3870, DELTA);
        assertEquals((double) statistics.get(StatsType.VARIANCE)[0].getResult(), 1354.2150, DELTA);
        assertEquals((double) statistics.get(StatsType.MAX)[0].getResult(), 1408.0, DELTA);

        id = "testpolygon.3";
        statistics = statisticsTool.getStatistics(id);
        extrema = (Extrema) statistics.get(StatsType.EXTREMA)[0];
        extremaResult = (double[]) extrema.getResult();
        assertEquals(extremaResult[1] - extremaResult[0], 127.0000, DELTA);
        assertEquals((double) statistics.get(StatsType.DEV_STD)[0].getResult(), 30.9412, DELTA);
        assertEquals((double) statistics.get(StatsType.MIN)[0].getResult(), 1173.0, DELTA);
        assertEquals((double) statistics.get(StatsType.MEAN)[0].getResult(), 1266.3876, DELTA);
        assertEquals((double) statistics.get(StatsType.VARIANCE)[0].getResult(), 957.3594, DELTA);
        assertEquals((double) statistics.get(StatsType.MAX)[0].getResult(), 1300.0, DELTA);

        // change the no data range and re-run
        noDataRange = RangeFactory.create(Float.valueOf(1300), false, Float.POSITIVE_INFINITY,
                true, true);

        // create the proper instance
        statisticsTool = new StatisticsTool(statsSet, coverage2D, bands, polygonList, noDataRange);

        // do analysis
        statisticsTool.run();

        // get the results
        id = "testpolygon.1";
        statistics = statisticsTool.getStatistics(id);
        extrema = (Extrema) statistics.get(StatsType.EXTREMA)[0];
        extremaResult = (double[]) extrema.getResult();

        assertEquals(extremaResult[1] - extremaResult[0], 45.0, DELTA);
        assertEquals((double) statistics.get(StatsType.DEV_STD)[0].getResult(), 11.7972, DELTA);
        assertEquals((double) statistics.get(StatsType.MIN)[0].getResult(), 1255.0, DELTA);
        assertEquals((double) statistics.get(StatsType.MEAN)[0].getResult(), 1283.1634, DELTA);
        assertEquals((double) statistics.get(StatsType.VARIANCE)[0].getResult(), 139.1754, DELTA);
        assertEquals((double) statistics.get(StatsType.MAX)[0].getResult(), 1300.0, DELTA);

        // change the no data range once more
        noDataRange = new MultiRange(
                new Range[] {
                        RangeFactory.create(Float.valueOf(0), true, Float.valueOf(1370), false,
                                true),
                        RangeFactory.create(Float.valueOf(1600), false, Float.POSITIVE_INFINITY,
                                true, true) });
        statisticsTool = new StatisticsTool(statsSet, coverage2D, bands, polygonList, noDataRange);
        // do analysis
        statisticsTool.run();

        // get the results
        statistics = statisticsTool.getStatistics(id);
        extrema = (Extrema) statistics.get(StatsType.EXTREMA)[0];
        extremaResult = (double[]) extrema.getResult();
        assertEquals(extremaResult[1] - extremaResult[0], 228.0, DELTA);
        assertEquals((double) statistics.get(StatsType.DEV_STD)[0].getResult(), 63.7335, DELTA);
        assertEquals((double) statistics.get(StatsType.MIN)[0].getResult(), 1370.0, DELTA);
        assertEquals((double) statistics.get(StatsType.MEAN)[0].getResult(), 1433.8979, DELTA);
        assertEquals((double) statistics.get(StatsType.VARIANCE)[0].getResult(), 4061.9665, DELTA);
        assertEquals((double) statistics.get(StatsType.MAX)[0].getResult(), 1598.0, DELTA);

        reader.dispose();
        coverage2D.dispose(true);
        image.flush();
    }
}
