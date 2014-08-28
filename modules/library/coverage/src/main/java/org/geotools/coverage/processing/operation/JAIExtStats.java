package org.geotools.coverage.processing.operation;

import it.geosolutions.jaiext.stats.Statistics;
import it.geosolutions.jaiext.stats.Statistics.StatsType;

import java.awt.Shape;
import java.awt.image.RenderedImage;
import java.awt.image.SampleModel;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Level;

import javax.media.jai.JAI;
import javax.media.jai.ParameterBlockJAI;
import javax.media.jai.RenderedOp;
import javax.media.jai.operator.ExtremaDescriptor;

import org.geotools.coverage.TypeMap;
import org.geotools.coverage.grid.GridCoverage2D;
import org.geotools.coverage.processing.BaseStatisticsOperationJAI;
import org.geotools.util.NumberRange;
import org.opengis.coverage.processing.OperationNotFoundException;
import org.opengis.parameter.ParameterValueGroup;
import org.opengis.referencing.crs.CoordinateReferenceSystem;
import org.opengis.referencing.operation.MathTransform;
import org.opengis.util.InternationalString;

/**
 * This operation simply wraps JAI-Ext Statistics operations described by {@link StatisticsDescriptor} inside a GeoTools operation in order to make it
 * spatial-aware. JAI-Ext Statistics supports filtering samples by ROI, and supports setting a "No Data" value range.
 * 
 * <p>
 * <strong>How to use this operation</strong> Here is a very simple example on how to use this operation to calculate the minimum and maximum, and a
 * histogram with 10 bins in the value range between 100 and 200 of the source coverage.
 * 
 * <code>
 * final OperationJAI op=new OperationJAI("JAI-EXT.stats");
 * ParameterValueGroup params = op.getParameters();
 * params.parameter("Source").setValue(coverage);
 * params.parameter("stats").setValue(new StatsType[] { StatsType.MIN, StatsType.MAX, StatsType.HISTOGRAM });
 * params.parameter("numBins").setValue(new int[] { 10 });
 * params.parameter("minBounds").setValue(new double[] { 100 });
 * params.parameter("maxBounds").setValue(new double[] { 200 });
 * coverage=(GridCoverage2D) op.doOperation(params,null);
 * Statistics[] minimumsPerBand = (Statistics[])coverage.getProperty(StatsType.MIN.name());
 * Statistics[] maximumsPerBand = (Statistics[])coverage.getProperty(StatsType.MAX.name());
 * Statistics[] histogramPerBand = (Statistics[])coverage.getProperty(StatsType.HISTOGRAM.name());
 * System.out.println("min band 0: " + (double)minimumsPerBand[0].getResult());
 * System.out.println("min band 1: " + (double)minimumsPerBand[1].getResult());
 * System.out.println("min band 2: " + (double)minimumsPerBand[2].getResult());
 * System.out.println("max band 0: " + (double)maximumsPerBand[0].getResult());
 * System.out.println("max band 1: " + (double)maximumsPerBand[1].getResult());
 * System.out.println("max band 2: " + (double)maximumsPerBand[2].getResult());
 * for (int band = 0; band < histogramPerBand.length; band++){
 *   double[] countsPerBin = (double[])histogramPerBand[band].getResult();
 *   for (int bin = 0; bin < countsPerBin.length; bin++){
 *     System.out.println("count for band " + band + " at bin " + bin + ": " + countsPerBin[bin]);
 *   } 
 * }
 * </code>
 * 
 * @author Rich Fecher
 * @since 13.0
 * 
 *
 *
 *
 * @source $URL$
 */
public class JAIExtStats extends BaseStatisticsOperationJAI {
    /** serialVersionUID */
    private static final long serialVersionUID = 556965726989407084L;

    /**
     * {@link String} key for getting the {@link it.geosolutions.jaiext.stats.Statistics} matrix.
     */
    public final static String GT_SYNTHETIC_PROPERTY_JAI_EXT_STATS = "jai_ext_stats";

    private static final int MAX_DEFAULT_NUM_BINS = 65536;

    /**
     * Default constructor for the {@link Histogram} operation.
     * 
     * @throws OperationNotFoundException
     */
    public JAIExtStats() throws OperationNotFoundException {
        super(getOperationDescriptor("Stats"));
    }

    /**
     * This operation MUST be performed on the geophysics data for this {@link GridCoverage2D}.
     * 
     * @param parameters {@link ParameterValueGroup} that describes this operation
     * @return always true.
     */
    protected boolean computeOnGeophysicsValues(ParameterValueGroup parameters) {
        return true;
    }

    /**
     * Prepare the {@link it.geosolutions.jaiext.stats.Statistics} property for this jai-ext stats operation.
     * <p>
     * See {@link it.geosolutions.jaiext.stats.StatisticsDescriptor} for more info.
     * 
     * @see OperationJAI#getProperties(RenderedImage, CoordinateReferenceSystem, InternationalString, MathTransform, GridCoverage2D[],
     *      org.geotools.coverage.processing.OperationJAI.Parameters),
     */
    protected Map<String, ?> getProperties(RenderedImage data, CoordinateReferenceSystem crs,
            InternationalString name, MathTransform toCRS, GridCoverage2D[] sources,
            Parameters parameters) {
        // /////////////////////////////////////////////////////////////////////
        //
        // If and only if data is a RenderedOp we prepare the properties for
        // stats as the output of the stats operation.
        //
        // /////////////////////////////////////////////////////////////////////
        if (data instanceof RenderedOp) {
            final RenderedOp result = (RenderedOp) data;

            // get the statistics
            final Statistics[][] statistics = (Statistics[][]) result
                    .getProperty(Statistics.STATS_PROPERTY);

            // return the map
            final Map<String, Object> synthProp = new HashMap<String, Object>();
            synthProp.put(GT_SYNTHETIC_PROPERTY_JAI_EXT_STATS, statistics);

            // also add an entry for each individual stat as a property for convenience
            int statsTypeIndex = parameters.parameters.indexOfParam("stats");
            if (statsTypeIndex >= 0 && statistics != null) {
                StatsType[] types = (StatsType[]) parameters.parameters.getParameters().get(
                        statsTypeIndex);
                for (int t = 0; t < types.length; t++) {
                    Statistics[] perBandStats = new Statistics[statistics.length];
                    for (int b = 0; b < statistics.length; b++) {
                        perBandStats[b] = statistics[b][t];
                    }
                    synthProp.put(types[t].name(), perBandStats);
                }
            }
            return Collections.unmodifiableMap(synthProp);

        }
        return super.getProperties(data, crs, name, toCRS, sources, parameters);
    }

    @Override
    protected ParameterBlockJAI prepareParameters(ParameterValueGroup parameters) {
        ParameterBlockJAI block = super.prepareParameters(parameters);
        Object o = parameters.parameter("roi").getValue();
        if (o != null && block.indexOfParam("roi") < 0) {
            block.setParameter("roi", o);
        }
        block.setParameter("useRoiAccessor", parameters.parameter("useRoiAccessor").getValue());
        block.setParameter("NoData", parameters.parameter("NoData").getValue());
        StatsType[] statsTypes = (StatsType[]) parameters.parameter("stats").getValue();
        if (statsTypes == null) {
            // stats type is required by JAI-ext, and it wouldn't be appropriate to make any assumptions or use a default
            throw new IllegalArgumentException(
                    "JAI-ext Statistics requires StatsType[] array with parameter name 'stats'");
        }
        double[] userMinBounds = (double[]) parameters.parameter("minBounds").getValue();
        double[] userMaxBounds = (double[]) parameters.parameter("maxBounds").getValue();
        int[] bands = (int[]) parameters.parameter("bands").getValue();
        int[] userNumBins = (int[]) parameters.parameter("numBins").getValue();
        double[] finalMinBounds = userMinBounds;
        double[] finalMaxBounds = userMaxBounds;
        int[] finalNumBins = userNumBins;
        // these are the required fields for complex stats, we can make a reasonable default based on the bounds of the sample model in the coverage,
        // if they are not set
        if (userMinBounds == null || userMaxBounds == null || userNumBins == null) {
            // check if its a complex stat, if it is it will require these
            // parameters
            boolean isSimpleStat = true;
            for (int i = 0; i < statsTypes.length; i++) {
                if (statsTypes[i].getStatsId() > 6) {
                    isSimpleStat = false;
                    break;
                }
            }
            if (!isSimpleStat) {
                finalMinBounds = new double[bands.length];
                finalMaxBounds = new double[bands.length];
                finalNumBins = new int[bands.length];
                // its a complex stat, we have to ensure minBounds, maxBounds, and
                // numBins are set because each is required
                final GridCoverage2D source = (GridCoverage2D) parameters.parameter(
                        operation.getSourceNames()[PRIMARY_SOURCE_INDEX]).getValue();
                SampleModel sampleModel = source.getRenderedImage().getSampleModel();
                for (int i = 0; i < bands.length; i++) {
                    int b = bands[i];
                    final NumberRange range = TypeMap.getRange(TypeMap.getSampleDimensionType(
                            sampleModel, b));
                    int bins;
                    double min;
                    double max;
                    if (userMinBounds == null) {
                        min = range.getMinimum(true);
                    } else {
                        min = userMinBounds[i];
                    }
                    if (userMaxBounds == null) {
                        max = range.getMaximum(true);
                    } else {
                        max = userMaxBounds[i];
                    }
                    if (userNumBins == null) {
                        if (Double.isInfinite(min) || Double.isInfinite(max) || Double.isNaN(min)
                                || Double.isNaN(max)) {
                            throw new IllegalArgumentException(
                                    "Attempt to infer number if bins has failed because the range for band "
                                            + b + " is infinite");
                        } else {
                            bins = (int) ((max - min) + 1);
                            if (bins > MAX_DEFAULT_NUM_BINS) {
                                if (LOGGER.isLoggable(Level.WARNING)) {
                                    LOGGER.warning("The range of values is too large (" + bins
                                            + "). Defaulting to " + MAX_DEFAULT_NUM_BINS + " bins.");
                                }
                                bins = MAX_DEFAULT_NUM_BINS;
                            }
                        }
                    } else {
                        bins = userNumBins[i];
                    }

                    finalMinBounds[i] = min;
                    finalMaxBounds[i] = max;
                    finalNumBins[i] = bins;
                }
            }
        }

        block.setParameter("minBounds", finalMinBounds);
        block.setParameter("maxBounds", finalMaxBounds);
        block.setParameter("numBins", finalNumBins);
        block.setParameter("stats", statsTypes);
        return block;
    }
}
